#pragma once


/** UNDER CONSTRUCTION:
-currently, only the blit-generation works - blep/blamp not yet implemented
-there's no window function applied yet (which boils down to ahving a recatngular window)

A generic implementation of the "BLEP" (= (B)and (L)imited st(EP)) technique for anti-aliasing
step discontinuities in a signal and/or its derivatives. Discontinuities in the signal itself show 
up as steps and produce spectra with a 6 dB/oct spectral roll-off, discontinuties in the first 
derivative show up as corners and produce spectra with a 12 dB/oct roll-off and so on. The basic 
idea is to take the difference between a bandlimited step (i.e. an integrated sinc function) and a 
naive step (this difference is called the residual) and add that residual to the output signal.

If a discontinuity happens between the previous sample n-1 and the current sample n, you should 
announce that event to the object by calling addImpulse/addStep/addRamp/etc. immediately before the
call to getSample at index n. The object will then prepare the appropriate correction signals and 
apply them in subsequent calls to getSample. In order to be able to correct past samples as well, a
delay is introduced. When announcing a discontinuity, you need to tell the object the fractional 
delay (how long ago is the instant of the discontinuity with respect to sample index n - if the 
discontinuity occured at n-frac, you should pass frac) and the size of the discontinuity. 


maybe rename to rsDiscontinuityBandLimiter - it also bandlimits impulses (i.e. implmennts BLITs)
....

References:


*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsStepBandLimiter
{

public:

  rsStepBandLimiter()
  {
    setLength(5);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the length of the underlying sinc function in terms of its numbers of zero crossings
  on the positive x-axis. With a value of 1, only the mainlobe is used, a value of 2 uses 
  additionally the first pair of sidelobes left and right of the mainlobe and so on. */
  void setLength(int newLength)
  {
    sincLength = newLength;
    bufferSize = rsNextPowerOfTwo(sincLength+1); 
    mask       = bufferSize - 1;
    updateTables();
    allocateBuffers();
  }

  /** Sets the number of samples that we take for each lobe of the sinc function in the lookup 
  tables. */
  void setTablePrecision(int newValue)
  {
    samplesPerLobe = newValue;
    updateTables();
  }

  /** Not yet used. Later these coefficients will be used to generate a window function defined by
  a sum over cosine terms.  */
  void setWindowCosineWeights(TTim* newWeights, int numWeights)
  {
    rsAssert(numWeights <= maxNumWindowCoeffs);
    for(int i = 0; i < numWeights; i++)
      windowCoeffs[i] = newWeights[i];
    for(int i = numWeights; i < maxNumWindowCoeffs; i++)
      windowCoeffs[i] = TTim(0);
  }

  //void setRelativeCutoff
  // 1: fs/2, <1: below fs/2

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the delay (in samples) that is introduced by the algorithm. May be used for delay 
  compensation. */
  int getDelay() { return sincLength; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Adds the residual for a bandlimited impulse into our correction buffer. Call this right 
  before calling getSample, if your naive generator will generate an impulse at the next sample. */
  void addImpulse(TTim delayFraction, TSig amplitude)
  {
    TTim eps = RS_EPS(TTim);
    if(delayFraction < eps || delayFraction >= TTim(1) || sincLength == 0 )
      return;
    // maybe we can get rid of this by having one extra sample in blitTbl, etc.
    // maybe when frac == 0.0 or frac == 1.0, one of the loops below should range from 0 to 
    // sincLength-1 and the other from 0 to sincLength+1 instead of both ranging from 0 to 
    // sincLength ...but that's a special case where the last samplemin the loop (at sincLength+1)
    // would evaluate to zero anyway, i think ...tests needed

    int i;
    TTim frac = delayFraction;
    int ic = sincLength;
    for(i = 0; i <  sincLength; i++) tempBuffer[ic+i] = blit(frac + i);
    for(i = 1; i <= sincLength; i++) tempBuffer[ic-i] = blit(frac - i);
    rsScale(tempBuffer, amplitude / rsSum(tempBuffer)); // sum of values should be "amplitude"
    //rsStemPlot(tempBuffer);
    // todo: optimize this - in each iteration of the loops, the "f" value in the called 
    // blit-function is the same - we don't need to recompute it in each iteration

    // apply correction to stored past samples:
    for(i = 0; i < sincLength; i++)
      delayline[wrap(bufIndex+i)] += tempBuffer[i];
    //rsStemPlot(delayline);

    // update correction to be applied to future samples:
    for(i = 0; i < sincLength; i++)
      corrector[wrap(bufIndex + i)] += tempBuffer[ic+i];
    corrector[bufIndex] -= amplitude;
    //rsStemPlot(corrector);
  }

  /** Adds the residual for a bandlimited step into our correction buffer. Call this right 
  before calling getSample, if your naive generator will generate a step at the next sample. */
  /*
  void addStep(TTim time, TSig amplitude)
  {
    for(int i = 0; i < delayLength; i++)
      corrector[wrap(bufIndex + i)] += amplitude * blepResidual(i, time);
  }
  */

  /** Adds the residual for a bandlimited ramp into our correction buffer. Call this right 
  before calling getSample, if your naive generator will generate a corner at the next sample. */
  /*
  void addRamp(TTim time, TSig amplitude)
  {
    for(int i = 0; i < delayLength; i++)
      corrector[wrap(bufIndex + i)] += amplitude * blampResidual(i, time);
  }
  */
  // maybe rename to addCorner

  // maybe have higher order bandlimited integrated impulses (qudaratic, cubic, quartic, etc.)
  // maybe store the intergrated impulses in a 2D array and unify these 3 functions (maybe keep the
  // separate functions as convenience functions anyway), like:
  /*
  void addDiscontinuity(Tim time, TSig amplitude, int order)
  {
    for(int i = 0; i < delayLength; i++)
      corrector[wrap(bufIndex + i)] += amplitude * residual(order, i, time);
  }
  void addImpulse(TTim time, TSig amplitude) { addDiscontinuity(time, amplitude, 0); }
  void addStep(   TTim time, TSig amplitude) { addDiscontinuity(time, amplitude, 1); }
  void addRamp(   TTim time, TSig amplitude) { addDiscontinuity(time, amplitude, 2); }
  */

  /** Produces one sample at a time. */
  inline TSig getSample(TSig in)
  {
    delayline[wrap(bufIndex+sincLength)] = in + corrector[bufIndex];
    corrector[bufIndex] = TSig(0);  // corrector at this position has been consumed
    TSig y = delayline[bufIndex];
    bufIndex = wrap(bufIndex + 1);
    return y;
  }

  /** Fills the delayline and correction buffer with all zeros and resets the buffer index. */
  void reset();



protected:

  inline TSig blit(TTim time)
  {
    return readTable(time, blitTbl);
    /*
    TTim p = (time + sincLength) * samplesPerLobe; // read position
    // + sincLength because our buffer index is shifted by that amount with respect to the time
    // origin - i.e. the center-index of blitTbl corresponds to time = 0

    int  i = (int) p;
    TTim f = p - i;
    return (1-f) * blitTbl[i] + f*blitTbl[i+1]; 
    */
  }

  inline TSig readTable(TTim time, const std::vector<TTim>& tbl)
  {
    TTim p = (time + sincLength) * samplesPerLobe; // read position
    // + sincLength because our buffer index is shifted by that amount with respect to the time
    // origin - i.e. the center-index of blitTbl corresponds to time = 0

    int  i = (int) p;
    TTim f = p - i;
    return (1-f) * tbl[i] + f*tbl[i+1]; 
  }

  


  /*
  inline TSig blitResidual(int i, TTim frac)
  {
    //int k = samplesPerLobe * i;
    return (1-frac) * blitTbl[i] + frac*blitTbl[i+1]; 
  }
  */

  /*
  inline TSig blepResidual(int i, TTim frac)
  {
    return (1-frac) * blepTbl[i] + frac*blepTbl[i+1]; 
    // preliminary - linear interpolation - later use hermite interpolation using the blit-values
    // for the derivative - but maybe precompute the polynomial coefficients
    // ...maybe, it should be the other way around (1-frac vs frac)...will dpend on which 
    // conventions we adopt for "frac" ....
  }

  inline TSig blampResidual(int i, TTim frac)
  {
    return (1-frac) * blampTbl[i] + frac*blampTbl[i+1]; 
  }
  */

  inline int wrap(int i)
  {
    return i & mask;
  }

  /** Fills our tables with blit, blep, blamp, etc. values. */
  void updateTables();

  /** Allocates the buffers for correction signal and delayed input signal. */
  void allocateBuffers();


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  int bufIndex = 0;  // index in circular buffers for delayline and correction
  int sincLength;    // number of zero-crossings to the right of y-axis
  int bufferSize;    // nextPowerOfTwo(sincLength+1); 
  int mask;          // bufferSize - 1

  int samplesPerLobe = 20;
  // Step-size to move from one sample to the next in the blit, blep, etc. buffers, i.e. the 
  // sample-rate at which the bandlimited functions are sampled/tabulated with respect to their 
  // zero corssings. A step-size of 1 means, they are sampled at the zero-crossings of the sinc, 
  // 2 means at zero-crossings and halfway in between (i.e. at the maxima of the underlying sine), 
  // etc.


  static const int maxNumWindowCoeffs = 5;
  TTim windowCoeffs[maxNumWindowCoeffs];  // not yet used


  std::vector<TTim> timeTbl, blitTbl, blitDrvTbl, blepTbl, blampTbl;
  // Tables for bandlimited impulse (windowed sinc), its derivative, first integral (bandlimited 
  // step) and second integral (bandlimited ramp)
  // blitDrvTbl not yet used - may not be needed (i was thinking to use it for Hermite 
  // interpolation of the blit table)
  // timeTbl may also not be used - and maybe we should use only the positive half of the table and
  // produce values for negative time instants by symmetry

  std::vector<TSig> corrector;  // buffer of correction samples to be added to future inputs
  std::vector<TSig> delayline;  // buffer of corrected, delayed input signal values
  std::vector<TSig> tempBuffer; // temporary storage of sampled blit/blep/blamp values

};