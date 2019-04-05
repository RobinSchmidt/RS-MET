#pragma once


/** UNDER CONSTRUCTION

A generic implementation of the "BLEP" (= (B)and (L)imited st(EP)) technique for anti-aliasing
step discontinuities in a signal and/or its derivatives. Discontinuities in the signal itself show 
up as steps and produce spectra with a 6 dB/oct spectral roll-off, discontinuties in the first 
derivative show up as corners and produce spectra with a 12 dB/oct roll-off and so on. The basic 
idea is to take the difference between a bandlimited step (i.e. an integrated sinc function) and a 
naive step (this difference is called the residual) and add that residual to the output signal.

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

  /** 
  
  */
  void setLength(int newLength)
  {
    sincLength = newLength;
    bufferSize = rsNextPowerOfTwo(sincLength+1); 

    bufferSize *= 2; // just for development purposes - remove later

    mask       = bufferSize - 1;
    updateTables();
    allocateBuffers();
  }

  void setTablePrecision(int newValue)
  {
    samplesPerLobe = newValue;
    updateTables();
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getDelay() { return sincLength; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Adds the residual for a bandlimited impulse into our correction buffer. Call this right 
  before calling getSample, if your naive generator will generate an impulse at the next sample. */
  void addImpulse(TTim delayFraction, TSig amplitude)
  {
    TTim eps = RS_EPS(TTim);
    if(delayFraction < eps)
      return;
    // maybe we can get rid of this by having one extra sample in blitTbl, etc.

    int i;
    TTim frac = TTim(1) - delayFraction;
    int ic = sincLength;
    for(i = 0; i <  sincLength; i++) tempBuffer[ic+i] = blit(frac + i);
    for(i = 1; i <= sincLength; i++) tempBuffer[ic-i] = blit(frac - i);
    // ...maybe scale contents of tempBuffer to normalize the mean
    //rsStemPlot(tempBuffer); // ok - that looks good

    // apply correction to stored past samples:
    //rsStemPlot(delayline);
    for(int i = 0; i < sincLength; i++)
      delayline[wrap(bufIndex - i)] += amplitude * tempBuffer[ic-i-1];
    //rsStemPlot(delayline);


    // update corrector to be applied to future samples:
    for(int i = 0; i < sincLength; i++)
      corrector[wrap(bufIndex + i)] += amplitude * tempBuffer[ic+i];
    //rsStemPlot(corrector);
    corrector[bufIndex] -= amplitude;
    //rsStemPlot(corrector);

    int dummy = 0;


    // this is still wrong! we must add the left side of the correcttion signal to the delayed 
    // inputs and store the right side of the correction signal into the buffer to be applied to
    // future samples...maybe make a version that only applies the correction to future samples 
    // (minblep) - maybe use a (windowed) impulse-, step- and ramp- response of an elliptic filter
    // to window step- and ramp-response, subtract the naive versions, apply the window and add the
    // naive versions back

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
    delayline[bufIndex] = in + corrector[bufIndex];
    corrector[bufIndex] = 0;  // clear corrector at this position - it has been consumed
    TSig y = delayline[wrap(bufIndex-sincLength)];
    bufIndex = wrap(bufIndex + 1);
    return y;
  }

  /** Fills the delayline and blit/blep/blamp/etc. buffers with all zeros. */
  void reset();



protected:

  inline TSig blit(TTim time)
  {
    TTim p = (time + sincLength) * samplesPerLobe; // read position
    // + sincLength because our buffer index is shifted by that amount with respect to the time
    // origin - i.e. the center-index of blitTbl corresponds to time = 0

    int  i = (int) p;
    TTim f = p - i;
    return (1-f) * blitTbl[i] + f*blitTbl[i+1]; 
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






  int bufIndex = 0;

  int samplesPerLobe = 20;
  // Step-size to move from one sample to the next in the blit, blep, etc. buffers, i.e. the 
  // sample-rate at which the bandlimited functions are sampled/tabulated with respect to their 
  // zero corssings. A step-size of 1 means, they are sampled at the zero-crossings of the sinc, 
  // 2 means at zero-crossings and halfway in between (i.e. at the maxima of the underlying sine), 
  // etc.

  int sincLength;    // number of zero-crossings to the right of y-axis
  //int delayLength;   // 2*sincLength+1 - wrong! it's the same 
  //int delaySize;     // nextPowerOf2(delayLength)

  int bufferSize;
  int mask;          // bufferSize - 1
  //int centerIndex;


  std::vector<TTim> timeTbl, blitTbl, blitDrvTbl, blepTbl, blampTbl;
  // Tables for bandlimited impulse (windowed sinc), its derivative, first integral (bandlimited 
  // step) and second integral (bandlimited ramp)

  //std::vector<TTim> blepBuf, blampBuf; // buffer scaled blep/blamp values to be added to output
  // maybe we can use a single buffer for blep and blamp - just accumulate them both - call it
  // residualBuffer or correctionBuffer or something

  std::vector<TSig> corrector;  // buffer of correction samples to be added to delayed input
  std::vector<TSig> delayline;  // buffer for delayed input signal values
  std::vector<TSig> tempBuffer;

};