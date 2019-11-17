#pragma once

/** Baseclass for table-based blit/blep/blamp applicators. When a naive sound-synthesis algorithm,
such as a sawtooth or triangle-oscillator, can supply the information about when exactly (i.e. the
subsample position) the discontinuity occurred and the size of the discontinuity (or these two
numbers can somehow be estimated from the signal), objects of subclasses of this can anti-alias the
naive signal as a post-processing step. There are different variants (linear- and minimum phase)
that share some common code, that's why we have this baseclass. */

template<class TSig, class TTim>
class rsTableBlep // maybe rename to rsTableBlitBlepBlamp
{

public:

  rsTableBlep();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the coefficients of the different cosine terms in the window function that is applied to
  the sinc function. We use a cosine-sum window of that kind:
  https://en.wikipedia.org/wiki/Window_function#Cosine-sum_windows
  but note that here, we don't use the convention of alternating signs as wikipedia does, i.e. the
  a_k that you pass should include the minus sign, if the k-th cosine term is subtracted. The
  maximum number of allowed weights is some constant defined here in the class (currently 5) */
  void setWindowCosineWeights(TTim* newWeights, int numWeights)
  {
    rsAssert(numWeights <= maxNumWindowCoeffs);
    for(int i = 0; i < numWeights; i++)
      windowCoeffs[i] = newWeights[i];
    for(int i = numWeights; i < maxNumWindowCoeffs; i++)
      windowCoeffs[i] = TTim(0);
    updateTables();
  }

  /** Sets the number of samples that we take for each lobe of the sinc function in the lookup
  tables. */
  void setTablePrecision(int newValue)
  {
    tablePrecision = newValue;
    updateTables();
  }

  //void setRelativeCutoff
  // 1: fs/2, <1: below fs/2

protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Computes the value of the window function at a given time instant (in samples). Used to
  generate the windowed sinc table. */
  inline TTim window(TTim t)
  {
    TTim w = windowCoeffs[0];
    for(int i = 1; i < maxNumWindowCoeffs; i++)
      w += windowCoeffs[i] * cos(2 * i * PI * t / blepLength);
    return w;
  }
  // todo: maybe have functions that directly evaluate the windowed-sinc, the integrated windowed
  // sinc the twice integrated windowed sinc etc. - i think, with this cosine-sum window, analytic
  // integrals should be possible

  /** Wrap-around function for the index into our delayline and corrector. */
  inline int wrap(int i)
  {
    return i & mask;
  }

  /** Should fill our tables with blit, blep, blamp, etc. values. */
  virtual void updateTables() = 0;

  /** Should allocates the buffers for correction signal and delayed input signal. */
  virtual void allocateBuffers() = 0;


  static const int maxNumWindowCoeffs = 5;
  TTim windowCoeffs[maxNumWindowCoeffs];
  // Coefficients for the cosine terms in the cosine-sum window that is used to window the sinc or
  // residual

  int tablePrecision = 20;
  // Step-size to move from one sample to the next in the blit, blep, etc. tables, i.e. the
  // sample-rate at which the bandlimited functions are sampled/tabulated with respect to their
  // zero crossings. A step-size of 1 means, they are sampled at the zero-crossings of the sinc,
  // 2 means at zero-crossings and halfway in between (i.e. at the maxima of the underlying sine),
  // etc. - roughly the number of samples per lobe (applies exactly only to the linear phase blep)

  std::vector<TTim> blitTbl, blepTbl, blampTbl;    // tables of blit/blep/blamp
  // Tables for bandlimited impulse (windowed-sinc/impulse-response), first integral (bandlimited
  // step) and second integral (bandlimited ramp). For blep and blamp, the residuals are tabulated,
  // i.e. the respective response minus the ideal input step/ramp. But tabulating a residual would
  // not make sense for the blit, because its continuous time version is the infinitely narrow delta
  // distribution, so the blit is stored directly and the residual must be computed on the fly in
  // preparForImpulse (which is just a single multiplication-and-subtraction in this case)

  int bufIndex = 0;             // index in circular buffer(s) for delayline and correction
  int bufferSize;               // a power of 2
  int mask;                     // = bufferSize - 1, used for efficient wrap-around
  int blepLength = 30;          // number of samples for the blep

  std::vector<TSig> corrector;  // buffer of correction samples to be added to future inputs

};

//=================================================================================================

/** UNDER CONSTRUCTION


A generic implementation of the "BLEP" (= (B)and (L)imited st(EP)) technique for anti-aliasing
step discontinuities in a signal and/or its derivatives. Discontinuities in the signal itself show
up as steps and produce spectra with a 6 dB/oct spectral roll-off, discontinuties in the first
derivative show up as corners and produce spectra with a 12 dB/oct roll-off and so on. The basic
idea is to take the difference between a bandlimited step (i.e. an integrated sinc function) and a
naive step (this difference is called the residual) and add that residual to the output signal.

If a discontinuity happens between the previous sample n-1 and the current sample n, you should
announce that event to the object by calling prepareForImpulse/prepareForStep/prepareForCorner
immediately before the call to getSample at index n. The object will then prepare the appropriate
correction signals and apply them in subsequent calls to getSample. In order to be able to correct
past samples as well, a delay is introduced. When announcing a discontinuity, you need to tell the
object the fractional delay (how long ago is the instant of the discontinuity with respect to
sample index n in upcoming call to getSample - if the discontinuity occured at continuous time
n-frac, you should pass frac) and the size of the discontinuity.




make also rsTableMinBlit/rsTableMinBlep/rsTableMinBlamp
rsPolyLinBlep/rsPolyLinBlamp

...all thes variations should have the same interface (as far as possible) such that client code
can conveniently interchange them - maybe factor out a common baseclass

maybe name the file DiscontinuityBandLimiters and have various classes with the various variants
in a single file

References:


*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsTableLinBlep : public rsTableBlep<TSig, TTim>
{

public:

  rsTableLinBlep();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the length of the blit/blep/blamp in samples. Should be an even number - if you pass an
  odd number, the next even number will be used for the actual length. */
  void setLength(int newLength)
  {
    this->blepLength = rsNextEvenNumber(newLength);
    this->halfLength = this->blepLength/2;
    this->bufferSize = rsNextPowerOfTwo(halfLength+1); // why +1?
    //bufferSize = 2*rsNextPowerOfTwo(halfLength+1); // why +1?
    this->mask       = this->bufferSize - 1;
    updateTables();
    allocateBuffers();
  }
  // maybe make virtual and override

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the delay (in samples) that is introduced by the algorithm. May be used for delay
  compensation. */
  int getDelay() const { return halfLength; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Adds the residual for a bandlimited impulse into our correction buffer. Call this right
  before calling getSample, if you are going pass the impulse itself in your upcoming call to
  getSample. */
  void prepareForImpulse(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    if(shouldReturnEarly(delayFraction))
      return;
    fillTmpBuffer(delayFraction, amplitude, this->blitTbl);
    rsScale(tempBuffer, amplitude / rsSum(tempBuffer)); // sum of values should be "amplitude"
    applyTmpBuffer();
    this->corrector[this->bufIndex] -= amplitude;
    //rsStemPlot(tempBuffer);
    //rsStemPlot(delayline);
    //rsStemPlot(corrector);
  }
  // maybe rename to prepareForSpike

  /** Adds the residual for a bandlimited step into our correction buffer. Call this right
  before calling getSample, if your naive generator will generate a step at the next sample. */
  void prepareForStep(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    if(shouldReturnEarly(delayFraction))
      return;

    fillTmpBuffer(delayFraction, amplitude, this->blepTbl);
    rsScale(tempBuffer, amplitude); // get rid
    // for the step, we may want to scale the tempBuffer such that the mean is 0.5*amplitude?

    //rsAssert(rsMax(tempBuffer) <= 0.7); // for debug

    applyTmpBuffer();


    // this is needed when not the blep-residual but the blep itself is tabulated:
    for(int i = 0; i < halfLength; i++)
      this->corrector[this->wrap(this->bufIndex + i)] -= amplitude;


    //rsAssert(rsMax(delayline) <= 0.7); // for debug
    //rsAssert(rsMax(corrector) <= 0.7); // for debug

    //rsStemPlot(tempBuffer);
    //rsStemPlot(delayline);
    //rsStemPlot(corrector);

  }

  void prepareForCorner(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    if(shouldReturnEarly(delayFraction))
      return;
    fillTmpBuffer(delayFraction, amplitude, this->blampTbl);
    rsScale(tempBuffer, amplitude); // get rid
    applyTmpBuffer();

    //rsStemPlot(delayline);
    //rsStemPlot(corrector);
  }
  // or maybe rename to prepareForRamp


  /** Produces one sample at a time. */
  inline TSig getSample(TSig in)
  {
    delayline[this->wrap(this->bufIndex+halfLength)] = in + this->corrector[this->bufIndex];
    this->corrector[this->bufIndex] = TSig(0);  // corrector at this position has been consumed
    TSig y = delayline[this->bufIndex];
    this->bufIndex = this->wrap(this->bufIndex + 1);
    return y;
  }

  /** Fills the delayline and correction buffer with all zeros and resets the buffer index. */
  void reset();


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Figures out, if one of our prepareForImpulse/Step/Ramp etc. functions should return early
  (i.e. when some special limiting cases occur such that nothing needs to be done and just doing
  it anyway may cause access violations)*/
  inline bool shouldReturnEarly(TTim /*delayFraction*/)
  {
    //return delayFraction == 1;

    return false;  // i guess, we can get rid of that function completely now

    //TTim eps = RS_EPS(TTim);
    //return delayFraction >= TTim(1) || halfLength == 0;
    //return delayFraction < eps || delayFraction >= TTim(1) || sincLength == 0;
    // maybe we can get rid of this by having one extra sample in blitTbl, etc.
    // maybe when frac == 0.0 or frac == 1.0, one of the loops below should range from 0 to
    // sincLength-1 and the other from 0 to sincLength+1 instead of both ranging from 0 to
    // sincLength ...but that's a special case where the last samplemin the loop (at sincLength+1)
    // would evaluate to zero anyway, i think ...tests needed
  }
  // get rid of that!

  /** Reads out the given table (one of blitTbl, blepTbal, blampTbl, etc.) at the given time
  instant (in samples). */
  inline TSig readTable(TTim time, const std::vector<TTim>& tbl)
  {
    TTim p = (time + halfLength) * this->tablePrecision; // read position
    // + sincLength because our buffer index is shifted by that amount with respect to the time
    // origin - i.e. the center-index of blitTbl corresponds to time = 0

    int  i = (int) p;


    // test - should be done only for the blep-table because at this position, the table itself
    // makes a jump in value, leading to spurious peaks in the output signal:
    //if( i == halfLength*tablePrecision - 1 )
    //  return 0;
    // ok - this seems to help - but is it really the right thing to do? actually correct would be
    // to linearly interpolate the blep itself (not the blep-residual) and the subtract the ideal
    // step in the discrete time domain - but this is a lot of additional computation - how can we
    // have the best of the two approaches?



    TTim f = p - i;
    rsAssert(i+1 < (int) tbl.size());
    return (1-f) * tbl[i] + f*tbl[i+1];
  }

  /** Fills our temporary buffer with sampled values from the given table  (one of blitTbl,
  blepTbl, blampTbl, etc.). */
  inline void fillTmpBuffer(TTim delayFraction, TSig amplitude, const std::vector<TTim>& tbl)
  {
    int i;
    TTim frac = delayFraction;
    int ic = halfLength;
    for(i = 0; i <  halfLength; i++) tempBuffer[ic+i] = readTable(frac + i, tbl);
    for(i = 1; i <= halfLength; i++) tempBuffer[ic-i] = readTable(frac - i, tbl);
    // todo: optimize this - in each iteration of the loops, the "f" value in the called
    // readTable function is the same - we don't need to recompute it in each iteration
    // also, use a single loop instead of two

    //rsStemPlot(tempBuffer);

    // todo: use the "amplitude" to scale values returned from readTable - avoid doing that in
    // prepareFor...
  }

  /** Applies the content of the temporary buffer to the delayline and corrector. */
  inline void applyTmpBuffer()
  {
    int i;
    int ic = halfLength;

    // apply correction to stored past samples:
    for(i = 0; i < halfLength; i++) this->delayline[this->wrap(this->bufIndex + i)] += tempBuffer[i];

    // update corrector to be applied to future samples:
    for(i = 0; i < halfLength; i++) this->corrector[this->wrap(this->bufIndex + i)] += tempBuffer[ic+i];

    // ...may be done in a single loop
  }


  virtual void updateTables() override;
  virtual void allocateBuffers() override;

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  int halfLength;    // number of zero-crossings to the right of y-axis

  std::vector<TSig> delayline;  // buffer of corrected, delayed input signal values
  std::vector<TSig> tempBuffer; // temporary storage of sampled blit/blep/blamp values
  // actually, the separate tempBuffer is only needed for the blit and in this case actually only
  // for the normalization of the sum to unity - otherwise, we could render directly into the
  // corrector and the delayline - do it like that in production code (which probably doesn't even
  // need the blit option)

};

//=================================================================================================

/** Implementation of a minimum phase bandlimited step "MinBLEP" correction based on the
step-response of an elliptic lowpass filter. We tabulate the oversampled step response of the
elltiptic filter and use that as correction table. It may also generate MinBLITs and MinBLAMPs
(bandlimited impulses and ramps).

Maybe let the user also select the filter type (elliptic, Butterworth, Bessel, etc. ...Butterworth
can actually be obtained as special case of an elliptic filter by setting both ripple parameters
zero)..."butterblep"....sounds moderately funny :-)
*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsTableMinBlep : public rsTableBlep<TSig, TTim>
// maybe rename to rsTableFilterBlep and make a version based on an actual minimum-phase
// transform (i think, i already have code for that somewhere)...this is perhaps mainly of academic
// interest
{


public:


  rsTableMinBlep();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setLength(int newLength)
  {
    this->blepLength = newLength;
    this->bufferSize = rsNextPowerOfTwo(this->blepLength+1); // why +1?
    this->mask       = this->bufferSize - 1;
    allocateBuffers();
    updateTables();
  }
  // hmmm...in the linear phase-version, the setLength function actually sets the half-length -
  // make them consistent - in the linear phase version, add 1 if the user passes an odd number

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Unlike the linear phase variants, MinBLEPs don't introduce a fixed delay. We provide this
  function anyway for interface compatibility and return always zero.

  ...There actually is some frequency dependent delay introduced though...maybe we should have a
  frequency-parameter and return the actual delay? ...and if so, should it be the phase-delay or
  the group-delay? with the current elliptic filter, it has a delay of around 3.4 samples, if you
  consider the time when it actually reaches 1 in the step response (look at the plot of the
  blepTbl and divide the time (68) by the tablePrecision). */
  //int getDelay() { return 0; }
  int getDelay() const { return 3; }
  // value 3 obtained by eyeballing - maybe we should use the peak of the blit (impulse response)
  // ...which happens to be located at 3*tablePrecision - so 3 seems indeed the most reasonable
  // value


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  void prepareForImpulse(TTim delayFraction, TSig amplitude)
  {
    fillCorrector(delayFraction, amplitude, this->blitTbl);
    this->corrector[this->bufIndex] -= amplitude;
  }
  // not yet tested

  void prepareForStep(TTim delayFraction, TSig amplitude)
  {
    fillCorrector(delayFraction, amplitude, this->blepTbl);
  }
  // this seems to work

  void prepareForCorner(TTim delayFraction, TSig amplitude)
  {
    fillCorrector(delayFraction, amplitude, this->blampTbl);
  }
  // peak amplitude of triangle wave looks too large but spectrum looks good


  /** Produces one sample at a time. */
  inline TSig getSample(TSig in)
  {
    TSig out = in + this->corrector[this->bufIndex]; // our "read-and-delete head" consumes and...
    this->corrector[this->bufIndex] = TSig(0);       // ...then clears the correction signal at this position
    this->bufIndex = this->wrap(this->bufIndex + 1); // ...and then moves on in our circular correction buffer
    return out;
  }

  /** Fills the correction buffer with all zeros and resets the buffer index. */
  void reset();


protected:

  /** Accumulates correction samples into our "corrector" buffer from the given table
  (blitTbl, blepTbl or blampTbl) with phase according to the fractional delay and scaled by the
  amplitude. */
  inline void fillCorrector(TTim delayFraction, TSig amplitude, const std::vector<TTim>& tbl)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");

    //if(delayFraction == 1.0)  return;
    // we get an access violation in this case -> fix that


    TTim ps = delayFraction * this->tablePrecision; // start read-position in table
    int  i  = (int) ps;                       // integer part of table read start position
    TTim f  = ps - i;                         // fractional part of table read start position
    TTim w0 = amplitude * (1-f);              // weight for left table entry (at i)
    TTim w1 = amplitude * f;                  // weight for right table entry (at i+1)
    for(int k = 0; k < this->blepLength; k++) {
      rsAssert(i+1 < (int) tbl.size());
      this->corrector[this->wrap(this->bufIndex+k)] += w0 * tbl[i] + w1 * tbl[i+1];
      i += this->tablePrecision;
    }
    // todo: re-order the table such that instead of i += tablePrecision we can also do i++
    // -> optimize memory access


    //rsPlotVector(corrector);
  }


  virtual void updateTables() override;
  virtual void allocateBuffers() override;

};

//=================================================================================================

/** Implementation of polynomial bandlimited steps and ramps (polyBLEPs and polyBLAMPS) with one
sample delay, i.e. it corrects one sample before and one sample after the occurence of the
step/corner (which itself occurs somewhere in between these two samples).

References:
(1) Perceptually informed synthesis of bandlimited classical waveforms using integrated
polynomial interpolation

*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsPolyBlep1
{

public:


  rsPolyBlep1() { reset(); }

  int getDelay() const { return 1; }

  /** Returns the current value of the "delayed" sample. This is the value to be returned in the
  next call to  getSample(). It's the past input sample with future and past correction already
  applied. The future-correction was applied in getSample, the past-correction was applied in
  prepareForStep, etc. */
  inline TSig getDelayed() const { return delayed; }

  /** Returns the current content of the "corrector" sample. This is the value that will be added
  to the incoming sample passed to getSample before writing it into the delayed sample buffer
  variable. It represents the future correction. The past correction is not separately available
  but directly accumulated into the "delayed" sample. */
  inline TSig getCorrector() const { return corrector; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  //void prepareForImpulse(TTim delayFraction, TSig amplitude);

  void prepareForStep(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    TTim d  = delayFraction;
    TSig a  = amplitude;
    TTim d2 = d*d;
    delayed   += a *   d2/2;
    corrector += a * (-d2/2 + d - 1./2);
    // formulas from (1), table I
  }

  void prepareForCorner(TTim delayFraction, TSig amplitude)
  {
    rsError("Not yet implemented");
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    // unfortunately, the papers don't give the formulas for this case -> re-derive their formulas
    // and along with them also the formula for this case
  }

  /** Produces one sample at a time. */
  inline TSig getSample(TSig in)
  {
    TSig out  = delayed;
    delayed   = in + corrector;
    corrector = TSig(0);
    return out;
  }

  /** Resets the state of the object. */
  inline void reset()
  {
    delayed   = TSig(0);
    corrector = TSig(0);
  }


protected:

  TSig delayed;    // delayed input sample
  TSig corrector;  // corrector to be applied to next incoming sample

};

//=================================================================================================

/** Two sample delay version - uses 3rd order polynomial approximation...


References:
(1) Perceptually informed synthesis of bandlimited classical waveforms using integrated
    polynomial interpolation
(2) Rounding corners with BLAMP

*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsPolyBlep2
{

public:

  rsPolyBlep2() { reset(); }

  int getDelay() const { return 2; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  //void prepareForImpulse(TTim delayFraction, TSig amplitude);

  void prepareForStep(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    TSig a  = amplitude;
    TTim d  = delayFraction;
    TTim d2 = d*d;
    TTim d3 = d2*d;
    TTim d4 = d2*d2;
    delayed[1]   += a * ( d4/24);
    delayed[0]   += a * (-d4/8  + d3/6 + d2/4 +   d/6 + 1./24);
    corrector[0] += a * ( d4/8  - d3/3        + 2*d/3 - 1./2);
    corrector[1] += a * (-d4/24 + d3/6 - d2/4 +   d/6 - 1./24);
    // formulas from (1), table VII
    // this is the B-Spline ...add Lagrange, too

    // maybe use Horners rule to evaluate the quartic...or maybe even call
    // rsPolynomial::evaluateQuartic and have static const coeff arrays defined here
  }

  void prepareForCorner(TTim delayFraction, TSig amplitude)
  {
    rsAssert(delayFraction >= TTim(0) && delayFraction <= TTim(1), "Delay out of range");
    TSig a  = amplitude;
    TTim d  = delayFraction;
    TTim d2 = d*d;
    TTim d3 = d2*d;
    TTim d4 = d2*d2;
    TTim d5 = d4*d;
    delayed[1]   += a * ( d5/120);
    delayed[0]   += a * (-d5/40  + d4/24 + d3/12 + d2/12 + d/24 + 1./120);
    corrector[0] += a * ( d5/40  - d4/12         + d2/3  - d/2  + 7./30 );
    corrector[1] += a * (-d5/120 + d4/24 - d3/12 + d2/12 - d/24 + 1./120);
    // formulas from (2), table 1
  }
  // todo: optimize: precompute values that are used multiple times like:
  // d4_24 = (1./24) * d4; d3_12 = ...

  // have different versions: prepareForStepLagrange, prepareForStepBSpline - maybe dispatch in
  // the unspecified ones according to a user option - client code may then choose to bypass
  // the dispatcher by directly calling the desired version but we still conform to the general
  // blep/blamp interface that the table-based versions also use

  /** Produces one sample at a time. */
  inline TSig getSample(TSig in)
  {
    TSig out     = delayed[1];
    delayed[1]   = delayed[0];
    delayed[0]   = in + corrector[0];
    corrector[0] = corrector[1];
    corrector[1] = TSig(0);
    return out;
  }

  /** Resets the state of the object. */
  inline void reset()
  {
    delayed[0]   = TSig(0);
    delayed[1]   = TSig(0);
    corrector[0] = TSig(0);
    corrector[1] = TSig(0);
  }

protected:


  TSig delayed[2];    // two delayed input samples
  TSig corrector[2];  // corrector to be applied to next two incoming samples
  // maybe use variables y1,y2,c1,c2


};

// -maybe try some other polynomial approximations - maybe based on directly matching function
//  values and derivative values to the ideal analytic bandlimited discontinuities (based on the
//  Si function)
// -try to make a 3-sample delay version with the next higher order B-spline and Lagrange
//  interpolators (it may require to re-derive the formulas for the 1- and 2-sample version to see
//  the general pattern)...maybe try even higher orders...4 samples delay is still a quite small
//  number...but maybe the cost of the polynomial evaluation may become too expensive...the whole
//  thing will scale with d^2 where d is the delay in samples (we need to compute 2*d correction
//  samples and the degree of the polynomial also increases linearly with d) ...that's an advantage
//  of the table-based versions - they scales linearly with d - so maybe teh general reccomendation
//  should be: for (very) small bleps, use polynomials, otherwise tables
// -as for (table-based) MinBlep vs LinBlep:
//  -i think, the LinBlep can be overall shorter because the window is only unilateral
//  -the MinBlep has a fixed (3 sample) delay, independent from kernel length
//  -the LinBlep has less overshoot (but maybe with other filter-designs (Bessel), we can get rid
//   of the overshoot in the MinBlep, too)
// -maybe write a wrapper class that lets client code select the type of the blep (min/lin/poly) at
//  runtime and also conforms to the general blep interface - can be used during development to
//  experiment with various variants and later be swapped for a fixed choice (so it doesn't have to
//  be efficient - it's just for experimentation to help select the best blep for the particular
//  application)

