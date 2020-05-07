#ifndef RAPT_ONEPOLEFILTER_H_INCLUDED
#define RAPT_ONEPOLEFILTER_H_INCLUDED


/** Baseclass for first order (1-pole/1-zero) filters. Maintains only the filter state and 
coefficients but no user parameters like cutoff, sampleRate, mode, etc. which makes it lean to use
in contexts where it is undesirable for each benign first-order filter to maintain the sampleRate 
and other things itself. It also provides the bare coefficient computation formulas as static 
member functions. Beware that the formulas as implemented here assume a positive sign convention 
for feedback coeffs, i.e. they compute coeffs for the difference equation:
y[n] = b0*x[n] + b1*x[n-1] + a1*y[n-1] */

template<class TSig, class TPar>
class rsFirstOrderFilterBase
{

public:

  // for convenience:
  typedef const TSig& CRSig;  // const reference to a signal value
  //typedef const TPar& CRPar;  // const reference to a parameter value
  // we need to pass signalvalues as const reference because otherwise, we get a compiler error 
  // when instantiating the template with rsFloat64x2 for the signal type an compile for 32 bit
  // ("formal parameter with requested alignment of 16 won't be aligned")
  // ...when we want to instantiate it with rsFloat64x2 for TPar, too, we'll probably have to pass
  // all parameters as const reference, too - that really sucks!


  //-----------------------------------------------------------------------------------------------
  /** \name Coefficient Computation */

  // functions can be used from outside - maybe move into a class "FilterDesignFormulas" or
  // something. The formulas here use a positive sign convention for feedback-coeffs...make the 
  // usage of positive or negative sign-convention consistent throughout the library...
  // or maybe support both conventions in the design formulas

  /** Trivial "bypass" coeffs. */
  template<class T>
  static inline void coeffsBypass(T* b0, T* b1, T* a1)
  {
    *b0 = 1;
    *b1 = 0;
    *a1 = 0;
  }

  /** Lowpass via impulse invariant transform (from dspguide) */
  template<class T>
  static inline void coeffsLowpassIIT(T w, T* b0, T* b1, T* a1)
  {
    *a1 = exp(-w); 
    *b0 = 1 - *a1;
    *b1 = 0;
    // w < 0:   ...
    // w = 0:   b0 = 0, a1 = 1 -> silence or DC
    // w = pi:  nothing special
    // w = inf: b0 = 1, a1 = 0 -> bypass
  }
  // https://en.wikipedia.org/wiki/Impulse_invariance

  /** Highpass via matched Z-transform (from dspguide) */
  template<class T>
  static inline void coeffsHighpassMZT(T w, T* b0, T* b1, T* a1)
  {
    *a1 = exp(-w);
    *b0 =  T(0.5)*(1 + *a1);
    *b1 = -(*b0);
    // w = 0:   a1 = 1, b0 = 1, b1 = -1 -> differentiate and integrate
    // w = pi:  
    // w = inf: 
  }
  // https://en.wikipedia.org/wiki/Matched_Z-transform_method

  /** Allpass via bilinear transform (from DAFX) */
  template<class T>
  static inline void coeffsAllpassBLT(T w, T* b0, T* b1, T* a1)
  {
    T t = tan(T(0.5)*w); // tan w/2
    *b0 = (t-1) / (t+1);
    *b1 = 1.0;
    *a1 = -(*b0);
    // w = 0:
    // w = pi:  t = inf -> b0, a1 = NaN
    // w = inf: 
    // tan(0) = 0, tan(pi/2) = inf, tan(pi) = 0
  }
  // https://en.wikipedia.org/wiki/Bilinear_transform
  // todo: catch the special case |t| == inf by setting coeffs to the limiting values
  // ...maybe make "safe" versions of these functions coeffsAllpassSafeBLT or something

  /** Lowpass via bilinear transform (from DAFX) */
  template<class T>
  static inline void coeffsLowpassBLT(T w, T* b0, T* b1, T* a1)
  {
    T t = tan(T(0.5)*w);

    rsAssert(rsIsFiniteNumber(t));
    // hmmm...maybe we should allow t == inf and set a1 = -1 in this case (that's the limiting 
    // value)


    *a1 = (1-t) / (1+t);
    *b0 = T(0.5)*(1 - *a1);
    *b1 = *b0;
    // w = 0:  t = 0 -> a1 = 1, b0,b1 = 0
    // w = pi: t = inf -> a1,b0,b1 = NaN
  }

  /** Highpass via bilinear transform (from DAFX) */
  template<class T>
  static inline void coeffsHighpassBLT(T w, T* b0, T* b1, T* a1)
  {
    T t = tan(T(0.5)*w);     // maybe re-use w
    *a1 = (1-t) / (1+t);
    *b0 = T(0.5)*(1 + *a1);
    *b1 = -(*b0);
  }

  /** Low shelving via bilinear transform, g is linear shelving gain. (from DAFX) */
  template<class T>
  static inline void coeffsLowShelfBLT(T w, T g, T* b0, T* b1, T* a1)
  {
    T t = tan(T(0.5)*w);  // re-use w
    t   = g >= T(1.0) ? (t-1)/(t+1) : (t-g)/(t+g);
    T c = T(0.5)*(g-1);   // re-use g
    c  += c*t;
    *b0 = 1 + c;
    *b1 = t + c;
    *a1 = -t;
  }

  template<class T>
  static inline void coeffsHighShelfBLT(T w, T g, T* b0, T* b1, T* a1)
  {
    T t = tan(T(0.5)*w);
    t   = g >= T(1.0) ? (t-T(1.0))/(t+T(1.0)) : (g*t-1)/(g*t+1);
    T c = T(0.5)*(g-1);
    c  -= c*t;
    *b0 = 1 + c;
    *b1 = t - c;
    *a1 = -t;
  }
  // maybe get rid of some of the T(1.0) - replace them by 1, if the compiler doesn't warn (but try
  // with highest warnign level)

  //// todo - but these need the sample-rate as input - can we reformulate the equations such that
  //// we don't need the sample-rate...if w corresponds to some frequency in Hz then pi corresponds 
  //// to the Nyquist frequency...
  //template<class T>
  //static void coeffsLowShelfNMM(T w, T g, T* b0, T* b1, T* a1);

  //template<class T>
  //static void coeffsHighShelfNMM(T w, T g, T* b0, T* b1, T* a1);


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the filter coefficients manually. */
  inline void setCoefficients(TPar newB0, TPar newB1, TPar newA1)
  {
    b0 = newB0;
    b1 = newB1;
    a1 = newA1;
  }

  /** Sets up the internal state variables for both channels. */
  inline void setInternalState(CRSig newX1, CRSig newY1)
  {
    x1 = newX1;
    y1 = newY1;
  }
  // rename to setState

  // todo:
  // void setup(int mode, TPar omega, TPar gain = TPar(1));


    //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  TPar getB0() const { return b0; }

  TPar getB1() const { return b1; }

  TPar getA1() const { return a1; }
  // maybe rename to getCoeffB0, etc...

  TSig getStateX() const { return x1; }

  TSig getStateY() const { return y1; }




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates a single filtered output-sample. */
  inline TSig getSample(CRSig in)
  {
    y1 = b0*in + b1*x1 + a1*y1;
    x1 = in;
    return y1;
  }

  /** Resets the internal buffers (for the \f$ x[n-1], y[n-1] \f$-samples) to zero. */
  inline void reset()
  {
    x1 = TSig(0);
    y1 = TSig(0);
  }

  /** Sets the internal state to the state to which the filter will converge when seeing a 
  constant input of value c for a sufficiently long time. */
  inline void setStateForConstInput(TSig c)
  {
    x1 = c; 
    y1 = (b0+b1)*c / (TPar(1)-a1);
  }

  /** When this filter is used bidirectionally, you can call this function between the forward and 
  the backward pass to set the internal states of the filter appropriately. It simulates running 
  the filter for a ringout phase using an all-zeros signal as input and then using the reversed 
  ringout tail to warm up the filter for the backward pass. It doesn't actually do this but instead
  uses closed form formulas for what the states should be. */
  void prepareForBackwardPass()
  {
    x1 = a1*y1 + b1*x1;
    y1 = (a1*b1 + b0)*x1 / (TPar(1) - a1*a1);
    // For a one-pole filter without a zero, the simplified formula for the y1 state would be:
    // y1 = (b0*y1*a1) / (1-a1*a1);
    // If a1 is close to 1 (which is typical for lowpass filters), the denominator 1-a1^2 will
    // suffer from precision loss due to cancellation - can this be avoided?
  }
  // maybe for image processing, it's better to assume that the missing pixels just repeat the last
  // value instead of going down to zero? todo: derive the equations for x1,y1 when we do not assume 
  // that the input drops to zero but instead goes to some arbitrary constant value - include the 
  // value as optional parameter which defaults to zero

  /** Like prepareForBackwardPass without parameter, but does not assume the input signal to go 
  down to zero in the tail but instead to settle to some arbitrary constant value c. */
  void prepareForBackwardPass(TSig c)
  {
    // use the simpler function above, if possible:
    if(c == TSig(0)) { 
      prepareForBackwardPass(); return; }

    // compute some intermediate values:
    TPar a2 = a1*a1;                              // a1^2
    TPar a3 = a1*a2;                              // a1^3
    TPar cx = b0*a1 + (a2-a1)*b1 - b0;            // coeff for x1 in update of y1
    TPar cc = b0*(b0*a1 + (a1+1)*b1) + b1*b1;     // coeff for c in update of y1

    // compute new state variables:
    x1 = b0*c + a1*y1 + b1*x1;                    // update x1
    y1 = (cc*c - cx*x1) / (a3-a2-a1+1);           // update y1 using updated x1
  }

  /** Applies the filter bidirectionally (once forward, once backward) to the input signal x and 
  stores the result in y. Both buffers are assumed to be of length N. Can be used in place, i.e. 
  x and y may point to the same buffer. We apply the forward pass first, then the backward pass 
  but it should make no difference, if it would be done the other way around (up to roundoff 
  errors). The optional parameters xL,xR (L,R stand for left,right), which default to zero, are 
  used for telling the filter, how the input signal x should be assumed to continue outside the 
  range of valid sample indices. We will assume that x[n] = xL for n < 0 and 
  x[n] = xR for n >= N. */
  void applyBidirectionally(TSig* x, TSig* y, int N, TSig xL = TSig(0), TSig xR = TSig(0))
  {
    // todo: have optional xL,xR parameters for x[n] for n < 0 and n >= N, both defaulting to 0

    // forward pass:
    //reset();                         // todo: generalize to setStateForConstInput(xL)
    setStateForConstInput(xL);
    for(int n = 0; n < N; n++)
      y[n] = getSample(x[n]);

    // backward pass:
    //prepareForBackwardPass();        // todo: pass xR as parameter
    prepareForBackwardPass(xR);
    for(int n = N-1; n >= 0; n--)
      y[n] = getSample(y[n]);
  }

  /** Applies the filter bidirectionally with a stride (i.e. index-distance between two successive 
  samples) that is not necessarrily unity. This may be useful for filtering along a particular 
  dimension in multidimensional arrays such as images. */
  void applyBidirectionally(TSig* x, TSig* y, int N, int stride, 
    TSig xL = TSig(0), TSig xR = TSig(0))
  {
    // forward pass:
    //reset();
    setStateForConstInput(xL);
    for(int n = 0; n < N; n++)
      y[n*stride] = getSample(x[n*stride]);

    // backward pass:
    //prepareForBackwardPass();
    prepareForBackwardPass(xR);
    for(int n = N-1; n >= 0; n--)
      y[n*stride] = getSample(y[n*stride]);
  }
  // generalize this to use optional xL,xR parameters, too
  // optimize: use n += stride and n -= stride in loop headers and get rid of the multiplications 
  // n*stride in loop bodies


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  // buffering:
  TSig x1 = 0, y1 = 0; // past input x[n-1] and output y[n-1]

  // filter coefficients:
  TPar b0 = 1, b1 = 0; // feedforward coeffs
  TPar a1 = 0;         // feedback coeff

};

//=================================================================================================
/** This is an implementation of a simple one-pole filter unit.

\todo 
-rename to FirstOrderFilter (it does not only have a pole but also a zero).
 or to rsFilter1p1z
-make it possible to set up time constants in terms of dB/sec 
*/

template<class TSig, class TPar>
class rsOnePoleFilter : public rsFirstOrderFilterBase<TSig, TPar>
{

public:



  /** This is an enumeration of the available filter modes. */
  enum modes
  {
    BYPASS = 0,
    LOWPASS_IIT,      // lowpass via impulse invariant transform
    HIGHPASS_MZT,     // highpass via matched-Z transform
    //ALLPASS_MZT
    LOWPASS_BLT,      // lowpass via bilinear transform
    HIGHPASS_BLT,     // highpass via bilinear transform
    ALLPASS_BLT,      // allpass via bilinear transform
    LOWSHELV_BLT,     // low shelving via bilinear transform
    HIGHSHELV_BLT,    // high shelving via bilinear transform

    LOWSHELV_NMM,     // low shelving via nyquist magnitude match
    HIGHSHELV_NMM,    // high shelving via nyquist magnitude match

  };
  // NMM maybe can also be called PMM for pointwise magnitude match

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsOnePoleFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Chooses the filter mode. See the enumeration for available modes. */
  void setMode(int newMode);

  /** Sets the cutoff-frequency for this filter. */
  void setCutoff(TPar newCutoff);

  /** This will set the time constant 'tau' for the case, when lowpass mode is chosen. This is
  the time, it takes for the impulse response to die away to 1/e = 0.368... or equivalently, the
  time it takes for the step response to raise to 1-1/e = 0.632... */
  void setLowpassTimeConstant(TPar newTimeConstant) 
  { setCutoff(TPar(1)/(TPar(2*PI)*newTimeConstant)); }

  /** Sets the gain factor for the shelving modes (this is not in decibels). */
  void setShelvingGain(TPar newGain);

  /** Sets the gain for the shelving modes in decibels. */
  void setShelvingGainInDecibels(TPar newGain);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline TPar getCutoff() const { return cutoff; }

  inline TPar getShelvingGain() const { return shelvingGain; }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the magnitude response of this filter at the given frqeuency (in Hz). */
  TPar getMagnitudeAt(TPar frequency);
  // hmmm...would be nice to have that in baseclass taking omega as input - rename this to
  // getMagnitudeAtHz

protected:

  /** \name Internal Functions */
  void calcCoeffs();  // calculates filter coefficients from filter parameters


  /** \name Data */

  // filter parameters:
  TPar cutoff;
  TPar shelvingGain;
  TPar freqToOmega; // = 2*PI/sampleRate, conversion factor from physical to digital frequency
  int  mode;

  //TPar sampleRate;
  //TPar sampleRateRec;  // reciprocal of the sampleRate 
  // instead of these two, keep only 2*PI/sampleRate, i.e. the freqToOmega factor
};

#endif
