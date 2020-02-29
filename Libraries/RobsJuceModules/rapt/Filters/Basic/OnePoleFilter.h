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
    *b0 =  0.5*(1 + *a1);
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
    T t = tan(0.5*w); // tan w/2
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
    T t = tan(0.5*w);

    rsAssert(rsIsFiniteNumber(t));
    // hmmm...maybe we should allow t == inf and set a1 = -1 in this case (that's the limiting 
    // value)


    *a1 = (1-t) / (1+t);
    *b0 = 0.5*(1 - *a1);
    *b1 = *b0;
    // w = 0:  t = 0 -> a1 = 1, b0,b1 = 0
    // w = pi: t = inf -> a1,b0,b1 = NaN
  }

  /** Highpass via bilinear transform (from DAFX) */
  template<class T>
  static inline void coeffsHighpassBLT(T w, T* b0, T* b1, T* a1)
  {
    T t = tan(0.5*w);     // maybe re-use w
    *a1 = (1-t) / (1+t);
    *b0 = 0.5*(1 + *a1);
    *b1 = -(*b0);
  }

  /** Low shelving via bilinear transform, g is linear shelving gain. (from DAFX) */
  template<class T>
  static inline void coeffsLowShelfBLT(T w, T g, T* b0, T* b1, T* a1)
  {
    T t = tan(0.5*w);  // re-use w
    t   = g >= 1.0 ? (t-1)/(t+1) : (t-g)/(t+g);
    T c = 0.5*(g-1);   // re-use g
    c  += c*t;
    *b0 = 1 + c;
    *b1 = t + c;
    *a1 = -t;
  }

  template<class T>
  static inline void coeffsHighShelfBLT(T w, T g, T* b0, T* b1, T* a1)
  {
    T t = tan(0.5*w);
    t   = g >= 1.0 ? (t-1.0)/(t+1.0) : (g*t-1)/(g*t+1);
    T c = 0.5*(g-1);
    c  -= c*t;
    *b0 = 1 + c;
    *b1 = t - c;
    *a1 = -t;
  }

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
  void setLowpassTimeConstant(TPar newTimeConstant) { setCutoff(1.0/(2*PI*newTimeConstant)); }

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
