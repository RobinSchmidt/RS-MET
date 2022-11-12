#ifndef RAPT_BIQUAD_H
#define RAPT_BIQUAD_H


// this file contains several classes for biquad filters, biquad cascades and related stuff


/** This class serves as baseclass for biquad filters which realize the difference equation:

\f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

it contains only the coefficients as members - past inputs/output samples are supposed to be
added in subclasses because different topolies require different state variables.

// todo: check, if the difference equation's sign convention is right (plusses for the recursive
part) */

template<class TCof> // coefficient type
class rsBiquad
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsBiquad();

  /** Destructor. */
  ~rsBiquad();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup: */

  /** Sets the filter coefficients to new values. */
  void setCoefficients(TCof newB0, TCof newB1, TCof newB2, TCof newA1, TCof newA2);



  //-----------------------------------------------------------------------------------------------
  /** \name Misc: */

  /** Initializes the biquad coefficients to b0=1.0, b1=b2=a1=a2=0.0 which is essentially a
  bypass 'filter'. */
  void initializeCoefficients();


protected:

  // direct form coefficients:
  TCof b0, b1, b2, a1, a2;

};


//=================================================================================================

/**

This class implements a biquad filter which realize the difference equation:

\f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

note the positive sign of the feedback part, as opposed to negative sign seen in many DSP
textbooks - this has been chosen to allow for better potential for parallelization.

*/

template<class TSig, class TCof>
class rsBiquadDF1 : public rsBiquad<TCof>
{

public:

  /** \name Construction/Destruction: */

  /** Constructor. */
  rsBiquadDF1();

  /** Destructor. */
  ~rsBiquadDF1();


  /** \name Audio processing */

  /** Calculates a single filtered output-sample in place. */
  RS_INLINE void getSampleInPlace(TSig &inOut);  // get rid

  /** Calculates a single filtered output-sample. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc: */

  /** Sets the buffers for the previous input and output samples to zero. */
  void reset();


protected:

  // internal state buffers:
  TSig x1, x2, y1, y2;

};

//-----------------------------------------------------------------------------------------------
// inline functions:

template<class TSig, class TCof>
RS_INLINE void rsBiquadDF1<TSig, TCof>::getSampleInPlace(TSig &inOut)
{
  TSig tmp = inOut;

  tmp = (this->b0 * tmp) + (this->b1 * x1 + this->b2 * x2) + (this->a1 * y1 + this->a2 * y2);
  // parentheses (may?) facilitate vectorization

  x2 = x1;
  x1 = inOut;
  y2 = y1;
  y1 = tmp;

  inOut = tmp;
}

template<class TSig, class TCof>
RS_INLINE TSig rsBiquadDF1<TSig, TCof>::getSample(TSig in)
{
  TSig tmp;

  tmp = (this->b0 * in) + (this->b1 * x1 + this->b2 * x2) + (this->a1 * y1 + this->a2 * y2);
  // parentheses (may?) facilitate vectorization

  x2 = x1;
  x1 = in;
  y2 = y1;
  y1 = tmp;

  return tmp;
}

//===============================================================================================

/** This class implements design formulas for the calculation of the coefficients of biquad filters
which realize the difference equation:

\f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

note the positive sign of the feedback part, as opposed to negative sign seen in many DSP
textbooks - this has been chosen to allow for better potential for parallelization.

\todo make this consistent with BiquadCascade and all other filters...
\todo move these functions into the BiquadClass
\todo get rid of the inlining, maybe pass w = 2*PI*f/fs instead of frequency and
          "oneOverSampleRate" */

class rsBiquadDesigner
{

public:

  /** \name Coefficient calculation: */

  /** Assigns the coefficients so as to realize a neutral filter. */
  template<class T>
  static RS_INLINE void makeBypassBiquad(T& b0, T& b1, T& b2, T& a1,
    T& a2);

  /** Calculates the cos(w), sin(w) with w=2*pi*frequency*oneOverSampleRate which are
  intermediate variables in many design formulas. */
  template<class T>
  static RS_INLINE void calculateSineAndCosine(T& sinResult, T& cosResult,
    const T& frequency, const T& oneOverSampleRate);

  /** Calculates the coefficients for a one-pole lowpass filter according to the formula from
  "The Scientist and Engineer's Guide to Digital Signal Processing". */
  template<class T>
  static RS_INLINE void calculateFirstOrderLowpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency);

  /** Calculates the coefficients for a one-pole lowpass filter by the bilinear transform
  method. */
  template<class T>
  static RS_INLINE void calculateFirstOrderLowpassCoeffsBilinear(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency);

  /** Calculates coefficients for a 1st order lowpass with a gain at the Nyquist frequency that
  matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
  template<class T>
  static void calculateFirstOrderLowpassCoeffsPrescribedNyquist(T &b0, T &b1,
    T &b2, T &a1, T &a2, const T &sampleRate, const T &frequency);

  /** Calculates the coefficients for a two-pole lowpass type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookLowpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q);

  /** Calculates the coefficients for a one-pole highpass filter according to the formula from
  "The Scientist and Engineer's Guide to Digital Signal Processing". */
  template<class T>
  static RS_INLINE void calculateFirstOrderHighpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency);

  /** Calculates the coefficients for a one-pole highpass filter by the bilinear transform
  method. */
  template<class T>
  static RS_INLINE void calculateFirstOrderHighpassCoeffsBilinear(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& oneOverSampleRate,
    const T& frequency);

  /** Calculates coefficients for a 1st order highpass with a gain at the Nyquist frequency that
  matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
  template<class T>
  static void calculateFirstOrderHighpassCoeffsPrescribedNyquist(T &b0, T &b1,
    T &b2, T &a1, T &a2, const T &sampleRate, const T &frequency);

  /** Calculates the coefficients for a two-pole highpass type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookHighpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q);

  /** Calculates the coefficients for a two-pole bandpass type of filter with constant skirt
  gain according to Robert Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookBandpassConstSkirtCoeffsViaQ(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& oneOverSampleRate,
    const T& frequency, const T& q);

  /** Calculates the coefficients for a two-pole bandpass type of filter with constant skirt
  gain according to Robert Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookBandpassConstSkirtCoeffsViaBandwidth(T& b0,
    T& b1, T& b2, T& a1, T& a2, const T& oneOverSampleRate,
    const T& frequency, const T& bandwidth);

  /** Calculates the coefficients for a two-pole bandreject type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookBandrejectCoeffsViaQ(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& oneOverSampleRate,
    const T& frequency, const T& q);

  /** Calculates the coefficients for a two-pole bandreject type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookBandrejectCoeffsViaBandwidth(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& oneOverSampleRate,
    const T& frequency, const T& bandwidth);

  /** Calculates the coefficients for a peaking equalizer type of filter according to Robert
  Bristow Johnson's Cookbook formula using the 'Q' parameter (as opposed to a bandwidth). */
  template<class T>
  static RS_INLINE void calculateCookbookPeakFilterCoeffsViaQ(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q, const T& gainFactor);

  /** Calculates the coefficients for a first order low shelving type of filter (according to the
  formula from DAFX?). */
  template<class T>
  static RS_INLINE void calculateFirstOrderLowShelvCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& gainFactor);

  /** Calculates the coefficients for a low shelving type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookLowShelvCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q, const T& gainFactor);

  /** Calculates coefficients for a 1st order low-shelving with a gain at the Nyquist frequency
  that matches the gain of the corresponding analog prototype filter at that Nyquist
  frequency. */
  template<class T>
  static void calculateFirstOrderLowShelvCoeffsPrescribedNyQuist(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& sampleRate, const T& frequency,
    const T& gainFactor);

  /** Calculates the coefficients for a first order high shelving type of filter (according to
  the formula from DAFX?). */
  template<class T>
  static RS_INLINE void calculateFirstOrderHighShelvCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& gainFactor);

  /** Calculates the coefficients for a high shelving type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookHighShelvCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q, const T& gainFactor);

  /** Calculates coefficients for a 1st order high-shelving with a gain at the Nyquist frequency
  that matches the gain of the corresponding analog prototype filter at that Nyquist
  frequency. */
  template<class T>
  static void calculateFirstOrderHighShelvCoeffsPrescribedNyQuist(T& b0, T& b1,
    T& b2, T& a1, T& a2, const T& sampleRate, const T& frequency,
    const T& gainFactor);

  /** Calculates the coefficients for a first order allpass filter (according to the formula from
  DAFX?). */
  template<class T>
  static RS_INLINE void calculateFirstOrderAllpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency);

  /** Calculates the coefficients for a second order allpass type of filter according to Robert
  Bristow Johnson's Cookbook formula. */
  template<class T>
  static RS_INLINE void calculateCookbookAllpassCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& q);

  /** Calculates coefficients for a biquad that represents a chain of a first order lowpass and a
  first order highpass filter. */
  template<class T>
  static RS_INLINE void calculateLowHighpassChainCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& highpassCutoff,
    const T& lowpassCutoff);

  /** Calculates the coefficients for a second order filter which is morphable between lowpass
  through peak to highpass according to some formulas derived by Robin Schmidt. */
  template<class T>
  static RS_INLINE void calculateLowPeakHighMorphCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, T& preGain, const T& sampleRate,
    const T& oneOverSampleRate, const T& frequency, const T& q,
    const T& morph, const bool &scaleGainForTwoStages);

  /** Calculates the coefficients for a peaking equalizer type of filter according to Sophocles
  Ofanidis' design procedure which matches the analog equalizers gain at the Nyquist
  frequency. */
  template<class T>
  static RS_INLINE void calculatePrescribedNyquistGainEqCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
    const T& bandwidthInOctaves, const T& gainFactor, const T& referenceGain);

  /** Calculates the coefficients for a one-pole filter where the pole is a number along the
  real axis in the z-plane (and should be less than 1 in magnitude for stable filters). */
  template<class T>
  static RS_INLINE void calculateOnePoleCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& pole);

  /** Calculates the coefficients for a one-zero filter where the pole is a number along the
  real axis in the z-plane. */
  template<class T>
  static RS_INLINE void calculateOneZeroCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& zero);

  /** Calculates the coefficients for a two-pole filter where the poles form a complex conjugate
  pair with the given radius and angle (the radius should be less than 1 in magnitude for stable
  filters). */
  template<class T>
  static RS_INLINE void calculatePolePairCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& poleRadius, const T& poleAngle);

  /** Calculates the coefficients for a two-zero filter where the zeros form a complex conjugate
  pair with the given radius and angle. */
  template<class T>
  static RS_INLINE void calculateZeroPairCoeffs(T& b0, T& b1, T& b2,
    T& a1, T& a2, const T& zeroRadius, const T& zeroAngle);


   /** \name Misc */

  /** Returns the magnitude response of a biquad filter with given coefficients at a given
  frequency (and sample-rate). */
  template<class T>
  static T getBiquadMagnitudeAt(const T& b0, const T& b1, const T& b2,
    const T& a1, const T& a2, const T& frequency, const T& sampleRate);

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
RS_INLINE void rsBiquadDesigner::makeBypassBiquad(T& b0, T& b1, T& b2, T& a1,
  T& a2)
{
  b0 = 1.0;
  b1 = 0.0;
  b2 = 0.0;
  a1 = 0.0;
  a2 = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateSineAndCosine(T& sinResult, T& cosResult,
  const T& frequency, const T& oneOverSampleRate)
{
  // calculate intermediate variables:
  T omega  = 2.0 * PI * frequency * oneOverSampleRate;
  rsSinCos(omega, &sinResult, &cosResult);
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderLowpassCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency)
{
  T omega = 2.0 * PI * frequency * oneOverSampleRate;
  T x     = exp(-omega);

  a1 = x;
  a2 = 0.0;
  b0 = 1.0-x;
  b1 = 0.0;
  b2 = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderLowpassCoeffsBilinear(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency)
{
  T omegaPreWarped = tan(PI * frequency * oneOverSampleRate);
  T pAnalog        = -omegaPreWarped;
  a1                    = (1.0+pAnalog)/(1.0-pAnalog);
  a2                    = 0.0;
  T g              = 0.5*rsSqrt(1.0+a1*a1-2.0*a1); // gain-factor for normalization at DC
  b0                    = g;
  b1                    = g;
  b2                    = 0.0;
}

/*
template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderLowpassCoeffsPrescribedNyquist(T &b0,
  T &b1, T &b2, T &a1, T &a2, const T &sampleRate,
  const T &frequency)
{
  T wc  = 2.0*PI*frequency/sampleRate;
  T Wc  = tan(wc/2.0);
  T Ws  = PI*sampleRate;
  T Wc2 = Wc*Wc;

  T Wca = 2*PI*frequency;                      // non-pre-warped analog cutoff frequency
  T k2  = 1.0 / (1.0 + ((Ws*Ws)/(Wca*Wca)) ); // gain of prototype at Nyquist-freq
  //k2  = 0; // test - should revert to standard bilinear design

  T A2  = 1.0 / (Wc2*(1.0-2.0*k2));
  T B2  = k2*A2;
  T g02 = 1.0;

  T A   = rsSqrt(A2);
  T B   = rsSqrt(B2);
  T g0  = rsSqrt(g02);      // == 1.0 -> optimize out
  T rD  = 1.0 / (1.0+A);  // reciprocal of denominator

  b0 =  rD * (g0+B);
  b1 =  rD * (g0-B);
  b2 =  0.0;
  a1 = -rD * (1 -A);
  a2 =  0.0;
}
*/

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookLowpassCoeffs(T& b0, T& b1,
  T& b2, T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
  const T& q)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(2*q);
  T a0Rec = 1/(1+alpha);

  a1 = 2*cosine   * a0Rec;
  a2 = (alpha-1)  * a0Rec;
  b1 = (1-cosine) * a0Rec;
  b0 = T(0.5)*b1;
  b2 = b0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderHighpassCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency)
{
  T omega = 2.0 * PI * frequency * oneOverSampleRate;
  T x     = exp(-omega);

  a1 = x;
  a2 = 0.0;
  b0 = 0.5*(1.0+x);
  b1 = -b0;
  b2 = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderHighpassCoeffsBilinear(T &b0,
  T &b1, T &b2, T &a1, T &a2, const T &oneOverSampleRate,
  const T &frequency)
{
  T omegaPreWarped = tan(PI * frequency * oneOverSampleRate);
  T pAnalog        = -omegaPreWarped;
  a1                    = (1.0+pAnalog)/(1.0-pAnalog);
  a2                    = 0.0;
  T g              = 0.5*rsSqrt(1.0+a1*a1+2.0*a1); // gain-factor for normalization at DC
  b0                    = g;
  b1                    = -g;
  b2                    = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookHighpassCoeffs(T& b0, T& b1,
  T& b2, T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
  const T& q)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(T(2)*q);
  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine    * a0Rec;
  a2 = (alpha-T(1))   * a0Rec;
  b1 = -(T(1)+cosine) * a0Rec;
  b0 = -T(0.5)*b1;
  b2 = b0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ(T &b0,
  T &b1, T &b2, T &a1, T &a2, const T &oneOverSampleRate,
  const T &frequency, const T &q)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(T(2)*q);
  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine   * a0Rec;
  a2 = (alpha-T(1))  * a0Rec;
  b1 = T(0);
  b0 = q*alpha       * a0Rec;
  b2 = -b0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaBandwidth(
  T &b0, T &b1, T &b2, T &a1, T &a2, const T &oneOverSampleRate,
  const T &frequency, const T &bandwidth)
{
  // to avoid divison by zero:
  T f = rsMax(T(0.0001), frequency);
  T b = rsMax(T(0.0001), bandwidth);

  T sine, cosine;
  calculateSineAndCosine(sine, cosine, f, oneOverSampleRate);

  // we need some if( sine == 0 )  ..check the limit
  //T alpha = sine * sinh(T(0.5)*log(T(2)) * b * T(2.0*PI)*f*oneOverSampleRate / sine);
  T alpha = sine * sinh( T(log(2.0)*PI) * b*f*oneOverSampleRate / sine);

  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine   * a0Rec;
  a2 = (alpha-T(1))  * a0Rec;
  b1 = T(0);
  b0 = T(0.5)*sine   * a0Rec;
  b2 = -b0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookBandrejectCoeffsViaQ(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency,
  const T &q)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(T(2)*q);
  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine  * a0Rec;
  a2 = (alpha-T(1)) * a0Rec;
  b0 = T(1)         * a0Rec;
  b1 = T(-2)*cosine * a0Rec;
  b2 = T(1)         * a0Rec;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookBandrejectCoeffsViaBandwidth(T &b0,
  T &b1, T &b2, T &a1, T &a2, const T &oneOverSampleRate,
  const T &frequency, const T &bandwidth)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  //T alpha = sine * sinh(T(0.5*log(2.0)) * bandwidth * T(2.0*PI)*frequency*oneOverSampleRate / sine);
  T alpha = sine * sinh(T(log(2.0)*PI) * bandwidth*frequency*oneOverSampleRate / sine);
  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine  * a0Rec;
  a2 = (alpha-T(1)) * a0Rec;
  b0 = T(1)         * a0Rec;
  b1 = T(-2)*cosine * a0Rec;
  b2 = T(1)         * a0Rec;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookPeakFilterCoeffsViaQ(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency,
  const T &q, const T &gainFactor)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(T(2)*q);
  T A     = gainFactor;
  T a0Rec = T(1)/(T(1)+alpha/A);

  a1 = T(2)*cosine        * a0Rec;
  a2 = ((alpha/A) - T(1)) * a0Rec;
  b0 = (T(1)+alpha*A)     * a0Rec;
  b1 = T(-2)*cosine       * a0Rec;
  b2 = (T(1)-alpha*A)     * a0Rec;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderLowShelvCoeffs(T& b0, T& b1,
  T& b2, T& a1, T& a2, const T& oneOverSampleRate, const T& frequency,
  const T& gainFactor)
{
  T A     = gainFactor;
  T omega = T(2) * PI * frequency * oneOverSampleRate;
  T x     = tan(T(0.5)*omega);
  T c     = T(0.5)*(A*A-T(1));  // A is the square root of the gain-factor -> V_0 = A^2
  T a     = T(0);

  if(A > T(1.0001)) // boost
  {
    a  = (x-T(1))/(x+T(1));
    b0 = T(1) + c*a + c;
    b1 = a   + c*a + c;
    b2 = T(0);
    a1 = -a;
    a2 = T(0);
  }
  else if(A < T(0.999)) // cut
  {
    a  = (x-A*A)/(x+A*A);
    b0 = T(1) + c*a + c;
    b1 = a   + c*a + c;
    b2 = T(0);
    a1 = -a;
    a2 = T(0);
  }
  else
  {
    b0 = T(1);
    b1 = T(0);
    b2 = T(0);
    a1 = T(0);
    a2 = T(0);
  }
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookLowShelvCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency,
  const T &q, const T &gainFactor)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T A     = gainFactor;
  T beta  = rsSqrt(A) / q;
  T a0Rec = T(1) / ((A+T(1)) + (A-T(1))*cosine + beta*sine);

  a1 = T(2) *     ((A-T(1)) + (A+T(1))*cosine) * a0Rec;
  a2 = -((A+1.0) + (A-T(1))*cosine - beta*sine) * a0Rec;
  b0 =        A * ((A+T(1)) - (A-T(1))*cosine + beta*sine) * a0Rec;
  b1 = T(2) * A * ((A-T(1)) - (A+T(1))*cosine) * a0Rec;
  b2 =        A * ((A+T(1)) - (A-T(1))*cosine - beta*sine) * a0Rec;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderHighShelvCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate,
  const T &frequency, const T &gainFactor)
{
  T A     = gainFactor;
  T omega = T(2) * PI * frequency * oneOverSampleRate;
  T x     = tan(T(0.5)*omega);
  T c     = T(0.5)*(A*A-T(1));  // A is the square root of the gain-factor -> V_0 = A^2
  T a     = T(0);

  if(A > T(1.0001)) // boost
  {
    a  = (x-T(1))/(x+T(1));
    b0 = T(1) - c*a + c;
    b1 = a    + c*a - c;
    b2 = T(0);
    a1 = -a;
    a2 = T(0);
  }
  else if(A < T(0.999)) // cut
  {
    a  = (A*A*x-T(1))/(A*A*x+T(1));
    b0 = T(1) - c*a + c;
    b1 = a    + c*a - c;
    b2 = T(0);
    a1 = -a;
    a2 = T(0);
  }
  else
  {
    b0 = T(1);
    b1 = T(0);
    b2 = T(0);
    a1 = T(0);
    a2 = T(0);
  }
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookHighShelvCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency,
  const T &q, const T &gainFactor)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T A     = gainFactor;
  T beta  = rsSqrt(A) / q;
  T a0Rec = T(1)  /  ((A+T(1)) - (A-T(1))*cosine + beta*sine);

  a1 = T(-2) *     ((A-T(1)) - (A+T(1))*cosine)             * a0Rec;
  a2 = -((A+T(1)) - (A-T(1))*cosine - beta*sine)            * a0Rec;
  b0 =         A * ((A+T(1)) + (A-T(1))*cosine + beta*sine) * a0Rec;
  b1 = T(-2) * A * ((A-T(1)) + (A+T(1))*cosine)             * a0Rec;
  b2 =         A * ((A+T(1)) + (A-T(1))*cosine - beta*sine) * a0Rec;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateFirstOrderAllpassCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency)
{
  T t = tan(PI*frequency*oneOverSampleRate);
  T x = (t-T(1)) / (t+T(1));

  b0 = x;
  b1 = T(1);
  b2 = T(0);
  a1 = -x;
  a2 = T(0);
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateCookbookAllpassCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &oneOverSampleRate, const T &frequency,
  const T &q)
{
  T sine, cosine;
  calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
  T alpha = sine/(T(2)*q);
  T a0Rec = T(1)/(T(1)+alpha);

  a1 = T(2)*cosine    * a0Rec;
  a2 = (alpha-T(1))   * a0Rec;
  b0 = (T(1)-alpha)   * a0Rec;
  b1 = (T(-2)*cosine) * a0Rec;
  b2 = (T(1)+alpha)   * a0Rec;
}

template<class T>
RS_INLINE void calculateLowHighpassChainCoeffs(T& b0, T& b1, T& b2,
  T& a1, T& a2, const T& oneOverSampleRate, const T& highpassCutoff,
  const T& lowpassCutoff)
{
  // variables for the coefficients of the 2 one-pole filters:
  T b0Lpf, b1Lpf, a1Lpf; // lowpass-coeffs
  T b0Hpf, b1Hpf, a1Hpf; // highpass-coeffs

  T x; // intermediate variable for calculation (x is the amount of
  // decay between adjacent samples for LPFs):

  // calculate lowpass coefficients (see dspguide for details):
  if(lowpassCutoff == T(0))
  {
    b0Lpf = T(1);
    b1Lpf = T(0);
    a1Lpf = T(0);
  }
  else
  {
    x = exp(T(-2.0 * PI) * lowpassCutoff * oneOverSampleRate);

    b0Lpf = T(1)-x;
    b1Lpf = T(0);
    a1Lpf = -x;
  }

  // calculate highpass coefficients (see dspguide for details):
  if(highpassCutoff == T(0))
  {
    b0Hpf = T(1);
    b1Hpf = T(0);
    a1Hpf = T(0);
  }
  else
  {
    x = exp(T(-2.0 * PI) * highpassCutoff * oneOverSampleRate);

    b0Hpf =  T(0.5)*(1+x);
    b1Hpf = T(-0.5)*(1+x);
    a1Hpf = -x;
  }

  // combine the 2 one-pole filters into one biquad filter:
  b0 = b0Lpf*b0Hpf;
  b1 = b0Lpf*b1Hpf + b1Lpf*b0Hpf;
  b2 = b1Lpf*b1Hpf;
  a1 = -(a1Lpf+a1Hpf);
  a2 = -(a1Lpf*a1Hpf);
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateLowPeakHighMorphCoeffs(T &b0, T &b1,
  T &b2, T &a1, T &a2, T& preGain, const T &sampleRate,
  const T &oneOverSampleRate, const T &frequency, const T &q,
  const T &morph, const bool &scaleGainForTwoStages)
{
  // calculate normalized radian frequencies (digital and pre-warped analog):
  T wd = 2.0 * PI * frequency * oneOverSampleRate; // optimize to one mul
  T wa = 2.0 * sampleRate * tan(0.5*wd);

  // calculate the position of the analog pole-pair:
  T tmp1 = -wa / (2.0*q);
  T tmp2;
  if(q >= 0.5)
    tmp2 = rsSqrt(-wa*wa*(0.25/(q*q) - 1.0));
  else
  {
    tmp2 = 0.0;
    RS_DEBUG_BREAK; // this design procedure assumes q >= 0.5
  }

  std::complex<T> pa = std::complex<T>(tmp1, tmp2);

  // calculate the position of the digital pole-pair:
  std::complex<T> z  = 0.5 * oneOverSampleRate * pa;
  std::complex<T> pd = (1.0+z) / (1.0-z);

  // convert poles to biquad feedback coefficients:
  a1 = 2.0 * pd.re;
  a2 = -(pd.re*pd.re + pd.im*pd.im);

  // calculate the (T) zero for a peak filter:
  T num  = (1.0-pd.re)*(1.0-pd.re) + pd.im*pd.im;
  T den  = (1.0+pd.re)*(1.0+pd.re) + pd.im*pd.im;
  T c    = rsSqrt(num / den);
  T z_pk = (1.0-c) / (1.0+c);

  // obtain the zero by interpolating the zero position between lowpass and peak or peak
  // and highpass:
  T m = morph;
  m = rsSign(m) * (m*m);
  //m = sign(m) * pow(abs(m), 1.5);
  T zd;
  //T a = -z_pk;
  if(m < -0.99)
    zd = -0.99;
  else if(m > 0.99)
    zd = +0.99;
  else
    zd = (m+z_pk) / (1.0+z_pk*m);

  // convert zeros to biquad feedforward coefficients:
  b0 = 1.0;
  b1 = -(zd+zd);
  b2 = zd*zd;

  // normalize gain:
  preGain = q / getBiquadMagnitudeAt(b0, b1, b2, a1, a2, frequency, sampleRate);
  if(scaleGainForTwoStages)
    preGain *= preGain;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(T& b0,
  T& b1, T& b2, T& a1, T& a2, const T& oneOverSampleRate,
  const T& frequency, const T& bandwidthInOctaves, const T& gainFactor,
  const T& referenceGain)
{
  if(fabs(rsAmp2dB(gainFactor/referenceGain)) < 0.001)
  {
    makeBypassBiquad(b0, b1, b2, a1, a2);
    return;
  }

  T G0 = referenceGain;               // reference gain at DC
  T G  = gainFactor;                  // boost/cut gain
  T fc = frequency;                   // center frequency in Hz
  T bw = bandwidthInOctaves;          // bandwidth in octaves

  // intermediate variables:
  T fLo = fc  / pow(T(2), T(0.5)*bw);  // lower bandedge frequency
  T fHi = fLo * pow(T(2), bw);         // upper bandedge frequency
  T GB  = rsSqrt(G);                   // gain at bandedge frequencies
  T w0  = 2*PI*fc*oneOverSampleRate;
  T wLo = 2*PI*fLo*oneOverSampleRate;
  T wHi = 2*PI*fHi*oneOverSampleRate;
  T Dw  = wHi-wLo;

  // the design procedure (ported from peq.m by Sophocles Orfanidis) - may be optimized:
  T F   = fabs(G*G   - GB*GB);
  T G00 = fabs(G*G   - G0*G0);
  T F00 = fabs(GB*GB - G0*G0);
  T num = G0*G0 * (w0*w0-PI*PI)*(w0*w0-PI*PI) + G*G * F00 * PI*PI * Dw*Dw / F;
  T den = (w0*w0-PI*PI)*(w0*w0-PI*PI) + F00 * PI*PI * Dw*Dw / F;
  T G1  = rsSqrt(num/den);

  //---------------------------------
  // inserted modification by Robin Schmidt: enforce InEq.37 to hold:
  {
    T Delta = 2*PI*(fHi-fLo);                              // bandwidth in 2*pi Hz
    T ASq   = Delta*Delta * (GB*GB-G0*G0) / (G*G-GB*GB);   // Eq.5

    // redefine bandedge gain to enforce InEq.37 to hold:
    GB = rsSqrt(G1*G);

    // re-calculate bandwidth from the new GB definition as the width where an analog equalizer
    // would have the redefined bandedge gain:
    Delta  = rsSqrt(ASq*(G*G-GB*GB) / (GB*GB-G0*G0));
    Dw     = Delta*oneOverSampleRate;

    // re-calculate affected intermediate variables:
    F   = fabs(G*G   - GB*GB);
    F00 = fabs(GB*GB - G0*G0);
    num = G0*G0 * (w0*w0-PI*PI)*(w0*w0-PI*PI) + G*G * F00 * PI*PI * Dw*Dw / F;
    den = (w0*w0-PI*PI)*(w0*w0-PI*PI) + F00 * PI*PI * Dw*Dw / F;
    G1  = rsSqrt(num/den);
  }
  //---------------------------------

  T G01 = fabs(G*G   - G0*G1);
  T G11 = fabs(G*G   - G1*G1);
  T F01 = fabs(GB*GB - G0*G1);
  T F11 = fabs(GB*GB - G1*G1);
  T ta  = tan(w0/2);
  T W2  = rsSqrt(G11 / G00) * ta*ta;
  T DW  = (1 + rsSqrt(F00 / F11) * W2) * tan(Dw/2);
  T C   = F11 * DW*DW - 2 * W2 * (F01 - rsSqrt(F00 * F11));
  T D   = 2 * W2 * (G01 - rsSqrt(G00 * G11));
  T A   = rsSqrt((C + D) / F);
  T B   = rsSqrt((G*G * C + GB*GB * D) / F);
  T s   = 1.0 / (1 + W2 + A);

  b0         =  (G1 + G0*W2 + B)   * s;
  b1         =  (-2*(G1 - G0*W2))  * s;
  b2         =  (G1 - B + G0*W2)   * s;
  a1         = -(-2*(1 - W2))      * s;
  a2         = -(1 + W2 - A)       * s;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateOnePoleCoeffs(T& b0, T& b1, T& b2,
  T& a1, T& a2, const T& pole)
{
  b0 = 1.0;
  b1 = 0.0;
  b2 = 0.0;
  a1 = pole;
  a2 = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateOneZeroCoeffs(T& b0, T& b1, T& b2,
  T& a1, T& a2, const T& zero)
{
  b0 = 1.0;
  b1 = -zero;
  b2 = 0.0;
  a1 = 0.0;
  a2 = 0.0;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculatePolePairCoeffs(T& b0, T& b1, T& b2,
  T& a1, T& a2, const T& poleRadius, const T& poleAngle)
{
  b0 = 1.0;
  b1 = 0.0;
  b2 = 0.0;
  a1 = 2*poleRadius*cos(poleAngle);
  a2 = -poleRadius*poleRadius;
}

template<class T>
RS_INLINE void rsBiquadDesigner::calculateZeroPairCoeffs(T& b0, T& b1, T& b2,
  T& a1, T& a2, const T& zeroRadius, const T& zeroAngle)
{
  b0 = 1.0;
  b1 = -2*zeroRadius*cos(zeroAngle);
  b2 = zeroRadius*zeroRadius;
  a1 = 0.0;
  a2 = 0.0;
}

#endif
