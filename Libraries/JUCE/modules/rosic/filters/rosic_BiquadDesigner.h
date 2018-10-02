#ifndef rosic_BiquadDesigner_h
#define rosic_BiquadDesigner_h

//// rosic-indcludes:
//#include "../math/rosic_ComplexFunctions.h"

namespace rosic
{

  /**

  This class implements design formulas for the calculation of the coefficients of biquad filters
  which realize the difference equation:

  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

  note the positive sign of the feedback part, as opposed to negative sign seen in many DSP
  textbooks - this has been chosen to allow for better potential for parallelization.

  \todo: make this consistent with BiquadCascade and all other filters...

  */

  class BiquadDesigner
  {

  public:

    /** \name Coefficient calculation: */

    /** Assigns the coefficients so as to realize a neutral filter. */
    static INLINE void makeBypassBiquad(double& b0, double& b1, double& b2, double& a1, double& a2);

    /** Calculates the cos(w), sin(w) with w=2*pi*frequency*oneOverSampleRate which are
    intermediate variables in many design formulas. */
    static INLINE void calculateSineAndCosine(double& sinResult, double& cosResult,
      const double& frequency, const double& oneOverSampleRate);

    /** Calculates the coefficients for a one-pole lowpass filter according to the formula from
    "The Scientist and Engineer's Guide to Digital Signal Processing". */
    static INLINE void calculateFirstOrderLowpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency);

    /** Calculates the coefficients for a one-pole lowpass filter by the bilinear transform
    method. */
    static INLINE void calculateFirstOrderLowpassCoeffsBilinear(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency);

    /** Calculates coefficients for a 1st order lowpass with a gain at the Nyquist frequency that
    matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
    static void calculateFirstOrderLowpassCoeffsPrescribedNyquist(double &b0, double &b1,
      double &b2, double &a1, double &a2, const double &sampleRate, const double &frequency);

    /** Calculates the coefficients for a two-pole lowpass type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookLowpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q);

    /** Calculates the coefficients for a one-pole highpass filter according to the formula from
    "The Scientist and Engineer's Guide to Digital Signal Processing". */
    static INLINE void calculateFirstOrderHighpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency);

    /** Calculates the coefficients for a one-pole highpass filter by the bilinear transform
    method. */
    static INLINE void calculateFirstOrderHighpassCoeffsBilinear(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& oneOverSampleRate,
      const double& frequency);

    /** Calculates coefficients for a 1st order highpass with a gain at the Nyquist frequency that
    matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
    static void calculateFirstOrderHighpassCoeffsPrescribedNyquist(double &b0, double &b1,
      double &b2, double &a1, double &a2, const double &sampleRate, const double &frequency);

    /** Calculates the coefficients for a two-pole highpass type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookHighpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q);

    /** Calculates the coefficients for a two-pole bandpass type of filter with constant skirt
    gain according to Robert Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookBandpassConstSkirtCoeffsViaQ(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& oneOverSampleRate,
      const double& frequency, const double& q);

    /** Calculates the coefficients for a two-pole bandpass type of filter with constant skirt
    gain according to Robert Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookBandpassConstSkirtCoeffsViaBandwidth(double& b0,
      double& b1, double& b2, double& a1, double& a2, const double& oneOverSampleRate,
      const double& frequency, const double& bandwidth);

    /** Calculates the coefficients for a two-pole bandreject type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookBandrejectCoeffsViaQ(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& oneOverSampleRate,
      const double& frequency, const double& q);

    /** Calculates the coefficients for a two-pole bandreject type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookBandrejectCoeffsViaBandwidth(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& oneOverSampleRate,
      const double& frequency, const double& bandwidth);

    /** Calculates the coefficients for a peaking equalizer type of filter according to Robert
    Bristow Johnson's Cookbook formula using the 'Q' parameter (as opposed to a bandwidth). */
    static INLINE void calculateCookbookPeakFilterCoeffsViaQ(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q, const double& gainFactor);

    /** Calculates the coefficients for a first order low shelving type of filter (according to the
    formula from DAFX?). */
    static INLINE void calculateFirstOrderLowShelvCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& gainFactor);

    /** Calculates the coefficients for a low shelving type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookLowShelvCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q, const double& gainFactor);

    /** Calculates coefficients for a 1st order low-shelving with a gain at the Nyquist frequency
    that matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
    static void calculateFirstOrderLowShelvCoeffsPrescribedNyQuist(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& sampleRate, const double& frequency,
      const double& gainFactor);

    /** Calculates the coefficients for a first order high shelving type of filter (according to
    the formula from DAFX?). */
    static INLINE void calculateFirstOrderHighShelvCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& gainFactor);

    /** Calculates the coefficients for a high shelving type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookHighShelvCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q, const double& gainFactor);

    /** Calculates coefficients for a 1st order high-shelving with a gain at the Nyquist frequency
    that matches the gain of the corresponding analog prototype filter at that Nyquist frequency. */
    static void calculateFirstOrderHighShelvCoeffsPrescribedNyQuist(double& b0, double& b1,
      double& b2, double& a1, double& a2, const double& sampleRate, const double& frequency,
      const double& gainFactor);

    /** Calculates the coefficients for a first order allpass filter (according to the formula from
    DAFX?). */
    static INLINE void calculateFirstOrderAllpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency);

    /** Calculates the coefficients for a second order allpass type of filter according to Robert
    Bristow Johnson's Cookbook formula. */
    static INLINE void calculateCookbookAllpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q);

    /** Calculates coefficients for a biquad that represents a chain of a first order lowpass and a
    first order highpass filter. */
    static INLINE void calculateLowHighpassChainCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& highpassCutoff,
      const double& lowpassCutoff);

    /** Calculates the coefficients for a second order filter which is morphable between lowpass
    through peak to highpass according to some formulas derived by Robin Schmidt. */
    static INLINE void calculateLowPeakHighMorphCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, double& preGain, const double& sampleRate,
      const double& oneOverSampleRate, const double& frequency, const double& q,
      const double& morph, const bool &scaleGainForTwoStages);

    /** Calculates the coefficients for a peaking equalizer type of filter according to Sophocles
    Ofanidis' design procedure which matches the analog equalizers gain at the Nyquist
    frequency. */
    static INLINE void calculatePrescribedNyquistGainEqCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& bandwidthInOctaves, const double& gainFactor, const double& referenceGain);

    /** Calculates the coefficients for a one-pole filter where the pole is a number along the
    real axis in the z-plane (and should be less than 1 in magnitude for stable filters). */
    static INLINE void calculateOnePoleCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& pole);

    /** Calculates the coefficients for a one-zero filter where the pole is a number along the
    real axis in the z-plane. */
    static INLINE void calculateOneZeroCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& zero);

    /** Calculates the coefficients for a two-pole filter where the poles form a complex conjugate
    pair with the given radius and angle (the radius should be less than 1 in magnitude for stable
    filters). */
    static INLINE void calculatePolePairCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& poleRadius, const double& poleAngle);

    /** Calculates the coefficients for a two-zero filter where the zeros form a complex conjugate
    pair with the given radius and angle. */
    static INLINE void calculateZeroPairCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& zeroRadius, const double& zeroAngle);


     /** \name Miscellaneous */

    /** Returns the magnitude response of a biquad filter with given coefficients at a given
    frequency (and sample-rate). */
    static double getBiquadMagnitudeAt(const double& b0, const double& b1, const double& b2,
      const double& a1, const double& a2, const double& frequency, const double& sampleRate);

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void BiquadDesigner::makeBypassBiquad(double& b0, double& b1, double& b2, double& a1,
    double& a2)
  {
    b0 = 1.0;
    b1 = 0.0;
    b2 = 0.0;
    a1 = 0.0;
    a2 = 0.0;
  }

  INLINE void BiquadDesigner::calculateSineAndCosine(double& sinResult, double& cosResult,
      const double& frequency, const double& oneOverSampleRate)
  {
    // calculate intermediate variables:
    double omega  = 2.0 * PI * frequency * oneOverSampleRate;
    sinCos(omega, &sinResult, &cosResult);
  }

  INLINE void BiquadDesigner::calculateFirstOrderLowpassCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency)
  {
    double omega = 2.0 * PI * frequency * oneOverSampleRate;
    double x     = exp(-omega);

    a1 = x;
    a2 = 0.0;
    b0 = 1.0-x;
    b1 = 0.0;
    b2 = 0.0;
  }

  INLINE void BiquadDesigner::calculateFirstOrderLowpassCoeffsBilinear(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate, const double &frequency)
  {
    double omegaPreWarped = tan(PI * frequency * oneOverSampleRate);
    double pAnalog        = -omegaPreWarped;
    a1                    = (1.0+pAnalog)/(1.0-pAnalog);
    a2                    = 0.0;
    double g              = 0.5*sqrt(1.0+a1*a1-2.0*a1); // gain-factor for normalization at DC
    b0                    = g;
    b1                    = g;
    b2                    = 0.0;
  }

  /*
  INLINE void BiquadDesigner::calculateFirstOrderLowpassCoeffsPrescribedNyquist(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &sampleRate, const double &frequency)
  {
    double wc  = 2.0*PI*frequency/sampleRate;
    double Wc  = tan(wc/2.0);
    double Ws  = PI*sampleRate;
    double Wc2 = Wc*Wc;

    double Wca = 2*PI*frequency;                      // non-pre-warped analog cutoff frequency
    double k2  = 1.0 / (1.0 + ((Ws*Ws)/(Wca*Wca)) ); // gain of prototype at Nyquist-freq
    //k2  = 0; // test - should revert to standard bilinear design

    double A2  = 1.0 / (Wc2*(1.0-2.0*k2));
    double B2  = k2*A2;
    double g02 = 1.0;

    double A   = sqrt(A2);
    double B   = sqrt(B2);
    double g0  = sqrt(g02);      // == 1.0 -> optimize out
    double rD  = 1.0 / (1.0+A);  // reciprocal of denominator

    b0 =  rD * (g0+B);
    b1 =  rD * (g0-B);
    b2 =  0.0;
    a1 = -rD * (1 -A);
    a2 =  0.0;
  }
  */

  INLINE void BiquadDesigner::calculateCookbookLowpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine   * a0Rec;
    a2 = (alpha-1.0)  * a0Rec;
    b1 = (1.0-cosine) * a0Rec;
    b0 = 0.5*b1;
    b2 = b0;
  }

  INLINE void BiquadDesigner::calculateFirstOrderHighpassCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency)
  {
    double omega = 2.0 * PI * frequency * oneOverSampleRate;
    double x     = exp(-omega);

    a1 = x;
    a2 = 0.0;
    b0 = 0.5*(1.0+x);
    b1 = -b0;
    b2 = 0.0;
  }

  INLINE void BiquadDesigner::calculateFirstOrderHighpassCoeffsBilinear(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate, const double &frequency)
  {
    double omegaPreWarped = tan(PI * frequency * oneOverSampleRate);
    double pAnalog        = -omegaPreWarped;
    a1                    = (1.0+pAnalog)/(1.0-pAnalog);
    a2                    = 0.0;
    double g              = 0.5*sqrt(1.0+a1*a1+2.0*a1); // gain-factor for normalization at DC
    b0                    = g;
    b1                    = -g;
    b2                    = 0.0;
  }

  INLINE void BiquadDesigner::calculateCookbookHighpassCoeffs(double& b0, double& b1, double& b2,
      double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
      const double& q)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine    * a0Rec;
    a2 = (alpha-1.0)   * a0Rec;
    b1 = -(1.0+cosine) * a0Rec;
    b0 = -0.5*b1;
    b2 = b0;
  }

  INLINE void BiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine   * a0Rec;
    a2 = (alpha-1.0)  * a0Rec;
    b1 = 0.0;
    b0 = q*alpha      * a0Rec;
    b2 = -b0;
  }

  INLINE void BiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaBandwidth(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &bandwidth)
  {
    // to avoid divison by zero:
    double f = rmax(0.0001, frequency);
    double b = rmax(0.0001, bandwidth);

    double sine, cosine;
    calculateSineAndCosine(sine, cosine, f, oneOverSampleRate);

    // we need some if( sine == 0 )  ..check the limit
    double alpha = sine * sinh( 0.5*log(2.0) * b * 2.0*PI*f*oneOverSampleRate / sine );
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine   * a0Rec;
    a2 = (alpha-1.0)  * a0Rec;
    b1 = 0.0;
    b0 = 0.5*sine     * a0Rec;
    b2 = -b0;
  }

  INLINE void BiquadDesigner::calculateCookbookBandrejectCoeffsViaQ(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine  * a0Rec;
    a2 = (alpha-1.0) * a0Rec;
    b0 = 1.0         * a0Rec;
    b1 = -2.0*cosine * a0Rec;
    b2 = 1.0         * a0Rec;
  }

  INLINE void BiquadDesigner::calculateCookbookBandrejectCoeffsViaBandwidth(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &bandwidth)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine * sinh( 0.5*log(2.0) * bandwidth * 2.0*PI*frequency*oneOverSampleRate / sine );
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine  * a0Rec;
    a2 = (alpha-1.0) * a0Rec;
    b0 = 1.0         * a0Rec;
    b1 = -2.0*cosine * a0Rec;
    b2 = 1.0         * a0Rec;
  }

  INLINE void BiquadDesigner::calculateCookbookPeakFilterCoeffsViaQ(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q, const double &gainFactor)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double A     = gainFactor;
    double a0Rec = 1.0/(1.0+alpha/A);

    a1 = 2.0*cosine        * a0Rec;
    a2 = ((alpha/A) - 1.0) * a0Rec;
    b0 = (1.0+alpha*A)     * a0Rec;
    b1 = -2.0*cosine       * a0Rec;
    b2 = (1.0-alpha*A)     * a0Rec;
  }

  INLINE void BiquadDesigner::calculateFirstOrderLowShelvCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& oneOverSampleRate, const double& frequency,
    const double& gainFactor)
  {
    double A     = gainFactor;
    double omega = 2.0 * PI * frequency * oneOverSampleRate;
    double x     = tan(0.5*omega);
    double c     = 0.5*(A*A-1.0);  // A is the square root of the gain-factor -> V_0 = A^2
    double a     = 0.0;

    if( A > 1.0001 ) // boost
    {
      a  = (x-1)/(x+1);
      b0 = 1.0 + c*a + c;
      b1 = a   + c*a + c;
      b2 = 0.0;
      a1 = -a;
      a2 = 0.0;
    }
    else if( A < 0.999 ) // cut
    {
      a  = (x-A*A)/(x+A*A);
      b0 = 1.0 + c*a + c;
      b1 = a   + c*a + c;
      b2 = 0.0;
      a1 = -a;
      a2 = 0.0;
    }
    else
    {
      b0 = 1.0;
      b1 = 0.0;
      b2 = 0.0;
      a1 = 0.0;
      a2 = 0.0;
    }
  }

  INLINE void BiquadDesigner::calculateCookbookLowShelvCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q, const double &gainFactor)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double A     = gainFactor;
    double beta  = sqrt(A) / q;
    double a0Rec = 1.0 / ( (A+1.0) + (A-1.0)*cosine + beta*sine);

    a1 = 2.0 *     ( (A-1.0) + (A+1.0)*cosine              ) * a0Rec;
    a2 = -         ( (A+1.0) + (A-1.0)*cosine - beta*sine  ) * a0Rec;
    b0 =       A * ( (A+1.0) - (A-1.0)*cosine + beta*sine  ) * a0Rec;
    b1 = 2.0 * A * ( (A-1.0) - (A+1.0)*cosine              ) * a0Rec;
    b2 =       A * ( (A+1.0) - (A-1.0)*cosine - beta*sine  ) * a0Rec;
  }

  INLINE void BiquadDesigner::calculateFirstOrderHighShelvCoeffs(double &b0, double &b1,
    double &b2, double &a1, double &a2, const double &oneOverSampleRate,
    const double &frequency, const double &gainFactor)
  {
    double A     = gainFactor;
    double omega = 2.0 * PI * frequency * oneOverSampleRate;
    double x     = tan(0.5*omega);
    double c     = 0.5*(A*A-1.0);  // A is the square root of the gain-factor -> V_0 = A^2
    double a     = 0.0;

    if( A > 1.0001 ) // boost
    {
      a  = (x-1)/(x+1);
      b0 = 1.0 - c*a + c;
      b1 = a   + c*a - c;
      b2 = 0.0;
      a1 = -a;
      a2 = 0.0;
    }
    else if( A < 0.999 ) // cut
    {
      a  = (A*A*x-1)/(A*A*x+1);
      b0 = 1.0 - c*a + c;
      b1 = a   + c*a - c;
      b2 = 0.0;
      a1 = -a;
      a2 = 0.0;
    }
    else
    {
      b0 = 1.0;
      b1 = 0.0;
      b2 = 0.0;
      a1 = 0.0;
      a2 = 0.0;
    }
  }

  INLINE void BiquadDesigner::calculateCookbookHighShelvCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q, const double &gainFactor)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double A     = gainFactor;
    double beta  = sqrt(A) / q;
    double a0Rec = 1.0  /  ( (A+1.0) - (A-1.0)*cosine + beta*sine);

    a1 = -2.0 *     ( (A-1.0) - (A+1.0)*cosine              ) * a0Rec;
    a2 = -          ( (A+1.0) - (A-1.0)*cosine - beta*sine  ) * a0Rec;
    b0 =        A * ( (A+1.0) + (A-1.0)*cosine + beta*sine  ) * a0Rec;
    b1 = -2.0 * A * ( (A-1.0) + (A+1.0)*cosine              ) * a0Rec;
    b2 =        A * ( (A+1.0) + (A-1.0)*cosine - beta*sine  ) * a0Rec;
  }

  INLINE void BiquadDesigner::calculateFirstOrderAllpassCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency)
  {
    double t = tan(PI*frequency*oneOverSampleRate);
    double x = (t-1.0) / (t+1.0);

    b0 = x;
    b1 = 1.0;
    b2 = 0.0;
    a1 = -x;
    a2 = 0.0;
  }

  INLINE void BiquadDesigner::calculateCookbookAllpassCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, const double &oneOverSampleRate, const double &frequency,
    const double &q)
  {
    double sine, cosine;
    calculateSineAndCosine(sine, cosine, frequency, oneOverSampleRate);
    double alpha = sine/(2.0*q);
    double a0Rec = 1.0/(1.0+alpha);

    a1 = 2.0*cosine    * a0Rec;
    a2 = (alpha-1.0)   * a0Rec;
    b0 = (1.0-alpha)   * a0Rec;
    b1 = (-2.0*cosine) * a0Rec;
    b2 = (1.0+alpha)   * a0Rec;
  }

  INLINE void calculateLowHighpassChainCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& oneOverSampleRate, const double& highpassCutoff,
    const double& lowpassCutoff)
  {
    // variables for the coefficients of the 2 one-pole filters:
    double b0Lpf, b1Lpf, a1Lpf; // lowpass-coeffs
    double b0Hpf, b1Hpf, a1Hpf; // highpass-coeffs

    double x; // intermediate variable for calculation (x is the amount of
    // decay between adjacent samples for LPFs):

    // calculate lowpass coefficients (see dspguide for details):
    if( lowpassCutoff == 0.0 )
    {
      b0Lpf = 1.0;
      b1Lpf = 0.0;
      a1Lpf = 0.0;
    }
    else
    {
      x = exp( -2.0 * PI * lowpassCutoff * oneOverSampleRate);

      b0Lpf = 1-x;
      b1Lpf = 0.0;
      a1Lpf = -x;
    }

    // calculate highpass coefficients (see dspguide for details):
    if( highpassCutoff == 0.0 )
    {
      b0Hpf = 1.0;
      b1Hpf = 0.0;
      a1Hpf = 0.0;
    }
    else
    {
      x = exp( -2.0 * PI * highpassCutoff * oneOverSampleRate);

      b0Hpf =  0.5*(1+x);
      b1Hpf = -0.5*(1+x);
      a1Hpf = -x;
    }

    // combine the 2 one-pole filters into one biquad filter:
    b0 = b0Lpf*b0Hpf;
    b1 = b0Lpf*b1Hpf + b1Lpf*b0Hpf;
    b2 = b1Lpf*b1Hpf;
    a1 = -( a1Lpf+a1Hpf );
    a2 = -( a1Lpf*a1Hpf );
  }

  INLINE void BiquadDesigner::calculateLowPeakHighMorphCoeffs(double &b0, double &b1, double &b2,
    double &a1, double &a2, double& preGain, const double &sampleRate,
    const double &oneOverSampleRate, const double &frequency, const double &q,
    const double &morph, const bool &scaleGainForTwoStages)
  {
    // calculate normalized radian frequencies (digital and pre-warped analog):
    double wd = 2.0 * PI * frequency * oneOverSampleRate; // optimize to one mul
    double wa = 2.0 * sampleRate * tan(0.5*wd);

    // calculate the position of the analog pole-pair:
    double tmp1 = -wa / (2.0*q);
    double tmp2;
    if( q >= 0.5 )
      tmp2 = sqrt( -wa*wa*(0.25/(q*q) - 1.0) );
    else
    {
      tmp2 = 0.0;
      DEBUG_BREAK; // this design procedure assumes q >= 0.5
    }
    Complex pa = Complex(tmp1, tmp2);

    // calculate the position of the digital pole-pair:
    Complex z  = 0.5 * oneOverSampleRate * pa;
    Complex pd = (1.0+z) / (1.0-z);

    // convert poles to biquad feedback coefficients:
    a1 = 2.0 * pd.re;
    a2 = -(pd.re*pd.re + pd.im*pd.im);

    // calculate the (double) zero for a peak filter:
    double num  = (1.0-pd.re)*(1.0-pd.re) + pd.im*pd.im;
    double den  = (1.0+pd.re)*(1.0+pd.re) + pd.im*pd.im;
    double c    = sqrt( num / den );
    double z_pk = (1.0-c) / (1.0+c);

    // obtain the zero by interpolating the zero position between lowpass and peak or peak
    // and highpass:
    double m = morph;
    m = sign(m) * (m*m);
    //m = sign(m) * pow(abs(m), 1.5);
    double zd;
    //double a = -z_pk;
    if( m < -0.99 )
      zd = -0.99;
    else if( m > 0.99 )
      zd = +0.99;
    else
      zd = (m+z_pk) / (1.0+z_pk*m);

    // convert zeros to biquad feedforward coefficients:
    b0 = 1.0;
    b1 = -(zd+zd);
    b2 = zd*zd;

    // normalize gain:
    preGain = q / getBiquadMagnitudeAt(b0, b1, b2, a1, a2, frequency, sampleRate);
    if( scaleGainForTwoStages )
      preGain *= preGain;
  }

  INLINE void BiquadDesigner::calculatePrescribedNyquistGainEqCoeffs(double& b0,
    double& b1, double& b2, double& a1, double& a2, const double& oneOverSampleRate,
    const double& frequency, const double& bandwidthInOctaves, const double& gainFactor,
    const double& referenceGain)
  {
    if( fabs( RAPT::rsAmpToDb(gainFactor/referenceGain) ) < 0.001 )
    {
      makeBypassBiquad(b0, b1, b2, a1, a2);
      return;
    }

    double G0 = referenceGain;               // reference gain at DC
    double G  = gainFactor;                  // boost/cut gain
    double fc = frequency;                   // center frequency in Hz
    double bw = bandwidthInOctaves;          // bandwidth in octaves

    // intermediate variables:
    double fLo = fc  / pow(2.0,0.5*bw);      // lower bandedge frequency
    double fHi = fLo * pow(2.0, bw);         // upper bandedge frequency
    double GB  = sqrt(G);                    // gain at bandedge frequencies
    double w0  = 2*PI*fc*oneOverSampleRate;
    double wLo = 2*PI*fLo*oneOverSampleRate;
    double wHi = 2*PI*fHi*oneOverSampleRate;
    double Dw  = wHi-wLo;

    // the design procedure (ported from peq.m by Sophocles Orfanidis) - may be optimized:
    double F   = fabs(G*G   - GB*GB);
    double G00 = fabs(G*G   - G0*G0);
    double F00 = fabs(GB*GB - G0*G0);
    double num = G0*G0 * (w0*w0-PI*PI)*(w0*w0-PI*PI) + G*G * F00 * PI*PI * Dw*Dw / F;
    double den = (w0*w0-PI*PI)*(w0*w0-PI*PI) + F00 * PI*PI * Dw*Dw / F;
    double G1  = sqrt(num/den);

    //---------------------------------
    // inserted modification by Robin Schmidt: enforce InEq.37 to hold:
    {
      double Delta = 2*PI*(fHi-fLo);                              // bandwidth in 2*pi Hz
      double ASq   = Delta*Delta * (GB*GB-G0*G0) / (G*G-GB*GB);   // Eq.5

      // redefine bandedge gain to enforce InEq.37 to hold:
      GB = sqrt(G1*G);

      // re-calculate bandwidth from the new GB definition as the width where an analog equalizer
      // would have the redefined bandedge gain:
      Delta  = sqrt( ASq*(G*G-GB*GB) / (GB*GB-G0*G0) );
      Dw     = Delta*oneOverSampleRate;

      // re-calculate affected intermediate variables:
      F   = fabs(G*G   - GB*GB);
      F00 = fabs(GB*GB - G0*G0);
      num = G0*G0 * (w0*w0-PI*PI)*(w0*w0-PI*PI) + G*G * F00 * PI*PI * Dw*Dw / F;
      den = (w0*w0-PI*PI)*(w0*w0-PI*PI) + F00 * PI*PI * Dw*Dw / F;
      G1  = sqrt(num/den);
    }
    //---------------------------------

    double G01 = fabs(G*G   - G0*G1);
    double G11 = fabs(G*G   - G1*G1);
    double F01 = fabs(GB*GB - G0*G1);
    double F11 = fabs(GB*GB - G1*G1);
    double ta  = tan(w0/2);
    double W2  = sqrt(G11 / G00) * ta*ta;
    double DW  = (1 + sqrt(F00 / F11) * W2) * tan(Dw/2);
    double C   = F11 * DW*DW - 2 * W2 * (F01 - sqrt(F00 * F11));
    double D   = 2 * W2 * (G01 - sqrt(G00 * G11));
    double A   = sqrt((C + D) / F);
    double B   = sqrt((G*G * C + GB*GB * D) / F);
    double s   = 1.0 / (1 + W2 + A);

    b0         =  (G1 + G0*W2 + B)   * s;
    b1         =  (-2*(G1 - G0*W2))  * s;
    b2         =  (G1 - B + G0*W2)   * s;
    a1         = -(-2*(1 - W2))      * s;
    a2         = -(1 + W2 - A)       * s;
  }

  INLINE void BiquadDesigner::calculateOnePoleCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& pole)
  {
    b0 = 1.0;
    b1 = 0.0;
    b2 = 0.0;
    a1 = pole;
    a2 = 0.0;
  }

  INLINE void BiquadDesigner::calculateOneZeroCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& zero)
  {
    b0 = 1.0;
    b1 = -zero;
    b2 = 0.0;
    a1 = 0.0;
    a2 = 0.0;
  }

  INLINE void BiquadDesigner::calculatePolePairCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& poleRadius, const double& poleAngle)
  {
    b0 = 1.0;
    b1 = 0.0;
    b2 = 0.0;
    a1 = 2*poleRadius*cos(poleAngle);
    a2 = -poleRadius*poleRadius;
  }

  INLINE void BiquadDesigner::calculateZeroPairCoeffs(double& b0, double& b1, double& b2,
    double& a1, double& a2, const double& zeroRadius, const double& zeroAngle)
  {
    b0 = 1.0;
    b1 = -2*zeroRadius*cos(zeroAngle);
    b2 = zeroRadius*zeroRadius;
    a1 = 0.0;
    a2 = 0.0;
  }


} // end namespace rosic

#endif // rosic_BiquadDesigner_h
