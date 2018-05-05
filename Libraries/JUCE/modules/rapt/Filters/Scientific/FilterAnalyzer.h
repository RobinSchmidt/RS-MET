#ifndef RAPT_FILTERANALYZER_H_INCLUDED
#define RAPT_FILTERANALYZER_H_INCLUDED


// \todo move these standalonefunctions into class rsFilterAnalyzer

/** Returns the magnitude of a one-pole filter of the form  
y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]  
at the normalized radian frequency w */
template<class T>
T onePoleMagnitudeAt(T b0, T b1, T a1, T w);

/** Returns the magnitude of a biquad filter of the form
y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-1] - a1*y[n-1] - a2*y[n-2]
at the normalized radian frequency w = 2*PI*f/fs, f: frequency, fs: samplerate */
template<class T>
T biquadMagnitudeAt(T b0, T b1, T b2, T a1, T a2, T w);

// \todo: check if these formulas really correspond to the difference equations with negative 
// sign convention on the feedback coefficients

/** Given a set of biquad coefficients, this function will determine whether or not the biquad 
is stable and minimum phase (all poles and zeros must be inside or on the unit circle). Strictly 
speaking, when there are poles on the unit circle, the system is only marginally stable and not 
really stable in a bounded-input/bounded-output (BIBO) sense. */
template<class T>
bool isBiquadStableAndMinimumPhase(T b0, T b1, T b2, T a1, T a2);

/** Returns the squared magnitude of an analog biquad filter with s-domain transfer function H(s) 
and magnitude-squared function G^2(w) given by:
        B0 + B1*s + B2*s^2             B0^2 + (B1^2-2*B0*B2)*w^2 + B2^2*w^4
H(s) = --------------------, G^2(w) = --------------------------------------
        A0 + A1*s + A2*s^2             A0^2 + (A1^2-2*A0*A2)*w^2 + A2^2*w^4
at the (radian) frequency w = 2*pi*f, f: frequency in Hz. */
template<class T>
T analogBiquadMagnitudeSquaredAt(T B0, T B1, T B2, T A0, T A1, T A2, T w);


/*
Analog biquad prototype transfer functions (from RBJ cookbook):

Pass-filters:
LPF: H(s) = 1     / (s^2 + s/Q + 1)
HPF: H(s) = s^2   / (s^2 + s/Q + 1)
BPF: H(s) = s     / (s^2 + s/Q + 1), constant skirt gain, peak gain = Q
BPF: H(s) = (s/Q) / (s^2 + s/Q + 1), constant 0 dB peak gain

Bandreject and allpass:
BRF: H(s) = (s^2 + 1) / (s^2 + s/Q + 1)
APF: H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)

Shelving filters (low, high, band),  A  = sqrt(10^(dBgain/20)):
LSF: H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)   / (A*s^2 + (sqrt(A)/Q)*s + 1)
HSF: H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1) / (s^2 + (sqrt(A)/Q)*s + A)
BSF: H(s) = (s^2 + s*(A/Q) + 1)             / (s^2 + s/(A*Q) + 1)
*/





/** This class implements a bunch of static functions to analyze various aspects of filters. The 
functions use the following conventions for naming the parameters: 'f' denotes a physical frequency 
in Hz, 'fs' denotes a sample-rate in Hz, 'w' a denotes normalized radian frequency (w = 2*pi*f/fs), 
'b' denote feedforward coefficients, 'a' denote feedback coefficients.

Feedback coefficients are assumed to be applied with negative sign such that, for example, a first 
order filter implements the difference equation: \f[ y[n] = b_0*x[n] + b_1*x[n-1] - a_1*y[n-1] \f]

Frequency-response functions generally need an array of normalized radian frequencies be passed (or 
alternatively an array of physical frequencies together with a sample-rate), and an array of the 
same length where the results are written into.

Some of them may optionally work accumulatively, which means that they either add or multiply the 
array contents with what is already there - this is useful for obtaining the overall frequency 
response for cascaded or parallel connections of filters.

\todo: split into 2 classes: AnalogFilterAnalyzer, DigitalFilterAnalyzer and get rid of the 
 Analog/Digital qualifiers in the member function names */

template<class T>
class rsFilterAnalyzer
{

  typedef std::complex<T> Complex; // preliminary

public:

  /** Enumeration of accumulation modes that are used by some functions to conveniently obtain 
  frequency responses of filters that are made from series- or parallel connections of simpler 
  filters. */
  enum accumulationModes
  {
    NO_ACCUMULATION = 0,
    MULTIPLICATIVE_ACCUMULATION,
    ADDITIVE_ACCUMULATION
  };

  //-----------------------------------------------------------------------------------------------
  // functions for analyzing filters in zero/pole/gain representation:

  /** Returns the complex frequency response of an "N"th analog filter with zeros, poles and gain 
  given in "z", "p", "k" at the radian frequency "w". */
  static Complex getAnalogFrequencyResponseAt(Complex* z, Complex* p, T k, int N, T w);

  /** Returns the magnitude response of an "N"th analog filter with zeros, poles and gain given in
  "z", "p", "k" at the radian frequency "w". */
  static T getAnalogMagnitudeResponseAt(Complex* z, Complex* p, T k, int N, T w);

  /** Writes the magnitude response of an "N"th analog filter with zeros, poles and gain given in 
  "z", "p", "k" at the radian frequencies given in "w" into the array "m". Arrays "w" and "m" 
  should be of length "numBins". */
  static void getAnalogMagnitudeResponse(Complex* z, Complex* p, T k, int N, T* w, T* m, 
    int numBins);

  /** Writes the phase response of an "N"th analog filter with zeros, poles and gain given in 
  "z", "p", "k" at the radian frequencies given in "w" into the array "phs". Arrays "w" and "phs" 
  should be of length "numBins". */
  static void getAnalogPhaseResponse(Complex* z, Complex* p, T k, int N, T* w, T* phs, 
    int numBins);

  // \todo: getAnalogPhaseDelay, getAnalogGroupDelay

  /** Given the zeros, poles and gain-constant (z,p,k) of an Nth order analog filter transfer 
  function, this function finds the radian frequency at which the filter's magnitude response takes 
  on the desired value passed in "magnitude". The caller should ensure that the filter actually 
  takes on the desired magnitude value somewhere, otherwise the retrun value will be meaningless. 
  You may pass an initial guess for the radian frequency. */
  static T findAnalogFrequencyWithMagnitude(Complex* z, Complex* p, T* k, int N, T magnitude, 
    T initialGuess = 1.0);

  //-----------------------------------------------------------------------------------------------
  // functions for analyzing (cascaded) biquad filters:

  /** Returns the magnitude response of a biquad filter that realizes the difference equation:
  \f[ y[n] = b_0*x[n] + b_1*x[n-1] + b_2*x[n-2] - a_1*y[n-1] - a_2*y[n-2] \f]
  with given coefficients at a given normalized radian frequency w. */
  static T getBiquadMagnitudeAt(const T b0, const T b1, const T b2, const T a1, const T a2, 
    const T w);

  // \todo getBiquadTransferFunctionAt, getBiquadFrequencyResponseAt, getBiquadPhaseResponseAt

  /** Writes the magnitude response of a biquad-cascade at the normalized radian frequencies given 
  in 'w' into the array 'mag'. */
  static void getBiquadMagnitudeResponse(const T b0, const T b1, const T b2, const T a1, 
    const T a2, T* w, T* mag, int numBins, bool inDecibels = false);

  /** Returns the value of the transfer-function of a biquad at the given value 'z'. */
  static Complex getBiquadTransferFunctionAt(const T b0, const T b1, const T b2, const T a1, 
    const T a2, const Complex z);

  /** Returns the value of the transfer-function of a biquad-cascade at the given value 'z'. */
  static Complex getBiquadCascadeTransferFunctionAt(T* b0, T* b1, T* b2, T* a1, T* a2, 
    int numBiquads, Complex z);

  /** Writes the complex frequency-response of a biquad-cascade at the normalized radian 
  frequencies given in 'w' into the array 'H'. */
  static void getBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, T* a2, int numBiquads, 
    T* w, Complex* H, int numBins, int accumulationMode = NO_ACCUMULATION);

  /** Multiplies a given frequency-response 'H' at the normalized radian frequencies 'w' with the 
  frequency-response of a biquad-cascade. This function is useful for accumulating 
  frequency-responses of filters that are connected in series. */
  static void multiplyWithBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
    int numBiquads, T* w, Complex* H, int numBins);

  /** Adds a given frequency-response 'H' at the normalized radian frequencies 'w' with the 
  frequency-response of a biquad-cascade. This function is useful for accumulating 
  frequency-responses of filters that are connected in parallel. */
  static void addWithBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
    int numBiquads, T* w, Complex* H, int numBins);

  /** Writes the magnitude response of a biquad-cascade at the normalized radian frequencies given 
  in 'w' into the array 'magnitudes'. */
  static void getBiquadCascadeMagnitudeResponse(T* b0, T* b1, T* b2, T* a1, T* a2, int numBiquads, 
    T* w, T* magnitudes, int numBins, bool inDecibels = false, bool accumulate = false);

  /** Writes the magnitude response of a biquad-cascade at the physical frequencies given in 
  'frequencies' into the array 'magnitudes'. */
  static void getBiquadCascadeMagnitudeResponse(T* b0, T* b1, T* b2, T* a1, T* a2, int numBiquads, 
    T* frequencies, T sampleRate, T* magnitudes, int numBins, bool inDecibels = false, 
    bool accumulate = false);


  /** Returns the value of the (unwrapped) phase-response of a digital biquad at the given 
  normalized radian frequency 'w'. The returned value will be in the range -2*PI...0. */
  static T getBiquadPhaseResponseAt(T b0, T b1, T b2, T a1, T a2, T w);

  /** Returns the value of the (unwrapped) phase-response of a digital biquad-cascade at the given 
  normalized radian frequency 'w'. */
  static T getBiquadCascadePhaseResponseAt(T* b0, T* b1, T* b2, T* a1, T* a2, int numBiquads, T w);

  /** Writes the (unwrapped) phase response of a digital biquad-cascade at the normalized radian 
  frequencies given in 'w' into the array 'phases'. */
  static void getBiquadCascadePhaseResponse(T* b0, T* b1, T* b2, T* a1, T* a2, int numBiquads, 
    T* w, T* phases, int numBins, bool accumulate = false);


  //-----------------------------------------------------------------------------------------------
  // functions for analyzing direct-form filters:


  //-----------------------------------------------------------------------------------------------
  // some general helper functions:

  /** Extracts the magnitudes from the complex-values in H and writed them int 'magnitudes'. */
  static void getMagnitudes(Complex* H, T* magnitudes, int length);

  /** Extracts the phases from the complex-values in H and writed them int 'phases'. */
  static void getPhases(Complex* H, T* phases, int length);

  /** Converts an array of values (presumably magnitudes) to decibels. Because this conversion may 
  lead to negative infinite (in case of zero amplitude) or undefined (in case of negative 
  amplitude) values, the amplitude-values may be clipped (from bottom) to some small positive value. 
  The default value of 0.0000000001 corresponds to -200 dB. */
  static void convertToDecibels(T* values, int length, 
    T clipLowAmplitudeAt = 0.0000000001);

  /** Clamps all values (in the array 'values') at frequencies higher than sampleRate/2 to the 
  given constant 'clampValue' - this is useful for 'clipping' irrelevant data which otherwise would 
  have been plotted. */
  static void clampValuesAboveNyquist(T* frequencies, T* values, int length, T sampleRate, 
    T clampValue);

};

#endif 
