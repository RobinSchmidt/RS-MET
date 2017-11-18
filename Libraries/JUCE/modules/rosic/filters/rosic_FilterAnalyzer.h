#ifndef rosic_FilterAnalyzer_h
#define rosic_FilterAnalyzer_h

namespace rosic
{

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

class rsFilterAnalyzer
{

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
  static Complex getAnalogFrequencyResponseAt(Complex* z, Complex* p, double k, int N, double w);

  /** Returns the magnitude response of an "N"th analog filter with zeros, poles and gain given in
  "z", "p", "k" at the radian frequency "w". */
  static double getAnalogMagnitudeResponseAt(Complex* z, Complex* p, double k, int N, double w);

  /** Writes the magnitude response of an "N"th analog filter with zeros, poles and gain given in 
  "z", "p", "k" at the radian frequencies given in "w" into the array "m". Arrays "w" and "m" 
  should be of length "numBins". */
  static void getAnalogMagnitudeResponse(Complex* z, Complex* p, double k, int N, double* w, 
    double* m, int numBins);

  /** Writes the phase response of an "N"th analog filter with zeros, poles and gain given in 
  "z", "p", "k" at the radian frequencies given in "w" into the array "phs". Arrays "w" and "phs" 
  should be of length "numBins". */
  static void getAnalogPhaseResponse(Complex* z, Complex* p, double k, int N, double* w, 
    double* phs, int numBins);

  // \todo: getAnalogPhaseDelay, getAnalogGroupDelay

  /** Given the zeros, poles and gain-constant (z,p,k) of an Nth order analog filter transfer 
  function, this function finds the radian frequency at which the filter's magnitude response takes 
  on the desired value passed in "magnitude". The caller should ensure that the filter actually 
  takes on the desired magnitude value somewhere, otherwise the retrun value will be meaningless. 
  You may pass an initial guess for the radian frequency. */
  static double findAnalogFrequencyWithMagnitude(Complex* z, Complex* p, double* k, int N, 
    double magnitude, double initialGuess = 1.0);

  //-----------------------------------------------------------------------------------------------
  // functions for analyzing (cascaded) biquad filters:

  /** Returns the magnitude response of a biquad filter that realizes the difference equation:
  \f[ y[n] = b_0*x[n] + b_1*x[n-1] + b_2*x[n-2] - a_1*y[n-1] - a_2*y[n-2] \f]
  with given coefficients at a given normalized radian frequency w. */
  static double getBiquadMagnitudeAt(const double b0, const double b1, const double b2, 
    const double a1, const double a2, const double w);

  // \todo getBiquadTransferFunctionAt, getBiquadFrequencyResponseAt, getBiquadPhaseResponseAt

  /** Writes the magnitude response of a biquad-cascade at the normalized radian frequencies given 
  in 'w' into the array 'mag'. */
  static void getBiquadMagnitudeResponse(const double b0, const double b1, const double b2, 
    const double a1, const double a2, double* w, double* mag, int numBins, 
    bool inDecibels = false);

  /** Returns the value of the transfer-function of a biquad at the given value 'z'. */
  static Complex getBiquadTransferFunctionAt(const double b0, const double b1, const double b2, 
    const double a1, const double a2, const Complex z);

  /** Returns the value of the transfer-function of a biquad-cascade at the given value 'z'. */
  static Complex getBiquadCascadeTransferFunctionAt(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, Complex z);

  /** Writes the complex frequency-response of a biquad-cascade at the normalized radian 
  frequencies given in 'w' into the array 'H'. */
  static void getBiquadCascadeFrequencyResponse(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, double* w, Complex* H, int numBins, 
    int accumulationMode = NO_ACCUMULATION);

  /** Multiplies a given frequency-response 'H' at the normalized radian frequencies 'w' with the 
  frequency-response of a biquad-cascade. This function is useful for accumulating 
  frequency-responses of filters that are connected in series. */
  static void multiplyWithBiquadCascadeFrequencyResponse(double* b0, double* b1, double* b2, 
    double* a1, double* a2, int numBiquads, double* w, Complex* H, int numBins);

  /** Adds a given frequency-response 'H' at the normalized radian frequencies 'w' with the 
  frequency-response of a biquad-cascade. This function is useful for accumulating 
  frequency-responses of filters that are connected in parallel. */
  static void addWithBiquadCascadeFrequencyResponse(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, double* w, Complex* H, int numBins);

  /** Writes the magnitude response of a biquad-cascade at the normalized radian frequencies given 
  in 'w' into the array 'magnitudes'. */
  static void getBiquadCascadeMagnitudeResponse(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, double* w, double* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false);

  /** Writes the magnitude response of a biquad-cascade at the physical frequencies given in 
  'frequencies' into the array 'magnitudes'. */
  static void getBiquadCascadeMagnitudeResponse(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, double* frequencies, double sampleRate, double* magnitudes, 
    int numBins, bool inDecibels = false, bool accumulate = false);


  /** Returns the value of the (unwrapped) phase-response of a digital biquad at the given 
  normalized radian frequency 'w'. The returned value will be in the range -2*PI...0. */
  static double getBiquadPhaseResponseAt(double b0, double b1, double b2, double a1, double a2, 
    double w);

  /** Returns the value of the (unwrapped) phase-response of a digital biquad-cascade at the given 
  normalized radian frequency 'w'. */
  static double getBiquadCascadePhaseResponseAt(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads, double w);

  /** Writes the (unwrapped) phase response of a digital biquad-cascade at the normalized radian 
  frequencies given in 'w' into the array 'phases'. */
  static void getBiquadCascadePhaseResponse(double* b0, double* b1, double* b2, double* a1, 
    double* a2, int numBiquads,  double* w, double* phases, int numBins, bool accumulate = false);


  //-----------------------------------------------------------------------------------------------
  // functions for analyzing direct-form filters:


  //-----------------------------------------------------------------------------------------------
  // some general helper functions:

  /** Extracts the magnitudes from the complex-values in H and writed them int 'magnitudes'. */
  static void getMagnitudes(Complex* H, double* magnitudes, int length);

  /** Extracts the phases from the complex-values in H and writed them int 'phases'. */
  static void getPhases(Complex* H, double* phases, int length);

  /** Converts an array of values (presumably magnitudes) to decibels. Because this conversion may 
  lead to negative infinite (in case of zero amplitude) or undefined (in case of negative 
  amplitude) values, the amplitude-values may be clipped (from bottom) to some small positive value. 
  The default value of 0.0000000001 corresponds to -200 dB. */
  static void convertToDecibels(double* values, int length, 
    double clipLowAmplitudeAt = 0.0000000001);

  /** Clamps all values (in the array 'values') at frequencies higher than sampleRate/2 to the 
  given constant 'clampValue' - this is useful for 'clipping' irrelevant data which otherwise would 
  have been plotted. */
  static void clampValuesAboveNyquist(double* frequencies, double* values, int length, 
    double sampleRate, double clampValue);

};

}

#endif 
