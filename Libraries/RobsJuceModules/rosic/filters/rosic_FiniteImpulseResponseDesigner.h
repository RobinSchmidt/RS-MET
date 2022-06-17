#ifndef rosic_FiniteImpulseResponseDesigner_h
#define rosic_FiniteImpulseResponseDesigner_h

namespace rosic
{

/** This class designs impulse responses that can be used for FIR filtering using the windowing 
method.

  \ todo: implement spline transition

*/

class FiniteImpulseResponseDesigner
{

public:

  /** Enumeration of the available filter modes. */
  enum modes
  {
    BYPASS = 0,
    LOWPASS,
    HIGHPASS,
    BANDPASS,
    BANDREJECT,
    DIFFERENTIATOR,
    HILBERT,
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FiniteImpulseResponseDesigner();

  /** Destructor. */
  ~FiniteImpulseResponseDesigner();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Selects the mode for the filter. */
  void setMode(int newMode);

  /** Sets the sample-rate on which the filter should operate. */
  void setSampleRate(double newSampleRate);

  /** Sets the characteristic frequency of the filter (if applicable). */
  void setFrequency(double newFrequency);

  /** Sets the bandwidth of the filter in octaves (if applicable). */
  void setBandwidth(double newBandwidth);

  /** Selects the lower cutoff-/corner frequency for bandpass- and bandreject-filters. */
  void setLowerFrequency(double newLowerFrequency);

  /** Selects the upper cutoff-/corner frequency for bandpass- and bandreject-filters. */
  void setUpperFrequency(double newUpperFrequency);

  /** Sets the window that is to be used for the design. @see WindowDesigner::windowTypes */
  void setWindowType(int newWindow);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the sample-rate. */
  double getSampleRate() const { return sampleRate; }

  /** Returns true if the currently selected mode supports a bandwidth parameter. */
  bool hasCurrentModeBandwidthParameter();

  /** Generates the impulse-response of desired length according to the settings and writes it
  into the passed array. */
  void getImpulseResponse(double* impulseResponse, int length);

  // void getMagnitudeResponse ...nah, this goes into FiniteImpulseResponseFilter

  //-----------------------------------------------------------------------------------------------
  // design functions:

  /** Writes a lowpass response with normalized radian cutoff frequency 'omega' into the passed 
  array. Currently, only odd lengths are supported. */
  static void getLowpassResponse(double* impulseResponse, int length, double omega, 
    int windowType = WindowDesigner::BLACKMAN);

  /** Writes a highpass response with normalized radian cutoff frequency 'omega' into the passed 
  array. Currently, only odd lengths are supported. */
  static void getHighpassResponse(double* impulseResponse, int length, double omega, 
    int windowType = WindowDesigner::BLACKMAN);

  /** Writes a bandpass response with normalized radian bandedge frequencies 'omegaLow', 
  'omegaHigh' into the passed array. Currently, only odd lengths are supported. */
  static void getBandpassResponse(double* impulseResponse, int length, double omegaLow, 
    double omegaHigh, int windowType = WindowDesigner::BLACKMAN);

  /** Writes a bandreject response with normalized radian bandedge frequencies 'omegaLow',
  'omegaHigh' into the passed array. Currently, only odd lengths are supported. */
  static void getBandrejectResponse(double* impulseResponse, int length, double omegaLow, double omegaHigh,
    int windowType = WindowDesigner::BLACKMAN);

  /** Writes a differentiator response into the passed array - it is recommended to use filters
  with even lengths because the differentiator is antisymmetric and odd lengths together with 
  antisymmetry (aka type III) enforce zeros at z=1 and z=-1 (DC and Nyquist-frequency), whereas
  even lengths together with antisymmetry (aka type IV) enforce a zero only at z=1 (DC). As the
  differentiator is highpass in nature, a zero at z=-1 is undesirable and leads to a bad 
  approximation of the magnitude response. On the other hand, an even length filter will return a
  signal that approximates the derivative at time instants that are shifted by an (additional) 
  half-sample from the original sampling instants. */
  static void getDifferentiatorResponse(double* impulseResponse, int length,
    int windowType = WindowDesigner::BLACKMAN);

  /** Writes a Hilbert transformer response into the passed array - this approximates a phase-shift
  of 90 degrees at all frequencies.
  Currently, only odd lengths are supported. */
  static void getHilbertTransformerResponse(double* impulseResponse, int length, int windowType = WindowDesigner::BLACKMAN);

  /** Performs spectral inversion of the filter by changing the signs of the coefficients (aka 
  impulse response) and adding one at the center of symmetry. The resulting filter is equivalent to
  a parallel connection of the original filter (with negative sign) and an appropriate compensation
  delay. This amounts to flipping the frequency response top-for-bottom, thus turning a lowpass 
  into a highpass with the same cutoff frqeuency, for example. For this approach to make sense, the
  impulse response must have left/right symmetry and the length must be an odd number samples 
  (otherwise there is no sample at the center of symmetry). */
  static void spectralInversion(double* impulseResponse, int length);

  /** Performs spectral reversal of the filter by changing the sign of the even coefficients. This is equivalent to modulating the
  impulse-response by a sinusoid at sampleRate/2 which shifts the whole frequency respose by sampleRate/2. This turns a lowpass into a
  highpass, where the cutoff-frequency of the higpass is sampleRate/2 minus the cutoff-frequency of the original lowpass. */
  static void spectralReversal(double* impulseResponse, int length);

  /** Normalizes the sum of the impulse-response samples to unity.  For a lowpass filter, this amounts to a normalization of the DC gain
  to unity. */
  static void normalizeSumToUnity(double* impulseResponse, int length);

  //=====================================================================================================================================

protected:

  /** Calculates the lower and upper characteristic frequencies according to frequency and bandwidth settings. */
  void calculateLowerAndUpperFrequency();

  double sampleRate;
  double frequency, bandwidth;
  double lowerFrequency, upperFrequency;
  int    mode, windowType;

};

}  // end namespace rosic

#endif 
