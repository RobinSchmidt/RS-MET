#ifndef rosic_InfiniteImpulseResponseDesigner_h
#define rosic_InfiniteImpulseResponseDesigner_h

namespace rosic
{

/** This class designs (calculates coefficients for) a high order digital infinite impulse response 
(IIR) filter a la Orfanidis. */

class rsInfiniteImpulseResponseDesigner
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
    LOW_SHELV,
    HIGH_SHELV,
    PEAK
  };

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsInfiniteImpulseResponseDesigner();

  /** Destructor. */
  ~rsInfiniteImpulseResponseDesigner();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample-rate on which the filter should operate. */
  void setSampleRate(double newSampleRate);

  /** Selects the mode for the filter. */
  void setMode(int newMode);

  /** Chooses one of the approximation methods as defined in 
  PrototypeDesigner::approximationMethods. */
  void setApproximationMethod(int newApproximationMethod);

  /** Selects the order of the prototype (lowpass-, or low-shelving) filter. The order of the final 
  filter will be either the same number (for lowpass, highpass, lo-shelving, high-shelving modes) 
  or twice as high (for bandpass, bandreject and peak-filter modes). */
  void setPrototypeOrder(int newOrder);

  /** Sets the characteristic frequency of the filter. For lowpass and highpass filters, this is 
  the cutoff frequency, for bandpass, bandreject and peak filters the center frequency and for 
  Butterworth low-shelving and high shelving filters the half-gain frequency (for the approximation 
  methods, the ripple settings define the gain at that frequency). It calls setLowerFrequency and
  setUpperFrequency (which you also may use yourself alternatively) */
  void setFrequency(double newFrequency);

  /** Sets the bandwidth of the filter in octaves. For Butterworth bandpass and bandreject filters 
  it is defined by the -3.01 points on both edges of the passband and for Butterworth peak filters 
  it is defined by half-gain frequencies. For other approximation methods, the ripple settings are 
  taken into account. It calls setLowerFrequency and setUpperFrequency (which you also may use 
  yourself alternatively) */
  void setBandwidth(double newBandwidth);

  /** Selects the lower cutoff-/corner frequency for bandpass, bandreject and peaking filters or 
  the one and only cutoff-/corner frequency for lowpass, low-shelving, highpass and high-shelving 
  modes. */
  void setLowerFrequency(double newLowerFrequency);

  /** Selects the upper cutoff-/corner frequency for bandpass, bandreject and peaking filters. For 
  other filter-modes, this is irrelevant. */
  void setUpperFrequency(double newUpperFrequency);

  /** Sets the gain for peak and shelving modes (in decibels). */
  void setGain(double newGain);

  /** Sets the ripple in the passband for pass filter designs in decibels and also the (in-band 
  and/or out-of-band) ripple in the shelving filter designs in terms of a percentage of the peak 
  gain. */
  void setRipple(double newPassbandRipple);

  /** Sets the rejection in the stopband for lowpass designs in decibels. */
  void setStopbandRejection(double newStopbandRejection);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the approximation method to be used
  @see enum PrototypeDesigner::approximationMethods. */
  int getApproximationMethod() { return prototypeDesigner.getApproximationMethod(); }

  /** Some filter modes will result in a final filter that has twice the order of the prototype 
  filter. This function return true when the selected mode is one of those. */
  bool doesModeDoubleTheOrder();

  /** Returns the order of the final filter - this may either be equal to the order of the 
  prototype filter or twice this value (depending on the selected mode). */
  int getFinalFilterOrder();

  /** Returns the required number of biquad stages - this may be equal to the getFinalOrder()/2 
  or (getFinalOrder()+1)/2 depending on the even- or oddness of the final order. */
  int getNumBiquadStages();

  /** Returns the passbandRipple in dB. */
  double getPassbandRipple() const { return prototypeDesigner.getPassbandRipple(); }

  /** Returns true if the currently selected mode supports a bandwidth parameter. */
  bool hasCurrentModeBandwidthParameter();

  /** Returns true if the currently selected mode supports a gain parameter. */
  bool hasCurrentModeGainParameter();

  /** Returns true if the currently selected mode/method combination supports a ripple 
  parameter. */
  bool hasCurrentModeRippleParameter();

  /** Returns true if the currently selected mode/method combination supports a rejection 
  parameter. */
  bool hasCurrentModeRejectionParameter();

  //-----------------------------------------------------------------------------------------------
  /** \name Coefficient Retrieval */

  /** Returns the z-domain poles and zeros. */
  void getPolesAndZeros(Complex* poles, Complex* zeros);

  /** Calculates and stores the coefficients for a biquad cascade which realizes the filter with 
  the desired specifications. */
  void getBiquadCascadeCoefficients(double *b0, double *b1, double *b2, double *a1, double *a2);

  /** Calculates and returns the coefficients for a direct-form implementation. */
  void getDirectFormCoefficients(double *b, double *a);

  //===============================================================================================

protected:

  /** Calculates the lower and upper characteristic frequencies according to frequency and 
  bandwidth settings. */
  void calculateLowerAndUpperFrequency();

  /** Normalizes the biquad stages such that each stage has unit magnitude response at some 
  properly chosen normalized radian frequency, taking into account all the ripple, method, mode 
  settings. wc is the center frequency which is only relevant for bandpass/bandreject and peak 
  modes. */
  void normalizeGain(double *b0, double *b1, double *b2, double *a1, double *a2, double wc, 
    int numBiquads);

  /// Ensures that all frequencies have meainingful
  //void frequenciesSanityCheck();

  double sampleRate;
  double frequency, bandwidth;
  double lowerFrequency, upperFrequency;
  double gain;                              // gain in decibels (for shelving- and peak modes)
  int    mode;
  int    prototypeOrder;

  // embedded objects:
  rsPrototypeDesigner prototypeDesigner;

  /*
  double getBiquadMagnitudeAt(double b0, double b1, double b2,
    double a1, double a2, double omega);*/
  /**< Calculates the magnitude-response of a digital biquad filter with coefficeints b0, b1, b2,
  a0, a1, a1 at the normalized radian frequency 'omega'. */

  /*
  void normalizeBiquadStages(double *b0, double *b1, double *b2,
    double *a1, double *a2, double omega, int numStages);*/
  /**< Normalizes the biquad stages described by the given coefficients in such a way that each
  stage has unit magnitude at the normalized radian frequency 'omega'. */

};

}

#endif
