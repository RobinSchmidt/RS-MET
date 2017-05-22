#ifndef rosic_EngineersFilter_h
#define rosic_EngineersFilter_h

//// rosic-indcludes:
//#include "rosic_BiquadCascade.h"
//#include "rosic_InfiniteImpulseResponseDesigner.h"

namespace rosic
{

  /** This class implements scientific filters such as Butterworth, Chebychev I/II, elliptic, 
  Bessel, etc. You can create the standard lowpass, highpass, bandpass, bandreject  versions of
  all the filters and also shelving and peaking versions. */

  class EngineersFilter : public BiquadCascadeStereo
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    EngineersFilter();   

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate on which the filter should operate. */
    void setSampleRate(double newSampleRate);

    /** Selects the mode for the filter (lowpass, highpass, etc.) as defined in 
    InfiniteImpulseResponseDesigner::modes */
    void setMode(int newMode);

    /** Chooses one of the approximation methods (Butterworth, elliptic, etc.) as defined in 
    PrototypeDesigner::approximationMethods. */
    void setApproximationMethod(int newApproximationMethod);

    /** Sets the characteristic frequency of the filter. For lowpass and highpass filters, this is 
    the cutoff frequency, for bandpass, bandreject and peak filters the center frequency and for 
    Butterworth low-shelving and high shelving filters the half-gain frequency (for the other
    approximation methods, the ripple settings define the gain at that frequency). */
    void setFrequency(double newFrequency);

    /** Selects the order of the prototype (lowpass-, or low-shelving) filter. The order of the 
    final filter will be either the same number (for lowpass, highpass, lo-shelving, 
    high-shelving modes) or twice as high (for bandpass, bandreject and peak-filter modes). */
    void setPrototypeOrder(int newOrder);

    /** Sets the bandwidth of the filter in octaves. For Butterworth bandpass and bandreject 
    filters it is defined by the -3.01 points on both edges of the passband and for Butterworth 
    peak filters it is defined by half-gain frequencies. For other approximation methods, the 
    ripple settings are taken into account. */
    void setBandwidth(double newBandwidth);

    /** Sets the gain for peak and shelving modes (in decibels). */
    void setGain(double newGain);

    /** Sets the ripple in the passband for pass filter designs in decibels and also the 
    (in-band and/or out-of-band) ripple in the shelving filter designs in terms of a percentage 
    of the peak gain. */
    void setRipple(double newPassbandRipple);

    /** Sets the rejection in the stopband for pass-filters in decibels. */
    void setStopbandRejection(double newStopbandRejection);

    /** Selects the lower cutoff-/corner frequency for bandpass, bandreject and peaking filters or 
    the one and only cutoff-/corner frequency for lowpass, low-shelving, highpass and high-shelving
    modes. */
    //void setLowerFrequency(double newLowerFrequency);

    /** Selects the upper cutoff-/corner frequency for bandpass, bandreject and peaking filters. 
    For other filter-modes, this is irrelevant. */
    //void setUpperFrequency(double newUpperFrequency);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the approximation method to be used 
    @see enum PrototypeDesigner::approximationMethods. */
    int getApproximationMethod() { return designer.getApproximationMethod(); }

    /** Returns true if the currently selected mode supports a bandwidth parameter. */
    bool hasCurrentModeBandwidthParameter() { return designer.hasCurrentModeBandwidthParameter(); }

    /** Returns true if the currently selected mode supports a gain parameter. */
    bool hasCurrentModeGainParameter() { return designer.hasCurrentModeGainParameter(); }

    /** Returns true if the currently selected mode supports a ripple parameter. */
    bool hasCurrentModeRippleParameter() { return designer.hasCurrentModeRippleParameter(); }

    /** Returns true if the currently selected mode supports a rejection parameter. */
    bool hasCurrentModeRejectionParameter() { return designer.hasCurrentModeRejectionParameter(); }

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array 
    "frequencies" (in Hz) and stores them in the array "magnitudes". Both arrays are assumed to be 
    "numBins" long. "inDecibels" indicates, if the frequency response should be returned in 
    decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be 
    multiplied with (or added to, when "inDecibels" is true) to the magnitudes which are already 
    there in the "magnitudes"-array. This is useful for calculating the magnitude response of 
    several biquad-cascades in series. */
    void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
      bool inDecibels = false, bool accumulate = false);

  protected:

    /** Triggers a re-calculation of the biquad coefficients. */
    void updateCoefficients();

    InfiniteImpulseResponseDesigner designer;

    double sampleRate;

  };

}

#endif 
