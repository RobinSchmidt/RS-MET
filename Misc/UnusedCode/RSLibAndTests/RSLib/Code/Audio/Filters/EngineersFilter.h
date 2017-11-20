#ifndef RS_ENGINEERSFILTER_H
#define RS_ENGINEERSFILTER_H

namespace RSLib
{

  /**

  This class implements classical filter design methods such as Butterworth, Chebychev, elliptic, 
  etc. and realizes the filters in the digital domain using a biquad cascade filter structure.

  */

  class RSLib_API rsEngineersFilter : public rsBiquadCascadeStereo
  {

  public:


    /** \name Construction/Destruction */

    /** Constructor. */
    rsEngineersFilter();   


    /** \name Setup */

    /** Sets the sample-rate on which the filter should operate. */
    void setSampleRate(double newSampleRate);

    /** Selects the mode for the filter defined in rsInfiniteImpulseResponseDesigner::modes. */
    void setMode(int newMode);

    /** Chooses one of the approximation methods as defined in 
    rsPrototypeDesigner::approximationMethods. */
    void setApproximationMethod(int newApproximationMethod);

    /** Sets the characteristic frequency of the filter. For lowpass and highpass filters, this is 
    the cutoff frequency, for bandpass, bandreject and peak filters the center frequency and for 
    Butterworth low-shelving and high shelving filters the half-gain frequency (for the 
    approximation methods, the ripple settings define the gain at that frequency). It calls 
    setLowerFrequency and setUpperFrequency (which you also may use yourself alternatively) */
    void setFrequency(double newFrequency);

    /** Selects the order of the prototype (lowpass-, or low-shelving) filter. The order of the 
    final filter will be either the same number (for lowpass, highpass, lo-shelving, high-shelving 
    modes) or twice as high (for bandpass, bandreject and peak-filter modes). */
    void setPrototypeOrder(int newOrder);

    /** Sets the bandwidth of the filter in octaves. For Butterworth bandpass and bandreject 
    filters it is defined by the -3.01 points on both edges of the passband and for Butterworth 
    peak filters it is defined by half-gain frequencies. For other approximation methods, the 
    ripple settings are taken into account. It calls setLowerFrequency and setUpperFrequency (which 
    you also may use yourself alternatively) */
    void setBandwidth(double newBandwidth);

    /** Sets the gain for peak and shelving modes (in decibels). */
    void setGain(double newGain);

    /** Sets the ripple in the passband for pass filter designs in decibels and also the (in-band 
    and/or out-of-band) ripple in the shelving filter designs in terms of a percentage of the peak 
    gain. */
    void setRipple(double newPassbandRipple);

    /** Sets the rejection in the stopband for lowpass designs in decibels. */
    void setStopbandRejection(double newStopbandRejection);

    /** Selects the lower cutoff-/corner frequency for bandpass, bandreject and peaking filters or 
    the one and only cutoff-/corner frequency for lowpass, low-shelving, highpass and 
    high-shelving modes. */
    //void setLowerFrequency(double newLowerFrequency);

    /** Selects the upper cutoff-/corner frequency for bandpass, bandreject and peaking filters. 
    For other filter-modes, this is irrelevant. */
    //void setUpperFrequency(double newUpperFrequency);


    /** \name Inquiry */

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
      // factor out into rsBiquadCascade

  protected:

    /** Triggers a re-calculation of the biquad coefficients. */
    void updateCoefficients();

    /** Make inherited method unavailable for client code, such that it doesn't accidently call it
    when it should actually call setPrototypeOrder. */
    void setOrder(int newOrder) { rsBiquadCascadeStereo::setOrder(newOrder); }

    rsInfiniteImpulseResponseDesigner designer;

    double sampleRate;



  };

}

#endif 
