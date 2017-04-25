#ifndef rosic_TwoPoleFilter_h
#define rosic_TwoPoleFilter_h

// rosic-indcludes:
#include "rosic_BiquadStereoDF1.h"
#include "rosic_BiquadStereoDF2.h"
#include "rosic_BiquadDesigner.h"

namespace rosic
{

  /**

  This class implements a TwoPole filter which can be used as band in a serial equalizer. For peak filter mode, the filter is designed by 
  Sophocles Orfanidis' design procedure to match the magnitude of an analog equalizer at the Nyquist frequency.

  */

  class TwoPoleFilter : public BiquadStereoDF2
  {

  public:

    /** Enumeration of the available filter modes. */
    enum modes
    {
      BYPASS = 0,
      PEAK,
      LOW_SHELF,
      HIGH_SHELF,
      LOWPASS6,
      LOWPASS12,
      HIGHPASS6,
      HIGHPASS12,
      BANDREJECT,
      BANDPASS,
      //LOWHIGH_PASS_CHAIN,

      REAL_POLE,
      REAL_ZERO,
      POLE_PAIR,
      ZERO_PAIR,

      NUM_FILTER_MODES
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    TwoPoleFilter();

    /** Destructor. */
    ~TwoPoleFilter();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate (in Hz) at which the filter runs. */
    void setSampleRate(double newSampleRate);

    /** Sets the filter mode as one of the values in enum modes. */
    void setMode(int newMode);

    /** Sets the center frequency in Hz. */
    void setFrequency(double newFrequency);

    /** Sets the boost/cut gain in dB. */
    void setGain(double newGain);

    /** Sets the bandwidth in octaves. */
    void setBandwidth(double newBandwidth);

    /** Sets the radius for the pole/zero filter types. */
    void setRadius(double newRadius);

    /** Sets up the frequency, gain and bandwidth at once. */
    void setParameters(int newMode, double newFrequency, double newGain, double newBandwidth, bool updateCoefficients = true);

    /** Sets the lower bandedge frequency by updating the bandwidth parameter (keeping the center frequency constant). */
    void setLowerBandedgeFrequency(double newLowerBandedgeFrequency);

    /** Sets the upper bandedge frequency by updating the bandwidth parameter (keeping the center frequency constant). */
    void setUpperBandedgeFrequency(double newUpperBandedgeFrequency);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Sets the filter mode as one of the values in enum modes. */
    int getMode() const { return mode; }

    /** Returns the center frequency in Hz. */
    double getFrequency() const { return frequency; }

    /** Returns the boost/cut gain in dB. */
    double getGain() const { return gain; }

    /** Returns the bandwidth in octaves. */
    double getBandwidth() const { return bandwidth; }

    /** Returns the magnitude of this equalizer band at a given frequency in Hz. */
    double getMagnitudeAt(double frequency) const 
    { return BiquadDesigner::getBiquadMagnitudeAt(b0, b1, b2, a1, a2, frequency, sampleRate); }

    /** Returns the lower bandedge frequency (defined at the half dB-gain point). */
    double getLowerBandedgeFrequency() const { return (frequency/pow(2.0, 0.5*bandwidth)); }

    /** Returns the upper bandedge frequency (defined at the half dB-gain point). */
    double getUpperBandedgeFrequency() const { return (frequency*pow(2.0, 0.5*bandwidth)); }

    /** Returns true if the current mode supports a bandwidth parameter, false otherwise. */
    bool doesModeSupportBandwidth() const;

    /** Returns true if the current mode supports a gain parameter, false otherwise. */
    bool doesModeSupportGain() const;

    /** Returns true if the current mode supports a frequency parameter, false otherwise. */
    bool doesModeSupportFrequency() const { return mode != BYPASS; }

    /** Returns true if the current mode supports a radius parameter, false otherwise. */
    bool doesModeSupportRadius() const;

    /** Converts a lower bandedge-frequency into the corresponding bandwidth (in octaves) for a given center frequency. */
    static double lowerBandedgeFrequencyToBandwdith(double lowerBandEdgefrequency, double centerFrequency);

    /** Converts an upper bandedge-frequency into the corresponding bandwidth (in octaves) for a given center frequency. */
    static double upperBandedgeFrequencyToBandwdith(double upperBandEdgefrequency, double centerFrequency);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Calculates the biquad coefficients from the specification. */
    void updateCoeffs();

    //=====================================================================================================================================

  protected:

    double frequency, gain, bandwidth, radius;
    double sampleRate;
    int    mode;

  };

} 

#endif 
