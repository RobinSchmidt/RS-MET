#ifndef rosic_DualTwoPoleFilter_h
#define rosic_DualTwoPoleFilter_h

//// rosic-indcludes:
//#include "rosic_TwoPoleFilter.h"

namespace rosic
{

  /**

  This class combines two TwoPoleFilter objects and allows continuous blend between serial and
  parallel connection.

  \todo: implement cutoff-scale and maybe gain-scale

  */

  class DualTwoPoleFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    DualTwoPoleFilter();

    /** Destructor. */
    ~DualTwoPoleFilter();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the blend factor between serial and parallel connection 
    (range 0...1, 0: fully serial, 1: fully parallel) */
    void setSerialParallelBlend(double newBlend);

    /** Sets the sample-rate (in Hz) at which the filter runs. */
    void setSampleRate(double newSampleRate)
    { 
      filter1.setSampleRate(newSampleRate);
      filter2.setSampleRate(newSampleRate);
    }

    /** Sets the mode for the first filter. @see: TwoPoleFilter::modes. */
    void setMode1(int newMode)
    { 
      filter1.setMode(newMode); 
    }

    /** Sets the mode for the second filter. @see: TwoPoleFilter::modes. */
    void setMode2(int newMode)
    { 
      filter2.setMode(newMode); 
    }

    /** Sets the frequency for the first filter in Hz. */
    void setFrequency1(double newFrequency)
    { 
      freq1 = newFrequency; 
      filter1.setFrequency(freqScale*freq1); 
    }

    /** Sets the frequency for the second filter in Hz. */
    void setFrequency2(double newFrequency)
    { 
      freq2 = newFrequency; 
      filter2.setFrequency(freqScale*freq2); 
    }

    /** Sets a scale factor to scale the frequencies of both filters. */
    void setFrequencyScale(double newScale)
    {
      freqScale = newScale; 
      filter1.setFrequency(freqScale*freq1);    
      filter2.setFrequency(freqScale*freq2);
    }

    /** Sets the boost/cut gain in dB for the first filter. */
    void setGain1(double newGain)
    { 
      gain1 = newGain; 
      filter1.setGain(gainScale*gain1); 
    }

    /** Sets the boost/cut gain in dB for the second filter. */
    void setGain2(double newGain)
    { 
      gain2 = newGain; 
      filter2.setGain(gainScale*gain2); 
    }

    /** Sets a scale factor to scale the gains of both filters. */
    void setGainScale(double newScale)
    {
      gainScale = newScale; 
      filter1.setGain(gainScale*gain1);    
      filter2.setGain(gainScale*gain2);
    }

    /** Sets the bandwidth in octaves for the first filter. */
    void setBandwidth1(double newBandwidth)
    { 
      bandwidth1 = newBandwidth; 
      filter1.setBandwidth(bandwidthScale*bandwidth1); 
    }

    /** Sets the bandwidth in octaves for the second filter. */
    void setBandwidth2(double newBandwidth)
    { 
      bandwidth2 = newBandwidth; 
      filter2.setBandwidth(bandwidthScale*bandwidth2); 
    }

    /** Sets a scale factor to scale the bandwidths of both filters. */
    void setBandwidthScale(double newScale)
    {
      bandwidthScale = newScale; 
      filter1.setBandwidth(bandwidthScale*bandwidth1);    
      filter2.setBandwidth(bandwidthScale*bandwidth2);
    }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the blend factor between serial and parallel connection 
    (range 0...1, 0: fully serial, 1: fully parallel) */
    double getSerialParallelBlend() const { return gp; }

    /** Returns the mode of the first filter @see TwoPoleFilter::modes. */
    int getMode1() const { return filter1.getMode(); }

    /** Returns the mode of the second filter @see TwoPoleFilter::modes. */
    int getMode2() const { return filter2.getMode(); }

    /** Returns the frequency of the first filter in Hz. */
    double getFrequency1() const { return freq1; }

    /** Returns the frequency of the second filter in Hz. */
    double getFrequency2() const { return freq2; }

    /** Returns the boost/cut gain of the first filter in dB. */
    double getGain1() const { return gain1; }

    /** Returns the boost/cut gain of the second filter in dB. */
    double getGain2() const { return gain2; }

    /** Returns the bandwidth of the first filter in octaves. */
    double getBandwidth1() const { return bandwidth1; }

    /** Returns the bandwidth of the second filter in octaves. */
    double getBandwidth2() const { return bandwidth2; }

    /** Returns the magnitude of this equalizer band at a given frequency in Hz. */
    //double getMagnitudeAt(double frequency) const
    //{ return BiquadDesigner::getBiquadMagnitudeAt(b0, b1, b2, a1, a2, frequency, sampleRate); }

    /** Returns true if the current mode of the first filter supports a bandwidth parameter, false 
    otherwise. */
    bool doesModeSupportBandwidth1() const { return filter1.doesModeSupportBandwidth(); }

    /** Returns true if the current mode of the second filter supports a bandwidth parameter, false 
    otherwise. */
    bool doesModeSupportBandwidth2() const { return filter2.doesModeSupportBandwidth(); }

    /** Returns true if the current mode of the first filter supports a gain parameter, false 
    otherwise. */
    bool doesModeSupportGain1() const { return filter1.doesModeSupportGain(); }

    /** Returns true if the current mode of the second filter supports a gain parameter, false 
    otherwise. */
    bool doesModeSupportGain2() const { return filter2.doesModeSupportGain(); }

    /** Returns true if the current mode of the first filter supports a frequency parameter, false 
    otherwise. */
    bool doesModeSupportFrequency1() const { return filter1.doesModeSupportFrequency(); }

    /** Returns true if the current mode of the second filter supports a frequency parameter, false 
    otherwise. */
    bool doesModeSupportFrequency2() const { return filter2.doesModeSupportFrequency(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Convenience function to be used when only mono output is desired. */
    INLINE double getSample(double x)
    {
      double dummy (0);
      getSampleFrameStereo(&x, &dummy);
      return x;
    }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Calculates the biquad coefficients from the specification. */
    void updateCoeffs() { filter1.updateCoeffs(); filter2.updateCoeffs(); }

    /** Resets the internal buffers of both filters. */
    void reset() { filter1.reset(); filter2.reset(); }

    //=============================================================================================

  protected:

    double freq1, freq2, freqScale, gain1, gain2, gainScale, bandwidth1, bandwidth2, bandwidthScale;

    TwoPoleFilter filter1, filter2;
    double gp;  // parallel gain

  };


  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void DualTwoPoleFilter::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double tmp1L = *inOutL;
    double tmp1R = *inOutR;

    filter1.getSampleFrameStereo(&tmp1L, &tmp1R);

    double tmp2L = (1.0-gp) * tmp1L + gp * (*inOutL);
    double tmp2R = (1.0-gp) * tmp1R + gp * (*inOutR);
      // maybe this can be optimized?

    filter2.getSampleFrameStereo(&tmp2L, &tmp2R);

    *inOutL = tmp2L + gp*tmp1L;
    *inOutR = tmp2R + gp*tmp1R;
  }


} // end namespace rosic

#endif // rosic_DualTwoPoleFilter_h
