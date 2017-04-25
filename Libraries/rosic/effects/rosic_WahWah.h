#ifndef rosic_WahWah_h
#define rosic_WahWah_h

// rosic-indcludes:
#include "rosic_ModulationEffect.h"
#include "../filters/rosic_TwoPoleFilter.h"

namespace rosic
{

  /**

  This class implements a wahwah effect with sinusoidal modulators for left and right channel
  separately.

  ToDo: use simple (bilinear) peak-filter coeff calculation

  */

  class WahWah : public ModulationEffect
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WahWah();   

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the depth of the LFO modulation in octaves. */
    void setDepth(double newDepth) { depth = newDepth; }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) { dryWet = newDryWet; }

    /** Sets the filter's mode @see TwoPoleFilter::modes. */
    void setFilterMode(int newMode);

    /** Sets the filter's nominal frequency in Hz. */
    void setFrequency(double newFrequency);

    /** Sets the boost/cut gain in dB. */
    void setGain(double newGain);

    /** Sets the bandwidth in octaves. */
    void setBandwidth(double newBandwidth);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns true if the current mode supports a bandwidth parameter, false otherwise. */
    bool doesModeSupportBandwidth() const { return filterL.doesModeSupportBandwidth();  }

    /** Returns true if the current mode supports a gain parameter, false otherwise. */
    bool doesModeSupportGain() const { return filterL.doesModeSupportGain(); }

    /** Returns true if the current mode supports a frequency parameter, false otherwise. */
    bool doesModeSupportFrequency() const { return filterL.doesModeSupportFrequency(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=============================================================================================

  protected:

    // embedded modules:
    TwoPoleFilter filterL, filterR;

    double pitch;
    double dryWet;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void WahWah::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double left, right;
    lfo.getSampleFrameStereo(&left, &right);
    double pitchL = pitch + 0.5*depth*left;
    double pitchR = pitch + 0.5*depth*right;
    //double pitchL = pitch + 0.5*depth*lfoL.getSample();
    //double pitchR = pitch + 0.5*depth*lfoR.getSample();

    filterL.setFrequency(pitchToFreq(pitchL));
    filterR.setFrequency(pitchToFreq(pitchR));

    double tmpL = filterL.getSample(*inOutL);
    double tmpR = filterR.getSample(*inOutR);

    *inOutL = (1.0-dryWet)*(*inOutL) + dryWet*tmpL;
    *inOutR = (1.0-dryWet)*(*inOutR) + dryWet*tmpR;
  }

} // end namespace rosic

#endif // #ifndef rosic_WahWah_h
