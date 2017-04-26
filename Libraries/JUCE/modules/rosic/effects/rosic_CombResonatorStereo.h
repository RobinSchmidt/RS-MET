#ifndef rosic_CombResonatorStereo_h
#define rosic_CombResonatorStereo_h

// rosic-indcludes:
#include "../filters/rosic_CombResonator.h"

namespace rosic
{

  /**

  This class implements a stereo comb resonator effect.

  */

  class CombResonatorStereo
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    CombResonatorStereo();   

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { combL.setDryWetRatio(newDryWet); combR.setDryWetRatio(newDryWet); }

    /** Sets the level for the wet signal. */
    void setLevel(double newLevel) { gain = dB2amp(newLevel); }

    /** Sets the filter's nominal frequency in Hz. */
    void setFrequency(double newFrequency) 
    { 
      pitch = freqToPitch(newFrequency); 
      combL.setFrequency( pitchToFreq(pitch+pitchOffsetL) );
      combR.setFrequency( pitchToFreq(pitch+pitchOffsetR) );
    }

    /** Sets a offset between the the pitches of the two combs (for left and right channel)
    in semitones. */
    void setDetune(double newDetune)
    {
      pitchOffsetL = +0.5*newDetune;
      pitchOffsetR = -0.5*newDetune;
      combL.setFrequency( pitchToFreq(pitch+pitchOffsetL) );
      combR.setFrequency( pitchToFreq(pitch+pitchOffsetR) );
    }

    /** Sets the panorama position for the 1st comb (between -1...+1). */
    void setPan1(double newPan1) { pan1 = 0.5*(newPan1+1.0); }

    /** Sets the panorama position for the 2nd comb (between -1...+1). */
    void setPan2(double newPan2) { pan2 = 0.5*(newPan2+1.0); }

    /** Sets the time for the comb to decay to -60 dB (in seconds). */
    void setDecayTime(double newDecayTime) 
    { combL.setDecayTime(newDecayTime); combR.setDecayTime(newDecayTime); }

    /** Scales the decay-time for the high-frequency band. */
    void setHighDecayScale(double newScale)
    { combL.setHighDecayScale(newScale); combR.setHighDecayScale(newScale); }

    /** Scales the decay-time for the low-frequency band. */
    void setLowDecayScale(double newScale)
    { combL.setLowDecayScale(newScale); combR.setLowDecayScale(newScale); }

    /** Sets the crossover-frequency between mid and high frequencies. */
    void setHighCrossoverFreq(double newFreq)
    { combL.setHighCrossoverFreq(newFreq); combR.setHighCrossoverFreq(newFreq); }

    /** Sets the crossover-frequency between low and mid frequencies. */
    void setLowCrossoverFreq(double newFreq)
    { combL.setLowCrossoverFreq(newFreq); combR.setLowCrossoverFreq(newFreq); }

    /** Switches the resonator into a mode where it produces only odd harmonics. */
    void setOddOnlyMode(bool shouldCreateOnlyOddHarmonics)
    { 
      combL.setOddOnlyMode(shouldCreateOnlyOddHarmonics); 
      combR.setOddOnlyMode(shouldCreateOnlyOddHarmonics); 
    }

    /** Sets the polarity of the wet signal negative (or positive, if false). */
    void setNegativePolarity(bool shouldBeNegative) 
    { combL.setNegativePolarity(shouldBeNegative); combR.setNegativePolarity(shouldBeNegative); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo ouput frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=============================================================================================

  protected:

    // embedded modules:
    CombResonator combL, combR;

    double pitch, pitchOffsetL, pitchOffsetR;
    double pan1, pan2, gain;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void CombResonatorStereo::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double tmp1 = gain * combL.getSample(*inOutL);
    double tmp2 = gain * combR.getSample(*inOutR);

    *inOutL = (1.0-pan1)*tmp1 + (1.0-pan2)*tmp2;
    *inOutR =      pan1 *tmp1 +      pan2 *tmp2;
  }

} // end namespace rosic

#endif // #ifndef rosic_CombResonatorStereo_h
