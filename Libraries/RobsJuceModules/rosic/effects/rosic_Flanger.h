#ifndef rosic_Flanger_h
#define rosic_Flanger_h

//// rosic-indcludes:
//#include "rosic_ModulationEffect.h"
//#include "../filters/rosic_CombFilter.h"

namespace rosic
{

  /**

  This class implements a flanger effect based on a comb-filter with sinusoidal modulation.

  todo (maybe): 
  -make LFO switchable to linear mode

  */

  class Flanger : public ModulationEffect
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Flanger();   

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the depth of the LFO modulation in octaves. */
    void setDepth(double newDepth) { depth = newDepth; }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { combL.setDryWetRatio(newDryWet); combR.setDryWetRatio(newDryWet); }

    /** Sets the filter's nominal frequency in Hz. */
    void setFrequency(double newFrequency) { pitch = RAPT::rsFreqToPitch(newFrequency); }

    /** Sets the feedback around the comb as raw factor. */
    void setFeedbackFactor(double newFeedbackFactor)
    { combL.setFeedbackFactor(newFeedbackFactor); combR.setFeedbackFactor(newFeedbackFactor); }

    /** Sets the feedback around the comb as in percent. */
    void setFeedbackInPercent(double newFeedback) { setFeedbackFactor(0.01*newFeedback); }

    /** Sets the polarity of the wet signal negative (or positive, if false). */
    void setNegativePolarity(bool shouldBeNegative) 
    { combL.setNegativePolarity(shouldBeNegative); combR.setNegativePolarity(shouldBeNegative); }

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
    CombFilter combL, combR;

    double pitch;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Flanger::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double left, right;
    lfo.getSampleFrameStereo(&left, &right);
    double pitchL = pitch + 0.5*depth*left;
    double pitchR = pitch + 0.5*depth*right;
    //double pitchL = pitch + 0.5*depth*lfoL.getSample();
    //double pitchR = pitch + 0.5*depth*lfoR.getSample();

    combL.setFrequency(RAPT::rsPitchToFreq(pitchL));
    combR.setFrequency(RAPT::rsPitchToFreq(pitchR));

    *inOutL = combL.getSample(*inOutL);
    *inOutR = combR.getSample(*inOutR);
  }

} // end namespace rosic

#endif // #ifndef rosic_Flanger_h
