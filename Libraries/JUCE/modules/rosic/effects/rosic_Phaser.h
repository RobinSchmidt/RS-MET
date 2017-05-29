#ifndef rosic_Phaser_h
#define rosic_Phaser_h

//// rosic-indcludes:
//#include "rosic_ModulationEffect.h"
//#include "../filters/rosic_AllpassChain.h"

namespace rosic
{

  /**

  This class implements a phaser effect based on a chain of first or second order allpass filters.

  \todo (maybe): 
  -polarity switch for wet signal (maybe for left and right separately)
  -make LFO switchable to linear mode

  */

  class Phaser : public ModulationEffect
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Phaser();   

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the depth of the LFO modulation in octaves. */
    void setDepth(double newDepth) { depth = newDepth; }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) { dryWet = newDryWet; }

    /** Sets the filter's mode @see AllpassChain::modes. */
    void setFilterMode(int newMode);

    /** Sets the filter's nominal frequency in Hz. */
    void setFrequency(double newFrequency);

    /** Sets the Q-factor for second order allpass mode. */
    void setQ(double newQ);

    /** Sets the feedback around the allpass chain as raw factor. */
    void setFeedbackFactor(double newFeedbackFactor);

    /** Sets the feedback around the comb as in percent. */
    void setFeedbackInPercent(double newFeedback) { setFeedbackFactor(0.01*newFeedback); }

    /** Sets the number of allpass stages to use. */
    void setNumStages(int newNumStages) 
    { 
      allpassL.setNumStages(newNumStages); 
      allpassR.setNumStages(newNumStages); 
    }

    /** Bypasses the effect entirely. */
    void setBypass(bool shouldBeBypassed) { bypass = shouldBeBypassed; }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns true if the current mode supports a Q parameter, false otherwise. */
    bool doesModeSupportQ() const { return allpassL.doesModeSupportQ(); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=====================================================================================================================================

  protected:

    // embedded modules:
    AllpassChain allpassL, allpassR;

    double pitch, feedbackFactor, yL, yR;
    double dryWet;
    bool   bypass;

  };

  //-------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void Phaser::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    if( bypass == true )
      return;

    double left, right;
    lfo.getSampleFrameStereo(&left, &right);
    double pitchL = pitch + 0.5*depth*left;
    double pitchR = pitch + 0.5*depth*right;
    //double pitchL = pitch + 0.5*depth*lfoL.getSample();
    //double pitchR = pitch + 0.5*depth*lfoR.getSample();

    allpassL.setFrequency(pitchToFreq(pitchL));
    allpassR.setFrequency(pitchToFreq(pitchR));

    yL = allpassL.getSample(*inOutL + feedbackFactor*yL);
    yR = allpassR.getSample(*inOutR + feedbackFactor*yR);

    *inOutL = (1.0-dryWet)*(*inOutL) + dryWet*yL;
    *inOutR = (1.0-dryWet)*(*inOutR) + dryWet*yR;
  }

}

#endif
