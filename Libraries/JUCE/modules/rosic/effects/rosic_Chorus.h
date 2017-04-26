#ifndef rosic_Chorus_h
#define rosic_Chorus_h

// rosic-indcludes:
#include "rosic_Vibrato.h"
//#include "../filters/rosic_OnePoleFilter.h"

namespace rosic
{

  /**

  This class implements a chorus effect effect by extending the Vibrato class by some more LFOs 
  and ouput taps.

  */

  class Chorus : public Vibrato
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - will allocate a delay-buffer with a given maximum number of samples delay. */
    Chorus(int bufferLengthToAllocate = 65536);

    /** Destructor */
    ~Chorus();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Switches one of the voices on or off. */
    void activateVoice(int voiceIndex, bool shouldBeActive) 
    { voiceIsActive[voiceIndex] = shouldBeActive; calculateGainCompensation(); }

    /** Scales the delaytime of one of the voices by the given factor. */
    void setVoiceDelayScale(int voiceIndex, double newScale) { delayScales[voiceIndex] = newScale; }

    /** Scales the modulation depth of one of the voices by the given factor. */
    void setVoiceDepthScale(int voiceIndex, double newScale) { depthScales[voiceIndex] = newScale; }

    /** Scales the amplitude of one of the voices by the given factor. */
    void setVoiceAmpScale(int voiceIndex, double newScale) 
    { ampScales[voiceIndex] = newScale; calculateGainCompensation();}

    /** Sets an overall factor for feedback around all (active) voices. This has to be refined 
    before becoming useful - currently it sounds like shit. */
    void setGlobalFeedback(double newFeedback) { globalFeedback = newFeedback; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    /** Calculates the member variable voiceGainCompensator from the number of active voices and 
    their amplitude scale factors. */
    void calculateGainCompensation();

    static const int numVoices = 4;
    double delayScales[numVoices];
    double depthScales[numVoices];
    double ampScales[numVoices];
    double voiceGainCompensator;
    double globalFeedback; 
    double yOldL, yOldR;
    bool   voiceIsActive[numVoices];

    //OnePoleFilter feedbackFilter;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Chorus::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    //feedIn(*inOutL+globalFeedback*yOldL, *inOutR+globalFeedback*yOldR);
    feedIn(*inOutL, *inOutR);

    double m1 = lfo.getSample();
    double m3 = -2*m1*m1+1;
    double m5 = -2*m3*m3+1;
    double m7 = -2*m5*m5+1;

    double yL = 0.0;                // accumulators for wet signal
    double yR = 0.0;
 
    double dL, dR;                  // instantaneous delays

    // accumulate the outputs of the voices:
    if( voiceIsActive[0] )
    {
      dL  = delayScales[0]*dA + depthScales[0]*d*m1;   
      dR  = delayScales[0]*dA - depthScales[0]*d*m1;  
      yL += ampScales[0]*getLeftOutputHermiteAt(dL);
      yR += ampScales[0]*getRightOutputHermiteAt(dR);
    }
    if( voiceIsActive[1] )
    {
      dL  = delayScales[1]*dA + depthScales[1]*0.5*d*m3;   
      dR  = delayScales[1]*dA - depthScales[1]*0.5*d*m3;   
      yL += ampScales[1]*getLeftOutputHermiteAt(dL);
      yR += ampScales[1]*getRightOutputHermiteAt(dR);
    }
    if( voiceIsActive[2] )
    {
      dL  = delayScales[2]*dA + depthScales[2]*0.25*d*m5;   
      dR  = delayScales[2]*dA - depthScales[2]*0.25*d*m5;   
      yL += ampScales[2]*getLeftOutputHermiteAt(dL);
      yR += ampScales[2]*getRightOutputHermiteAt(dR);
    }
    if( voiceIsActive[3] )
    {
      dL  = delayScales[3]*dA + depthScales[3]*0.125*d*m7;   
      dR  = delayScales[3]*dA - depthScales[3]*0.125*d*m7;   
      yL += ampScales[3]*getLeftOutputHermiteAt(dL);
      yR += ampScales[3]*getRightOutputHermiteAt(dR);
    }

    // compensate for the loudness gain due to the adding of multiple voices:
    yR *= voiceGainCompensator;
    yL *= voiceGainCompensator;

    //feedbackFilter.getSampleFrameStereo(&yL, &yR, &yOldL, &yOldR);
    //yOldL = clip(yL, -1.0, 1.0);
    //yOldR = clip(yR, -1.0, 1.0);

    *inOutL = (1.0-dryWetRatio)*(*inOutL) + dryWetRatio*yL;
    *inOutR = (1.0-dryWetRatio)*(*inOutR) + dryWetRatio*yR;

    incrementWritePointer();
  }

} // end namespace rosic

#endif // #ifndef rosic_Chorus_h
