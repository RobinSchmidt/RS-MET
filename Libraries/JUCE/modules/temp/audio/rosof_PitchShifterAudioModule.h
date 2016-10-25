#ifndef rosof_PitchShifterAudioModule_h
#define rosof_PitchShifterAudioModule_h

#include "../../../rosic/effects/rosic_PitchShifterGrainAdaptive.h"
using namespace rosic;

#include "../rosof_AudioModule.h"

namespace rosof
{

  /**

  This class wraps rosic::PitchShifter into a rosof::AudioModule to facilitate its use as plugIn.

  */

  class PitchShifterAudioModule : public AudioModule
  {

    friend class PitchShifterModuleEditor;

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    PitchShifterAudioModule(CriticalSection *newPlugInLock, rosic::PitchShifterGrainAdaptive *pitchShifterToWrap);   

    //-------------------------------------------------------------------------------------------------------------------------------------
    // automation and state management:

    /** Creates the static parameters for this module (i.e. parameters that are not created dynamically and are thus always there). */
    virtual void createStaticParameters();

    /** Restores the state of this module from an XmlElement (which was presumably previously created via getStateAsXml). */
    //virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean);

    /** Converts a state which might possibly be from an older version to the current patch-format. */
    virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState); 

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    virtual void setSampleRate(double newSampleRate) 
    { wrappedPitchShifter->setSampleRate(newSampleRate); }

    virtual void setBeatsPerMinute(double newBpm) 
    { wrappedPitchShifter->setBeatsPerMinute(newBpm); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
    { wrappedPitchShifter->getSampleFrameStereo(inOutL, inOutR); }

    virtual void processBlockStereo(float *left, float *right, int numSamples)
    { 
      double dL, dR;
      for(int n=0; n<numSamples; n++)
      {
        dL = (double) left[n];
        dR = (double) right[n];
        getSampleFrameStereo(&dL, &dR);
        left[n]  = (float) dL;
        right[n] = (float) dR;
      }
    }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // event processing:

    virtual void reset() 
    { wrappedPitchShifter->reset(); }

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    rosic::PitchShifterGrainAdaptive *wrappedPitchShifter; 

  };

}

#endif 
