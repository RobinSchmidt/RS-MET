#ifndef rosof_PolyphonicInstrumentAudioModule_h
#define rosof_PolyphonicInstrumentAudioModule_h

#include "../../rosic/infrastructure/rosic_PolyphonicInstrument.h"
using namespace rosic;

#include "rosof_AudioModule.h"

namespace rosof
{

  /**

  This class wraps rosic::PolyphonicInstrument into a rosof::AudioModule to facilitate its use as plugIn.

  */

  class PolyphonicInstrumentAudioModule : public AudioModule
  {

    friend class PolyphonicInstrumentModuleEditor;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock, rosic::PolyphonicInstrument *instrumentToWrap);   

    //---------------------------------------------------------------------------------------------
    // automation and state management:

    /** Overrides the parameterChanged() method of the indirect AutomationListener base class in 
    order to respond to interesting automation events. */
    virtual void parameterChanged(rojue::Parameter* parameterThatHasChanged);

    /** Returns the state (i.e. the settings of all relevant parameters) in form of an 
    XmlElement. */
    //virtual XmlElement* getStateAsXml(XmlElement* xmlElementToStartFrom = NULL);

    /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
    //virtual void setStateFromXml(const XmlElement& xmlState);

    /** Overriden to call allNotesOff before restoring the state. */
    virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean);

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    virtual void setSampleRate(double newSampleRate) 
    { 
      if( underlyingRosicInstrument != NULL )
        underlyingRosicInstrument->setSampleRate(newSampleRate); 
    }

    virtual void setBeatsPerMinute(double newBpm)
    {
      if( underlyingRosicInstrument != NULL )
        underlyingRosicInstrument->setBeatsPerMinute(newBpm);
    }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    virtual void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR)
    { 
      if( underlyingRosicInstrument != NULL )
        underlyingRosicInstrument->getSampleFrameStereo(outL, outR); 
    }

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Fills the array of automatable parameters. */
    void initializeAutomatableParameters();

  };

}

#endif
