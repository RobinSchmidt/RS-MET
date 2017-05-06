#ifndef jura_EchoLab_h
#define jura_EchoLab_h

/** This class wraps rosic::EchoLabDelayLine into a rosof::AudioModule. */

class EchoLabDelayLineAudioModule : public AudioModule, public ChangeBroadcaster
{

public:


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  EchoLabDelayLineAudioModule(CriticalSection *newPlugInLock, 
    rosic::EchoLabDelayLine *echoLabDelayLineToWrap);

  //-----------------------------------------------------------------------------------------------
  // setup:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  virtual EqualizerAudioModule* getInputEqualizerModule()    const { return inputEqualizerModule; }
  virtual EqualizerAudioModule* getFeedbackEqualizerModule() const { return feedbackEqualizerModule; }


  //-----------------------------------------------------------------------------------------------
  // audio processing:

  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR);
  //virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

  //-----------------------------------------------------------------------------------------------
  // others:

  //virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
  //  bool markAsClean);
    // temporarily overriden

  virtual void reset();

protected:

  void initializeAutomatableParameters();

  rosic::EchoLabDelayLine *wrappedEchoLabDelayLine;

  EqualizerAudioModule    *inputEqualizerModule;
  EqualizerAudioModule    *feedbackEqualizerModule;

  friend class EchoLabDelayLineModuleEditor;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

#endif 
