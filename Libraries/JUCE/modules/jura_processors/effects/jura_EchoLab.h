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

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
    {
      wrappedEchoLabDelayLine->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n],
        &inOutBuffer[0][n], &inOutBuffer[1][n], true);
    }
  }

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

class EchoLabDelayLineModuleEditor : public AudioModuleEditor, public ChangeBroadcaster, 
  public AudioModuleDeletionWatcher
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  EchoLabDelayLineModuleEditor(CriticalSection *newPlugInLock, 
    EchoLabDelayLineAudioModule* newEchoLabDelayLineAudioModule);

  //-----------------------------------------------------------------------------------------------
  // setup:


  virtual void setDelayLineModuleToEdit(EchoLabDelayLineAudioModule* newEchoLabDelayLineModuleToEdit);
  virtual void setHueOffsetForFilterEditors(float hueOffset);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual callback function inherited from AudioModuleDeletionWatcher 
  in order to invalidate our pointer-member "delayLineModuleToEdit" when the delayline will be 
  deleted. */
  virtual void audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted);

  virtual void paint(Graphics &g);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  virtual void updateWidgetVisibility();


  EchoLabDelayLineAudioModule* delayLineModuleToEdit;

  EqualizerModuleEditor *inputEqualizerEditor, *feedbackEqualizerEditor;

  juce::Rectangle<int> middleRectangle;
  RSlider *timeSlider, *gainSlider, *feedbackSlider, *panSlider;
  RButton *pingPongButton, *muteButton, *soloButton, *flushButton;

  
  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class EchoLabAudioModule : public AudioModule
{

  friend class EchoLabModuleEditor;
  friend class EchoLabPlotEditor;

public:

  EchoLabAudioModule(CriticalSection *newPlugInLock, rosic::EchoLab *delayDesignerToWrap);   

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean);
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //virtual XmlElement EchoLabAudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState);

  // delegations to rosic::EchoLab withh added thread-safety and possible additional actions:
  virtual void setSampleRate(double newSampleRate);
  virtual void setBeatsPerMinute(double newBpm);
  int addDelayLine(double newDelayTime, double newGainFactor);

  bool removeDelayLine(int index);

  void removeAllDelayLines();

  /** Returns a pointer to the delayline audiomodule with given index. */
  EchoLabDelayLineAudioModule* getDelayLineModule(int index) const;

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR);
  virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

  virtual void reset();

protected:

  virtual void initializeAutomatableParameters();

  //CriticalSection wrappedEchoLabLock;

  rosic::EchoLab              *wrappedEchoLab;

  juce::Array<EchoLabDelayLineAudioModule*> delayLineModules;

  //EchoLabDelayLineAudioModule               *selectedDelayLineModule;

  //EqualizerAudioModule *inputEqualizerModule;
  //EqualizerAudioModule *feedbackEqualizerModule;
  //...nah -> goes to EchoLabDelayLineAudioModule

  //juce::Array<EqualizerAudioModule*> inputEqualizerModules;
  //juce::Array<EqualizerAudioModule*> feedbackEqualizerModules;

  //juce::Array<EchoLabDelayLineAudioModule*> echoLabDelayLineModules;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif 
