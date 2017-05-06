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

#endif 
