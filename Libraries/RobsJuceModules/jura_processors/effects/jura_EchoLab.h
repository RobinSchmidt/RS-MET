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
  virtual EqualizerAudioModule* getFeedbackEqualizerModule() const 
  { return feedbackEqualizerModule; }

  //-----------------------------------------------------------------------------------------------
  // audio processing:

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

  rosic::EchoLabDelayLine* wrappedEchoLabDelayLine;
  EqualizerAudioModule* inputEqualizerModule;
  EqualizerAudioModule* feedbackEqualizerModule;

  friend class EchoLabDelayLineModuleEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EchoLabDelayLineAudioModule)
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

  virtual void setDelayLineModuleToEdit(
    EchoLabDelayLineAudioModule* newEchoLabDelayLineModuleToEdit);
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

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EchoLabDelayLineModuleEditor)
};

//=================================================================================================

class EchoLabAudioModule : public AudioModule
{

  friend class EchoLabModuleEditor;
  friend class EchoLabPlotEditor;

public:

  EchoLabAudioModule(CriticalSection *newPlugInLock, rosic::EchoLab *echoLabToWrap = nullptr);   

  virtual ~EchoLabAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void setBeatsPerMinute(double newBpm) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void reset() override;

  // delegations to rosic::EchoLab withh added thread-safety and possible additional actions:
  int addDelayLine(double newDelayTime, double newGainFactor);
  bool removeDelayLine(int index);
  void removeAllDelayLines();

  /** Returns a pointer to the delayline audiomodule with given index. */
  EchoLabDelayLineAudioModule* getDelayLineModule(int index) const;

protected:

  virtual void createParameters();

  rosic::EchoLab *wrappedEchoLab;
  bool wrappedEchoLabIsOwned = false;

  juce::Array<EchoLabDelayLineAudioModule*> delayLineModules; // use std::vector


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class plots a schematic impulse-response of a rosic::EchoLab object and allows for editing 
parameters like the delay-time and gains. */

class EchoLabPlotEditor	: virtual public rsDataPlot, 
  virtual public rsPlotEditor, public ParameterObserver,
  public ChangeBroadcaster, public AudioModuleDeletionWatcher
{

  /** Enumeration of the handles that can be grabbed and dragged by the mouse.  */
  enum dragHandles
  {
    NONE = 0,
    TIME_AND_GAIN
  };

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  EchoLabPlotEditor(CriticalSection *newPlugInLock, EchoLabAudioModule* newEchoLabModuleToEdit);

  /** Destructor. */
  virtual ~EchoLabPlotEditor(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** This function should be used to pass a pointer to an EchoLabDelayLineModuleEditor object. 
  This pointer will be dereferenced in order to make sure that the editor always shows the 
  delayline that is currently selected here. */
  virtual void setDelayLineModuleEditor(EchoLabDelayLineModuleEditor *delayLineEditorToUse);

  /** Marks one of the delaylines as selected and returns true when selection was successful. */
  virtual bool selectDelayLine(int indexToSelect);

  /** De-selects the currently selected delayline (if any). */
  virtual void deSelectDelayLine();

  /** Removes one of the delaylines and returns true when removal was successful. */
  virtual bool removeDelayLine(int indexToRemove);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the index of the currently selected delayline (or -1, if none is selected). */
  //virtual int getSelectedIndex() { return selectedIndex; }

  /** Returns the index of the delayline the is represented by a dot at the given pixel position. 
  It will return -1 if there isn't any delayline-handle at the position in question. */
  virtual int getIndexAtPixelPosition(int x, int y);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual callback function inherited from AudioModuleDeletionWatcher in 
  order to invalidate our pointer-member "selectedDelayLine" when the delayline will be deleted. */
  virtual void audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted);

  /** Implements the purely virtual callback function inherited from ParameterObserver in order to 
  update the plot when one of our observed parameters has chenged. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  /** Implements the purely virtual callback function inherited from ParameterObserver - but has 
  nothing to do. */
  virtual void parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted) { }

  /** Overrides mouseMove in order to update the cursor according to what is under the mouse. */
  virtual void mouseMove(const MouseEvent &e);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop up 
  when the right button is clicked for MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides mouseUp to reset the currentDragHandle to NONE. */
  virtual void mouseUp(const MouseEvent& e);

  /** Overrides mouseWheelMove to adjust the bandwidth on wheel moves. */
  virtual void mouseWheelMove(const MouseEvent& e, const MouseWheelDetails& wheel);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the frequency response plot. */
  virtual void updatePlot();

  /** Refreshes the state internal state delayDesignerForPlot to bring it in sync with the state of 
  echoLabToEdit and updates the plot accordingly. This should be called called whenever the 
  echoLabToEdit has been modified  without accessing it through this object, for example, on 
  preset-change. */
  virtual void refreshCompletely();

protected:

  /** Returns the handle for mouse grab/drag under the specified position (in pixels) as one of the 
  values in enum dragHandles. */
  virtual int getDragHandleAt(int x, int y);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handles. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = nullptr, 
    XmlElement *targetSVG = nullptr);

  /** De-registers this object as ParameterObserver from the observed Parameters. */
  virtual void deRegisterFromObservedParameters();

  CriticalSection       *plugInLock;              // mutex to access the edited AudioModule object 
  EchoLabAudioModule    *echoLabModuleToEdit;

  EchoLabDelayLineModuleEditor *delayLineEditor;

  int selectedIndex;                              // index of currently selected delayline (-1 if none)
  EchoLabDelayLineAudioModule *selectedDelayLine; // pointer to currently selected delayline (NULL if none)

  int currentlyDraggedHandle;  // kind of the handle that is currently being dragged

                               // magnitude response display stuff:
  int    numSamplesInPlot;
  double *timeAxis, *impulseResponse;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class EchoLabModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  EchoLabModuleEditor(CriticalSection *newPlugInLock, EchoLabAudioModule* newEchoLabAudioModule);

  //-----------------------------------------------------------------------------------------------
  // setup:

  virtual void initializeColourScheme();
  virtual void updateSubEditorColourSchemes();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom);
  virtual void rButtonClicked(RButton  *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox  *rComboBoxThatHasChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  EchoLabAudioModule *echoLabModuleToEdit;

  EchoLabPlotEditor         *delayPlotEditor;
  rsPlotZoomer *delayPlotZoomer;

  RSlider           *dryWetSlider, *wetLevelSlider;
  RButton           *snapToTimeGridButton, *delaySyncButton;
  RTimeGridComboBox *timeGridComboBox;

  EchoLabDelayLineModuleEditor *delayLineModuleEditor;

  /*
  //Rectangle leftRectangle1, leftRectangle2, leftRectangle3, eqBandParamRectangle;
  Rectangle middleRectangle, feedbackEqBandParamRectangle, outputEqBandParamRectangle;


  RLabel *delayModLabel, *ampModLabel, *globalLabel;

  RSlider *dryWetSlider, *wetLevelSlider, *timeSlider, *gainSlider, *panSlider, *feedbackSlider, 
  *delayModCycleSlider, *delayModDepthSlider, *delayModPhaseLeftSlider, 
  *delayModPhaseRightSlider, *ampModCycleSlider, *ampModDepthSlider, *ampModPhaseLeftSlider, 
  *ampModPhaseRightSlider;

  RButton *showFeedbackButton, *snapToTimeGridButton, *delaySyncButton, *delayModSyncButton, 
  *ampModSyncButton, *delayModTriggerButton, *ampModTriggerButton, *muteButton, *soloButton, 
  *pingPongButton, *flushButton;
  */

  juce_UseDebuggingNewOperator;
};

#endif 
