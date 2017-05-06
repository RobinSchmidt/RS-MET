#ifndef jura_SamplePlayer_h
#define jura_SamplePlayer_h

/** This class wraps rosic::SamplePlayer into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin. */

class SamplePlayerAudioModule : public AudioModule
{

  friend class SamplePlayerEditorContextMenu;
  friend class BasicSamplePlayerModuleEditor;
  friend class SamplePlayerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SamplePlayerAudioModule(CriticalSection *newPlugInLock, 
    rosic::SamplePlayer *newSamplePlayerToWrap);

  //---------------------------------------------------------------------------------------------
  // automation:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedSamplePlayer->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    jassertfalse; // not yet implemented
  }

  //---------------------------------------------------------------------------------------------
  // others:

  /** Sets up the root-key from the loop-length and number of cycles in loop. */
  virtual void setRootKeyFromLoop();

  /** Loads a new sample file and updates the ranges of the start-/loop-/etc. parameters
  accordingly. */
  virtual bool setSampleFromFile(const juce::File &fileToLoadFrom);


protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::SamplePlayer *wrappedSamplePlayer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class is intended to be used as a graphical editor for a waveform.

-deadlock in the waveform display when destructor is called while resized has not yet returned
-to reproduce it, create an instance inside some GUI, call SamplePlayerEditorDisplay::resized in 
the GUIs resized method (and nothing thereafter), open the GUI in JUCEPluginHost
-perhaps we should implement drawComponent instead of overriding paint */

class SamplePlayerEditorDisplay	: virtual public WaveformDisplay, public ChangeBroadcaster
  //, virtual public InteractiveCoordinateSystem
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SamplePlayerEditorDisplay(AudioFileBuffer *newBufferToUse = NULL);   

  /** Destructor. */
  virtual ~SamplePlayerEditorDisplay(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /*
  virtual CoordinateSystemRangeOld getMaximumMeaningfulRange(double relativeMarginLeft = 10.0, 
  double relativeMarginRight  = 10.0, double relativeMarginTop  = 10.0, 
  double relativeMarginBottom = 10.0);
  */
  /**< Returns an Rectangle wich encloses the curve. Optionally, a margin can
  be specified in percent. */

  /** Sets the rosic::SamplePlayer object which is to be edited. */
  virtual void setSamplePlayerToEdit(rosic::SamplePlayer *newPlayerToEdit);

  /** Passes a new audiofile to be used here. Returns flase when there is some error in opening 
  the file */
  virtual bool setAudioFileToUse(const juce::File &newFileToUse);

  /** Sets the rosic::SampleBuffer object which is to be edited. */
  //virtual void setSampleBufferToEdit(rosic::SampleBuffer *newBufferToEdit);

  /** Sets the Pointer to the rosic::SamplePlaybackParameters which hold the playback info (loop 
  settings etc.) */
  //virtual void setSamplePlaybackParametersToEdit(
  //  rosic::SamplePlaybackParameters *newParametersToEdit);

  /** Sets the actual waveform-data and informs about success. */
  //virtual bool setWaveform(double** newWaveformData, int newNumSampleFrames, int newNumChannels);

  /** Sets the waveform-data from an AudioSampleBuffer and informs about success. */
  //virtual bool setWaveform(const AudioSampleBuffer& newWaveformBuffer);

  /** Causes the internal peakData (inherited from WaveformDisplay) to be recalculated from the
  waveform contained in the rosic::SampleModulator object which is being edited. */
  virtual void updatePlot(bool resetZoomFactor = false);

  /** Overrides the paint-function of the base-classes. */
  virtual void paint(Graphics &g);

  /** Overides the getStateAsXml()-method from the CoordinateSystem base-class. */
  virtual XmlElement* getStateAsXml(
    const juce::String& stateName = juce::String("SamplePlayerEditorDisplayState")) const;

  /** Overwrites the setStateFromXml()-method from the CoordinateSystem base-class. */
  virtual bool setStateFromXml(const XmlElement &xmlState);


protected:

  enum mousableObjects
  {
    NO_OBJECT = 0,
    START_LOCATOR,
    END_LOCATOR,
    LOOP_START_LOCATOR,
    LOOP_END_LOCATOR
  };

  int   locatorBeingDragged; // index of the locator which is being dragged 
                             // (see enum above)

                             // some overrides for mouse-events:
  virtual void mouseDown(const MouseEvent &e);
  virtual void mouseDrag(const MouseEvent &e);
  virtual void mouseMove(const MouseEvent &e);
  virtual void mouseUp  (const MouseEvent &e);

  /** Returns an index of the object which is currently under the mouse-cursor. 
  @see mousableObjects */
  virtual int whatIsUnderTheMouseCursor(const MouseEvent &e);

  rosic::SamplePlayer *samplePlayerToEdit;

  jura::AudioFileBuffer buffer;
  // we have a by-value member-variable here to pass its address to the inherited by-reference
  // member (see constructor implementation)

  //rosic::SampleBuffer *bufferToEdit;
  /**< Pointer to the actual rosic::SampleBuffer object which is being edited. */

  //rosic::SamplePlaybackParameters *playbackParametersToEdit;
  /**< Pointer to the rosic::SamplePlaybackParameters which hold the playback info (loop 
  settings etc.). */

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/**

This is a class for a context menu that can be opened to show and edit a comprehensive set of
parameters for a SamplesPlayer. It is used by SamplePlayerEditor by clicking on the 
'More' button.

*/

class SamplePlayerEditorContextMenu : public ColourSchemeComponent,  public ChangeBroadcaster 
  //public ButtonListener //, public RSliderListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  SamplePlayerEditorContextMenu(SamplePlayerAudioModule* newSamplePlayerModuleToEdit, 
    Component* componentToAttachTo);  
  virtual ~SamplePlayerEditorContextMenu();  

  //---------------------------------------------------------------------------------------------
  // callbacks:

  //virtual void rButtonClicked(RButton *buttonThatWasClicked);
  //virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void resized();
  //virtual void updateWidgetsAccordingToState();

  //---------------------------------------------------------------------------------------------
  // widgets:

  RButton    *closeButton;
  RButton    *loopButton, *muteButton, *soloButton, *phaseRandomizeButton;
  RSlider    *levelSlider, *levelByKeySlider, *levelByVelSlider, *panSlider, *midSideSlider,
    *tuneSlider, *tuneByKeySlider, *tuneByVelSlider, *rootKeySlider, 
    *startSlider, *startByVelSlider, *loopStartSlider, *loopLengthSlider, 
    *lowpassSlider, *highpassSlider,
    *phaseSeedSlider;
  RTextField *ampHeadline, *tuningHeadline, *timeHeadline, *filterHeadline, *miscHeadline;

  //=============================================================================================
  juce_UseDebuggingNewOperator;

protected:

  /** Pointer to the actual SamplePlayerAudioModule object which is being edited. */
  SamplePlayerAudioModule* samplePlayerModuleToEdit;

};



//===============================================================================================

/** This class implements an editor for the SamplePlayerAudio with waveform display.

\todo: 
-move the File-management stuff into the audio module - see LFO
-we must also maintain a state -> don't derive from AudioFileManager, embedd one instead  */

class SamplePlayerModuleEditor : virtual public AudioModuleEditor,  public ChangeBroadcaster, 
  public AudioFileManager, public RSliderListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SamplePlayerModuleEditor(CriticalSection *newPlugInLock, 
    SamplePlayerAudioModule* newSamplePlayerModuleToEdit);

  /** Destructor. */
  virtual ~SamplePlayerModuleEditor();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::SampleBufferPlayer object which is to be edited. 
  Make sure to call this function again with a NULL-pointer when the object get deleted for some 
  reason. */
  //virtual void setSamplePlayerToEdit(rosic::SamplePlayer* newSamplePlayerToEdit);

  /** Sets the juce::Label in which the descriptions for the widgets will appear. */
  //virtual void setDescriptionField(RLabel* newDescriptionField);

  /** Load a new sample. */
  //virtual bool setSampleFromFile(const File &fileToLoadFrom);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Implements the purely virtual changeListenerCallback()-method of the ChangeListener 
  base-class. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Implements the purely virtual rSliderValueChanged()-method of the RSliderListener 
  base-class. */
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);

  /** Overrides paint() */
  virtual void paint(Graphics &g);

  /** Overrides resized() */
  virtual void resized();

  /** Overrides the purely virtual AudioFileManager::setAudioData. */
  virtual bool setAudioData(AudioSampleBuffer* newBuffer, const juce::File& underlyingFile, 
    bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // others:

  /** Updates the sliders, buttons, etc. according to the state of the rosic::StereoOscillator 
  object which is being edited. */
  virtual void updateWidgetsAccordingToState(bool updateSampleDisplayAlso);

  /** Updates the sliders, buttons, etc. according to the state of the rosic::StereoOscillator 
  object which is being edited. */
  virtual void updateWidgetsAccordingToState();

  /** Opens the dialog for loading a new sample. */
  //virtual bool openSampleLoadingDialog();

  //---------------------------------------------------------------------------------------------
  // public data members:

  SamplePlayerEditorDisplay     *sampleDisplay;
  CoordinateSystemZoomer        *sampleDisplayZoomer;
  SamplePlayerEditorContextMenu *contextMenu;

  RTextField *fileLabel, *sampleFileNameLabel, *formatLabel, *formatInfoLabel; 

  RButton *sampleLoadButton, *samplePlusButton, *sampleMinusButton, *loopButton, *muteButton, 
    *soloButton, *phaseRandomizeButton, *moreButton, *fromLoopButton, *loopSnapButton, 
    *loopLengthLockButton, *autoNumCyclesButton;

  RSlider *levelSlider, *tuneSlider, *lowpassSlider, *highpassSlider, *rootKeySlider, 
    *rootDetuneSlider, *startSlider, *startByVelSlider, *endSlider, *loopStartSlider, 
    *loopEndSlider, *loopLengthSlider, *loopCrossfadeTimeSlider, *loopCrossfadeShapeSlider, 
    *loopNumCyclesSlider;

protected:

  //virtual void updateLoopWidgets

  /** Updates the waveform display plot. */   
  //virtual void updatePlot();

  /** Pointer to the actual SampleBufferPlayer object which is being edited. */
  //rosic::SamplePlayer* samplePlayerToEdit;
  SamplePlayerAudioModule* samplePlayerModuleToEdit;

  // some rectangles to define functional groups:
  juce::Rectangle<int> fileRectangle, loopRectangle;

  juce_UseDebuggingNewOperator;
};


#endif 
