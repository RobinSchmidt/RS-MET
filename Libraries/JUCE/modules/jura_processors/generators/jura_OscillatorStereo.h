#ifndef jura_OscillatorStereo_h
#define jura_OscillatorStereo_h

/** This class wraps rosic::OscillatorStereo into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin.  */

class OscillatorStereoAudioModule : public AudioModule
{

  friend class OscillatorStereoEditor;
  friend class OscillatorStereoEditorContextMenu;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  OscillatorStereoAudioModule(CriticalSection *newPlugInLock, 
    rosic::OscillatorStereo *newOscillatorStereoToWrap);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  //virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedOscillatorStereo->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
      wrappedOscillatorStereo->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void createParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::OscillatorStereo *wrappedOscillatorStereo;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This is a class for a context menu that can be opened to show and edit a comprehensive set of
parameters for a stereo-oscillator. It is used by OscillatorStereoEditor by clicking on the 'More' 
button. */

class OscillatorStereoEditorContextMenu : public ColourSchemeComponent,  public ChangeBroadcaster, 
  public RButtonListener, public RSliderListener//, public ComponentMovementWatcher,
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  OscillatorStereoEditorContextMenu(OscillatorStereoAudioModule* newOscillatorModuleToEdit, 
    Component* componentToAttachTo);  

  /** Destructor. */
  virtual ~OscillatorStereoEditorContextMenu();  

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Implements the purely virtual rSliderValueChanged()-method of the RSliderListener 
  base-class. */
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);

  /** Overrides the resized-method of the RobsEditorBase base-class. */
  virtual void resized();

  /** Updates the sliders, buttons, etc. according to the state of the rosic::StereoOscillator 
  object which is being edited. */
  //virtual void updateWidgetsAccordingToState();
  // may be obsolete

  /** Implements the purely virtual method of the ComponentMovementWatcher baseclass. */
  //virtual void componentMovedOrResized(bool wasMoved, bool wasResized);

  /** Implements the purely virtual method of the ComponentMovementWatcher baseclass. */
  //virtual void componentPeerChanged();

  /** Overrides paint. */
  //virtual void paint(Graphics &g);

  RTextField *ampHeadline, *tuningHeadline, *timeHeadline, *magSpectrumHeadline, 
    *phaseSpectrumHeadline;

  RSlider *levelSlider, *levelByKeySlider, *levelByVelSlider, *panSlider, *midSideSlider, 
    *startPhaseSlider, *fullWavePhaseWarpSlider, *halfWavePhaseWarpSlider, 
    *combHarmonicSlider, *combAmountSlider, *detuneHzSlider, 
    *stereoDetuneSlider, *stereoDetuneHzSlider, *spectralContrastSlider, *spectralSlopeSlider, 
    *highestHarmonicSlider, *lowestHarmonicSlider, *evenOddSlider, *evenOddPhaseShiftSlider, 
    *phaseScaleSlider, *phaseShiftSlider, *stereoPhaseShiftSlider, 
    *evenOddStereoPhaseShiftSlider;

  TuningSlider *tuneSlider;

  RButton *invertButton, *reverseButton, *closeButton;

  // other ideas: stereo-shift, invert left and right separately

protected:

  virtual void createWidgets();

  /** Pointer to the actual OscillatorStereoAudioModule object which is being edited. */
  OscillatorStereoAudioModule* oscillatorModuleToEdit;

  juce_UseDebuggingNewOperator;
};

//===============================================================================================

/** Editor for the StereoOscillator */

class OscillatorStereoEditor : virtual public AudioModuleEditor, public AudioFileManager,
  public RSliderListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  OscillatorStereoEditor(CriticalSection *newPlugInLock, 
    OscillatorStereoAudioModule* newOscillatorStereoAudioModule);  

  /** Destructor. */
  virtual ~OscillatorStereoEditor();  

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Changes the text of the oscillators on/off button. */
  //virtual void setOnOffButtonText(const juce::String &newText);

  /** Sets the waveform from an audiofile. */
  //virtual bool setWaveformFromFile(const File &fileToLoadFrom);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void mouseDown(const MouseEvent &e);
  virtual void resized();

  /** Overrides the purely virtual AudioFileManager::setAudioData. */
  virtual bool setAudioData(AudioSampleBuffer* newBuffer, const juce::File& underlyingFile, 
    bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // others:

  /** Updates the sliders, buttons, etc. according to the state of the rosic::StereoOscillator 
  object which is being edited. */
  virtual void updateWidgetsAccordingToState();


  //virtual bool openWaveformLoadingDialog();

  //---------------------------------------------------------------------------------------------
  // public data members:

  WaveformDisplayOld                *waveformDisplay;
  CoordinateSystemOld               *emptyDisplay;
  OscillatorStereoEditorContextMenu *contextMenu;
  //Viewport                          *contextMenuViewport;

  RTextField   *waveFileLabel;
  RButton      *waveLoadButton, *wavePlusButton, *waveMinusButton, *moreButton;
  RSlider      *levelSlider, *pitchModulationSlider;
  TuningSlider *tuneSlider;

protected:

  /** Updates the waveform display plot. */
  virtual void updatePlot();

  /** Updates the visibility of the widgets according to the on/off state of the oscillator. */
  virtual void updateWidgetVisibility();

  /** Pointer to the actual StereoOscillator object which is being edited. */
  rosic::OscillatorStereo* oscillatorToEdit;

  int      numSamplesInPlot;
  double*  waveformBuffer;
  double** waveformPointers;

  juce_UseDebuggingNewOperator;
};



#endif 
