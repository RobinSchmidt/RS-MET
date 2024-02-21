#ifndef jura_OscillatorStereo_h
#define jura_OscillatorStereo_h
// rename file to WaveOscillator.h
// 

/** This class wraps rosic::OscillatorStereo into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin.  

todo:
-make parameters modulatable (use rsModulatableSlider/Parameter), where possible (and test it)
-check, how Straightliner can handle additional modulation
-add freq-offset parameter (as target for LFO)
-make a separate gui editor and use the current only for straightliner - remove the "Mod" slider
 from tune (only straightliner uses this
-maybe, if we have just one single WaveOsc, we may put the additional parameters from the context
 menu directly on the main editor - good opportunity to test the multiple editor types feature
 ->in the DualOsc case, we need the smaller gui
-i think, i need to factor out baseclass WaveOscEditorBase that contains what they all have in 
 common
-figure out, why setSate...is called so many times on start up
-make a DualWaveOscillator with interactions...maybe in rosic


-let the osc add its output to what comes in (done - do this for all sources and instruments)   

*/

class WaveOscModule : public AudioModuleWithMidiIn
{

  friend class WaveOscEditor;
  friend class WaveOscEditorContextMenu;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WaveOscModule(CriticalSection *newPlugInLock, 
    rosic::OscillatorStereo *oscToWrap = nullptr);

  virtual ~WaveOscModule();

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  //virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual AudioModuleEditor* createEditor(int type) override;
  // The type parameter is a provision to later enable different types of editors, i.e. more or
  // less compact or detailed ones, with different layouts, etc. It's not yet used - and should
  // probably not be an int but some sort of enum class anyway - Hmm - but different kinds of 
  // AudioModules may have different sets of editor-types - som maybe a regular enum that can 
  // decay to an integer is better here.

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    double tmpL, tmpR;
    for(int n = 0; n < numSamples; n++)
    {
      wrappedOsc->getSampleFrameStereo(&tmpL, &tmpR);
      inOutBuffer[0][n] += tmpL; 
      inOutBuffer[1][n] += tmpR;
    }
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    //wrappedOscillatorStereo->getSampleFrameStereo(left, right);
    double tmpL, tmpR;
    wrappedOsc->getSampleFrameStereo(&tmpL, &tmpR);
    *left  += tmpL;
    *right += tmpR;
  }
  // try to optimize the temporaries away - maybe the osc itself can accumulate its output into
  // the passed pointers

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedOsc->setSampleRate(newSampleRate);
  }

  virtual void reset() override
  {
    wrappedOsc->reset();
  }

  virtual void noteOn(int noteNumber, int velocity) override
  {
    wrappedOsc->setKeyAndVel(noteNumber, velocity);
    wrappedOsc->setFrequencyNominal(RAPT::rsPitchToFreq(noteNumber));
    wrappedOsc->reset();
     // preliminary - use tuning table - maybe baseclass AudioModuleWithMidiIn should maintain
     // a pointer to a TuningManager object - and maybe a subclass AudioModulePolyphonic can
     // also maintain a VoiceManager ...soo many managers...meta,smoothing,modulation,... and now
     // these :-)
  }


protected:

  /** Fills the array of automatable parameters. */
  virtual void createParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::OscillatorStereo *wrappedOsc;

  bool wrappedOscIsOwned = false;
  rosic::MipMappedWaveTableStereo *waveTable = nullptr; 
  // only needed when the wrapped osc is owned - in this case, we also need an owned wavetable


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(WaveOscModule)
};

//=================================================================================================

/** This is a class for a context menu that can be opened to show and edit a comprehensive set of
parameters for a stereo-oscillator. It is used by WaveOscEditor by clicking on the 'More' 
button. */

class WaveOscEditorContextMenu : public ColourSchemeComponent,  public ChangeBroadcaster, 
  public RButtonListener, public RSliderListener//, public ComponentMovementWatcher,
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WaveOscEditorContextMenu(WaveOscModule* newOscillatorModuleToEdit, 
    Component* componentToAttachTo);  

  /** Destructor. */
  virtual ~WaveOscEditorContextMenu();  

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

  /** Pointer to the actual WaveOscModule object which is being edited. */
  WaveOscModule* oscillatorModuleToEdit;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(WaveOscEditorContextMenu)
};

//===============================================================================================

/** Editor for the WaveOsc */

class WaveOscEditor : public SampleBasedAudioModuleEditor, public RSliderListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WaveOscEditor(CriticalSection *newPlugInLock, WaveOscModule* newWaveOscModule);  

  /** Destructor. */
  virtual ~WaveOscEditor();  

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

  rsWaveformPlot *waveformDisplay;
  rsPlot         *emptyDisplay;
  WaveOscEditorContextMenu *contextMenu;
  //Viewport *contextMenuViewport;

  RButton      *moreButton;
  RSlider      *levelSlider, *pitchModulationSlider;
  TuningSlider *tuneSlider;
    // remove the pitch-modulation slider - let it remain visible only for straightliner

protected:

  virtual void createWidgets();

  /** Updates the waveform display plot. */
  virtual void updatePlot();

  /** Updates the visibility of the widgets according to the on/off state of the oscillator. */
  virtual void updateWidgetVisibility();

  /** Pointer to the actual StereoOscillator object which is being edited. */
  rosic::OscillatorStereo* oscillatorToEdit;

  // get rid of this - it doesn't belog here:
  int      numSamplesInPlot;
  double*  waveformBuffer;
  double** waveformPointers;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(WaveOscEditor)
};



#endif 
