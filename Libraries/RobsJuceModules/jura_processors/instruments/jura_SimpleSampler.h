#ifndef jura_SimpleSampler_h
#define jura_SimpleSampler_h

/** This class wraps rosic::SimpleSampler into a rosof::AudioModule to facilitate its use as 
plugIn.  */

class SimpleSamplerAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class SimpleSamplerModuleEditor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SimpleSamplerAudioModule(CriticalSection *newPlugInLock);

  virtual ~SimpleSamplerAudioModule();

  AudioModuleEditor* createEditor(int type) override;


  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedSimpleSampler->setSampleRate(newSampleRate);
  }

  virtual void reset() override
  {
    wrappedSimpleSampler->resetAllVoices();
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedSimpleSampler->getSampleFrameStereo(inOutL, inOutR);
  }


protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule        *samplePlayerModule;
  MultiModeFilterAudioModule     *filterModule;
  BreakpointModulatorAudioModule *pitchEnvModule, *filterEnvModule, *ampEnvModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::SimpleSampler *wrappedSimpleSampler;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class SimpleSamplerModuleEditor : public PolyphonicInstrumentEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SimpleSamplerModuleEditor(CriticalSection *newPlugInLock, 
    SimpleSamplerAudioModule* newSimpleSamplerAudioModule);

  /** Destructor. */
  //virtual ~SimpleSamplerModuleEditor();

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  //virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Overrides changeListenerCallback() in order to start or stop the timers in the embedded
  module editors when the user selects another tab.*/
  virtual void changeListenerCallback(ChangeBroadcaster* objectThatHasChanged);

  //---------------------------------------------------------------------------------------------
  // others:

  /** Overrides the resized-method. */
  virtual void resized();

protected:

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  SimpleSamplerAudioModule *simpleSamplerAudioModule;

  // the sub-editors:
  SamplePlayerModuleEditor    *samplePlayerEditor;
  MultiModeFilterModuleEditor *filterEditor;
  BreakpointModulatorEditor   *pitchEnvEditor, *filterEnvEditor, *ampEnvEditor;

  // the tabber to choose which of the embedded module editors is shown:
  TabbedButtonBar* moduleEditorTabButtonBar;


  // the file-browser stuff:
  WildcardFileFilter*   fileFilter;
  //FilePreviewComponent* filePreviewer;
  //FileBrowserComponent* browser;

  juce_UseDebuggingNewOperator;
};

#endif 
