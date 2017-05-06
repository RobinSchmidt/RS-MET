#ifndef jura_KeyShotAudioModule_h
#define jura_KeyShotAudioModule_h

class KeyShotAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class KeyShotModuleEditor;

public:

  KeyShotAudioModule(CriticalSection *newPlugInLock, rosic::KeyShot *keyShotToWrap);

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate)
  {
    wrappedKeyShot->setSampleRate(newSampleRate);
  }

  virtual void reset()
  {
    wrappedKeyShot->resetAllVoices();
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedKeyShot->getSampleFrameStereo(inOutL, inOutR);
  }

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule        *samplePlayerModule;
  BreakpointModulatorAudioModule *ampEnvModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::KeyShot *wrappedKeyShot;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class KeyShotModuleEditor : public PolyphonicInstrumentEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  KeyShotModuleEditor(CriticalSection *newPlugInLock, KeyShotAudioModule* newKeyShotAudioModule);

  /** Destructor. */
  //virtual ~KeyShotModuleEditor();

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
  KeyShotAudioModule *keyShotAudioModule;

  // the sub-editors:
  //SamplePlayerModuleEditor  *samplePlayerEditor;
  //DualWaveformDisplay       *waveformDisplay;
  WaveformDisplayOld        *waveformDisplay;
  //BreakpointModulatorEditor *ampEnvEditor;

  // the tabber to choose which of the embedded module editors is shown:
  //TabbedButtonBar* moduleEditorTabButtonBar;

  // the file-browser stuff:
  WildcardFileFilter*   fileFilter;
  //FilePreviewComponent* filePreviewer;
  //FileBrowserComponent* browser;

  juce_UseDebuggingNewOperator;
};

#endif 
