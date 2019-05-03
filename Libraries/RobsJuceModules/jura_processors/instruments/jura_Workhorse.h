#ifndef jura_Workhorse_h
#define jura_Workhorse_h

class WorkhorseAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class WorkhorseModuleEditor;

public:

  WorkhorseAudioModule(CriticalSection *newPlugInLock);

  virtual ~WorkhorseAudioModule();

  AudioModuleEditor* createEditor(int type) override;



  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate) override
  {
    if(wrappedWorkhorse != NULL)
      wrappedWorkhorse->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedWorkhorse->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset() override
  {
    if(wrappedWorkhorse != NULL)
      wrappedWorkhorse->resetAllVoices();
  }

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule *samplePlayerTopLeftModule, *samplePlayerTopRightModule,
    *samplePlayerBottomLeftModule, *samplePlayerBottomRightModule;
  MultiModeFilterAudioModule     *filterModule;
  BreakpointModulatorAudioModule *pitchEnvModule, *filterEnvModule, *ampEnvModule;
  VectorMixerAudioModule *vectorMixerModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::Workhorse *wrappedWorkhorse;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class WorkhorseModuleEditor : public PolyphonicInstrumentEditor
{

public:

  WorkhorseModuleEditor(CriticalSection *newPlugInLock, 
    WorkhorseAudioModule* newWorkhorseAudioModule);

  // callbacks:
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();
  virtual void paint(Graphics &g);

  /** Overrides changeListenerCallback() in order to start or stop the timers in the embedded
  module editors when the user selects another tab. ...??? */
  virtual void changeListenerCallback(ChangeBroadcaster* objectThatHasChanged);


protected:

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** Makes all the sub-editors for the various modules invisible. */
  virtual void makeSubEditorsInvisible();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  WorkhorseAudioModule *workhorseAudioModule;

  // the sub-editors:
  //CoordinateSystem* vectorMixerXYPad;
  VectorMixerModuleEditor* vectorMixerPad;
  SamplePlayerModuleEditor *samplePlayerTopLeftEditor, *samplePlayerTopRightEditor, 
    *samplePlayerBottomLeftEditor, *samplePlayerBottomRightEditor;
  MultiModeFilterModuleEditor *filterEditor;
  BreakpointModulatorEditor *pitchEnvEditor, *filterEnvEditor, *ampEnvEditor;
  // LowFreqOscEditors, EqualizerEditor, PoleZeroModelEditor, DelayEditor
  // the XY pad

  // ToDo: replace the buttons with mini-preview components (showing only a small plot, 
  // preset-management stuff and maybe an edit button (although we could also use the whole 
  // components like buttons)
  PlotPreviewButton *samplePlayerTopLeftButton, *samplePlayerTopRightButton, 
    *samplePlayerBottomLeftButton, *samplePlayerBottomRightButton, *lowFreqOscXButton, 
    *lowFreqOscYButton, *filterButton, *pitchEnvButton, *filterEnvButton, *ampEnvButton, 
    *masterFilterButton, *masterFilterEnvButton, *masterAmpEnvButton, *equalizerButton, 
    *delayButton, *reverbButton;

  juce_UseDebuggingNewOperator;
};


#endif 
