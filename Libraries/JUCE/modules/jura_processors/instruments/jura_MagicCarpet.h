#ifndef jura_MagicCarpet_h
#define jura_MagicCarpet_h

class DelayPhaserAudioModule : public AudioModule
{
  friend class DelayPhaserModuleEditor;
public:
  DelayPhaserAudioModule(CriticalSection *newPlugInLock, rosic::DelayPhaser *newDelayPhaserToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedDelayPhaser->getSampleFrameStereo(inOutL, inOutR);
  }
  juce_UseDebuggingNewOperator;
protected:
  virtual void initializeAutomatableParameters();
  rosic::DelayPhaser      *wrappedDelayPhaser;
  PhaserAudioModule       *phaser1Module, *phaser2Module;
  PingPongEchoAudioModule *delayModule;
};

//=================================================================================================

class MagicCarpetAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class MagicCarpetModuleEditor;

public:


  MagicCarpetAudioModule(CriticalSection *newPlugInLock);

  virtual ~MagicCarpetAudioModule();

  AudioModuleEditor* createEditor(int type) override;



  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedMagicCarpet->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedMagicCarpet->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset() override
  {
    wrappedMagicCarpet->resetAllVoices();
  }

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  VectorSamplePlayerAudioModule  *oscSectionModule;
  BreakpointModulatorAudioModule *filterEnvModule, *ampEnvModule;
  FourPoleFilterAudioModule      *filterModule;
  EqualizerAudioModule           *equalizerModule;
  DelayPhaserAudioModule         *delayPhaserModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::MagicCarpet *wrappedMagicCarpet;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Subclass of FourPoleFilterModuleEditor to customize the arrangement of widgets in resized. */

class MagicCarpetFilterEditor : public FourPoleFilterModuleEditor
{
public:
  MagicCarpetFilterEditor(CriticalSection *newPlugInLock, 
    FourPoleFilterAudioModule* newFourPoleFilterAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  RSlider *frequencyByKeySlider, *frequencyByVelSlider, *gainByKeySlider, *gainByVelSlider;
};

/** Subclass of PhaserModuleEditor to customize the arrangement of widgets in 
resized. */
class PhaserModuleEditorCompact : public PhaserModuleEditor
{
public:
  PhaserModuleEditorCompact(CriticalSection *newPlugInLock, PhaserAudioModule* newPhaserAudioModule);
  virtual void resized();
  RButton *onOffButton, *secondOrderButton;
  juce_UseDebuggingNewOperator;
};

/** Subclass of PingPongEchoModuleEditor to customize the arrangement of widgets in 
resized. */
class PingPongEchoModuleEditorCompact : public PingPongEchoModuleEditor
{
public:
  PingPongEchoModuleEditorCompact(CriticalSection *newPlugInLock, 
    PingPongEchoAudioModule* newPingPongEchoAudioModule);
  virtual void resized();
  RButton *onOffButton;
  juce_UseDebuggingNewOperator;
};

/** Editor for the DelayPhaser effect. */
class DelayPhaserModuleEditor : public AudioModuleEditor
{
public:
  DelayPhaserModuleEditor(CriticalSection *newPlugInLock, 
    DelayPhaserAudioModule* newDelayPhaserAudioModule);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  RTextField *feedbackLabel;
  RSlider    *dryWetSlider, *feedback1Slider, *feedback2Slider, *feedback3Slider;
  PhaserModuleEditorCompact       *phaser1Editor, *phaser2Editor;
  PingPongEchoModuleEditorCompact *delayEditor;
};

//=================================================================================================

/** Editor for the Magic Carpet synthesizer. */

class MagicCarpetModuleEditor : public PolyphonicInstrumentEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MagicCarpetModuleEditor(CriticalSection *newPlugInLock, 
    MagicCarpetAudioModule* newMagicCarpetAudioModule);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Overrides changeListenerCallback() in order to start or stop the timers in the embedded
  module editors when the user selects another tab.*/
  virtual void changeListenerCallback(ChangeBroadcaster* objectThatHasChanged);

  virtual void resized();

protected:

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  MagicCarpetAudioModule *magicCarpetAudioModule;

  // the sub-editors:
  VectorSamplePlayerEditor      *oscSectionEditor;
  //EnvelopeGeneratorModuleEditor *filterEnvEditor, *ampEnvEditor;
  BreakpointModulatorEditorCompact *filterEnvEditor, *ampEnvEditor;
  MagicCarpetFilterEditor       *filterEditor;
  EqualizerModuleEditor         *equalizerEditor;
  DelayPhaserModuleEditor       *delayPhaserEditor;

  juce::Rectangle<int> performanceRect, effectRect;

  juce_UseDebuggingNewOperator;
};

#endif 
