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


  MagicCarpetAudioModule(CriticalSection *newPlugInLock, rosic::MagicCarpet *magicCarpetToWrap);

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate)
  {
    wrappedMagicCarpet->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedMagicCarpet->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset()
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


#endif 
