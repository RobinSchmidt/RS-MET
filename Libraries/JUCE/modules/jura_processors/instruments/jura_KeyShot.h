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


#endif 
