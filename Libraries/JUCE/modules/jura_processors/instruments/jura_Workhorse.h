#ifndef jura_Workhorse_h
#define jura_Workhorse_h

class WorkhorseAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class WorkhorseModuleEditor;

public:

  WorkhorseAudioModule(CriticalSection *newPlugInLock, rosic::Workhorse *workhorseToWrap);

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedWorkhorse != NULL)
      wrappedWorkhorse->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedWorkhorse->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset()
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



#endif 
