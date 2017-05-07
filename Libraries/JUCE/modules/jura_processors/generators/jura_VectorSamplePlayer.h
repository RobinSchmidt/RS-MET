#ifndef jura_VectorSamplePlayer_h
#define jura_VectorSamplePlayer_h


class VectorSamplePlayerAudioModule : public AudioModule
{

  friend class VectorSamplePlayerEditor;

public:


  VectorSamplePlayerAudioModule(CriticalSection *newPlugInLock, 
    rosic::VectorSamplePlayer *vectorSamplePlayerToWrap);

  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedVectorSamplePlayer != NULL)
      wrappedVectorSamplePlayer->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedVectorSamplePlayer->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset()
  {
    if(wrappedVectorSamplePlayer != NULL)
      wrappedVectorSamplePlayer->reset();
  }

protected:

  // we maintain wrappped versions (into jura::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule *samplePlayerTopLeftModule, *samplePlayerTopRightModule,
    *samplePlayerBottomLeftModule, *samplePlayerBottomRightModule;
  VectorMixerAudioModule *vectorMixerModule;
  LowFrequencyOscillatorAudioModule *xLfoModule, *yLfoModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::VectorSamplePlayer *wrappedVectorSamplePlayer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================



#endif 
