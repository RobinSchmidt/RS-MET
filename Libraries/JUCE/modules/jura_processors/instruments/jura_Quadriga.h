#ifndef jura_Quadriga_h
#define jura_Quadriga_h

class QuadrigaAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class QuadrigaModuleEditor;

public:

  QuadrigaAudioModule(CriticalSection *newPlugInLock, rosic::Quadriga *quadrigaToWrap);

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset()
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->resetAllVoices();
  }

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  QuadrigenAudioModule *quadrigenModule;
  QuadrifexAudioModule *quadrifexModule;
  EqualizerAudioModule *equalizerModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::Quadriga *wrappedQuadriga;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif 
