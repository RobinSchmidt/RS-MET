#ifndef jura_Moduluxury_h
#define jura_Moduluxury_h

/** This is supposed to be come a luxurious modulation plugin  - with freely adjustable modulations
for filters, phasers, flangers, panning, volume, etc. - maybe with routable modulation 
generators... */

class ModuluxuryAudioModule : public AudioModule
{

  friend class ModuluxuryModuleEditor;

public:

  ModuluxuryAudioModule(CriticalSection *newPlugInLock, 
    rosic::Moduluxury *moduluxuryToWrap = nullptr);

  virtual ~ModuluxuryAudioModule();

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedModuluxury->setSampleRate((float)newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    //  wrappedModuluxury->getSampleFrameStereo(inOutL, inOutR); 
  }

  virtual void reset()
  {
    //  wrappedModuluxury->reset(); 
  }


protected:

  void initializeAutomatableParameters();

  rosic::Moduluxury *wrappedModuluxury;
  bool wrappedModuluxuryIsOwned = false;

  juce_UseDebuggingNewOperator;
};

#endif 
