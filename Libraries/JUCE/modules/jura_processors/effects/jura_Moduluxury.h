#ifndef jura_Moduluxury_h
#define jura_Moduluxury_h


class ModuluxuryAudioModule : public AudioModule
{

  friend class ModuluxuryModuleEditor;

public:

  ModuluxuryAudioModule(CriticalSection *newPlugInLock, rosic::Moduluxury *moduluxuryToWrap);

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

  juce_UseDebuggingNewOperator;
};

//=================================================================================================



#endif 
