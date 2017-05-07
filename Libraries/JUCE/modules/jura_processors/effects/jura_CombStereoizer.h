#ifndef jura_CombStereoizer_h
#define jura_CombStereoizer_h

class CombStereoizerAudioModule : public AudioModule
{

  friend class CombStereoizerModuleEditor;

public:

  CombStereoizerAudioModule(CriticalSection *newPlugInLock, rosic::CombStereoizer *stereoizerToWrap);

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedCombStereoizer->setSampleRate((float)newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedCombStereoizer->getSampleFrameStereo(inOutL, inOutR);
  }

protected:

  void initializeAutomatableParameters();

  rosic::CombStereoizer *wrappedCombStereoizer;
  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif 
