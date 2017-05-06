#ifndef jura_MultiModeFilter_h
#define jura_MultiModeFilter_h

class MultiModeFilterAudioModule : public AudioModule
{

  friend class MultiModeFilterModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultiModeFilterAudioModule(CriticalSection *newPlugInLock, 
    rosic::MultiModeFilter *newMultiModeFilterToWrap);

  //---------------------------------------------------------------------------------------------
  // overrides:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedMultiModeFilter->getSampleFrameStereo(inOutL, inOutR);
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::MultiModeFilter *wrappedMultiModeFilter;


  juce_UseDebuggingNewOperator;
};

#endif 
