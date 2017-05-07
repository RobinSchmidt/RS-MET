#ifndef jura_LowFrequencyOscillator_h
#define jura_LowFrequencyOscillator_h


class LowFrequencyOscillatorAudioModule : public AudioModule
{

  friend class LowFrequencyOscillatorEditor;
  friend class LowFrequencyOscillatorEditorCompact;

public:

  LowFrequencyOscillatorAudioModule(CriticalSection *newPlugInLock, 
    rosic::LowFrequencyOscillator *newLowFrequencyOscillatorToWrap);


  virtual void parameterChanged(Parameter* parameterThatHasChanged);  // to be deprecated

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);

  virtual void setBeatsPerMinute(double newBpm)
  {
    wrappedLowFrequencyOscillator->setBeatsPerMinute(newBpm);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = wrappedLowFrequencyOscillator->getSample();
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::LowFrequencyOscillator *wrappedLowFrequencyOscillator;

  // child modules:
  WaveTableAudioModule *waveTableModule;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================





#endif 
