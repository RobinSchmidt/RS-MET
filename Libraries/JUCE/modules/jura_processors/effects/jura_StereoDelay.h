#ifndef jura_StereoDelay_h
#define jura_StereoDelay_h

class StereoDelayAudioModule : public AudioModule
{

  friend class StereoDelayModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  StereoDelayAudioModule(CriticalSection *newPlugInLock, rosic::StereoDelay *stereoDelayToWrap);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedStereoDelay->setSampleRate(newSampleRate);
  }

  virtual void setBeatsPerMinute(double newBpm)
  {
    wrappedStereoDelay->setBeatsPerMinute(newBpm);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    wrappedStereoDelay->getSampleFrameStereo(inOutL, inOutL);
  }

//---------------------------------------------------------------------------------------------
// event processing:

  virtual void reset()
  {
    wrappedStereoDelay->reset();
  }


protected:

  void initializeAutomatableParameters();

  rosic::StereoDelay *wrappedStereoDelay;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================




#endif 
