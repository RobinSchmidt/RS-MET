#ifndef jura_VectorMixer_h
#define jura_VectorMixer_h

class VectorMixerAudioModule : public AudioModule
{

  friend class VectorMixerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  VectorMixerAudioModule(CriticalSection *newPlugInLock, 
    rosic::VectorMixer *newVectorMixerToWrap);

  //---------------------------------------------------------------------------------------------
  // automation:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  /*
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);
  */

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }


protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::VectorMixer *wrappedVectorMixer;

  juce_UseDebuggingNewOperator;
};

#endif 
