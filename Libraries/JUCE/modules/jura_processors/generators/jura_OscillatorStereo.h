#ifndef jura_OscillatorStereo_h
#define jura_OscillatorStereo_h

/** This class wraps rosic::OscillatorStereo into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin.  */

class OscillatorStereoAudioModule : public AudioModule
{

  friend class OscillatorStereoEditor;
  friend class OscillatorStereoEditorContextMenu;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  OscillatorStereoAudioModule(CriticalSection *newPlugInLock, 
    rosic::OscillatorStereo *newOscillatorStereoToWrap);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedOscillatorStereo->getSampleFrameStereo(inOutL, inOutR);
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::OscillatorStereo *wrappedOscillatorStereo;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

#endif 
