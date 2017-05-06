#ifndef jura_SamplePlayer_h
#define jura_SamplePlayer_h

/** This class wraps rosic::SamplePlayer into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin. */

class SamplePlayerAudioModule : public AudioModule
{

  friend class SamplePlayerEditorContextMenu;
  friend class BasicSamplePlayerModuleEditor;
  friend class SamplePlayerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SamplePlayerAudioModule(CriticalSection *newPlugInLock, 
    rosic::SamplePlayer *newSamplePlayerToWrap);

  //---------------------------------------------------------------------------------------------
  // automation:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedSamplePlayer->getSampleFrameStereo(inOutL, inOutR);
  }

  //---------------------------------------------------------------------------------------------
  // others:

  /** Sets up the root-key from the loop-length and number of cycles in loop. */
  virtual void setRootKeyFromLoop();

  /** Loads a new sample file and updates the ranges of the start-/loop-/etc. parameters
  accordingly. */
  virtual bool setSampleFromFile(const juce::File &fileToLoadFrom);


protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::SamplePlayer *wrappedSamplePlayer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif 
