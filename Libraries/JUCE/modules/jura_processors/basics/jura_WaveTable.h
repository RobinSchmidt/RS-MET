#ifndef jura_WaveTable_h
#define jura_WaveTable_h

//===============================================================================================
// class StandardWaveformRendererAudioModule:

class StandardWaveformRendererAudioModule : public AudioModule
{
public:
  StandardWaveformRendererAudioModule(CriticalSection *newPlugInLock,
    rosic::StandardWaveformRenderer *newStandardWaveformRendererToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }
  juce_UseDebuggingNewOperator;
protected:
  virtual void initializeAutomatableParameters();
  rosic::StandardWaveformRenderer *wrappedStandardWaveformRenderer;
};

//===============================================================================================
// class WaveformBufferAudioModule:

class WaveformBufferAudioModule : public AudioModule, public AudioFileManager
{
public:
  WaveformBufferAudioModule(CriticalSection *newPlugInLock, 
    rosic::WaveformBuffer *newWaveformBufferToWrap);
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }
  virtual bool loadFile(const juce::File& fileToLoad);
  virtual bool saveToFile(const juce::File& fileToSaveTo);
  virtual bool setAudioData(AudioSampleBuffer* newBuffer, const juce::File& underlyingFile,
    bool markAsClean);
  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);
  juce_UseDebuggingNewOperator;
protected:
  rosic::WaveformBuffer *wrappedWaveformBuffer;
};


//===============================================================================================
// class WaveformRendererAudioModule:

class WaveformRendererAudioModule : public AudioModule
{
  friend class WaveformRendererEditor;
public:
  WaveformRendererAudioModule(CriticalSection *newPlugInLock, 
    rosic::WaveformRenderer *newWaveformRendererToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }
  juce_UseDebuggingNewOperator;
protected:
  virtual void initializeAutomatableParameters();
  rosic::WaveformRenderer *wrappedWaveformRenderer;

  StandardWaveformRendererAudioModule *standardRendererModule;
  WaveformBufferAudioModule           *waveformBufferModule;

};

//===============================================================================================
// class WaveTableAudioModule:

/**

This class wraps rosic::WaveTable into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin.

*/

class WaveTableAudioModule : public AudioModule //, public ChangeListener
{

  friend class WaveTableModuleEditorPopUp;
  friend class WaveTableModuleEditorCompact;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WaveTableAudioModule(CriticalSection *newPlugInLock, rosic::WaveTable *newWaveTableToWrap);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  /** We override this function inherited from StateManager to trigger an update of the waveform
  buffers whenever some parameter changes that should affect the wavefrom buffer (this does not
  happen automatically because changes in the WaveformRenderer do not go through its embedding
  WaveTable object). */
  virtual void markStateAsDirty();

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  // maintain them only temporarily for the transition of the Quadrifex-presets:
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);
  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. This does not really do anything but is needed to make
  the class non-abstract. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }

//=============================================================================================
  juce_UseDebuggingNewOperator;

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  // wrapped rosic object:
  rosic::WaveTable *wrappedWaveTable;

  // child modules:
  WaveformRendererAudioModule *rendererModule;

  //FileManager audioFileManager;  // keep only temporarily

};

#endif 
