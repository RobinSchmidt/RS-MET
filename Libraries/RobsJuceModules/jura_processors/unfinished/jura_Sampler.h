#pragma once

/** A sampler with functionality roughly based on the sfz specification */

class JUCE_API SamplerModule : public jura::AudioModuleWithMidiIn
{

public:

  SamplerModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr);



  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor(int type) override;

  void setSampleRate(double newSampleRate) override;
  void setGain(double newGain);


  int getNumActiveLayers() const { return engine.getNumActiveLayers(); }


  void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;


  // Midi Handling:
  void noteOn(int key, int vel) override;
  void noteOff(int key) override;
  //void handleMidiMessage(MidiMessage message) override;

  // Audio Processing:
  void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  void processStereoFrame(double *left, double *right) override;

  void reset() override;


protected:

  /** Creates the parameters of the sampler that sit directly in the xml file. They have nothing to
  do with the opcodes defined in the sfz files and they also do not alter the sfz settings that the
  engine works with. Anything that is controlled by these parameters is either post-processing step
  of the engine's output or controls a global behavior that is not specified in the sfz 
  specification such as the resampling quality, the selection whether opcodes should work 
  accumulatively or overridingly, polyphony, a global gain, etc.. */
  virtual void createParameters();

  /** Sets up the member variables that define where the app expects sfz-files, samples, etc. 
  Called in constructor. */
  virtual void setupDirectories();

  /** Checks, whether or not an .sfz file with the given path exists. The path is supposed to be 
  relative with respect to our "sfzRootDir" member. */
  virtual bool doesSfzFileExist(const juce::String& path);

  void setBusMode(bool shouldAccumulate);

  // Shorthands for convenience:
  //using Engine = rosic::rsSamplerEngine;  // old
  using Engine     = rosic::Sampler::rsSamplerEngine2;   // new
  using ReturnCode = rosic::Sampler::rsReturnCode;
  using Event      = rosic::Sampler::rsMusicalEvent<float>;

  Engine engine;
  //juce::File sfzFile;

  // under construction:
  juce::String sfzRootDir;
  //juce::String sampleRootDir;




  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};

//=================================================================================================

/** Editor for SamplerAudioModule */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, public jura::FileManager, 
  public Timer
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  // Overrides                                   // overriden from...
  bool loadFile(const juce::File& f) override;   //   FileManager
  bool saveToFile(const juce::File& f) override; //   FileManager
  void timerCallback() override;                 //   Timer
  void resized() override;                       //   AudioModuleEditor

protected:

  virtual void createWidgets();


  jura::RTextField *instrumentLabel;
  //jura::RTextField *numLayersLabel, *numLayersField; // numLayersField may become obsolete
  jura::MeteringDisplay *layersMeter;

  /*
  FileSelectionBox *sfzFileLoader;
  RTextField *numLayersLabel, *numLayersOfLabel, ;
  RDraggableNumber *maxNumLayersSlider;

  RTextField *cpuLoadLabel, *cpuLoadField, *ramLoadLabel, *ramLoadField;  
  // todo: diskLoad - maybe ram should also show the total occupation (by all apps) and the 
  // remaining available
  */







  SamplerModule* samplerModule = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};
// maybe it should derive from juce::Timer to periodically update the levelMeter and the 
// numLayersField