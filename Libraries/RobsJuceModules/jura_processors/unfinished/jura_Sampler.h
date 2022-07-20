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

/** A subclass of jura::FileManager that deals specifically with .sfz files and has a pointer to
a jura::SamplerModule which it maintains in sync with the loaded file...tbc...  */

class JUCE_API SfzFileManager : public jura::FileManager
{

public:

  SfzFileManager(SamplerModule *moduleToSetup) : samplerModule(moduleToSetup) {}
  bool loadFile(const juce::File& fileToLoad) override;
  bool saveToFile(const juce::File& fileToSaveTo) override;

protected:

  SamplerModule *samplerModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzFileManager)
};


//=================================================================================================

/** Editor for SamplerAudioModule */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, public jura::FileManager, 
  public Timer
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  virtual ~SamplerEditor();

  // Overrides                                   // overriden from...
  bool loadFile(const juce::File& f) override;   //   FileManager
  bool saveToFile(const juce::File& f) override; //   FileManager
  void timerCallback() override;                 //   Timer
  void resized() override;                       //   AudioModuleEditor

protected:

  virtual void createWidgets();


  SamplerModule* samplerModule = nullptr;

  jura::RTextField *instrumentLabel;
  jura::MeteringDisplayWithText *layersMeter;

  // SFZ text editor:
  jura::SfzFileManager *sfzFileManager;
  jura::FileSelectionBox *sfzFileLoader;
  juce::CodeDocument sfzDoc;           // Declare doc before the editor because the editor holds a 
  juce::CodeEditorComponent sfzEditor; // reference to it (-> order of construction/destruction)
  //juce::CodeTokeniser sfzTokenizer;
  // ToDo: implement this, see:
  //   https://docs.juce.com/master/classCodeTokeniser.html
  // The CodeTokeniser class is abstract. We need to make a subclass rsSfzTokenizer. Maybe someone
  // else already did that? Check open-source SFZ sampler projects. When we have that, we need to 
  // pass a pointer to our tokenizer to the constructor of the editor


  /*

  RTextField *numLayersLabel, *numLayersOfLabel, ;
  RDraggableNumber *maxNumLayersSlider;

  RTextField *cpuLoadLabel, *cpuLoadField, *ramLoadLabel, *ramLoadField;  
  // todo: diskLoad - maybe ram should also show the total occupation (by all apps) and the 
  // remaining available
  */









  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};
// maybe it should derive from juce::Timer to periodically update the levelMeter and the 
// numLayersField