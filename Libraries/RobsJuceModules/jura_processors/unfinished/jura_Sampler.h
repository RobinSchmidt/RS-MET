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

  virtual void createParameters();

  // Shorthands for convenience:
  using Engine = rosic::rsSamplerEngine;
  using ReturnCode = rosic::rsSamplerEngine::ReturnCode;
  using Event = rosic::rsMusicalEvent<float>;

  rosic::rsSamplerEngine engine;
  juce::File sfzFile;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};

//=================================================================================================

/** Editor for SamplerAudioModule */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, public jura::FileManager
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  bool loadFile(const juce::File& fileToLoad) override;
  bool saveToFile(const juce::File& fileToSaveTo) override;

  void resized() override;

protected:

  virtual void createWidgets();

  RTextField *instrumentLabel;
  FileSelectionBox *sfzFileLoader;

  RTextField *numLayersLabel, *numLayersField;
  RDraggableNumber *maxNumLayersSlider;

  RTextField *cpuLoadLabel, *cpuLoadField, *ramLoadLabel, *ramLoadField;  
  // todo: diskLoad - maybe ram should also show the total occupation (by all apps) and the 
  // remaining available


  SamplerModule* samplerModule = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};
// maybe it should derive from juce::Timer to periodically update the levelMeter and the 
// numLayersField