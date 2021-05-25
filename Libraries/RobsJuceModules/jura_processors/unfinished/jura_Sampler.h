#pragma once

/** A sampler with functionality roughly based on the sfz specification */

class JUCE_API SamplerModule : public jura::AudioModule
{

public:

  SamplerModule(CriticalSection *lockToUse);



  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor(int type) override;

  void setSampleRate(double newSampleRate) override;
  void setGain(double newGain);

  void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  void processStereoFrame(double *left, double *right) override;

  void reset() override;


protected:

  virtual void createParameters();

  rosic::rsSamplerEngine engine;
  juce::File sfzFile;

  using ReturnCode = rosic::rsSamplerEngine::ReturnCode;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};

//=================================================================================================

/** Editor for SamplerAudioModule */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  virtual void resized() override;

protected:

  SamplerModule* samplerModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};