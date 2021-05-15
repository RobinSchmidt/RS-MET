#pragma once

// preliminary - eventually, this should go into the rosic module:
//#include "../../rs_testing/Prototypes/SamplerEngine.h"

/** A sampler with functionality roughly based on the sfz specification */

class JUCE_API SamplerModule : public jura::AudioModule
{

public:

  SamplerModule(CriticalSection *lockToUse);



  // overriden from AudioModule baseclass:
  virtual AudioModuleEditor* createEditor(int type) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;


protected:

  //rsSamplerEngine* core;

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