#ifndef jura_DebugModule_h
#define jura_DebugModule_h

/** A module for framework debugging puposes */

class JUCE_API DebugAudioModule : public jura::AudioModule
{

public:

  DebugAudioModule(CriticalSection *lockToUse);
    
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DebugAudioModule)
};

#endif 