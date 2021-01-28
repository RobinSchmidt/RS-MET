// This infrastucture is not used anywhere, the files are not even included by
// jura_framework.h...maybe get rid of them...
// or maybe move code formAudioModulePoly here

class JUCE_API PolyAudioModule : public jura::AudioModule // maybe derive from ModulatableAudioModule
{

public:

  PolyAudioModule(CriticalSection *lockToUse) : AudioModule(lockToUse) {}

  // override audio processing function to process all voices (but don't add them up - there may be
  // other poly modules down the signal chain, that need them separately)

protected:

  // VoiceManager* voiceManager;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyAudioModule)
};
