
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
