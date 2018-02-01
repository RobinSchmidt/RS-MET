
class JUCE_API DualFilterAudioModule : public jura::AudioModule
{

public:

  DualFilterAudioModule(CriticalSection *lockToUse, rosic::rsDualFilter* coreToUse);

protected:

  rosic::rsDualFilter* filterCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DualFilterAudioModule)
};

//=================================================================================================

class JUCE_API DualFilterEditor : public jura::AudioModuleEditor
{

public:

  DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit);

protected:

  DualFilterAudioModule* filterModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DualFilterEditor)
};