
class JUCE_API QuadSourceAudioModule : public jura::AudioModule
{

public:

  QuadSourceAudioModule(CriticalSection *lockToUse, rosic::rsQuadSource* coreToUse);


protected:

  rosic::rsQuadSource* sourceCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceAudioModule)
};

//=================================================================================================

class JUCE_API QuadSourceEditor : public jura::AudioModuleEditor
{

public:

protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceEditor)
};
