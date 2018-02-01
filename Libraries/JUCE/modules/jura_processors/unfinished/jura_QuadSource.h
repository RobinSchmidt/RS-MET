
class JUCE_API QuadSourceAudioModule : public jura::AudioModule // maybe rename to QuadPolySource
{

public:

  QuadSourceAudioModule(CriticalSection *lockToUse/*, rosic::rsQuadSource* coreToUse*/);

  AudioModuleEditor* createEditor() override;


protected:

  PolySlotAudioModule *topLeftModule, *topRightModule, *bottomLeftModule, *bottomRightModule;
  //VectorMixerAudioModule *vectorMixerModule;

  //rosic::rsQuadSource* sourceCore; 
    // Perhaps it's a better idea to not delegate the dsp to a rsQuadSource object. Doing it 
    // directly here will give better reusability (it can be used independently from rsQuadSource)
    // The same consideration applies to the DualFilterAudioModule

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceAudioModule)
};

//=================================================================================================

class JUCE_API QuadSourceEditor : public jura::AudioModuleEditor
{

public:

  QuadSourceEditor(CriticalSection* lockToUse, QuadSourceAudioModule* sourceToEdit);

protected:

  PolySlotEditor *topLeftEditor, *topRightEditor, *bottomLeftEditor, *bottomRightEditor;
  QuadSourceAudioModule* sourceModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceEditor)
};
