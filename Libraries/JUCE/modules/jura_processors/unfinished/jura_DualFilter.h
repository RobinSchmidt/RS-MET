
class JUCE_API DualFilterAudioModule : public jura::AudioModule // maybe rename to DualPolyFilter
{

public:

  DualFilterAudioModule(CriticalSection *lockToUse/*, rosic::rsDualFilter* coreToUse*/);

  AudioModuleEditor* createEditor(int type) override;

  void setModuleFactory(AudioModuleFactory* newFactory);

protected:

  PolySlotAudioModule *leftModule, *rightModule;
  PolyAudioModule *vectorMixerModule; // preliminary
  //VectorMixerAudioModule *vectorMixerModule;

  //rosic::rsDualFilter* filterCore;

  friend class DualFilterEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DualFilterAudioModule)
};

//=================================================================================================

class JUCE_API DualFilterEditor : public jura::AudioModuleEditor
{

public:

  DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit);

  virtual void resized() override;

protected:

  PolySlotEditor *leftEditor, *rightEditor;
  AudioModuleEditor *vectorMixerEditor;

  DualFilterAudioModule* filterModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DualFilterEditor)
};