
class JUCE_API QuadSourceAudioModule : public jura::AudioModule // maybe rename to QuadPolySource
{

public:

  QuadSourceAudioModule(CriticalSection *lockToUse/*, rosic::rsQuadSource* coreToUse*/);

  AudioModuleEditor* createEditor(int type) override;

  void setModuleFactory(AudioModuleFactory* newFactory);


protected:

  PolySlotAudioModule *topLeftModule, *topRightModule, *bottomLeftModule, *bottomRightModule;
  PolyAudioModule *vectorMixerModule; // preliminary
  //VectorMixerAudioModule *vectorMixerModule;

  //rosic::rsQuadSource* sourceCore; 
    // Perhaps it's a better idea to not delegate the dsp to a rsQuadSource object. Doing it 
    // directly here will give better reusability (it can be used independently from rsQuadSource)
    // The same consideration applies to the DualFilterAudioModule

  friend class QuadSourceEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceAudioModule)
};

//=================================================================================================

class JUCE_API QuadSourceEditor : public jura::AudioModuleEditor
{

public:

  QuadSourceEditor(CriticalSection* lockToUse, QuadSourceAudioModule* sourceToEdit);

  virtual void resized() override;

protected:

  PolySlotEditor *topLeftEditor, *topRightEditor, *bottomLeftEditor, *bottomRightEditor;
  AudioModuleEditor *vectorMixerEditor;
  QuadSourceAudioModule* sourceModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(QuadSourceEditor)
};
