
/** A polyphonic modulator section. The user can plug in an arbitrary number of modulator 
modules. */

class JUCE_API PolyModulatorsAudioModule : public jura::PolySlotAudioModule
{

public:

  PolyModulatorsAudioModule(CriticalSection *lockToUse);

  AudioModuleEditor* createEditor() override;

protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyModulatorsAudioModule)
};

//=================================================================================================

class JUCE_API PolyModulatorsEditor : public jura::AudioModuleEditor
{

public:

  PolyModulatorsEditor(CriticalSection* lockToUse, PolyModulatorsAudioModule* modulatorsToEdit);

protected:

  PolyModulatorsAudioModule* modulatorsModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyModulatorsEditor)
};