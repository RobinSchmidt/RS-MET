/** An AudioModule that allows a PolyAudioModule to be plugged in. Users can select which type of
module they want to use in this slot. */

class JUCE_API PolySlotAudioModule : public jura::PolyAudioModule
{

public:

  PolySlotAudioModule(CriticalSection *lockToUse) : PolyAudioModule(lockToUse) {}

  AudioModuleEditor* createEditor() override;

  // override get/setState functions to transparently delegate to the slotInsert (maybe adding some
  // additional infor, if necessarry)

  // audio-processing functions should delegate to the slotInsert, too

protected:

  PolyAudioModule* slotInsert = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolySlotAudioModule)
};

//=================================================================================================

class JUCE_API PolySlotEditor : public jura::AudioModuleEditor
{

public:

  PolySlotEditor(CriticalSection* lockToUse, PolySlotAudioModule* slotToEdit);

protected:


  PolySlotAudioModule* slotModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolySlotEditor)
};