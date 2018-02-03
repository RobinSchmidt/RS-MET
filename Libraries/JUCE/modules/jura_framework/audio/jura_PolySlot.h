/** An AudioModule that allows a PolyAudioModule to be plugged in. Users can select which type of
module they want to use in this slot. */

class JUCE_API PolySlotAudioModule : public jura::PolyAudioModule
{

public:

  PolySlotAudioModule(CriticalSection *lockToUse) : PolyAudioModule(lockToUse) {}

  AudioModuleEditor* createEditor() override;

  /** Sets the factory object that is used to create the slot insert modules. */
  void setModuleFactory(AudioModuleFactory* newFactory) { moduleFactory = newFactory; }

  // override get/setState functions to transparently delegate to the slotInsert (maybe adding some
  // additional infor, if necessarry)

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    if(slotInsert)
      slotInsert->processBlock(inOutBuffer, numChannels, numSamples);
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    if(slotInsert)
      slotInsert->processStereoFrame(left, right);
  }

protected:

  PolyAudioModule* slotInsert = nullptr;  // wrapped/inserted module
  AudioModuleFactory* moduleFactory;      // to create and plug in a new module

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolySlotAudioModule)
};

//=================================================================================================

class JUCE_API PolySlotEditor : public jura::AudioModuleEditor
{

public:

  PolySlotEditor(CriticalSection* lockToUse, PolySlotAudioModule* slotToEdit);

protected:


  PolySlotAudioModule* slotModule;      // slot module that wraps around the slot insert
  AudioModuleEditor* slotInsertEditor;  // editor for plugged in module
  AudioModuleSelector* moduleSelector;  // widget to plug in a new module into this slot

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolySlotEditor)
};