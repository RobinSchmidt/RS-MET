
class JUCE_API rsMidiControllerModulatorModulePoly : public ModulatorModulePoly
{

public:

  rsMidiControllerModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    setModulationSourceName("MidiController");
  }

  virtual void setMidiController(int controllerNumber, float controllerValue) override
  {
    this->controllerNumber = (char)controllerNumber;
    this->controllerValue  = (double)controllerValue;
  }

  double renderModulation(int voiceIndex) override
  {
    return controllerValue;
  }
  // Do we need to override this or does the infrastructure already make sure that the value 
  // produced by renderVoiceModulation is used in case we don't override this? -> try it.
  // Although overriding may nevertheless be advisable as optimization...

  double renderVoiceModulation(int voiceIndex) override
  {
    return controllerValue;
  }

  // we may have to override both renderModuation and renderVoiceModulation

protected:

  double controllerValue  = 0.5;  // use the middle value by default
  char   controllerNumber = 74;   // CC#74 is typically used folr filter cutoff

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsMidiControllerModulatorModulePoly)
};
// hmm - how are we supposed to distinguish between the controller numbers? maybe the user needs
// to create a controller object and somehow set the controller number from the gui? maybe, when
// midi controller is selected as modluation source when creating a connection, the user is asked 
// to enter the controller number? ...try it with acid-devil - modualte Cutoff, stepLength, etc by 
// midi-controller...but AcidDevil is actually monophonic - however, it should work for mono and 
// poly modules the same way, from a user perspective - only under the hood the handling may be 
// different. midi controllers apply to the whole instrument anyway and not to a particular voice,
// as notes do. ...but the consept of voice is used only internally - incoming midi data has 
// actually no concept of voices




















/** A polyphonic modulator section. The user can plug in an arbitrary number of modulator 
modules. */

/*
class JUCE_API PolyModulatorsAudioModule : public jura::PolySlotAudioModule
{

public:

  PolyModulatorsAudioModule(CriticalSection *lockToUse);

  AudioModuleEditor* createEditor(int type = 0) override;

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
*/