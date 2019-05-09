#pragma once

class JUCE_API BlepOscArrayModule : public jura::AudioModule
{

public:

  BlepOscArrayModule(CriticalSection *lockToUse/*, rosic::rsDualFilter* coreToUse*/);

  AudioModuleEditor* createEditor(int type) override;


protected:


  rosic::rsOscArrayPolyBlep1 core;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayModule)
};

//=================================================================================================

class JUCE_API BlepOscArrayEditor : public jura::AudioModuleEditor
{

public:

  BlepOscArrayEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit);

  virtual void resized() override;

protected:

  BlepOscArrayModule* module;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayEditor)
};