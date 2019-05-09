#pragma once

/** The FooBarModule class serves as a template for creating new AudioModule subclasses. */

class JUCE_API FooBarModule : public jura::AudioModule
{

public:

  FooBarModule(CriticalSection *lockToUse);

  AudioModuleEditor* createEditor(int type) override;

protected:

  rosic::rsFooBar core;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FooBarModule)
};

//=================================================================================================

/** Template for creating editors for new AudioModule subclasses. */

class JUCE_API FooBarEditor : public jura::AudioModuleEditor
{

public:

  FooBarEditor(FooBarModule* oscArrayToEdit);

  virtual void resized() override;

protected:

  FooBarModule* fooBarModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FooBarEditor)
};