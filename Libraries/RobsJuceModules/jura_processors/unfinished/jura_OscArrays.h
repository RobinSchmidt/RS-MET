#pragma once




// maybe rename to AnlogOscArray, MorphOscArray (when we later allow waveshape-morphing)
class JUCE_API BlepOscArrayModule : public jura::AudioModule
{

public:

  BlepOscArrayModule(CriticalSection *lockToUse);

  //AudioModuleEditor* createEditor(int type) override;
  // override this later to create a custom editor - for the moment, we use the generic editor

protected:

  void createParameters();

  RAPT::rsRatioGenerator<double> ratioGenerator;
  rosic::rsOscArrayPolyBlep1 oscArrayCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayModule)
};

//=================================================================================================

class JUCE_API BlepOscArrayEditor : public jura::AudioModuleEditor
{

public:

  BlepOscArrayEditor(BlepOscArrayModule* oscArrayToEdit);

  virtual void resized() override;

protected:

  BlepOscArrayModule* oscArrayModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BlepOscArrayEditor)
};