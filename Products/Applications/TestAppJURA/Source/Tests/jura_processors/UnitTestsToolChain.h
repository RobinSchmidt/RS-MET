#ifndef jura_UnitTestsToolChain_h
#define jura_UnitTestsToolChain_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

class JUCE_API UnitTestToolChain : public juce::UnitTest
{

public:


  UnitTestToolChain() : juce::UnitTest("ToolChain", "Processors") {}

  void runTest() override;


protected:

  // Helper functions:
  bool isInDefaultState(const jura::AudioModule* m);



  // called from runTest:
  void runTestVoiceManager();
  void runTestEqualizer();
  void runTestWaveOscillator();
  void runTestQuadrifex();
  void runTestEditorCreation();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestToolChain)
};

#endif