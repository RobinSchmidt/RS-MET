#ifndef jura_UnitTestsToolChain_h
#define jura_UnitTestsToolChain_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

class JUCE_API UnitTestToolChain : public juce::UnitTest
{

public:


  UnitTestToolChain() : juce::UnitTest("ToolChain", "Processors") {}

  void runTest() override;


protected:

  // called from runTest:
  void runTestVoiceManager();
  void runTestWaveOscillator();
  void runTestQuadrifex();
  void runTestEditorCreation();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestToolChain)
};

#endif