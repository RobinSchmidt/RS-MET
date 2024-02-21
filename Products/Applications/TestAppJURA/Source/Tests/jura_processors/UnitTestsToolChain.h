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

  /** Checks if all parameters of the module have their default values. */
  bool isInDefaultState(const jura::AudioModule* m);

  /** Randomizes all the parameters of the given jura::AudioModule. */
  void randomizeParameters(jura::AudioModule* m, int seed = 0);



  // called from runTest:
  void runTestVoiceManager();
  void runTestEqualizer();
  void runTestMultiAnalyzer();
  void runTestWaveOscillator();
  void runTestQuadrifex();
  void runTestEditorCreation(int seed);


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestToolChain)
};

#endif