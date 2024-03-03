#ifndef jura_UnitTestsToolChain_h
#define jura_UnitTestsToolChain_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

namespace jura
{

/**


*/

class JUCE_API UnitTestToolChain : public juce::UnitTest
{

public:


  UnitTestToolChain() : juce::UnitTest("ToolChain", "Processors") {}

  void runTest() override;


protected:


  //-----------------------------------------------------------------------------------------------
  // \name Helper functions

  /** Checks if all parameters of the module have their default values. */
  bool isInDefaultState(const jura::AudioModule* m);

  /** Resets all the parameters of the given jura::AudioModule to their default values. */
  void resetParameters(jura::AudioModule* m);

  /** Randomizes all the parameters of the given jura::AudioModule. */
  void randomizeParameters(jura::AudioModule* m, int seed = 0);

  /** Generates a mock mouse event t´hat can be used to test GUI stuff. */
  juce::MouseEvent getMockMouseEvent(float mouseX = 0.f, float mouseY = 0.f);

  //-----------------------------------------------------------------------------------------------
  // \name Tests

  // Tests of the infrastructure:
  void runTestVoiceManager();
  void runTestEditorCreation(int seed);
  void runTestStateRecall(int seed);
  // Some tests use some randomization internally. They can be called with a seed parameter for the 
  // PRNG.

  // Tests of individual modules:
  void runTestEqualizer();
  void runTestMultiAnalyzer();
  void runTestQuadrifex();
  void runTestStraightliner();
  void runTestWaveOscillator();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestToolChain)
};


}

#endif