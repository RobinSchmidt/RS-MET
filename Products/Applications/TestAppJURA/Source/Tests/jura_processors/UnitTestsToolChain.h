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
  /** \name State Verification

       Functions to verify if a given module or editor is in some expected state. */

  /** Checks if all parameters of the module have their default values. */
  bool isInDefaultState(const jura::AudioModule* m);

  /** Checks if the slots of the passed ToolChain are of the correct type. When you have a pointer
  to a ToolChain object "toolChain", then you can call the function like:  

    doSlotsContain(toolChain, { "Equalizer", "FuncShaper", "None" } );

  when you expect an Equalizer in the 1st slot, etc.  */
  bool doSlotsContain(const jura::ToolChain* toolChain, 
    const std::vector<juce::String>& typeNames);

  /** The ToolChain object and its ToolChainEditor have some parallel arrays (of the plugged in 
  modules, their editors, their selectors, etc.). This function checks, if they are all 
  consistent with one another. */
  bool areArraysConsistent(jura::ToolChainEditor* toolChainEditor);
  // ToDo: Use a const pointer - it currently doesn't compile when doing so. Figure out why and 
  // fix it. It should be possible to work with a const pointer here.


  //-----------------------------------------------------------------------------------------------
  /** \name State Manipulation 
  
      Functions to manipulate the state of audio modules or their editors. This includes setting up
      parameters in certain ways like randomizing or resetting them. */

  /** Resets all the parameters of the given jura::AudioModule to their default values. */
  void resetParameters(jura::AudioModule* m);

  /** Randomizes all the parameters of the given jura::AudioModule. You can pass a seed for the
  pseudo random number generator that is used for the randomization. You can also select whether
  or not the randomization shall be applied recursively to the child modules if the given module 
  has any. */
  void randomizeParameters(jura::AudioModule* m, int randomSeed = 0, bool recursively = true);

  void randomizeAudioModule(jura::AudioModule* m, int seed = 0, bool recursively = true);
  // Not yet implemented.
  // Should randomize the state of the given AudioModule. This includes randomization of the 
  // parameters, loading random samples or other sorts of files (ColorMaps, Tunings, etc.), 
  // randomizing strings, etc. What exactly is done depends on the kind of module - we dispatch
  // depending on the kind


  //-----------------------------------------------------------------------------------------------
  /** \name Module Tests

      Function to perform a specific test on a given module or editor. */

  /** Tests the state recall for the given module by setting up the parameters randomly (using the
  given randomSeed for the PRNG), then retrieving the state-xml, then resetting all parameters to 
  defaults, then recalling the state, then rertrieving the the new, recalled state and the 
  checking the two xmls for equivalence.  */
  bool testStateRecall(jura::AudioModule* m, int randomSeed = 0);


  //-----------------------------------------------------------------------------------------------
  /** \name Element Retrieval

      Functions to retrieve certain elements of audio modules or their editors. These "elements" 
      can be widgets, parameters, sub-modules, etc. */

  /** Returns a vector of pointer to all the widgets in the given editor that have no parameter 
  assigned to them. A widget without an assigned parameter, i.e. an "orphaned" widget, may 
  indicate a bug. It may not always be a bug, though. Some widgets may legitimately have no 
  parameter assigned but these are more the exception rather than the rule. */
  std::vector<jura::RWidget*> getWidgetsWithoutParameter(jura::AudioModuleEditor* editor);

  /** Filters out only widgets of a specific type from the given widgets array. Can be used like:

  std::vector<jura::RSlider*> sliders = filterWidgets<RSlider>(widgets);

  to filter out only the sliders from a given array of widgets, for example.  */
  template<class WidgetType>
  std::vector<WidgetType*> filterWidgets(const std::vector<jura::RWidget*>& widgets);


  //-----------------------------------------------------------------------------------------------
  /** \name Mocking

      Functions for mocking GUI events. */

  /** Generates a mock mouseDown event that can be used to test GUI stuff. */
  juce::MouseEvent getMockMouseDownEvent(float mouseX = 0.f, float mouseY = 0.f, 
    juce::Component* eventComp = nullptr, juce::Component* originatorComp = nullptr);
  // What does the distinction between eventComp and originatorComp mean?


  //-----------------------------------------------------------------------------------------------
  /** \name Unit Tests

      These functions implement the actual unit tests which are called from the main runTest() 
      function. Some tests use some randomization internally. They can be called with a seed 
      parameter for the PRNG.  */

  // Tests of the infrastructure - ToDo: document them all:
  void runTestVoiceManager();

  /** Tests insertion, removal, replacement and swapping of modules. */
  void runTestSlotInsertRemoveEtc();

  // Tests that run through all modules available in ToolChain except some specified ignored ones:
  void runTestEditorCreation(int randomSeed, const juce::StringArray& ignoreList);
  void runTestStateRecall(int randomSeed, const juce::StringArray& ignoreList);
  //  Maybe these tests should take an array of strings of modules to ignore

  // Tests of individual modules:
  void runTestEqualizer();
  void runTestFuncShaper();
  void runTestMultiAnalyzer();
  void runTestQuadrifex();
  void runTestStraightliner();
  void runTestWaveOscillator();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestToolChain)
};


}

#endif