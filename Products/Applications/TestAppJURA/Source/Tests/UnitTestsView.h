#ifndef jura_UnitTestsView_h
#define jura_UnitTestsView_h  

#include "jura_framework/UnitTestsParameter.h"
#include "jura_framework/UnitTestsModulation.h"
#include "jura_framework/UnitTestsMisc.h"

#include "jura_processors/UnitTestsToolChain.h"
#include "jura_processors/UnitTestsSampler.h"


/** A component to perform unit tests for jura classes, print results, etc. */

class JUCE_API UnitTestsView : public jura::Editor, public juce::UnitTestRunner, 
  public jura::RButtonListener
{

public:

  enum class testIndices  // rename to TestIndex
  {
    ALL = 1,
    PARAMETERS,
    MODULATION,
    MISC,
    TOOL_CHAIN,
    SFZ_SAMPLER,


    NUM_UNIT_TESTS
  };
  // todo: use an enum class

  UnitTestsView();

  /** Runs the tests with given index (see enum testIndices). */
  void runTest(testIndices testIndex);

  virtual void resized() override;
  virtual void rButtonClicked(jura::RButton* button) override;
  virtual void resultsUpdated() override;

protected:

  /** Returns true, if the test with given index should be included in the next run. This 
  depends on what the user hase selected via the testSelectorBox. */
  bool includeTest(testIndices testIndex);

  void createWidgets();

  // widgets:
  jura::RTextField *runTestsLabel;
  jura::RComboBox  *testSelectorBox;
  jura::RClickButtonNotifyOnMouseUp *runButton;
  jura::RTextEditor *testResultView;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsView)
};

#endif