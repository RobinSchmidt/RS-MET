#ifndef jura_UnitTestsView_h
#define jura_UnitTestsView_h  

#include "jura_framework/UnitTestsParameter.h"

/** A component to perform unit tests for jura classes, print results, etc. */

class JUCE_API UnitTestsView : public jura::Editor, public juce::UnitTestRunner, 
  public jura::RButtonListener
{

public:

  enum testIndices
  {
    ALL = 1,
    PARAMETERS,

    NUM_UNIT_TESTS
  };

  UnitTestsView();

  /** Runs the tests with given index (see enum testIndices). */
  void runTest(int testIndex);

  virtual void resized() override;
  virtual void rButtonClicked(jura::RButton* button) override;

protected:



  void createWidgets();

  // widgets:
  jura::RTextField *runTestsLabel;
  jura::RComboBox  *testSelectorBox;
  jura::RClickButtonNotifyOnMouseUp *runButton;
  jura::RTextEditor *testResultView;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsView)
};

#endif