#ifndef jura_UnitTestsView_h
#define jura_UnitTestsView_h  

//#include "../../JuceLibraryCode/JuceHeader.h"
#include "jura_framework/UnitTestsParameter.h"

/** A component to perform unit tests for jura classes, print results, etc. */

class JUCE_API UnitTestsView : public jura::Editor
{

public:

  enum testIndices
  {
    ALL = 1,
    PARAMETERS,

    NUM_UNIT_TESTS
  };

  UnitTestsView();

  virtual void resized() override;

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