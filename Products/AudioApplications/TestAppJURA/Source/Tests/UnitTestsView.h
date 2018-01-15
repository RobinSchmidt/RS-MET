#ifndef jura_UnitTestView_h
#define jura_UnitTestView_h  

#include "../../JuceLibraryCode/JuceHeader.h"

/** A component to perform unit tests for jura classes, print results, etc. */

class JUCE_API UnitTestsView : public jura::Editor
{

public:

  UnitTestsView();


protected:

  void createWidgets();

  // we need a Label that says "Run Tests:", a combobox to select which tests (Parameters, Widgets,
  // .., All) and a click-button "Run"....and an output text window that informs about the results

  jura::RTextField *runTestsLabel;
  jura::RTextEditor *testResultView;
  jura::RClickButtonNotifyOnMouseUp *runButton;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsView)
};

#endif