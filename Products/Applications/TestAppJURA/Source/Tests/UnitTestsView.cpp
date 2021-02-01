#include "UnitTestsView.h"
using namespace jura;

UnitTestsView::UnitTestsView()
{
  createWidgets();
}

void UnitTestsView::runTest(int testIndex)
{
  juce::Array<UnitTest*> tests;

  if(includeTest(PARAMETERS)) tests.add(new UnitTestParameter);
  if(includeTest(MODULATION)) tests.add(new UnitTestModulation);
  if(includeTest(TOOL_CHAIN)) tests.add(new UnitTestToolChain);

  //beginTest();
  runTests(tests);



  for(int i = 0; i < tests.size(); i++)
    delete tests[i];
  tests.clear();
}

void UnitTestsView::resized()
{
  int m  = 4;   // margin
  int x  = m;
  int y  = m;

  runTestsLabel  ->setBounds(x, y,  76, 20); x +=  76+m;
  testSelectorBox->setBounds(x, y, 120, 20); x += 120+m;
  runButton      ->setBounds(x, y,  40, 20); x +=  40+m;

  x = 0;
  y = runButton->getBottom() + m;
  int w = getWidth();
  int h = getHeight();
  testResultView->setBounds(x, y, w, h);
}

void UnitTestsView::rButtonClicked(RButton* b)
{
  if(b == runButton)
    runTest(ALL); // todo: inquire combobox, which test to run
}

void UnitTestsView::resultsUpdated()
{
  // this is a "shlemiel-the-painter" algorithm - try to avoid that later
  int numResults = getNumResults();
  juce::String str;
  for(int i = 0; i < numResults; i++)
  {
    const TestResult* result = getResult(i);
    str += result->unitTestName; 
    //str += result->subcategoryName;
    str += ": Passes: " + String(result->passes);
    str += ", Fails: "  + String(result->failures);
    //str += ", Messages: " + String(result->messages());
    str += "\n";
  }
  testResultView->setText(str);
}

bool UnitTestsView::includeTest(int testIndex)
{
  juce::String selection = testSelectorBox->getSelectedItemText();
  if(selection == "All Tests")                             return true;
  if(selection == "Parameters" && testIndex == PARAMETERS) return true;
  if(selection == "Modulation" && testIndex == MODULATION) return true;
  if(selection == "ToolChain"  && testIndex == TOOL_CHAIN) return true;
  return false;
}

void UnitTestsView::createWidgets()
{
  runTestsLabel  = new RTextField("Select Test:");
  addWidget(runTestsLabel);

  testSelectorBox = new RComboBox;
  testSelectorBox->addItem(ALL,        "All Tests");
  testSelectorBox->addItem(PARAMETERS, "Parameters");
  testSelectorBox->addItem(MODULATION, "Modulation");
  testSelectorBox->addItem(MODULATION, "ToolChain");
  //testSelectorBox->addItem(WIDGETS, "Widgets");
  testSelectorBox->selectItemFromText("All Tests", false);
  addWidget(testSelectorBox);

  runButton = new RClickButtonNotifyOnMouseUp("Run");
  runButton->addRButtonListener(this);
  addWidget(runButton);

  testResultView = new RTextEditor;
  testResultView->setMultiLine(true);
  addWidget(testResultView);

  // todo: add descriptions
}

