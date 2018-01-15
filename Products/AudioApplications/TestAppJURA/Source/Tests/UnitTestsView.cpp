#include "UnitTestsView.h"
using namespace jura;

UnitTestsView::UnitTestsView()
{
  createWidgets();
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

void UnitTestsView::createWidgets()
{
  runTestsLabel  = new RTextField("Select Test:");
  addWidget(runTestsLabel);

  testSelectorBox = new RComboBox;
  testSelectorBox->addItem(ALL,        "All Tests");
  testSelectorBox->addItem(PARAMETERS, "Parameters");
  addWidget(testSelectorBox);

  runButton = new RClickButtonNotifyOnMouseUp("Run");
  addWidget(runButton);

  testResultView = new RTextEditor;
  addWidget(testResultView);

  // todo: add descriptions
}

