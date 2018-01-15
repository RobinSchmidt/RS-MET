#include "UnitTestsView.h"
using namespace jura;

UnitTestsView::UnitTestsView()
{
  createWidgets();
}

void UnitTestsView::createWidgets()
{
  runTestsLabel  = new RTextField("Run Tests:");
  addWidget(runTestsLabel);




  //RTextEditor *testResultView;
  //RClickButtonNotifyOnMouseUp *runButton;

}

