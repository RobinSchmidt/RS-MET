#include "UnitTestsParameter.h"

void UnitTestParameter::runTest()
{
  runTestParameter();
  runTestSmoothableParameter();
  runTestMetaControlledParameter();
  runTestModulatableParameter();
}

void UnitTestParameter::runTestParameter()
{
  beginTest("Parameter");

  jura::Parameter p("TestParam");

  p.setValue(0.25, false, false);
  expectEquals(p.getValue(), 0.25, "Parameter::getValue failed");

  // test, if callbacks and notifications work correctly


  int dummy = 0;
}

void UnitTestParameter::runTestSmoothableParameter()
{

}

void UnitTestParameter::runTestMetaControlledParameter()
{

}

void UnitTestParameter::runTestModulatableParameter()
{

}