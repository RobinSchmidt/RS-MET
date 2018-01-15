#include "UnitTestsParameter.h"

void UnitTestParameter::runTest()
{
  runTestParameter();
  runTestSmoothableParameter();
  runTestMetaControlledParameter();
  runTestModulatableParameter();
}

void UnitTestParameter::callbackTargetDouble(double value)
{
  lastCallbackValue = value;
  numCallbacksReceived++;
}

void UnitTestParameter::runTestParameter()
{
  beginTest("Parameter");

  numCallbacksReceived = 0;
  numNotificationsReceived = 0; 

  jura::Parameter p("TestParam", 0.0, 10.0, 5.0);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  expectEquals(numCallbacksReceived, 1,   "");

  p.setValue(2.5, false, false);
  expectEquals(p.getValue(),           2.5,  "Parameter::getValue failed");
  expectEquals(p.getNormalizedValue(), 0.25, "Parameter::getNormalizedValue failed");
  expectEquals(numCallbacksReceived,   1, "");

  p.setValue(7.5, true, true);
  expectEquals(p.getValue(),           7.5,  "Parameter::getValue failed");
  expectEquals(p.getNormalizedValue(), 0.75, "Parameter::getNormalizedValue failed");
  expectEquals(numCallbacksReceived,   2,   "");
  expectEquals(lastCallbackValue,      7.5, "");

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