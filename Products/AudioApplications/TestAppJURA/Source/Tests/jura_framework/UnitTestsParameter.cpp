#include "UnitTestsParameter.h"
using namespace juce;
using namespace jura;

void UnitTestParameter::runTest()
{
  runTestParameter();
  runTestSmoothableParameter();
  runTestMetaControlledParameter();
  runTestModulatableParameter();
}

void UnitTestParameter::parameterChanged(Parameter* parameterThatHasChanged)
{
  numNotificationsReceived++;
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
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  expectEquals(numCallbacksReceived,   1);

  p.setValue(2.5, false, false);
  expectEquals(p.getValue(),             2.5);
  expectEquals(p.getNormalizedValue(),   0.25);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(numNotificationsReceived, 0);

  p.setValue(7.5, true, true);
  expectEquals(p.getValue(),             7.5);
  expectEquals(lastCallbackValue,        7.5);
  expectEquals(p.getNormalizedValue(),   0.75);
  expectEquals(numCallbacksReceived,     2);
  expectEquals(numNotificationsReceived, 1);

  p.setNormalizedValue(0.25, true, true);
  expectEquals(p.getValue(),             2.5);
  expectEquals(lastCallbackValue,        2.5);
  expectEquals(p.getNormalizedValue(),   0.25);
  expectEquals(numCallbacksReceived,     3);
  expectEquals(numNotificationsReceived, 2);

  // test custom mapper, value quantization



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