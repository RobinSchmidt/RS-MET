#include "UnitTestsParameter.h"
using namespace juce;
using namespace jura;

/** A ModulationSource subclass that just outputs a constant value. To be used for testing 
ModulatableParameter. */
class JUCE_API ConstantModulationSource : public ModulationSource
{
public:
  void setModulationValue(double newValue) { modValue = newValue; }
  virtual void updateModulationValue() override 
  {
    // nothing to do here, modValue is just a constant
  }
};

UnitTestParameter::UnitTestParameter() 
  : juce::UnitTest("Parameter", "Control"), modManager(&lock)
{
  smoothingManager.setMutexLock(&lock);
}

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
  resetCounters();

  jura::Parameter p("TestParam", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  expectEquals(numCallbacksReceived, 1);
  expectEquals(lastCallbackValue,      5.0);
  expectEquals(p.getValue(),           5.0);
  expectEquals(p.getNormalizedValue(), 0.5);

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
  beginTest("ModulatableParameter");
  resetCounters();

  ConstantModulationSource modSource;
  modManager.registerModulationSource(&modSource);
  modSource.setModulationValue(1.0);

  jura::ModulatableParameter p("TestParam", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  p.setSmoothingManager(&smoothingManager);
  p.setMetaParameterManager(&metaManager);
  p.setModulationManager(&modManager);
  modManager.registerModulationTarget(&p);
  expectEquals(numCallbacksReceived, 1);
  expectEquals(lastCallbackValue,      5.0);
  expectEquals(p.getValue(),           5.0);
  expectEquals(p.getNormalizedValue(), 0.5);

  // establish a modulation connection:
  modManager.addConnection(&modSource, &p);
  jura::ModulationConnection* modCon = modManager.getConnectionBetween(&modSource, &p);
  modCon->setDepthRangeAndValue(-10.0, +10.0, 1.0);

  // applying modulations should result in a callback to be called with value 6, the results 
  // returned by get(Normalized)Value should not be affected:
  modManager.applyModulations();
  expectEquals(numCallbacksReceived,   2);
  expectEquals(lastCallbackValue,      6.0);
  expectEquals(p.getValue(),           5.0);
  expectEquals(p.getNormalizedValue(), 0.5);



  int dummy = 0;
}

void UnitTestParameter::resetCounters()
{
  numCallbacksReceived = 0;
  numNotificationsReceived = 0; 
}