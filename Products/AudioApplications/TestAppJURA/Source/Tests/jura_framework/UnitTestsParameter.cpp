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

void UnitTestParameter::resetCounters()
{
  numCallbacksReceived = 0;
  numNotificationsReceived = 0; 
}

void UnitTestParameter::runTestParameter()
{
  beginTest("Parameter");
  resetCounters();

  jura::Parameter p("Parameter", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  expectEquals(numNotificationsReceived, 0);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p.getValue(),             5.0);
  expectEquals(p.getNormalizedValue(),   0.5);

  testParameter(&p);
}

void UnitTestParameter::runTestSmoothableParameter()
{
  beginTest("SmoothableParameter");
  resetCounters();

  jura::rsSmoothableParameter p("SmoothableParameter", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  p.setSmoothingManager(&smoothingManager);
  expectEquals(numNotificationsReceived, 0);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p.getValue(),             5.0);
  expectEquals(p.getNormalizedValue(),   0.5);

  testParameter(  &p);
  testSmoothable( &p);
}

void UnitTestParameter::runTestMetaControlledParameter()
{
  beginTest("MetaControlledParameter");
  resetCounters();

  jura::MetaControlledParameter p("MetaControlledParameter", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  p.setSmoothingManager(&smoothingManager);
  p.setMetaParameterManager(&metaManager);
  expectEquals(numNotificationsReceived, 0);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p.getValue(),             5.0);
  expectEquals(p.getNormalizedValue(),   0.5);

  testParameter(  &p);
  testSmoothable( &p);
  testMetaControl(&p);
}

void UnitTestParameter::runTestModulatableParameter()
{
  beginTest("ModulatableParameter");
  resetCounters();

  jura::ModulatableParameter p("ModulatableParameter", 0.0, 10.0, 5.0);
  p.registerParameterObserver(this);
  p.setValueChangeCallback<UnitTestParameter>(this, &UnitTestParameter::callbackTargetDouble);
  p.setSmoothingManager(&smoothingManager);
  p.setMetaParameterManager(&metaManager);
  p.setModulationManager(&modManager);
  expectEquals(numNotificationsReceived, 0);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p.getValue(),             5.0);
  expectEquals(p.getNormalizedValue(),   0.5);

  testParameter(  &p);
  testSmoothable( &p);
  testMetaControl(&p);
  testModulation( &p);
}

void UnitTestParameter::testParameter(jura::Parameter* p)
{
  resetCounters();

  p->setValue(2.5, false, false);
  expectEquals(p->getValue(),             2.5);
  expectEquals(p->getNormalizedValue(),   0.25);
  expectEquals(numCallbacksReceived,      0);
  expectEquals(numNotificationsReceived,  0);

  p->setValue(7.5, true, true);
  expectEquals(p->getValue(),            7.5);
  expectEquals(lastCallbackValue,        7.5);
  expectEquals(p->getNormalizedValue(),  0.75);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(numNotificationsReceived, 1);

  p->setNormalizedValue(0.25, true, true);
  expectEquals(p->getValue(),            2.5);
  expectEquals(lastCallbackValue,        2.5);
  expectEquals(p->getNormalizedValue(),  0.25);
  expectEquals(numCallbacksReceived,     2);
  expectEquals(numNotificationsReceived, 2);

  // test custom mapper, value quantization, values outside range, calling multiple times with
  // the same value, state recall

  int dummy = 0;
}

void UnitTestParameter::testSmoothable(jura::rsSmoothableParameter* p)
{

}

void UnitTestParameter::testMetaControl(jura::MetaControlledParameter* p)
{

}

void UnitTestParameter::testModulation(jura::ModulatableParameter* p)
{
  resetCounters();

  modManager.registerModulationTarget(p);

  ConstantModulationSource modSource;
  modManager.registerModulationSource(&modSource);
  modSource.setModulationValue(1.0);

  // establish a modulation connection:
  modManager.addConnection(&modSource, p);
  jura::ModulationConnection* modCon = modManager.getConnectionBetween(&modSource, p);
  modCon->setDepthRangeAndValue(-10.0, +10.0, 1.0);

  p->setNormalizedValue(0.5, true, true);
  expectEquals(p->getNormalizedValue(),  0.5);
  expectEquals(p->getValue(),            5.0);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(numNotificationsReceived, 1);

  // applying modulations should result in a callback to be called with value 6, the results 
  // returned by get(Normalized)Value should not be affected:
  modManager.applyModulations();
  expectEquals(lastCallbackValue,        6.0);
  expectEquals(p->getValue(),            5.0);
  expectEquals(p->getNormalizedValue(),  0.5);
  expectEquals(numCallbacksReceived,     2);  
  expectEquals(numNotificationsReceived, 1);

  modManager.deRegisterModulationTarget(p);
  modManager.deRegisterModulationSource(&modSource);
}

