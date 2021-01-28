#include "UnitTestsParameter.h"
using namespace juce;
using namespace jura;

/** A ModulationSource subclass that just outputs a constant value. To be used for testing 
ModulatableParameter. */
class JUCE_API ConstantModulationSource : public ModulationSource
{
public:
  void setModulationValue(double newValue) { modValue = newValue; }

  virtual double renderModulation() override
  {
    return modValue;
  }

  //virtual void updateModulationValue() override 
  //{
  //  // nothing to do here, modValue is just a constant
  //}
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

  testParameter( &p);
  testSmoothable(&p);
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
  resetCounters();
  smoothingManager.setSampleRate(1000.0);

  expectEquals(p->getSmoothingTime(), 0.0); // 0 should be the default setting
  p->setSmoothingTime(2.0);                 // 2 milliseconds

  // init:
  p->setNormalizedValue(0.0, false, false);
  expectEquals(p->getValue(),            0.0);
  expectEquals(p->getNormalizedValue(),  0.0);
  expectEquals(numCallbacksReceived,     0);
  expectEquals(numNotificationsReceived, 0);

  // set a target value:
  p->setNormalizedValue(0.5, true, true);
  expectEquals(p->getValue(),            5.0);
  expectEquals(p->getNormalizedValue(),  0.5);
  expectEquals(numCallbacksReceived,     0);   // callback should be defered to smoother update
  expectEquals(numNotificationsReceived, 1);

  // perform smoothing:
  int i = doSmoothingUntilDone();
  expectEquals(p->getValue(),            5.0);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p->getNormalizedValue(),  0.5);
  expectEquals(numCallbacksReceived,     i);
  expectEquals(numNotificationsReceived, 2);   // receives pre- and post-smoothing notification

  p->setSmoothingTime(0.0); // reset to not thwart subsequent tests
}

void UnitTestParameter::testMetaControl(jura::MetaControlledParameter* p)
{
  resetCounters();

  // test remapping:

  // set up a nonmonotonic mapping function:
  rsMetaParameterMapper* mapper = p->getMetaMapper(); // (0,0),(1,1)
  mapper->addNode(0.5, 1.0);                          // (0,0),(0.5,1),(1,1)
  mapper->moveNode(2, 1.0, 0.0);                      // (0,0),(0.5,1),(1,0)

  // init:
  p->setNormalizedValue(0.0, false, false);
  expectEquals(p->getValue(),            0.0);
  expectEquals(p->getNormalizedValue(),  0.0);
  expectEquals(numCallbacksReceived,     0);
  expectEquals(numNotificationsReceived, 0);

  p->setNormalizedValue(0.25, true, true);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p->getValue(),            5.0);
  expectEquals(p->getNormalizedValue(),  0.25);
  expectEquals(numCallbacksReceived,     1);
  expectEquals(numNotificationsReceived, 1);

  p->setNormalizedValue(0.5, true, true);
  expectEquals(lastCallbackValue,        10.0);
  expectEquals(p->getValue(),            10.0);
  expectEquals(p->getNormalizedValue(),   0.5);
  expectEquals(numCallbacksReceived,      2);
  expectEquals(numNotificationsReceived,  2);

  p->setNormalizedValue(0.75, true, true);
  expectEquals(lastCallbackValue,        5.0);
  expectEquals(p->getValue(),            5.0);
  expectEquals(p->getNormalizedValue(),  0.75);
  expectEquals(numCallbacksReceived,     3);
  expectEquals(numNotificationsReceived, 3);

  // maybe somehow adding/moving/removing nodes should also trigger callbacks

  // test mapping with smoothing:
  p->setSmoothingTime(2.0);
  p->setNormalizedValue(0.25, true, true);
  expectEquals(p->getValue(),            5.0);
  expectEquals(p->getNormalizedValue(),  0.25);
  expectEquals(numNotificationsReceived, 4);
  expectEquals(numCallbacksReceived,     3);
  int i = 3;  // number of callbacks
  i += doSmoothingUntilDone();
  expectEquals(i, 4);         // 0.25 and 0.75 both map to 5, so we reach it immediately
  expectEquals(numCallbacksReceived,     i);
  expectEquals(numNotificationsReceived, 5);  // because of post-smoothing notification
  expectEquals(lastCallbackValue,        5.0);

  p->setNormalizedValue(0.5, true, true);
  expectEquals(p->getValue(),           10.0);
  expectEquals(p->getNormalizedValue(),  0.5);
  expectEquals(numNotificationsReceived, 6);
  expectEquals(numCallbacksReceived,     i);
  i += doSmoothingUntilDone();
  expectEquals(i, 16);  
  expectEquals(numCallbacksReceived,     i);
  expectEquals(lastCallbackValue,       10.0);

  p->setSmoothingTime(0.0); 

  // todo: test meta-attachment, cross-coupling


  mapper->initToIdentity(); // to not thwart subsequent tests
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

  // init:
  p->setNormalizedValue(0.0, false, false); 
   // maybe get rid of this, but then we will hit a jassert because setting the normalizedValue in 
   // the next call will immediately return, normalized and unnormalized are out of sync because
   // of the resetting of the mapping function which didn't trigger an update of dependent 
   // variables

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

  modCon->setDepth(2.0);
  modManager.applyModulations();
  expectEquals(lastCallbackValue,    7.0);
  expectEquals(numCallbacksReceived, 3); 

  // todo: test different mod-modes, limiting of the final modulated value, how it works together
  // with smoothing


  modManager.deRegisterModulationTarget(p);
  modManager.deRegisterModulationSource(&modSource);
}

int UnitTestParameter::doSmoothingUntilDone()
{
  int i = 0; // iteration counter
  while(smoothingManager.needsSmoothing())
  {
    smoothingManager.updateSmoothedValuesNoLock();
    i++;
  }
  return i;
}
