#ifndef jura_UnitTestsParameter_h
#define jura_UnitTestsParameter_h  

#include "../../../JuceLibraryCode/JuceHeader.h"


/** Unit tests for jura::Parameter and its various subclasses. */

class JUCE_API UnitTestParameter : public juce::UnitTest, public jura::ParameterObserver
{

public:

  UnitTestParameter();

  void runTest() override;

  /** Overriden to increment our notification counter. */
  virtual void parameterChanged(jura::Parameter* parameterThatHasChanged) override;

  /** Used as target function for the callback in parameter objects. */
  void callbackTargetDouble(double value);

  /** Resets counters for recived notifications and callbacks. */
  void resetCounters();

protected:

  // called from runTest:
  void runTestParameter();
  void runTestSmoothableParameter();
  void runTestMetaControlledParameter();
  void runTestModulatableParameter();

  // functions that perform tests on the passed Parameter pointers:
  void testParameter(  jura::Parameter* p);
  void testSmoothable( jura::rsSmoothableParameter* p);
  void testMetaControl(jura::MetaControlledParameter* p);
  void testModulation( jura::ModulatableParameter* p);

  /** Performs smoothing iterations until target value has been reached and returns the number of
  iterations that were needed. */
  int doSmoothingUntilDone();

  // bookkeeping data:
  double lastCallbackValue;
  int numCallbacksReceived;
  int numNotificationsReceived;

  // managers:
  juce::CriticalSection lock;
  jura::rsSmoothingManager   smoothingManager;
  jura::MetaParameterManager metaManager;
  jura::ModulationManager    modManager;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestParameter)
};

#endif