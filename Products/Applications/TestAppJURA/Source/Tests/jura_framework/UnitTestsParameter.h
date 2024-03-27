#ifndef jura_UnitTestsParameter_h
#define jura_UnitTestsParameter_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

namespace jura
{

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

  // Called from runTest:
  void runTestBasicParameter();  // Maybe rename to runTestBasicParameter
  void runTestSmoothableParameter();
  void runTestMetaControlledParameter();
  void runTestModulatableParameter();
  void runTestParameterMapping();


  // Functions that perform tests on the passed Parameter pointers:
  void testCallbacks_0_10(jura::Parameter* p);        // _0_10 is the assumed range of p
  void testCallbacksSmoothable_0_10(jura::rsSmoothableParameter* p);
  void testCallbacksMetaControlled_0_10(jura::MetaControlledParameter* p);
  void testCallbacksModulated_0_10(jura::ModulatableParameter* p);

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

}

#endif