#ifndef jura_UnitTestsParameter_h
#define jura_UnitTestsParameter_h  

#include "../../../JuceLibraryCode/JuceHeader.h"


/** Unit tests for jura::Parameter and its various subclasses. */

class JUCE_API UnitTestParameter : public juce::UnitTest, public jura::ParameterObserver
{

public:

  UnitTestParameter() : juce::UnitTest("Parameter", "Control") {}

  virtual void runTest() override;
  virtual void parameterChanged(jura::Parameter* parameterThatHasChanged) override;

  /** Used as target function for the callback in parameter objects. */
  void callbackTargetDouble(double value);

protected:

  void runTestParameter();
  void runTestSmoothableParameter();
  void runTestMetaControlledParameter();
  void runTestModulatableParameter();

  double lastCallbackValue;
  int numCallbacksReceived;
  int numNotificationsReceived;

  juce::CriticalSection lock;

  // managers:
  jura::rsSmoothingManager   smoothingManager;
  jura::MetaParameterManager metaManager;
  //jura::ModulationManager    modManager(&lock); // doesn't work


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestParameter)
};


#endif