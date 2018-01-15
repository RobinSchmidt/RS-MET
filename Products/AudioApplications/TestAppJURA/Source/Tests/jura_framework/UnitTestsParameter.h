#ifndef jura_UnitTestsParameter_h
#define jura_UnitTestsParameter_h  

#include "../../../JuceLibraryCode/JuceHeader.h"


/** Unit tests for jura::Parameter and its various subclasses. */

class JUCE_API UnitTestParameter : public juce::UnitTest
{

public:

  UnitTestParameter() : juce::UnitTest("Parameter", "Control") {}

  virtual void runTest() override;

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
  //jura::rsSmoothingManager   smoothingManager;
  //jura::MetaParameterManager metaManager;
  //jura::ModulationManager    modManager(&lock); // doesn't work



};


#endif