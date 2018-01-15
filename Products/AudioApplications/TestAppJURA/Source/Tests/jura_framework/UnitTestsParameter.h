#ifndef jura_UnitTestsParameter_h
#define jura_UnitTestsParameter_h  

#include "../../../JuceLibraryCode/JuceHeader.h"


/** Unit tests for jura::Parameter and its various subclasses. */

class JUCE_API UnitTestParameter : public juce::UnitTest
{

public:

  UnitTestParameter() : juce::UnitTest("Parameter", "Control") {}

  virtual void runTest() override;

protected:


};


#endif