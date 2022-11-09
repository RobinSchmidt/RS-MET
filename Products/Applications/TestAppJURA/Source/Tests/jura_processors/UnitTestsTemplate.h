#ifndef jura_UnitTestsTemplate_h
#define jura_UnitTestsTemplate_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

/** A template meant for copy-and-paste to create a new unit test. */

class JUCE_API UnitTestsTemplate : public juce::UnitTest
{

public:


  UnitTestsTemplate() : juce::UnitTest("Template", "Processors") {}

  void runTest() override;


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsTemplate)
};

#endif