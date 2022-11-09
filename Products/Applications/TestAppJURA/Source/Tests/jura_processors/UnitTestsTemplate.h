#ifndef jura_UnitTestsTemplate_h
#define jura_UnitTestsTemplate_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

/** A template meant for copy-and-paste to create a new unit test. */

class JUCE_API jura_UnitTestsTemplate_h : public juce::UnitTest
{

public:


  jura_UnitTestsTemplate_h() : juce::UnitTest("Tamplate", "Processors") {}

  void runTest() override;


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(jura_UnitTestsTemplate_h)
};

#endif