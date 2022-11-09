#ifndef jura_UnitTestsSampler_h
#define jura_UnitTestsSampler_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

/** A Sampler meant for copy-and-paste to create a new unit test. */

class JUCE_API jura_UnitTestsSampler_h : public juce::UnitTest
{

public:


  jura_UnitTestsSampler_h() : juce::UnitTest("Tamplate", "Processors") {}

  void runTest() override;


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(jura_UnitTestsSampler_h)
};

#endif