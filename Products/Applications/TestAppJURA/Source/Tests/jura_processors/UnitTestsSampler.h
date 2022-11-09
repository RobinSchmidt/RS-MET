#ifndef jura_UnitTestsSampler_h
#define jura_UnitTestsSampler_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

/** A Sampler meant for copy-and-paste to create a new unit test. */

class JUCE_API UnitTestsSampler : public juce::UnitTest
{

public:


  UnitTestsSampler() : juce::UnitTest("Sampler", "Processors") {}

  void runTest() override;


protected:

  void testSamplerAudioModule();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsSampler)
};

#endif