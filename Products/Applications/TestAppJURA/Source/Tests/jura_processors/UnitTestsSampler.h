#ifndef jura_UnitTestsSampler_h
#define jura_UnitTestsSampler_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

/**  */

class JUCE_API UnitTestsSampler : public juce::UnitTest
{

public:


  UnitTestsSampler() : juce::UnitTest("Sampler", "Processors") {}

  void runTest() override;


protected:

  void testSamplerAudioModule();


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsSampler)
};




class JUCE_API SamplerModuleTest : public jura::SamplerModule
{

public:

  using jura::SamplerModule::SamplerModule; // inherit constructors

  jura::AudioModuleEditor* createEditor(int type) override;

};


class JUCE_API SamplerEditorTest : public jura::SamplerEditor
{

public:

  using jura::SamplerEditor::SamplerEditor; // inherit constructors

  bool testTreeViewNodeSelection();


};



#endif