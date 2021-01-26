#ifndef jura_UnitTestsModulation_h
#define jura_UnitTestsModulation_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

class JUCE_API UnitTestModulation : public juce::UnitTest
{

public:

  UnitTestModulation();

  virtual void runTest() override;

protected:

  // called from runTest:
  void runTestPolyModulation();

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestModulation)
};


#endif