#ifndef jura_UnitTestsMisc_h
#define jura_UnitTestsMisc_h  

#include "../../../JuceLibraryCode/JuceHeader.h"


/** Unit tests for miscellanneous smaller classes and functions of jura. */

class JUCE_API UnitTestMisc : public juce::UnitTest
{

public:

  UnitTestMisc();

  void runTest() override;

protected:

  // called from runTest:
  void runTestColor();
  void runTestFileManager();



  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestMisc)
};

#endif