#ifndef jura_UnitTestsModulation_h
#define jura_UnitTestsModulation_h  

#include "../../../JuceLibraryCode/JuceHeader.h"

class JUCE_API UnitTestModulation : public juce::UnitTest
{

public:

  UnitTestModulation();

  virtual void runTest() override;

protected:

  void reset();


  // called from runTest:
  void runTestPolyToMono();
  void runTestPolyToPoly();



  juce::CriticalSection lock;
  jura::rsVoiceManager voiceMan;
  jura::ModulationManagerPoly modMan;
  std::vector<double> voiceBuffer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestModulation)
};


#endif