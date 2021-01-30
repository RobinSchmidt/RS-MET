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

  /** Convenience function to add a modulation connection to our modMan. */
  void addConnection(jura::ModulationSource* source, jura::ModulationTarget* target, double depth)
  {
    jura::ModulationConnection* c = new jura::ModulationConnection(source, target, nullptr);
    c->setDepth(depth);
    modMan.addConnection(c);
  }

  // called from runTest:
  void runTestMonoToMono();
  void runTestMonoToPoly();
  void runTestPolyToMono();
  void runTestPolyToPoly();



  juce::CriticalSection lock;
  jura::rsVoiceManager voiceMan;
  jura::ModulationManagerPoly modMan;
  std::vector<double> voiceBuffer;

  // what about meta-control and smoothing?

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestModulation)
};


#endif