#include "UnitTestsModulation.h"
//using namespace juce;
//using namespace jura;



class JUCE_API TestPolyModulator : public jura::ModulatorModulePoly
{
  using jura::ModulatorModulePoly::ModulatorModulePoly;  // inherit constructors


  double getModulatorOutputSample(int voiceIndex) override
  {
    return 1; // preliminary
  }

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyModulator)
};

class JUCE_API TestPolyAudioModule : public jura::AudioModulePoly
{
  using jura::AudioModulePoly::AudioModulePoly;  // inherit constructors

  void processStereoFrameVoice(double* left, double* right, int voice) override
  {

  }


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyAudioModule)
};

//=================================================================================================

UnitTestModulation::UnitTestModulation() : juce::UnitTest("Modulation", "Control")
{

}

void UnitTestModulation::runTest()
{
  runTestPolyModulation();
}

void UnitTestModulation::runTestPolyModulation()
{
  // Create the mutex lock and the managers:
  juce::CriticalSection lock;
  jura::rsVoiceManager voiceMan;
  jura::ModulationManagerPoly modMan(&lock);
  modMan.setVoiceManager(&voiceMan);

  // Create and register modulation sources:
  static const int numSources = 5;
  TestPolyModulator* sources[numSources];
  for(int i = 0; i < numSources; i++) {
    sources[i] = new TestPolyModulator(&lock, nullptr, &modMan, &voiceMan);
    modMan.registerModulationSource(sources[i]); }

  // Create the target module and register its parameters as modulation targets:
  TestPolyAudioModule targetModule(&lock, nullptr, &modMan, &voiceMan);



  int dummy = 0;
}


