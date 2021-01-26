#include "UnitTestsModulation.h"

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
public:

  using jura::AudioModulePoly::AudioModulePoly;  // inherit constructors

  void init()
  {
    jura::rsVoiceManager* vm = getVoiceManager();
    int numVoices = vm->getNumVoices();
    voiceFreqs.resize(numVoices);
    voiceAmps.resize(numVoices);
    createParameters();
  }

  void createParameters()
  {
   // ScopedLock scopedLock(*lock);
    typedef jura::ModulatableParameterPoly Param;
    Param* p;

    p = new Param("Frequency", 20.0, 20000.0, 1000.0, Param::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallbackPoly([this](double v, int i) { setFrequency(v, i); });

    p = new Param("Amplitude", -1.0, +1.0, +1.0, Param::LINEAR);
    addObservedParameter(p);
    p->setValueChangeCallbackPoly([this](double v, int i) { setAmplitude(v, i); });
  }


  void processStereoFrameVoice(double* left, double* right, int voice) override
  {
    int dummy = 0;
  }


  void setFrequency(double newFreq, int voice) 
  { 
    voiceFreqs[voice] = newFreq;
  }

  void setAmplitude(double newAmp, int voice) 
  { 
    voiceAmps[voice] = newAmp;
  }

  //double monoFreq, monoAmp;  
  std::vector<double> voiceFreqs, voiceAmps;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyAudioModule)
};

//=================================================================================================

UnitTestModulation::UnitTestModulation() : juce::UnitTest("Modulation", "Control")
{

}

void UnitTestModulation::runTest()
{
  UnitTest::beginTest("Modulation");
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
  /*
  static const int numSources = 5;
  TestPolyModulator* sources[numSources];
  for(int i = 0; i < numSources; i++) {
    sources[i] = new TestPolyModulator(&lock, nullptr, &modMan, &voiceMan);
    modMan.registerModulationSource(sources[i]); }
  */

  int iVal;

  // Create modulation sources:
  jura::rsConstantOneModulatorModulePoly  constant( &lock, nullptr, &modMan, &voiceMan);
  jura::rsNotePitchModulatorModulePoly    notePitch(&lock, nullptr, &modMan, &voiceMan);
  jura::rsNoteVelocityModulatorModulePoly noteVel(  &lock, nullptr, &modMan, &voiceMan);

  // Register modulation sources:
  modMan.registerModulationSource(&constant);
  modMan.registerModulationSource(&notePitch);
  modMan.registerModulationSource(&noteVel);
  iVal = (int) modMan.getAvailableModulationSources().size();
  expectEquals(iVal, 3, "Failed to register sources in modulation manager");

  // Create the target module and register its parameters as modulation targets:
  TestPolyAudioModule targetModule(&lock, nullptr, &modMan, &voiceMan);
  targetModule.init();
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");


  int dummy = 0;
}


