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
    *left  = voiceFreqs[voice];
    *right = voiceAmps[voice];
  }

  void handleMidiMessage(MidiMessage message) override
  {
    voiceManager->handleMidiMessage(message);
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
  // Create the basic infrastructure:
  juce::CriticalSection lock;
  jura::rsVoiceManager voiceMan;
  jura::ModulationManagerPoly modMan(&lock);
  modMan.setVoiceManager(&voiceMan);
  std::vector<double> voiceBuffer;
  voiceBuffer.resize(2*voiceMan.getMaxNumVoices());
  voiceMan.setVoiceSignalBuffer(&voiceBuffer[0]);
  voiceMan.setKillMode(jura::rsVoiceManager::KillMode::immediately);

  // Temporary values for the tests:
  int iVal;
  double dVal1, dVal2;

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
  targetModule.setVoiceSignalBuffer(&voiceBuffer[0]);
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");

  // Create and add connections:
  const std::vector<jura::ModulationTarget*>& targets = modMan.getAvailableModulationTargets();
  auto addConnection = [&](jura::ModulationSource* source, jura::ModulationTarget* target, 
    double depth)  // convenience function
  {
    jura::ModulationConnection* c = new jura::ModulationConnection(source, target, nullptr);
    c->setDepth(depth);
    modMan.addConnection(c);
  };
  addConnection(&notePitch, targets[0], 1.0);
  addConnection(&noteVel,   targets[1], 0.5);
  iVal = (int) modMan.getModulationConnections().size();
  expectEquals(iVal, 2, "Failed add connection to modulation manager");


  // Let the target module process an audio frame:
  targetModule.processStereoFrameVoice(&dVal1, &dVal2, 0);
  expectEquals(dVal1, 0.0);
  expectEquals(dVal2, 0.0);

  // Trigger a note-on - this should cause the 0th voice to start playing
  using Msg = juce::MidiMessage;
  int    key1 = 69;
  float  vel  = 0.5f;
  double velR = round(vel * 127.0) / 127.0;  // midi-byte roundtrip messes the velocity up
  targetModule.handleMidiMessage(Msg::noteOn(1, key1, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 1);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + key1    );  // default freq of 1kHz plus the note-number
  expectEquals(dVal2,    1.0 + 0.5*velR);  // default amp of 1 plus depth*vel2

  // Add a 2nd connection to the 1st parameter:
  addConnection(&constant, targets[0], 1.0);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + key1 + 1.0); // the + 1.0 is the addition from the 2nd connection
  expectEquals(dVal2,    1.0 + 0.5*velR  ); // this should be the same value as before

  // Trigger a second note:
  int key2 = 50;
  targetModule.handleMidiMessage(Msg::noteOn(1, key2, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 2000.0 + key1 + key2 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + 0.5*velR)  );

  // ...and a third:
  int key3 = 100;
  targetModule.handleMidiMessage(Msg::noteOn(1, key3, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 3);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 3000.0 + key1 + key2 + key3 + 3.0);
  expectEquals(dVal2,  3*( 1.0 + 0.5*velR)  );

  // Release the first note:
  targetModule.handleMidiMessage(Msg::noteOn(1, key1, 0.f));  
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 2000.0 + key2 + key3 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + 0.5*velR)  );

  // Trigger another note:
  int key4 = 20;
  targetModule.handleMidiMessage(Msg::noteOn(1, key4, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 3);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 3000.0 + key2 + key3 + key4 + 3.0);
  expectEquals(dVal2,  3*( 1.0 + 0.5*velR)  );

  // Release 3rd note:
  targetModule.handleMidiMessage(Msg::noteOn(1, key3, 0.f));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 2000.0 + key2 + key4 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + 0.5*velR)  );

  // Remove the connection between the constant 1 and the frequency parameter:
  modMan.removeConnectionsWith(&constant);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 2000.0 + key2 + key4);

  // Change the value of the frequency parameter:
  targetModule.getParameterByName("Frequency")->setValue(500.0, false, true);
  modMan.applyModulationsNoLock();
  targetModule.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + key2 + key4);





  // What if we trigger a note again before a corresponding note-off was received? i think, best
  // would be to just do nothing, if the note is being held and just turn it back on when it was
  // releasing ....but that's a matter of voice management




  int dummy = 0;
}

