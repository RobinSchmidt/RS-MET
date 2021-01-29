#include "UnitTestsModulation.h"

// Define Some subclasses of AudioModule that are used in the tests. They have two parameters:
// Frequency and Amplitude which are also used as output signals for the per-sample callback, so we
// can figure out what their current values are.


// todo: a monophonic modulator

//=================================================================================================

/** A polyphonic modulator. */

class JUCE_API TestPolyModulator : public jura::ModulatorModulePoly
{

public:

  using jura::ModulatorModulePoly::ModulatorModulePoly;  // inherit constructors


  double renderModulation() override
  {
    return value;
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    return value;
  }

  // do we need to override the monophonic renderModulation too?

  double value = 1.0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyModulator)
};


//=================================================================================================

/** A monophonic audio module with parameters of class ModulatebleParameter. */

class JUCE_API TestMonoAudioModule : public jura::ModulatableAudioModule
{

public:

  TestMonoAudioModule(juce::CriticalSection* lockToUse, 
    jura::MetaParameterManager* metaManagerToUse = nullptr,
    jura::ModulationManager* modManagerToUse = nullptr) :
    ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse) 
  {
    createParameters();
  }

  void processStereoFrame(double* left, double* right) override
  {
    *left  = freq;
    *right = amp;
  }

  void setFrequency(double newFreq) 
  { 
    freq = newFreq; 
  }

  void setAmplitude(double newAmp)  
  { 
    amp  = newAmp; 
  }




private:

  void createParameters()
  {
    // ScopedLock scopedLock(*lock);
    typedef TestMonoAudioModule TMAM;
    typedef jura::ModulatableParameter Param;
    Param* p;

    p = new Param("Frequency", 20.0, 20000.0, 1000.0, Param::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallback<TMAM>(this, &TMAM::setFrequency);

    p = new Param("Amplitude", -1.0, +1.0, +1.0, Param::LINEAR);
    addObservedParameter(p);
    p->setValueChangeCallback<TMAM>(this, &TMAM::setAmplitude);
  }


  double freq = 0, amp = 0;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestMonoAudioModule)
};


//=================================================================================================

/** A monophonic audio module with parameters of class ModulatebleParameterPoly. */

class JUCE_API TestPolyAudioModule : public jura::AudioModulePoly
{
public:

  TestPolyAudioModule(juce::CriticalSection* lockToUse,
    jura::MetaParameterManager* metaManagerToUse = nullptr,
    jura::ModulationManager* modManagerToUse = nullptr)
    : AudioModulePoly(lockToUse, metaManagerToUse, modManagerToUse)
  {
    createParameters();
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

  void allocateVoiceResources() override
  {
    int numVoices = 0;
    jura::rsVoiceManager* vm = getVoiceManager();
    if(vm)
      numVoices = vm->getNumVoices();
    voiceFreqs.resize(numVoices);
    voiceAmps.resize(numVoices);
  }


private:


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

  //double monoFreq, monoAmp;  
  std::vector<double> voiceFreqs, voiceAmps;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyAudioModule)
};

//=================================================================================================

UnitTestModulation::UnitTestModulation() 
  : juce::UnitTest("Modulation", "Control"), modMan(&lock)
{
  // Set up the basic infrastructure:
  modMan.setVoiceManager(&voiceMan);
  voiceBuffer.resize(2*voiceMan.getMaxNumVoices());
  voiceMan.setVoiceSignalBuffer(&voiceBuffer[0]);
}

void UnitTestModulation::reset()
{
  modMan.deRegisterAllSources();
  modMan.deRegisterAllTargets();
  //modMan.removeAllConnections();  // should happen automatically in either of the 2 calls above

  voiceMan.reset();
}

void UnitTestModulation::runTest()
{
  UnitTest::beginTest("Modulation");
  runTestPolyToPoly();
  runTestPolyToMono();

}


void UnitTestModulation::runTestPolyToMono()
{
  reset();

  int iVal;
  double dVal1, dVal2;

  // Create and register modulation sources:
  TestPolyModulator mod1(&lock, nullptr, &modMan); mod1.setVoiceManager(&voiceMan);
  TestPolyModulator mod2(&lock, nullptr, &modMan); mod2.setVoiceManager(&voiceMan);
  TestPolyModulator mod3(&lock, nullptr, &modMan); mod3.setVoiceManager(&voiceMan);
  modMan.registerModulationSource(&mod1);
  modMan.registerModulationSource(&mod2);
  modMan.registerModulationSource(&mod3);
  iVal = (int) modMan.getAvailableModulationSources().size();
  expectEquals(iVal, 3, "Failed to register sources in modulation manager");
  // hmm...maybe we should use the velocity and freq/pitch as modulators

  // Create the target module and register its parameters as modulation targets:
  TestMonoAudioModule gen(&lock, nullptr, &modMan); // "gen" for "generator"
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");

  // Create and add connections:
  const std::vector<jura::ModulationTarget*>& targets = modMan.getAvailableModulationTargets();
  double depth1 = 1.0;
  double depth2 = 0.5;
  addConnection(&mod1, targets[0], depth1);
  addConnection(&mod2, targets[1], depth2);
  iVal = (int) modMan.getModulationConnections().size();
  expectEquals(iVal, 2, "Failed add connection to modulation manager");

  // Let the target module process an audio frame - this should produce unmodulated values:
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0); 
  expectEquals(dVal2,    1.0);

  // Set up modulator values and produce a frame - this should produce modulated values:
  mod1.value = 3.0;
  mod2.value = 5.0;
  modMan.applyModulationsNoLock();
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + depth1 * mod1.value); 
  expectEquals(dVal2,    1.0 + depth2 * mod2.value);
  // this fails! 
  // -TestPolyModulator::renderModulation is called 3 times (-> ok)
  // -setFreqiuencvy is called, but with 1000 - 
  //  affectedTargets[0].modulatedValue is still 1000



  // Trigger a note-on on the voice manager (the gen module itself does not have a midi event 
  // handler because it's just a ModulatableAudioModule with no midi in), call applyModulations 
  // and produce a sample again - this should produce a value with the modulations applied:
  using Msg = juce::MidiMessage;
  int    key1 = 69;
  float  vel  = 0.5f;
  double velQ = voiceMan.quantize7BitUnsigned(vel);       // midi roundtrip quantizes velocity
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel));


  // todo: try changing the parameter(s)


  int dummy = 0;
}


void UnitTestModulation::runTestPolyToPoly()
{
  reset();

  voiceMan.setKillMode(jura::rsVoiceManager::KillMode::immediately);

  // Temporary values for the tests:
  int iVal;
  double dVal1, dVal2;

  // Create modulation sources:
  jura::rsConstantOneModulatorModulePoly  constant( &lock, nullptr, &modMan);
  jura::rsNotePitchModulatorModulePoly    notePitch(&lock, nullptr, &modMan);
  jura::rsNoteVelocityModulatorModulePoly noteVel(  &lock, nullptr, &modMan);
  constant.setVoiceManager(&voiceMan);
  notePitch.setVoiceManager(&voiceMan);
  noteVel.setVoiceManager(&voiceMan);

  // Register modulation sources:
  modMan.registerModulationSource(&constant);
  modMan.registerModulationSource(&notePitch);
  modMan.registerModulationSource(&noteVel);
  iVal = (int) modMan.getAvailableModulationSources().size();
  expectEquals(iVal, 3, "Failed to register sources in modulation manager");

  // Create the target module and register its parameters as modulation targets:
  TestPolyAudioModule targetModule(&lock, nullptr, &modMan);
  targetModule.setVoiceManager(&voiceMan);
  //targetModule.init();
  targetModule.setVoiceSignalBuffer(&voiceBuffer[0]);
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");

  // Create and add connections:
  const std::vector<jura::ModulationTarget*>& targets = modMan.getAvailableModulationTargets();
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


