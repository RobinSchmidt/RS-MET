#include "UnitTestsModulation.h"

// Define Some subclasses of AudioModule that are used in the tests. They have two parameters:
// Frequency and Amplitude which are also used as output signals for the per-sample callback, so we
// can figure out what their current values are.


// todo: a monophonic modulator


//=================================================================================================

/** A monophonic modulator. */

class JUCE_API TestMonoModulator : public jura::ModulatorModuleMono
{

public:

  using jura::ModulatorModuleMono::ModulatorModuleMono;  // inherit constructors


  double renderModulation() override
  {
    numMonoRenders++;
    return value;
  }


  double value = 1.0;
  int numMonoRenders = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestMonoModulator)
};

//=================================================================================================

/** A polyphonic modulator. */

class JUCE_API TestPolyModulator : public jura::ModulatorModulePoly
{

public:

  using jura::ModulatorModulePoly::ModulatorModulePoly;  // inherit constructors


  double renderModulation() override
  {
    numMonoRenders++;
    return value;
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    numPolyRenders++;
    return value;
  }

  // do we need to override the monophonic renderModulation too?

  double value = 1.0;

  int numMonoRenders = 0;
  int numPolyRenders = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPolyModulator)
};
// todo: count, how many times the render callbacks have been called and check in the test, if 
// that matches the expectation


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
    numSetFreqCalls++;
  }

  void setAmplitude(double newAmp)  
  { 
    amp  = newAmp;
    numSetAmpCalls++;
  }



  int numSetFreqCalls = 0, numSetAmpCalls = 0;

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

    // We don't want to count the calls that happen during construction. The counters get 
    // incremented during setValueChangeCallback because this setter actually also invokes the 
    // callback once, to make sure, everything is in sync after construction. So we decrement them
    // again:
    numSetFreqCalls--;
    numSetAmpCalls--;
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
    numSetFreqCalls++;
  }

  void setAmplitude(double newAmp, int voice) 
  { 
    voiceAmps[voice] = newAmp;
    numSetAmpCalls++;
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

  int numSetFreqCalls = 0, numSetAmpCalls = 0;


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

    //numSetFreqCalls--;
    //numSetAmpCalls--;
    // nope - no need to decrement here - the poly callback is not invoked when its wired up 
    // because that would make no sense - we don't even know for which voice we should call it
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
  runTestMonoToPoly();
  runTestPolyToMono();
  runTestPolyToPoly();

  // todo: make a test that combines all sorts of modulators and targets - both mono/poly and wire
  // everything to everything

  // The modulation depths in the connections should be between 0 and 1 because they are clipped at
  // these values by default because depth is actually a Parameter object that has a default range
  // of 0..1.
}


void UnitTestModulation::runTestMonoToPoly()
{
  reset();

  int iVal;
  double dVal1, dVal2;  // values
  double dTgt1, dTgt2;  // targets

  // Create and register modulation sources:
  TestMonoModulator mod1(&lock, nullptr, &modMan);
  TestMonoModulator mod2(&lock, nullptr, &modMan);
  modMan.registerModulationSource(&mod1);
  modMan.registerModulationSource(&mod2);
  iVal = (int) modMan.getAvailableModulationSources().size();
  expectEquals(iVal, 2, "Failed to register sources in modulation manager");

  // Create the target module and register its parameters as modulation targets:
  TestPolyAudioModule gen(&lock, nullptr, &modMan);
  gen.setVoiceManager(&voiceMan);
  gen.setVoiceSignalBuffer(&voiceBuffer[0]);
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

  // Let the target module process an audio frame. This should produce zeros because we have not 
  // yet triggered a note. Polyphonic sources produce nonzero output only when voices are active 
  // (by default, but subclasses can override this behavior):
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 0.0); 
  expectEquals(dVal2, 0.0);

  // Now do the same thing but call applyModulationsNoLock before. We still expect zero outputs and
  // we also don't expect the setFrequency/Amplitude modulation callback to get called because when
  // the output is zero anyway it would be pointless to modulate any parameters. However, we do 
  // render the modulation signals...because they could be wired to monophonic targets as well..
  // ...hmm...does that make sense?
  mod1.value = 3.0;
  mod2.value = 5.0;
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 1);
  expectEquals(mod2.numMonoRenders, 1);
  expectEquals(gen.numSetFreqCalls, 0);
  expectEquals(gen.numSetAmpCalls,  0);
  expectEquals(dVal1, 0.0); 
  expectEquals(dVal2, 0.0);

  // Trigger a note and do the same. Now we actually expect the setFreq callbacks to be called and 
  // a nonzero output to be produced:
  using Msg = juce::MidiMessage;
  int    key1  = 69;
  float  vel1  = 1.0f;
  double vel1Q = voiceMan.quantize7BitUnsigned(vel1);       // midi roundtrip quantizes velocity
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel1));
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 2);
  expectEquals(mod2.numMonoRenders, 2);
  expectEquals(gen.numSetFreqCalls, 1);
  expectEquals(gen.numSetAmpCalls,  1);
  gen.processStereoFrame(&dVal1, &dVal2);
  dTgt1 = 1000.0 + depth1 * mod1.value;
  dTgt2 =    1.0 + depth2 * mod2.value;
  expectEquals(dVal1, 1000.0 + depth1 * mod1.value); 
  expectEquals(dVal2,    1.0 + depth2 * mod2.value);

  // Turn the nore off again. After noteOff, we should again get a zero output and setFreq should
  // not get called:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, 0.f));
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 3);
  expectEquals(mod2.numMonoRenders, 3);
  expectEquals(gen.numSetFreqCalls, 1);
  expectEquals(gen.numSetAmpCalls,  1);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 0.0); 
  expectEquals(dVal2, 0.0);

  // Create and register modulators for note-key and note-vel and connect them to the parameters:
  jura::rsNotePitchModulatorModulePoly    notePitch(&lock, nullptr, &modMan);
  jura::rsNoteVelocityModulatorModulePoly noteVel(  &lock, nullptr, &modMan);
  notePitch.setVoiceManager(&voiceMan);
  noteVel.setVoiceManager(&voiceMan);
  modMan.registerModulationSource(&notePitch);
  modMan.registerModulationSource(&noteVel);
  double depthKey = 0.25;
  double depthVel = 0.75;
  addConnection(&notePitch, targets[0], depthKey);
  addConnection(&noteVel,   targets[1], depthVel);
  iVal = (int) modMan.getModulationConnections().size();
  expectEquals(iVal, 4, "Failed add connection to modulation manager");

  // Trigger the note again and check, if now the additional note-based modulators are applied:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel1));
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 4);
  expectEquals(mod2.numMonoRenders, 4);
  expectEquals(gen.numSetFreqCalls, 2);
  expectEquals(gen.numSetAmpCalls,  2);
  gen.processStereoFrame(&dVal1, &dVal2);
  dTgt1 = 1000.0 + depth1 * mod1.value + depthKey * key1;
  dTgt2 =    1.0 + depth2 * mod2.value + depthVel * vel1Q;
  expectEquals(dVal1, dTgt1); 
  expectEquals(dVal2, dTgt2);

  // Trigger a 2nd, 3rd, 4th, 5th note, release the 2nd and 5th and 3rd. That means that key1 and 
  // key4 are still active after that sequence:
  int   key2 = 20,  key3 = 30,   key4 = 40,   key5 = 50;
  float vel2 = 0.2, vel3 = 0.3f, vel4 = 0.4f, vel5 = 0.5f;
  double vel4Q = voiceMan.quantize7BitUnsigned(vel4);  
  voiceMan.handleMidiMessage(Msg::noteOn(1, key2, vel2));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key3, vel3));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key4, vel4));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key5, vel5));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key2, 0.f));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key5, 0.f));
  voiceMan.handleMidiMessage(Msg::noteOn(1, key3, 0.f));
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 5);
  expectEquals(mod2.numMonoRenders, 5);
  expectEquals(gen.numSetFreqCalls, 4);  // incremented by two
  expectEquals(gen.numSetAmpCalls,  4); 
  gen.processStereoFrame(&dVal1, &dVal2);
  dTgt1  = 1000.0 + depth1 * mod1.value + depthKey * key1;  // contribution of voice with key1
  dTgt2  =    1.0 + depth2 * mod2.value + depthVel * vel1Q;
  dTgt1 += 1000.0 + depth1 * mod1.value + depthKey * key4;  // contribution of voice with key4
  dTgt2 +=    1.0 + depth2 * mod2.value + depthVel * vel4Q;
  expectEquals(dVal1, dTgt1); 
  expectEquals(dVal2, dTgt2);

  // Change parameter value - set freq from 1000 to 500:
  //gen.getParameterByName("Frequency")->setValue(500.0, false, false); // or false, true?
  gen.getParameterByName("Frequency")->setValue(500.0, false, true);
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 6);
  expectEquals(mod2.numMonoRenders, 6);
  expectEquals(gen.numSetFreqCalls, 6);
  expectEquals(gen.numSetAmpCalls,  6);
  dTgt1  = 500.0 + depth1 * mod1.value + depthKey * key1;
  dTgt2  =   1.0 + depth2 * mod2.value + depthVel * vel1Q;
  dTgt1 += 500.0 + depth1 * mod1.value + depthKey * key4;
  dTgt2 +=   1.0 + depth2 * mod2.value + depthVel * vel4Q;
  expectEquals(dVal1, dTgt1);
  expectEquals(dVal2, dTgt2);
  // This test still fails. Maybe the unmodulated value is not getting updated in setValue?

  int dummy = 0;
}


void UnitTestModulation::runTestPolyToMono()
{
  reset();

  int iVal;
  double dVal1, dVal2;

  // Create and register modulation sources:
  TestPolyModulator mod1(&lock, nullptr, &modMan); mod1.setVoiceManager(&voiceMan);
  TestPolyModulator mod2(&lock, nullptr, &modMan); mod2.setVoiceManager(&voiceMan);
  TestPolyModulator mod3(&lock, nullptr, &modMan); mod3.setVoiceManager(&voiceMan); // not yet used
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
  expectEquals(mod1.numMonoRenders, 1);
  expectEquals(mod1.numPolyRenders, 0);
  expectEquals(mod2.numMonoRenders, 1);
  expectEquals(mod2.numPolyRenders, 0);
  expectEquals(gen.numSetFreqCalls, 1);
  expectEquals(gen.numSetAmpCalls,  1);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + depth1 * mod1.value); 
  expectEquals(dVal2,    1.0 + depth2 * mod2.value);

  // Trigger a note-on on the voice manager (the gen module itself does not have a midi event 
  // handler because it's just a ModulatableAudioModule with no midi in), call applyModulations 
  // and produce a sample again - this should produce a value with the modulations applied:
  using Msg = juce::MidiMessage;
  int    key1 = 69;
  float  vel  = 1.0f;
  double velQ = voiceMan.quantize7BitUnsigned(vel);       // midi roundtrip quantizes velocity
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel));
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 1);
  expectEquals(mod1.numPolyRenders, 1);
  expectEquals(mod2.numMonoRenders, 1);
  expectEquals(mod2.numPolyRenders, 1);
  expectEquals(gen.numSetFreqCalls, 2); // fails! it got called 3 times. there's an extra
  expectEquals(gen.numSetAmpCalls,  2); // invocation of the callback somewhere!
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + depth1 * mod1.value); 
  expectEquals(dVal2,    1.0 + depth2 * mod2.value);

  // todo: 
  // -try changing the parameter(s) - for example, set the frequency to 500
  // -connect the a note-based modulator (pitch or vel) and check, if the applied note-based 
  //  modulation corresponds to the mostrecently triggered (via noteOn) note that has not yet 
  //  been released (via noteOff)

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


