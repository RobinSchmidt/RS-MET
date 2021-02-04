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
  // Do we need to override the monophonic renderModulation too? What happens if we don't? Maybe 
  // make a TestPolyModulator2 that overrides only the polyphonic version and test, how it behaves.
  // I think, when it is connected to a monophonic parameter, we should expect it to produce a zero 
  // signal in cases when there's no active voice and the modulation signal of the newest active 
  // voice when at least one voice is active.

  double renderVoiceModulation(int voiceIndex) override
  {
    numPolyRenders++;
    return value;
  }

  void allocateVoiceModResources() override {}



  double value = 1.0;

  int numMonoRenders = 0;
  int numPolyRenders = 0;  // maybe rename to numVoiceRenders

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
    const jura::rsVoiceManager* vm = getVoiceManager();
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
  //UnitTest::beginTest("Modulation");
  voiceMan.setKillMode(jura::rsVoiceManager::KillMode::immediately);
  voiceMan.setMaxNumVoices(16);
  voiceMan.setNumVoices(8);

  runTestMonoToMono();
  runTestMonoToPoly();
  runTestPolyToMono();
  runTestPolyToPoly();

  // ToDo: make a test that combines all sorts of modulators and targets - both mono/poly and wire
  // everything to everything

  // When adding new tests, keep in mind:
  // The modulation depths in the connections should be between 0 and 1 because they are clipped at
  // these values by default because depth is actually a Parameter object that has a default range
  // of 0..1.
}

void UnitTestModulation::runTestMonoToMono()
{
  UnitTest::beginTest("Modulation: mono to mono");
  reset();

  int iVal;
  double dVal1, dVal2;  // values

  // Create and register modulation sources:
  TestMonoModulator mod1(&lock, nullptr, &modMan);
  TestMonoModulator mod2(&lock, nullptr, &modMan);
  modMan.registerModulationSource(&mod1);
  modMan.registerModulationSource(&mod2);
  iVal = (int) modMan.getAvailableModulationSources().size();
  expectEquals(iVal, 2, "Failed to register sources in modulation manager");

  // Create the target module and register its parameters as modulation targets:
  TestMonoAudioModule gen(&lock, nullptr, &modMan); // "gen" for "generator"
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");

  // Create and add connections:
  const std::vector<jura::ModulationTarget*>& targets = modMan.getAvailableModulationTargets();
  double depth1 = 0.75;
  double depth2 = 0.5;
  addConnection(&mod1, targets[0], depth1);
  addConnection(&mod2, targets[1], depth2);
  iVal = (int) modMan.getModulationConnections().size();
  expectEquals(iVal, 2, "Failed add connection to modulation manager");

  // Let the target module process an audio frame. This should produce the unmodulated values.
  // Monophonic sources always potentially produce output, regardless whether a note is active or 
  // not. Of course, a monophonic synth can choose to produce no output unless a note is active 
  // but that has nothing to do with how the framework should behave. Monophonic modulators that 
  // are subclasses of ModulatorModuleMono actually can and typically do respond to midi messages
  // even though there is no voice manager involved, because they are already inheriting 
  // noteOn/Off callbacks from AudioModuleWithMidiIn which is a baseclass of ModulatorModuleMono.
  // This makes sense because modulators are typically things that we want to trigger via midi even
  // in absence of a polyphony framework. 
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0); 
  expectEquals(dVal2,    1.0);

  // Apply modulations and render a sample again. This time, we should get modulated outputs:
  
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 1);
  expectEquals(mod2.numMonoRenders, 1);
  expectEquals(gen.numSetFreqCalls, 1);
  expectEquals(gen.numSetAmpCalls,  1);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectEquals(dVal1, 1000.0 + depth1 * mod1.value); 
  expectEquals(dVal2,    1.0 + depth2 * mod2.value);

  // Test changing the frequency parameter:
  jura::Parameter* freqParam = gen.getParameterByName("Frequency");
  freqParam->setValue(500.0, false, false);
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 2);
  expectEquals(mod2.numMonoRenders, 2);
  expectEquals(gen.numSetFreqCalls, 2);
  expectEquals(gen.numSetAmpCalls,  2);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 500.0 + depth1 * mod1.value, 1.e-13); 
  expectEquals(dVal2,   1.0 + depth2 * mod2.value);
  // We get roundoff error here: the unmodulated value after calling setValue is something like 
  // 499.9999.... It's because setValue lets the value makes a roundtrip through the normalized 
  // value. In MetaControlledParameter::setValue, we do:
  //   double y = mapper->unmap(newValue);
  //   double x = metaMapper.unmap(y);
  // This is somehow to keep the value in sync with the host automation, also taking into account
  // the user-defined mapping between normalized meta-parameters. I don't know exactly anymore
  // why that is necessarry and how it works but i guess, it's ok. And it's definitely not the
  // fault of the modulation system. It has to do with the automation/metaparameter system.
  // The roundoff error does not seem to happen for the polyphonic parameter. That seems to be
  // because we use other modulations there - the unmodulated value is still wrong but the 
  // modulated value apparently happens to come out right again. 

  // Test, if the engine is also updated when there are no connections:
  modMan.removeAllConnections();
  expectEquals(gen.numSetAmpCalls, 3);      // removing a connection will also trigger a callback, 
  expectEquals(gen.numSetAmpCalls, 3);      // so we do not get stuck at the modulated value
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 500.0 , 1.e-13); 
  expectEquals(dVal2,   1.0);
  freqParam->setValue(1000.0, false, true); // simulates a change from the GUI
  expectEquals(gen.numSetFreqCalls, 4);
  modMan.applyModulationsNoLock();          // should have no effect on numSetFreq/AmpCalls
  expectEquals(mod1.numMonoRenders, 3);
  expectEquals(mod2.numMonoRenders, 3);
  expectEquals(gen.numSetFreqCalls, 4);
  expectEquals(gen.numSetAmpCalls,  3);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 1000.0, 1.e-12); 
  expectEquals(dVal2,    1.0);
  // ToDo: maybe provide facilities that suppress the rendering of the disconnected modulation 
  // sources. In ModulationManager(Poly), we always loop over all available sources to let them 
  // render their outputs, so the mod1/2.numMonoRenders are here actually expected to go up by 
  // one. But these calls could be optimized away in a way similar as we do with the 
  // affctedTargets array. In ToolChain, we would not expect these useless render calls in this 
  // configuration because it skips the call to applyModulationsNoLock entirely, when there are no 
  // connections at all. But here, we call it. But we really should optimize away also the 
  // rendering calls to individual disconnected sources.


  // What about smoothing? I think, the smoothing callback should set the unmodulatedValue. ToDo:
  // test smoothing here, too. Set up a smoothing manager "smoothMan" and repeat some of the tests 
  // above

  int dummy = 0;
}

void UnitTestModulation::runTestMonoToPoly()
{
  UnitTest::beginTest("Modulation: mono to poly");
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
  int   key2 = 20,   key3 = 30,   key4 = 40,   key5 = 50;
  float vel2 = 0.2f, vel3 = 0.3f, vel4 = 0.4f, vel5 = 0.5f;
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
  jura::Parameter* freqParam = gen.getParameterByName("Frequency");
  freqParam->setValue(500.0, false, false);
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders, 6);
  expectEquals(mod2.numMonoRenders, 6);
  expectEquals(gen.numSetFreqCalls, 6);
  expectEquals(gen.numSetAmpCalls,  6);
  gen.processStereoFrame(&dVal1, &dVal2);
  dTgt1  = 500.0 + depth1 * mod1.value + depthKey * key1;
  dTgt2  =   1.0 + depth2 * mod2.value + depthVel * vel1Q;
  dTgt1 += 500.0 + depth1 * mod1.value + depthKey * key4;
  dTgt2 +=   1.0 + depth2 * mod2.value + depthVel * vel4Q;
  expectEquals(dVal1, dTgt1);
  expectEquals(dVal2, dTgt2);

  // Test, if the engine is also updated when there are no connections:
  modMan.removeAllConnections();
  expectEquals(gen.numSetFreqCalls, 8);              // should have increased by 2 because 
  expectEquals(gen.numSetAmpCalls,  8);              // 2 voices are still playing
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 1000.0, 1.e-12);  // 2 notes, no mods, freq is 500
  expectEquals(dVal2,    2.0);
  freqParam->setValue(1000.0, false, true);          // simulates setting from gui
  expectEquals(gen.numSetFreqCalls, 10);             // gui should trigger callback
  modMan.applyModulationsNoLock();
  expectEquals(mod1.numMonoRenders,  7);
  expectEquals(mod2.numMonoRenders,  7);
  expectEquals(gen.numSetFreqCalls, 10);             // updateModulations should not trigger any
  expectEquals(gen.numSetAmpCalls,   8);             // more callbacks
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 2000.0, 1.e-12);  // 2 notes, no mods, freq is 1000
  expectEquals(dVal2, 2.0);

  // OK - i think the problem is that when a parameter has no modulators connected to it, it's per 
  // voice versions of the value are not updated because it's not in the affectedTargets array 
  // anymore. What should we do about this? Maybe when the last modulator is disconnected from a 
  // poly-target, we should call the per-voice callback for all available voices with the 
  // unmodulated value? but that would be expensive. Or maybe it's sufficient to call it for the 
  // active voices? But then, what when later an additional voice is triggered? Will it then have 
  // the wrong setting? But maybe it's practical to call the callback for available voices when the
  // called back function just sets a "dirty" flag and defers coefficient recomputations to the 
  // next getSample call which is a good pattern anyway because it consolidates multiple setter 
  // calls into a single update computation and makes the updates thread-safe.

  // So, there's a general problem with polyphonic parameters that are not connected to any 
  // modulator: we have no infrastructure in place that calls its per-voice callbacks in this case.
  // The valueChangeCallbackPoly gets called only in doVoiceModulationUpdate which is invoked from
  // applyModulationsNoLock, iff the parameter is connected to at least one modulation source. When 
  // it's disconnected, we somehow need to directly call it from elsewhere...maybe noteOn of the
  // AudioModulePoly class and/or setValue of ModulatableParameterPoly class. But that may, again,
  // lead to redundant callbacks in the case of connected parameters. But noteOn events are rare
  // so a redundant callback seems acceptable. setParameter occurs from the gui and from automation
  // redundant callbacks from the gui may also be rare enough to be accpetable but automation 
  // events could potentially be dense. maybe the best thing to do is to keep track of the number
  // of connected sources and invoke the callback from setValue only if the parameter is 
  // disconnected. maybe the same should be done also for monophonic parameters?

  // ....OK - this seems to be solved - the next probelm is when a new note is triggered:

  // At the moment key1 and key4 are still held. We now trigger key2 again. In this case, gen must 
  // also receive the noteOn, so we call its noteOn as well. Formerly, with modulation connections 
  // in place this wasn't necessary, because the parameter callbacks were invoked by the modulation
  // system anyway. But without connections, that doesn't happen, so we must call noteOn. In the 
  // context of ToolChain, we can assume that this is always called by ToolChains midi-handler, 
  // here we need to do it ourselves:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key2, 1.f));
  gen.noteOn(key2, 127);
  expectEquals(gen.numSetFreqCalls, 11);    // gen.noteOn should not trigger a callback
  expectEquals(gen.numSetAmpCalls,   9);    // for the newly allocated voice
  modMan.applyModulationsNoLock();          // should trigger no callbacks
  expectEquals(gen.numSetFreqCalls, 11);
  expectEquals(gen.numSetAmpCalls,   9);
  expectEquals(mod1.numMonoRenders,  8);    // ..but it does (still) trigger rendering..
  expectEquals(mod2.numMonoRenders,  8);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 3000.0, 1.e-12);  // 3 notes, no mods, freq is 1000

  // It should also work, when we do it like this:
  int voice = voiceMan.noteOnReturnVoice(key4, 127);
  gen.noteOnForVoice(key4, 127, voice);
  expectEquals(gen.numSetFreqCalls, 12);
  expectEquals(gen.numSetAmpCalls,  10);
  gen.processStereoFrame(&dVal1, &dVal2);
  expectWithinAbsoluteError(dVal1, 4000.0, 1.e-12);  // 4 notes, no mods, freq is 1000




  int dummy = 0;
}

void UnitTestModulation::runTestPolyToMono()
{
  UnitTest::beginTest("Modulation: poly to mono");
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
  UnitTest::beginTest("Modulation: poly to poly");
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
  TestPolyAudioModule gen(&lock, nullptr, &modMan);
  gen.setVoiceManager(&voiceMan);
  gen.setVoiceSignalBuffer(&voiceBuffer[0]);
  iVal = (int) modMan.getAvailableModulationTargets().size();
  expectEquals(iVal, 2, "Failed to register parameters in modulation manager");

  // Create and add connections:
  const std::vector<jura::ModulationTarget*>& targets = modMan.getAvailableModulationTargets();
  double depthVel = 0.5;
  addConnection(&notePitch, targets[0], 1.0);
  addConnection(&noteVel,   targets[1], depthVel);
  iVal = (int) modMan.getModulationConnections().size();
  expectEquals(iVal, 2, "Failed add connection to modulation manager");


  // Let the target module process an audio frame:
  gen.processStereoFrameVoice(&dVal1, &dVal2, 0);
  expectEquals(dVal1, 0.0);
  expectEquals(dVal2, 0.0);

  // Define a convenience function for applying the modulations and generating a sample frame:
  auto process = [&]()
  {
    modMan.applyModulationsNoLock();
    gen.processStereoFrame(&dVal1, &dVal2);
  };


  // Trigger a note-on - this should cause the 0th voice to start playing
  using Msg = juce::MidiMessage;
  int    key1 = 69;
  float  vel  = 0.5f;
  double velR = round(vel * 127.0) / 127.0;  // midi-byte roundtrip messes the velocity up
  gen.handleMidiMessage(Msg::noteOn(1, key1, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 1);
  process();
  expectEquals(dVal1, 1000.0 + key1    );       // default freq of 1kHz plus the note-number
  expectEquals(dVal2,    1.0 + depthVel*velR);  // default amp of 1 plus depth*vel2

  // Add a 2nd connection to the constant 1 the frequency parameter. It should immediately take 
  // effect:
  addConnection(&constant, targets[0], 1.0);
  process();
  expectEquals(dVal1, 1000.0 + key1 + 1.0);      // + 1.0 is the addition from the 2nd connection
  expectEquals(dVal2,    1.0 + depthVel*velR  ); // this should be the same value as before

  // Trigger a second note:
  int key2 = 50;
  gen.handleMidiMessage(Msg::noteOn(1, key2, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  process();
  expectEquals(dVal1, 2000.0 + key1 + key2 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + depthVel*velR)  );

  // ...and a third:
  int key3 = 100;
  gen.handleMidiMessage(Msg::noteOn(1, key3, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 3);
  process();
  expectEquals(dVal1, 3000.0 + key1 + key2 + key3 + 3.0);
  expectEquals(dVal2,  3*( 1.0 + depthVel*velR)  );

  // Release the first note:
  gen.handleMidiMessage(Msg::noteOn(1, key1, 0.f));  
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  process();
  expectEquals(dVal1, 2000.0 + key2 + key3 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + depthVel*velR)  );

  // Trigger another note:
  int key4 = 20;
  gen.handleMidiMessage(Msg::noteOn(1, key4, vel));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 3);
  process();
  expectEquals(dVal1, 3000.0 + key2 + key3 + key4 + 3.0);
  expectEquals(dVal2,  3*( 1.0 + depthVel*velR)  );

  // Release 3rd note:
  gen.handleMidiMessage(Msg::noteOn(1, key3, 0.f));
  iVal = voiceMan.getNumActiveVoices();
  expectEquals(iVal, 2);
  process();
  expectEquals(dVal1, 2000.0 + key2 + key4 + 2.0);
  expectEquals(dVal2,  2*( 1.0 + depthVel*velR)  );

  // Remove the connection between the constant 1 and the frequency parameter:
  modMan.removeConnectionsWith(&constant);
  process();
  expectEquals(dVal1, 2000.0 + key2 + key4);

  // Change the value of the frequency parameter:
  gen.getParameterByName("Frequency")->setValue(500.0, false, false);
  process();
  expectEquals(dVal1, 1000.0 + key2 + key4);

  // Switch the generator into monophonic mode, trigger another note with key3 (2 voices are still
  // active, playing key2 and key4, so this note will use the 3rd active voice). We expect that an
  // output corresponding to key3 only:
  gen.setMonophonic(true);
  gen.handleMidiMessage(Msg::noteOn(1, key3, vel));
  process();
  expectEquals(dVal1, 500.0 + key3);
  expectEquals(dVal2,   1.0 + depthVel*velR);

  // Switch back to polyphonic mode again:
  gen.setMonophonic(false);
  process();
  expectEquals(dVal1, 1500.0 + key2 + key4 + key3);
  expectEquals(dVal2,  3*(1.0 + depthVel*velR));

  // Connect the constant 1 modulator to the frequency again:
  addConnection(&constant, targets[0], 1.0);
  process();
  expectEquals(dVal1, 1500.0 + key2 + key4 + key3 + 3.0);
  expectEquals(dVal2,  3*(1.0 + depthVel*velR));

  // Switch the constant modulator into mono mode. This should make no difference:
  constant.setMonophonic(true);
  process();
  expectEquals(dVal1, 1500.0 + key2 + key4 + key3 + 3.0);

  // Switch the notePitch modulator into mono mode, too. This means that instead of seeing
  // key2 + key4 + key3 as modulation signal, the module should see 3*key3:
  notePitch.setMonophonic(true);
  process();
  expectEquals(dVal1, 1500.0 + 3*key3 + 3.0);

  // Now switch the generator also into mono-mode. We expect an output corresponding only to the
  // latest actiive note:
  gen.setMonophonic(true);
  process();
  expectEquals(dVal1, 500.0 + key3 + 1.0);

  // ToDo: check also the number of randering and callback calls in each test

  // What if we trigger a note again before a corresponding note-off was received? i think, best
  // would be to just do nothing, if the note is being held and just turn it back on when it was
  // releasing ....but that's a matter of voice management

  int dummy = 0;
}


