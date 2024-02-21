#include "UnitTestsToolChain.h"
using namespace juce;
using namespace jura;


void UnitTestToolChain::runTest()
{
  runTestVoiceManager();
  runTestWaveOscillator();

  runTestQuadrifex();
  // This test is currently called last because it creates an actual jura::ToolChain object which 
  // in turn instantiates all modules once in populateModuleFactory - which is annyoing during 
  // debugging because certain initialization functions for ToolChain's built in AudioModules will
  // get called more often than one would expect in the tests, i.e. the breakpoints will trigger
  // more often than the actual running test justifies. Putting this test to the end fixes this.
  // ToDo: Factor out the creation of a ToolChain object from the Quadrifex test. It should be a 
  // test in its own right. Then, it shouldn't matter where we put the test for Quadrifex.
}


void UnitTestToolChain::runTestVoiceManager()
{
  beginTest("rsVoiceManager");

  rsVoiceManager voiceMan;
  std::vector<double> voiceBuffer(2*voiceMan.getMaxNumVoices());
  voiceMan.setVoiceSignalBuffer(&voiceBuffer[0]);

  // Test voice killing:
  using KM = rsVoiceManager::KillMode;
  double killThresh = 0.2;
  double killTime   = 0.1;
  double sampleRate = 100;
  voiceMan.setKillMode(KM::afterSilence);
  voiceMan.setSampleRate(sampleRate);
  voiceMan.setKillThreshold(killThresh);
  voiceMan.setKillTime(killTime);
  int killSamples = voiceMan.getKillTimeSamples();

  using Msg = juce::MidiMessage;
  int    key1  = 10;
  float  vel1  = 1.0f;
  double vel1Q = voiceMan.quantize7BitUnsigned(vel1);
  int    active, releasing;

  // Trigger a voice, then simulate a generated signal for the voice above the kill-threshhold by
  // assigning that value to the first two slots of the voiceBuffer:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel1));
  active    = voiceMan.getNumActiveVoices();
  releasing = voiceMan.getNumReleasingVoices();
  expectEquals(active,    1);
  expectEquals(releasing, 0);
  voiceBuffer[0] = voiceBuffer[1] = 1.1*killThresh;

  // Release the voice just triggered. It should nevertheless remain active (in release state),
  // because the signals in the corresponding voiceBuffer slots are above the kill threshold:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, 0.f));
  active    = voiceMan.getNumActiveVoices();
  releasing = (int)voiceMan.getNumReleasingVoices();
  expectEquals(active,    1);
  expectEquals(releasing, 1);
  for(int n = 0; n < 2*killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();               // should find no voice to kill...
    active    = voiceMan.getNumActiveVoices();          // ...so this should remain at 1
    releasing = voiceMan.getNumReleasingVoices();       // ...and this should be 1, too
    expectEquals(active,    1);
    expectEquals(releasing, 1); }

  // Now set the content of the voiceBuffer equal to the kill threshold. The voice should still 
  // remain active for "killSamples" samples and then turn off:
  voiceBuffer[0] = voiceBuffer[1] = killThresh;
  for(int n = 1; n < killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();               // should find no voice to kill...
    active    = voiceMan.getNumActiveVoices();          // ...so this should remain at 1
    releasing = voiceMan.getNumReleasingVoices();       // ...and this should be 1, too
    expectEquals(active,    1);
    expectEquals(releasing, 1); }
  voiceMan.findAndKillFinishedVoices();               // should find and kill the voice
  active    = voiceMan.getNumActiveVoices();          // ...so this should go to 0
  releasing = voiceMan.getNumReleasingVoices();       // ...and this too
  expectEquals(active,    0);
  expectEquals(releasing, 0);

  // Maybe later test other kill modes. Currently the only other mode is the trivial 
  // "immediately" mode, but if we later come up with more modes, test for them should go here. But
  // i really can't think of any other meaningful modes at the moment.

  // Test voice stealing:
  using SM = rsVoiceManager::StealMode;
  voiceMan.setNumVoices(4);
  voiceMan.setStealMode(SM::oldest);
  voiceMan.setKillMode(KM::immediately);
  voiceMan.reset();
  int key2 = 20, key3 = 30, key4 = 40, key5 = 50, key6 = 60;
  int voice;
  voice = voiceMan.noteOnReturnVoice(key1, 100); expectEquals(voice, 0);
  voice = voiceMan.noteOnReturnVoice(key2, 100); expectEquals(voice, 1);
  voice = voiceMan.noteOnReturnVoice(key3, 100); expectEquals(voice, 2);
  voice = voiceMan.noteOnReturnVoice(key4, 100); expectEquals(voice, 3);

  // Now all 4 voices are used up and stealing should take place:
  voice = voiceMan.noteOnReturnVoice(key5, 100); expectEquals(voice, 0);
  voice = voiceMan.noteOnReturnVoice(key6, 100); expectEquals(voice, 1); // fails!

  // Releae voice 3 and trigger another note - for the new note, the juts freed voice should be 
  // used again:
  voice = voiceMan.noteOffReturnVoice(key4);      expectEquals(voice, 3);
  voice = voiceMan.noteOnReturnVoice( key1, 100); expectEquals(voice, 3);


  voiceMan.setNumVoices(2);
  voiceMan.setStealMode(SM::oldest);
  voiceMan.setKillMode(KM::afterSilence);
  voiceMan.reset();

  // noteOn for key1:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    1);
  expectEquals(voiceMan.getNumReleasingVoices(), 0);

  // noteOff for key1 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    1);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOn for key2:
  voice = voiceMan.noteOnReturnVoice(key2, 100); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key2 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key2,   0); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key3, should steal voice 0, so it leaves release mode:
  voice = voiceMan.noteOnReturnVoice(key3, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key3 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key3,   0); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key1 (again) - should steal voice 1:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key1 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key1 (again) - should steal voice 0:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key1 - should release voice 0 and 1
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 0);
  //expectEquals(voice, voiceMan.manyVoices);  // fails - returns the last voice that was released
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // todo: check that both voices are released in a single noteOff in cases where they both play 
  // the same note

  // Trigger 3 different notes, then release them and loop through the killSamples, feeding a value
  // below the threshold, so it should eventually kill the voices:
  voiceMan.reset();
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key2, 100);
  voice = voiceMan.noteOnReturnVoice(key3, 100);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key2,   0);
  voice = voiceMan.noteOnReturnVoice(key3,   0);

  double* vb = &voiceBuffer[0];
  vb[0] = vb[1] = vb[2] = vb[3] = killThresh;
  bool ok = true;
  for(int n = 1; n < killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();
    active    = voiceMan.getNumActiveVoices();
    releasing = voiceMan.getNumReleasingVoices();
    ok &= active == 2 && releasing == 2; 
    jassert(ok);  }
  expect(ok);
  voiceMan.findAndKillFinishedVoices();         // should find and kill the 2 voices
  active    = voiceMan.getNumActiveVoices();    // ...so this should go to 0
  releasing = voiceMan.getNumReleasingVoices(); // ...and this too
  expectEquals(active,    0);
  expectEquals(releasing, 0);

  int dummy = 0;

  /*
  voiceMan.reset();
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  */


  // now the releasingVoices array is overfull (size == 3 where numVoices is only 2) - todo:
  // handle releasingVoices in the same way as activeVoices (resize to maxNumVoices and keep a
  // numReleasingVoices variable


  // ToDo: test voice stealing in the various modes, voice retriggering, etc.

}
// actually, that test belongs into the framework tests - maybe make files UnitTestsAudio/Midi
// and put this test there


void UnitTestToolChain::runTestQuadrifex()
{
  // This code parallels what is done in createPluginFilter() in the ToolChain project. This is the
  // factory function that creates the juce::AudioProcessor object in the juce framework:
  int numMetaParams = 10;  // Number of exposed automatable (meta) parameters
  jura::ToolChain*      dummy = nullptr;
  juce::AudioProcessor* proc  = createPluginWithMidi(dummy, numMetaParams);

  // Try to cast it into a jura::AudioPlugIn. This class is the glue between juce::AudioProcessor
  // and jura::AudioModule:
  jura::AudioPlugin* plug  = dynamic_cast<jura::AudioPlugin*>(proc);
  expect(plug != nullptr);

  // Extract the wrapped jura::AudioModule that is wrapped into the jura::AudioPlugin and cast it
  // into a pointer to jura::ToolChain:
  jura::AudioModule* mod   = plug->wrappedAudioModule;  // Maybe use a getter
  jura::ToolChain*   tlChn = dynamic_cast<jura::ToolChain*>(mod);
  expect(tlChn != nullptr);

  // Check that initially, there is one module of type "None" in the slot 1 with index 0:
  expect(tlChn->getNumModules() == 1);
  mod = tlChn->getModuleAt(0);
  expect(mod != nullptr);
  jura::DummyModule* dum = dynamic_cast<jura::DummyModule*>(mod);
  expect(dum != nullptr);

  // The tests up to here could perhaps be factored out into a function runTestToolChainCreation.
  // They are not specific to Quadrifex. Maybe we could also use soem sort of factory function that
  // we can call in a one-line like:
  //   jura::ToolChain* tlChn = createToolChain(numMetaParams);
  // because that could be useful for other ToolChain tests.


  // Let the ToolChain module insert a Quadrifex into slot 2 with index 1:
  bool ok = tlChn->addModule("Quadrifex");
  expect(ok);
  mod = tlChn->getModuleAt(1);
  expect(mod != nullptr);
  jura::QuadrifexAudioModule* qfx = dynamic_cast<jura::QuadrifexAudioModule*>(mod);
  expect(qfx != nullptr);


  // This old test code previously triggered access violations - this has been fixed. It occurred
  // in the constructor of rosic:: FrequencyShifter. There is now unit test in TestsRosicAndRapt in
  // place that makes sure that this doesn't happen again:

  // Try to create a rosic::FrequencyShifter object. We seem to have an access violation in its
  // constructor, specifically in the line 
  //   halfbandFilter2.setApproximationMethod(...)
  // ToDo: maybe move that to the rosic unit tests:
  //rosic::FrequencyShifter freqShifter;
  // it happens in  rsEngineersFilter<TSig, TPar>::updateCoefficients(bool resetState) in the line
  // rsBiquadCascade<TSig, TPar>::initBiquadCoeffs(); and I think TSig=rsfloat64x2, TPar=double.

  // Let the Quadrifex load FrequencyShifterStereoModule - this is where we have an access 
  // violation:
  //using QFX = rosic::Quadrifex;
  //qfx->setEffectAlgorithm(0, QFX::FREQUENCY_SHIFTER);
  // YES! We successfully trigger it here! Now we can figure out what is going wrong....
  // It seems to happen in the constructor of rosic::FrequencyShifter, so above, we just create one
  // outside of the Quadrifex setting.


  // Clean up memory:
  delete proc;
}

void UnitTestToolChain::runTestWaveOscillator()
{
  CriticalSection lock;                   // Mocks the pluginLock.
  jura::WaveOscModule wvOsc1(&lock);

  // Helper function to check if the given oscillator is in default state, i.e. all parameters are
  // at the default values:
  auto isInDefaultState = [](const jura::WaveOscModule* osc)
  {
    bool ok = true;

    ok &= osc->getParameterByName("StereoPhaseShift")->getValue() == 0.0;
    // ...more checks to come....
  
    return ok;
  };

  // Check that all parameters have their intial/default values:
  expect(isInDefaultState(&wvOsc1));

  // Create an editor object for the osc:
  jura::AudioModuleEditor* ed1 = wvOsc1.createEditor(0);

  // Check that the right kind of editor was created:
  expect( dynamic_cast<jura::WaveOscEditor*>(ed1) );

  // Check that creating the editor didn't mess with the state, i.e. didn't change any parameters:
  expect(isInDefaultState(&wvOsc1));
  // THIS TRIGGERS! That was the goal of writing this unit test. Now the debugging can begin...
  // Creating an editor for an AudioModule should never alter the state of that module!
  // To debug this, we need a breakpoint in:
  //   rosic::MipMappedWaveTableStereo::setStereoPhaseShift
  // Creating the editor calls it with  newPhaseShift == 1. Why? This should not happen. Also, it 
  // seems to get called excessively often during startup - why? One call would be expected but 
  // there are many more. Ah - I think, the additional calls are coming from creating a ToolChain
  // object in another test -> change order of the tests -> yep! that fixes it!



  int dummy = 0;


  // ToDo: 
  // -Test creating a jura::WaveOscModule that creates the underlying DSP core of the oscillator 
  //  itself and one that wraps an existing oscillator.
}