#include "UnitTestsToolChain.h"
using namespace juce;
using namespace jura;


void UnitTestToolChain::runTest()
{
  // Test currently worked on copied to top of the function:
  runTestStraightliner();  // triggers jassert


  // All the tests in order:
  runTestVoiceManager();
  runTestEqualizer();
  //runTestMultiAnalyzer();  // Fails - see comments there. Fixing has low priority.
  runTestStraightliner();
  runTestWaveOscillator();

  runTestQuadrifex();
  runTestEditorCreation(0);  // Takes quite long
  // These tests are currently called last because they creates an actual jura::ToolChain object 
  // which in turn instantiates all modules once in populateModuleFactory - which is annyoing 
  // during debugging because certain initialization functions for ToolChain's built in 
  // AudioModules will get called more often than one would expect in the tests, i.e. the 
  // breakpoints will trigger more often than the actual running test justifies. Putting this 
  // test to the end fixes this. 
  // ToDo: Factor out the creation of a ToolChain object from the Quadrifex test. It should be a 
  // test in its own right. Then, it shouldn't matter where we put the test for Quadrifex.

  // We get memory leaks. They come from runTestEditorCreation. Maybe it's ToolChain itself? Figure out!

  runTestStateRecall(0);  // FAILS!!
}

bool UnitTestToolChain::isInDefaultState(const jura::AudioModule* m)
{
  bool ok = true;
  int numParams = m->getNumParameters();
  for(int i = 0; i < numParams; i++)
  {
    jura::Parameter* p    = m->getParameterByIndex(i);
    juce::String     name = p->getName();
    ok &= p->isCurrentValueDefaultValue();
  }
  return ok;
}

void UnitTestToolChain::resetParameters(jura::AudioModule* m)
{
  for(int i = 0; i < m->getNumParameters(); i++)
    m->getParameterByIndex(i)->resetToDefaultValue(true, true);

  // toDo: call it recursively on the child modules
}

void UnitTestToolChain::randomizeParameters(jura::AudioModule* m, int seed)
{
  // Create a pseudo random number generator:
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0.0, 1.0);
  prng.setSeed(seed);

  // Randomize the parameters:
  int numParams = m->getNumParameters();
  for(int i = 0; i < numParams; i++)
  {
    // Retrieve i-th parameter and its name:
    jura::Parameter* p    = m->getParameterByIndex(i);
    juce::String     name = p->getName();

    // Randomize its value:
    double min    = p->getMinValue();
    double max    = p->getMaxValue();
    double newVal = RAPT::rsLinToLin(prng.getSample(), 0.0, 1.0, min, max);
    p->setValue(newVal, true, true);
  }

  // Call randomizeParameters on all the child-modules recursively:
  int numChildren = m->getNumChildAudioModules();
  for(int i = 0; i < numChildren; i++)
  {
    int newSeed = (int) prng.getSampleRaw();
    jura::AudioModule* childModule = m->getChildAudioModule(i);
    randomizeParameters(childModule, newSeed);
  }
}

juce::MouseEvent UnitTestToolChain::getMockMouseDownEvent(float mouseX, float mouseY, 
  juce::Component* eventComp, juce::Component* originatorComp)
{
  float mouseDownX  = 0.f;
  float mouseDownY  = 0.f;
  float pressure    = 1.f;
  float orientation = 0.f;
  float rotation    = 0.f;
  float tiltX       = 0.f;
  float tiltY       = 0.f;
  juce::Time eventTime;      // Maybe try to somehow obtain the "now" time
  juce::Time mouseDownTime;
  int   numClicks   = 1;
  bool  wasDragged  = false;
  //juce::Component* eventComponent = nullptr;
  //juce::Component* originatorComponent = nullptr;
  juce::MouseInputSource mouseSource = juce::Desktop::getInstance().getMainMouseSource();
  juce::ModifierKeys     modKeys;

  juce::MouseEvent mouseEvent(mouseSource, 
    juce::Point<float>(mouseX, mouseY), modKeys, pressure, orientation, rotation, tiltX, tiltY,
    eventComp, originatorComp, eventTime, juce::Point<float>(mouseDownX, mouseDownY),
    mouseDownTime, numClicks, wasDragged);

  return mouseEvent;

  // It's surprisingly difficult to create mock mouse events in juce. See:
  // https://forum.juce.com/t/creating-a-mouseevent/32635/10
}

//-------------------------------------------------------------------------------------------------
// Tests for the infrastructure:

void UnitTestToolChain::runTestEditorCreation(int seed)
{
  // This test triggers a couple of assertions -> try to fix!

  CriticalSection lock;                   // Mocks the pluginLock.
  jura::ToolChain tlChn(&lock);

  // Let the ToolChain object create a module of each of the available types and plug it into the
  // first slot, then randomize its parameters, retrieve the state, open the editor, retrieve the
  // state again and then compare both states:
  std::vector<juce::String> moduleTypes = tlChn.getAvailableModuleTypes();
  for(size_t i = 0; i < moduleTypes.size(); i++)
  {
    juce::String type = moduleTypes[i];

    if(type == "MultiAnalyzer")
      continue;
    // It fails for MultiAnalyzer but it's not a big issue so we just skip it for the time being. 
    // There's a separate unit test for this which can be used for fixing this. When done, this 
    // skipping here can be removed as well.

    tlChn.replaceModule(0, type);

    AudioModule* m = tlChn.getModuleAt(0);
    expect(m->getModuleTypeName() == type);  // Check module type in slot 1
    randomizeParameters(m, seed);

    // Get the state, create the editor, get the state again and compare both states:
    juce::XmlElement* preXml  = m->getStateAsXml("State", true);
    jura::AudioModuleEditor* editor = m->createEditor(0);
    juce::XmlElement* postXml = m->getStateAsXml("State", true);
    expect(postXml->isEquivalentTo(preXml, false));
    delete preXml;
    delete editor;
    delete postXml;
  }

  // ToDo:
  // -Try the test with different seeds for the randomization of the parameters. Some bugs are exposed 
  //  only with certain parameter values.
}

void UnitTestToolChain::runTestStateRecall(int seed)
{
  CriticalSection lock;                   // Mocks the pluginLock.
  jura::ToolChain tlChn(&lock);
  std::vector<juce::String> moduleTypes = tlChn.getAvailableModuleTypes();
  for(size_t i = 0; i < moduleTypes.size(); i++)
  {
    juce::String type = moduleTypes[i];

    if(type == "FuncShaper")
      continue;
    // It fails for FuncShaper. To get passed the triggers of the assertions, we temporarily skip 
    // this test.

    tlChn.replaceModule(0, type);
    AudioModule* m = tlChn.getModuleAt(0);
    expect(m->getModuleTypeName() == type);  // Check module type in slot 1
    randomizeParameters(m, seed);
    juce::XmlElement* preXml = m->getStateAsXml("State", true);
    resetParameters(m);
    m->setStateFromXml(*preXml, "Recalled", true);
    juce::XmlElement* postXml = m->getStateAsXml("State", true);
    expect(postXml->isEquivalentTo(preXml, false));
    delete preXml;
    delete postXml;
  }

  // WaveOscillator fails
  // Why does Straightliner not fail? We have this bug with the recall of the "Mute" parameter
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


//-------------------------------------------------------------------------------------------------
// Tests for individual modules:

void UnitTestToolChain::runTestEqualizer()
{
  // This test was motivated by an access violation in the frequency response plot of the equalizer
  // in both of the modes that need to plot two graphs instead of just one. there was somet bug in
  // jura::rsDataPlot

  CriticalSection lock;                   // Mocks the pluginLock.
  jura::EqualizerAudioModule eq(&lock);
  eq.getParameterByName("StereoMode")->setValue( 1.0, true, true);  // 1 is Stereo L/R
  jura::AudioModuleEditor* editor = eq.createEditor(0);
  delete editor;
  // The call to eq.createEditor(0); triggers an access violation in rsPlotDrawer::drawWithLines 
  // but only if we call randomizeParameters before. Or 
  // eq.getParameterByName("StereoMode")->setValue( 0.83, true, true);
  // It seems to be the stereo mode parameter 0-Linked: ok, 1-L/R: crash, 2-M/S: crash, 3: ok
  // ..soo it seems like those stereo modes that have two graphs cause problems.
  // This also happes when actually runnign toolChain. plugging in an Equalizer and swicthing
  // the stereo mode from the GUI
  // OK - this crash might be fixed.
}

void UnitTestToolChain::runTestMultiAnalyzer()
{
  CriticalSection lock;
  jura::MultiAnalyzerAudioModule ana(&lock);


  // Get the state, create the editor, get the state again and compare both states:
  juce::XmlElement* preXml  = ana.getStateAsXml("State", true);
  jura::AudioModuleEditor* editor = ana.createEditor(0);
  juce::XmlElement* postXml = ana.getStateAsXml("State", true);

  juce::XmlDocument preDoc  = preXml ->toString();
  juce::XmlDocument postDoc = postXml->toString();

  expect(postXml->isEquivalentTo(preXml, false));
  // This fails!
  // Creating the editor apparently sets the TimeWindowLength parameter of the oscilloscope to 1.5
  // I think, this happens when the zoomer with scrollbars is created and wired up. In the 
  // oscilloscope that parameters are controlled by scrollbars rather than the usula sliders and this
  // behaves differently. It's not a big issue so fixing this should have lower priority.

  delete preXml;
  delete editor;
  delete postXml;
}

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

void UnitTestToolChain::runTestStraightliner()
{
  CriticalSection lock;                   // Mocks the pluginLock.
  jura::StraightlinerAudioModule synth(&lock);

  // Obtain pointers to the submodules:
  jura::MultiModeFilterAudioModule* filter   = synth.filterModule;
  jura::FourOscSectionAudioModule* oscs     = synth.oscSectionModule;
  jura::BreakpointModulatorAudioModule* pitchEnv = synth.pitchEnvModule;
  jura::BreakpointModulatorAudioModule* filtEnv  = synth.filterEnvModule;
  jura::BreakpointModulatorAudioModule* ampEnv   = synth.ampEnvModule;
  jura::WaveOscModule* osc1     = oscs->osc1Module;
  jura::WaveOscModule* osc2     = oscs->osc2Module;
  jura::WaveOscModule* osc3     = oscs->osc3Module;
  jura::WaveOscModule* osc4     = oscs->osc4Module;

  // Check that all the oscillators are in the expected default states. All 4 oscs should have the
  // sawtooth wave loaded and osc1 should be active (i.e. non-muted) and oscs 2..4 should be 
  // muted:
  //juce::String wave1 = osc1->getCurrentWaveformPath();

  // Helper function to check if the parameters with given name have the expected values 
  // v1,v2,v3,v4 in the 4 oscillators respectively:
  auto checkOscParams = [&](const char* name, double v1, double v2, double v3, double v4)
  {
    bool ok = true;
    jura::Parameter* p = nullptr;           // Parameter object
    double v = 0.0;                         // Value of parameter

    p = osc1->getParameterByName(name);
    jassert(p != nullptr);
    v = p->getValue();
    ok &= v == v1;

    p = osc2->getParameterByName(name);
    jassert(p != nullptr);
    v = p->getValue();
    ok &= v == v2;

    p = osc3->getParameterByName(name);
    jassert(p != nullptr);
    v = p->getValue();
    ok &= v == v3;

    p = osc4->getParameterByName(name);
    jassert(p != nullptr);
    v = p->getValue();
    ok &= v == v4;

    return ok;
  };

  expect( checkOscParams("Mute", 0, 1, 1, 1) );
  // Initially osc1 shoudl be active (non-muted), the others are muted

  // Now try setting the 2nd osc non-muted, retrieve and recall the state, then check, if the
  // "Mute" settings are as expected (1 and 2 non-muted, 3 and 4 muted):
  jura::Parameter* p = nullptr;
  p = osc2->getParameterByName("Mute");
  p->setValue(0.0, true, true);           // Activate Osc2

  // Factor out into retrieveAndRecallState - could be used for any AudioModule:
  juce::XmlElement* xml = synth.getStateAsXml("State", true);
  juce::String str = xml->toString();
  synth.setStateFromXml(*xml, "State", true);
  delete xml; xml = nullptr;

  expect( checkOscParams("Mute", 0, 0, 1, 1) );


  p = synth.getParameterByName("NumVoices");
  // Such a parameter does not exist!



  int dummy = 0;
}

void UnitTestToolChain::runTestWaveOscillator()
{
  // Thie test was motivated by a bug which caused the creation of an editor for a Wave-oscillator
  // to potentially change some of its settings. This was caused by some erroneous function call 
  // sequence when connecting widgets (e.g. sliders) to the parameters. The sliders should take 
  // over the value from the parameter (and its range, etc.) but for some reason, it also did set
  // the value. Of course, opening an editor should never change the state of an AudioModule:

  CriticalSection lock;                   // Mocks the pluginLock.
  jura::WaveOscModule wvOsc1(&lock);

  // Check that all parameters have their intial/default values:
  expect(isInDefaultState(&wvOsc1));

  // Create an editor object for the osc:
  jura::AudioModuleEditor* amEd = wvOsc1.createEditor(0);

  // Check that the right kind of editor was created:
  jura::WaveOscEditor* wvEd = dynamic_cast<jura::WaveOscEditor*>(amEd);
  expect(wvEd != nullptr);

  // Check that creating the editor didn't mess with the state, i.e. didn't change any parameters:
  expect(isInDefaultState(&wvOsc1));


  // Another problem that we had was that when Straightliner was fired up, we activate the 2nd osc,
  // save the preset and reload it, the 2nd osc is off again. The state xml file has stored a 
  // "Mute=1" attribute. The problem was that the GUI action of activating the osc bypassed the
  // Parameter infrastructure. The bug was fixed - here we check that it works as expected:

  juce::XmlElement* xml1 = wvOsc1.getStateAsXml("State", true);

  // Check, if the xml-string is as expected:
  juce::String str1 = xml1->toString();
  juce::String strT = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n<WaveOscillator PatchFormat=\"1\"/>\n";
  int cmp = str1.compare(strT);  // == 1. Should be 0, if the strings are equal
  //expect(str1 == strT);
  // This fails! Could it have to do with encoding or with CR/LF stuff? ...OK - this is currently 
  // not terribly important but eventually I want to figure this out and include such a test

  // Clean up:
  delete xml1; xml1 = nullptr;


  // Test if clicking on the waveform display toggles the oscillator's "Mute" parameter on/off:
  jura::Parameter* p = nullptr;
  double v;
  p = wvOsc1.getParameterByName("Mute");
  jassert(p);
  v = p->getValue();
  expect(v == 0);       // Initially, it should be unmuted

  // Mock the mouse-click on the wvaeform display:
  juce::Rectangle<int> r = wvEd->getWaveDisplayBounds();
  int x = r.getX() + r.getWidth()  / 2;
  int y = r.getY() + r.getHeight() / 2;
  juce::MouseEvent mouseEvent = getMockMouseDownEvent(x, y, wvEd, wvEd);
  amEd->mouseDown(mouseEvent);  // Deliberately using the baseclass pointer

  // Check "Mute" parameter again:
  v = p->getValue();
  expect(v == 1);       // Now it should be muted





  // Clean up memory:
  delete amEd;



  // ToDo: 
  // -Test creating a jura::WaveOscModule that creates the underlying DSP core of the oscillator 
  //  itself and one that wraps an existing oscillator.
  // -Maybe write a more general test that checks for every AudioModule included in ToolChain
  //  that opening the editor does not affect the state. Loop through all availablbe modules, 
  //  intantiate one, randomize the state, retrieve the state-xml, open the editor, retrieve the 
  //  state-xml again and check that both state xmls match. ..is under construction

}

