#include "UnitTestsToolChain.h"
using namespace juce;
using namespace jura;


void UnitTestToolChain::runTest()
{
  // For certain unit tests that loop over all available module, we know that the test will fail
  // for certain modules but this is not problematic, if these are modules that should not yet be
  // included in the next release anyway. So we make some ignore lists. The goal is that they 
  // should get shorter over time
  juce::StringArray ignoreListForStateRecall = { "FuncShaper", "MultiBandEffect", "SamplePlayer" };
  juce::StringArray ignoreListForEditorCreation = { "MultiAnalyzer" };


  // Test currently worked on copied to top of the function:
  runTestEchoLab();

  
  /*
  // All the tests in order:
  runTestVoiceManager();
  runTestSlotInsertRemoveEtc();
  runTestStateRecall(0, ignoreListForStateRecall);
  runTestEqualizer();
  runTestFuncShaper();
  //runTestMultiAnalyzer();  // Fails - see comments there. Fixing has low priority.
  runTestStraightliner();
  runTestWaveOscillator();
  runTestEchoLab();
  runTestQuadrifex();
  runTestEditorCreation(0, ignoreListForEditorCreation);  // Takes quite long
  */
  // These tests are currently called last because they creates an actual jura::ToolChain object 
  // which in turn instantiates all modules once in populateModuleFactory - which is annyoing 
  // during debugging because certain initialization functions for ToolChain's built in 
  // AudioModules will get called more often than one would expect in the tests, i.e. the 
  // breakpoints will trigger more often than the actual running test justifies. Putting this 
  // test to the end fixes this. 
  // ToDo: Factor out the creation of a ToolChain object from the Quadrifex test. It should be a 
  // test in its own right. Then, it shouldn't matter where we put the test for Quadrifex.


  // 


  // Notes:
  //
  // -FuncShaper is actually included in the upcoming release but we have to exclude it from the 
  //  state recall test because of an idiosynchrasy it has: It has the aMin, aMax, bMin, ...etc. 
  //  parameters which have infinite range. That trips up the state recall unit test in a way that 
  //  is understood and expected (see comments in runtestFuncShaper for details). If that's the 
  //  only problem, it should be fine because that state recall test produces states that are 
  //  invalid for for FuncShaper in its randomization of the parameters which are sanitized in the
  //  attempted recall - which causes the test to fail because the sanitized state doesn't match 
  //  the radomized state.
  //
  // We get memory leaks. They come from runTestEditorCreation. Maybe it's ToolChain itself? 
  // Figure out! When coomenting out the line
  //
  //   jura::AudioModuleEditor* editor = m->createEditor(0);
  //
  // the memleak goes away. Apparently, one of the editors causes it. Figure out which one it is!
  // Maybe comment out some sections in ToolChainAudioModule::populateModuleFactory. OK - it looks
  // like EchoLab is the offender.
  //
  // ToDo:
  // -Make a test for EchoLab - it should create the editor. We want to fix the memory leak that it
  //  causes.
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

bool UnitTestToolChain::doSlotsContain(const jura::ToolChain* toolChain, 
  const std::vector<juce::String>& typeNames)
{
  bool ok = toolChain->getNumModules() == (int) typeNames.size();
  if(!ok)
    return false;
  for(int i = 0; i < (int)typeNames.size(); i++)
    ok &= toolChain->isModuleOfType(i, typeNames[i]);
  return ok;
}

bool UnitTestToolChain::areArraysConsistent(jura::ToolChainEditor* editor)
{
  bool ok = true;
  const jura::ToolChain* tlChn = dynamic_cast<const jura::ToolChain*>(editor->getModuleToEdit());
  expect(tlChn != nullptr);

  // Check if array sizes match:
  int numModules = tlChn->getNumModules();
  ok &= editor->editors.size()   == numModules;
  ok &= editor->selectors.size() == numModules;
  if(!ok)
    return false;

  // Check if array contents match (plus some additional stuff):
  for(int i = 0; i < numModules; i++)
  {
    // Check if i-th editor corresponds to i-th module:
    jura::AudioModule*       m  = tlChn->modules[i];           // i-th module
    jura::AudioModuleEditor* e  = editor->getEditorForSlot(i); // i-th editor (call may create it)
    const jura::AudioModule* em = e->getModuleToEdit();        // module edited by i-th editor
    ok &= m == em;

    // Check if i-th selector corresponds to the type of i-th module:
    jura::AudioModuleSelector* s = editor->selectors[i];       // i-th selector (combo box)
    juce::String selectorText = s->getSelectedItemText();
    juce::String typeName   = m->getModuleTypeName();
    ok &= selectorText == typeName;

    // Check, if the display name in the module is correct. It should be Slot + "i-" + "TypeName"
    juce::String moduleName = m->getModuleName();
    ok &= moduleName == "Slot" + juce::String(i+1) + "-" + typeName;

    // Check if display name in the editor is correct:
    juce::String editorName = e->getHeadlineString();
    ok &= editorName == moduleName;
    int dummy = 0;
  }

  return ok;

  // Notes:
  // -The function getEditorForSlot(i) may either return an existing editor or create one if the 
  //  editor for the respective slot did not already exist.
  // -A "selector" is a combo-box by which the type of AudioModule can be selected
}


void UnitTestToolChain::resetParameters(jura::AudioModule* m)
{
  for(int i = 0; i < m->getNumParameters(); i++)
    m->getParameterByIndex(i)->resetToDefaultValue(true, true);

  // toDo: call it recursively on the child modules
}

void UnitTestToolChain::randomizeParameters(jura::AudioModule* m, int seed, bool recursively)
{
  // Create a pseudo random number generator:
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0.0, 1.0);
  prng.setSeed(seed);

  // Randomize the parameters:
  int numParams = m->getNumParameters();
  for(int i = 0; i < numParams; i++)
  {
    jura::Parameter* p    = m->getParameterByIndex(i);   // Retrieve i-th parameter
    juce::String     name = p->getName();                // Name for inspection in the debugger


    // OLD:
    double min    = p->getMinValue();
    double max    = p->getMaxValue();
    double newVal = RAPT::rsLinToLin(prng.getSample(), 0.0, 1.0, min, max);
    p->setValue(newVal, true, true);   // Randomize its value


    // NEW:
    //p->setNormalizedValue(prng.getSample(), true, true); // Randomize its value
    // Using this new version breaks the state recall test. I think, this breakage could have 
    // something to do with with parameter quantization. Comparing some pre/post state XMLs with
    // this new version of the code, some parameters are quantized in one of the states but not in
    // the other.
  }

  // Call randomizeParameters on all the child-modules recursively:
  if(recursively == true)
  {
    int numChildren = m->getNumChildAudioModules();
    for(int i = 0; i < numChildren; i++)
    {
      int newSeed = (int)prng.getSampleRaw();
      jura::AudioModule* childModule = m->getChildAudioModule(i);
      randomizeParameters(childModule, newSeed, true);
    }
  }
}

void UnitTestToolChain::randomizeAudioModule(jura::AudioModule* m, int seed, bool recursively)
{
  // Randomize the parameters:
  randomizeParameters(m, seed);
  // ToDo: give this a bool parameter "recurseChildModules" and pass false to it from here because
  // we already call ourselves recursively here which will then also randomize the parameters of 
  // the child modules

  // Randomize other aspects - this is individually different for the different types of modules,
  // so we dispatch based on the module type:
  // ....


  // Create a pseudo random number generator:
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0.0, 1.0);
  prng.setSeed(seed);





  // Call randomizeAudioModule on all the child-modules recursively:
  if(recursively == true)
  {
    int numChildren = m->getNumChildAudioModules();
    for(int i = 0; i < numChildren; i++)
    {
      int newSeed = (int)prng.getSampleRaw();
      jura::AudioModule* childModule = m->getChildAudioModule(i);
      randomizeAudioModule(childModule, newSeed, true);
    }
  }
  // This code is repetitive. The last section of randomizeParameters looks almost the same. Try
  // to refactor to get rid of the duplication! Maybe we could have a function that could be called
  // like:
  //
  //   callRecursively(&randomizeAudioModule, m, prng.getSampleRaw));
  //
  // The problem is probably that we cannot easily pass a pointer to a member function, i.e.
  // something like &randomizeAudioModule will likely not compile. We may need some sort of member
  // function callback. Or maybe use std::function for the parameter and wrap the call into a 
  // lambda here.
  // Or maybe at least try to reduce the amount of repetitive code. Get rid of the variables
  // newSeed, childModule, numChildren.
}

bool UnitTestToolChain::testStateRecall(jura::AudioModule* m, int seed)
{
  randomizeParameters(m, seed);    // ToDo: use randomizeAudioModule
  juce::XmlElement* preXml = m->getStateAsXml("State", true);
  resetParameters(m);              // ToDo: use resetAudioModule
  m->setStateFromXml(*preXml, "Recalled", true);
  juce::XmlElement* postXml = m->getStateAsXml("State", true);

  juce::String preStr  = preXml->toString();
  juce::String postStr = postXml->toString();

  bool ok = postXml->isEquivalentTo(preXml, false);
  delete preXml;
  delete postXml;
  return ok;
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

std::vector<jura::RWidget*> UnitTestToolChain::getWidgetsWithoutParameter(
  jura::AudioModuleEditor* editor)
{
  std::vector<jura::RWidget*> orphans;
  for(size_t i = 0; i < editor->widgets.size(); i++)
  {
    RWidget* w = editor->widgets[i];
    if(!w->hasAssignedParameter())
      orphans.push_back(w);
  }
  return orphans;
}

template<class WidgetType>
std::vector<WidgetType*> UnitTestToolChain::filterWidgets(const std::vector<jura::RWidget*>& widgets)
{
  std::vector<WidgetType*> matches;
  for(size_t i = 0; i < widgets.size(); i++)
  {
    WidgetType* casted = dynamic_cast<WidgetType*>(widgets[i]);
    if(casted != nullptr)
      matches.push_back(casted);
  }
  return matches;
}


//-------------------------------------------------------------------------------------------------
// Tests for the infrastructure:

void UnitTestToolChain::runTestEditorCreation(int seed, const juce::StringArray& ignoreList)
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

    if(ignoreList.contains(type))
      continue;

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

void UnitTestToolChain::runTestStateRecall(int seed, const juce::StringArray& ignoreList)
{
  CriticalSection lock;                   // Mocks the pluginLock.
  jura::ToolChain tlChn(&lock);
  std::vector<juce::String> moduleTypes = tlChn.getAvailableModuleTypes();
  for(size_t i = 0; i < moduleTypes.size(); i++)
  {
    juce::String type = moduleTypes[i];

    if(ignoreList.contains(type))
      continue;

    tlChn.replaceModule(0, type);
    AudioModule* m = tlChn.getModuleAt(0);
    expect(m->getModuleTypeName() == type);  // Check module type in slot 1

    expect(testStateRecall(m, seed));
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



void UnitTestToolChain::runTestSlotInsertRemoveEtc()
{
  CriticalSection lock;                   // Mocks the pluginLock.
  jura::ToolChain tlChn(&lock);

  jura::ToolChainEditor* editor = dynamic_cast<jura::ToolChainEditor*> (tlChn.createEditor(0));
  expect(editor != nullptr);

  // Some shorthands for the names of the different mdoule types that we intend to add:
  //std::vector<juce::String> moduleTypes = tlChn.getAvailableModuleTypes();

  juce::String eq  = "Equalizer";
  juce::String fs  = "FuncShaper";
  juce::String ldr = "Ladder";
  juce::String scp = "Scope";
  juce::String non = "None";

  // Insert some modules into the slots via the editor and check if after each insertion, the slots
  // are filled with the expected modules. Replacing the "None" module at the end with some actually
  // useful module will also have the side effect of selecting/activating the just inserted module.

  expect(doSlotsContain(&tlChn, { non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 0); 

  editor->replaceModule(0, eq);
  expect(doSlotsContain(&tlChn, { eq, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 0); 

  editor->replaceModule(1, fs);
  expect(doSlotsContain(&tlChn, { eq, fs, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 1); 

  editor->replaceModule(2, ldr);
  expect(doSlotsContain(&tlChn, { eq, fs, ldr, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 2); 

  editor->replaceModule(3, scp);
  expect(doSlotsContain(&tlChn, { eq, fs, ldr, scp, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 3); 


  // Test swapping modules. If one of the swapped modules is the active one, the new active one 
  // will also switch to the other. This is the behavior we wnat when bubbling up or down the 
  // selectors - we want the selected module to remain selected, even when it's no at a different
  // position.

  editor->swapModules(1, 3);
  expect(doSlotsContain(&tlChn, { eq, scp, ldr, fs, non }));
  expect(areArraysConsistent(editor));  // FAILS!
  expect(tlChn.activeSlot == 1);

  // Now do the swapping via mock GUI events:
  auto moveUp = [&]()  // Mock click on "Up" button
  {
    RButton* btn = editor->moveUpButton;
    juce::MouseEvent mouseEvent = getMockMouseDownEvent(8.f, 8.f, btn, btn);
    btn->mouseDown(mouseEvent);
  };
  auto moveDown = [&]()  // Mock click on "Up" button
  {
    RButton* btn = editor->moveDownButton;
    juce::MouseEvent mouseEvent = getMockMouseDownEvent(8.f, 8.f, btn, btn);
    btn->mouseDown(mouseEvent);
  };

  moveUp();    // Moves Scope from slot 1 to slot 0
  expect(doSlotsContain(&tlChn, { scp, eq, ldr, fs, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 0);

  moveUp();    // Does nothing because the active slot is on top
  expect(doSlotsContain(&tlChn, { scp, eq, ldr, fs, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 0);

  moveDown();  // Moves the Scope back into slot 1
  expect(doSlotsContain(&tlChn, { eq, scp, ldr, fs, non }));
  expect(areArraysConsistent(editor));
  expect(tlChn.activeSlot == 1);


  // Now check the behavior when "None" modules are moved to/from the bottom
  expect(tlChn.getNumModules() == 5);
  tlChn.activeSlot = 3;  // 2nd to last, FuncShaper
  moveDown();
  expect(tlChn.getNumModules() == 6);
  tlChn.activeSlot = 4; 
  moveUp();
  expect(tlChn.getNumModules() == 5);
  tlChn.activeSlot = 3;





  delete editor;  
  // We need to delete the editor before ToolChain gets out of scope - otherwise we trigger an
  // access violation. But this should be expected behavior - we expect anyway that an AudioModule
  // will always outlive its editor. The exception is when the editor is also an 
  // AudioModuleDeletionWatcher. This mechanism can be used to un-attach editors from modules that
  // may somehow get deleted before their editor. This is useful for re-use existing editors like
  // it is done in EchoLabDelayLineModuleEditor.


  // ToDo:
  // -Test what happens when we move a "None" module to the bottom and when we swap the "None" 
  //  module at the bottom with the module before it. Also check what happens when we move the
  //  second-to-bottom module down to the last slot. In any case, there should alway be exactly 
  //  one "None" moduel at the bottom. Check also when a none module comes from above and hits
  //  the secodn to last slot. 
  // -Maybe instead of calling replaceModule, mock the GUI actions that would in practice trigger
  //  these calls.
  // -I think, swapping two modules will potentially lead to problems when the state is saved and 
  //  recalled when the swapped modules are modulators. The evaluation order may be different. I 
  //  think we may need to also change their order in the list of modulators. Maybe set up some 
  //  cross-feedback modulation between two modulator modules, then swap them, then check, if the 
  //  output signal is still the differen. Maybe two LFO that cross-feedback modulate their 
  //  frequencies. In such a configuration, the evaluation order matters because of the feedback 
  //  delay
}

//-------------------------------------------------------------------------------------------------
// Tests for individual modules:

void UnitTestToolChain::runTestEqualizer()
{
  // This test was motivated by an access violation in the frequency response plot of the equalizer
  // in both of the modes that need to plot two graphs instead of just one. there was somet bug in
  // jura::rsDataPlot

  CriticalSection lock;                   // Mocks the pluginLock.

  // Tests for the Equalizer in isolation:
  {
    jura::EqualizerAudioModule eq(&lock);
    eq.getParameterByName("StereoMode")->setValue(1.0, true, true);  // 1 is Stereo L/R
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


  // A test using a pointer:
  {
    jura::EqualizerAudioModule* eq = new jura::EqualizerAudioModule(&lock);
    jura::AudioModuleEditor* editor = eq->createEditor(0);

    delete editor;
    // Needed to avoid the access violation below.

    delete eq;  
    // This trigger the access violation that triggers a debug break unless we delete the editor 
    // first. And the call stack in the debugger is *very* strange. There are functions in it that 
    // do not seem to be called by outer functions. From the stack, it looks like 
    // FileManager::setActiveFileIfInList is called from EqualizerPlotEditor::updatePlot. 
    // But it isn't called there. It calls equalizerModuleToEdit->getMagnitudeResponse. Very weird! 
    //
    // Now I'm not totally sure but I think within ToolChain in it's usual context, the AudioModule
    // will always outlive its editor, so a situation in which a module is deleted before the 
    // corresponding editor is deleted will not occur in practice. Or will it? We actually do have
    // some infrastructure in place to handle such situations - namely 
    // jura::AudioModuleDeletionWatcher. 
  }



  
  /*
  // Tests for the Equalizer in ToolChain:
  {
    jura::ToolChain tlChn(&lock);
    jura::ToolChainEditor* editor = dynamic_cast<jura::ToolChainEditor*> (tlChn.createEditor(0));
    expect(editor != nullptr);

    editor->replaceModule(0, "Equalizer"); 

  }
  */
  
}

void UnitTestToolChain::runTestFuncShaper()
{
  CriticalSection lock;                   // Mocks the pluginLock.

  // The state recall unit test for FuncShaper triggered a jassert. When the state is recalled, it 
  // want to set up the aMin, aMax, bMin, etc. to NaN. The stored xml element actually contains
  // these NaNs. The min/max parameters get set to NaN in randomizeParameters.


  jura::FuncShaperAudioModule fnSh(&lock);

  randomizeParameters(&fnSh, 0);

  jura::Parameter*  p = fnSh.getParameterByName("aMin");
  // Value of p is NaN. It happened in the randomization.


  juce::XmlElement* xml = fnSh.getStateAsXml("", false);
  juce::String      str = xml->toString();
  delete xml;
  // That xml string does only contain the offending NaNs if we didi the randomization


  //expect(testStateRecall(&fnSh)); 
  // This test would fail because the randomization of the parameters does does not respect the
  // aMin <= a <= aMax etc. conditions for a valid FuncShaper state. Upon recall, these randomized
  // and invalid values are sanitized leading to a different state and therefore to a failure in
  // the unit test. One solution to this could be to treat different kinds of modules differently 
  // in the randomization by dispatching to different randomization functions based on 
  // dynamic_cast. With this, we could also randomize more than just the parameters - for example,
  // the sample loaded or a string parameter. Yeah...that sound attractive.
  



  int dummy = 0;
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
  // oscilloscope that parameters are controlled by scrollbars rather than the usual sliders and 
  // this behaves differently. It's not a big issue so fixing this should have lower priority.

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


  jura::AudioModuleEditor* synthEditor = synth.createEditor(0);
  std::vector<jura::RWidget*> widgets = getWidgetsWithoutParameter(synthEditor);
  // Whoa! It has 10 widgets without parameter. Some of them are labels and text fields - which is
  // OK.

  // Filter out only the sliders:
  std::vector<jura::RSlider*> sliders = filterWidgets<RSlider>(widgets);
  expect(sliders.size() == 0); 
 
  // Filter out only the buttons:
  std::vector<jura::RButton*> buttons = filterWidgets<RButton>(widgets);
  expect(buttons.size() == 5); 
  // 5 buttons without an assigned Parameter: Tuning Load/Plus/Minus, the (invisible) "Setup" 
  // button and the (also invisible) weblink. Try to get rid of the latter two

  // We have additionally the "Glide" button orphaned - that's a bug!



  delete synthEditor;

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


void UnitTestToolChain::runTestEchoLab()
{
  CriticalSection lock;                   // Mocks the pluginLock.

  // We want to trigger the memory leak in the editor:
  jura::EchoLabAudioModule el(&lock);
  jura::AudioModuleEditor* editor = el.createEditor(0);
  delete editor;
  // It seems to be the call to setSize() at the bottom of the constructor of the editor. When 
  // commenting it out, the leak goes away. setSize triggers a call to resized. Commenting out 
  // some code there lets us also get rid of the leak. It's the call to
  //
  //   delayLineModuleEditor->setBounds() 
  //
  // at the bottom of resized -> check
  //
  //   EchoLabDelayLineModuleEditor::resized()
  //
  // it's in the calls to 
  //
  //  inputEqualizerEditor   ->setBounds(x,                 y, w+2, h);
  //  feedbackEqualizerEditor->setBounds(x+w+middleWidth-2, y, w+2, h);
  //
  // check
  //
  //   EqualizerModuleEditor::resized()
  //
  // ...

}


/*


ToDo:

-Write a unit test that: 
 -Randomizes all parameters while remembering their new values
 -Retrieves the state as xml
 -Randomizes the parameters again
 -Recalls the state from the xml
 -Checks the recalled paramters against the remembered ones
 -Retrieves the state again
 -Compares the 2nd xml to the 1st
 That unit test should expose the not recalled NumVoices parameter - and maybe more.


*/

