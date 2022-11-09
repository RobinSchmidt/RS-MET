#include "UnitTestsSampler.h"
//using namespace juce;
//using namespace jura;

void UnitTestsSampler::runTest()
{
  testSamplerAudioModule();
}


void UnitTestsSampler::testSamplerAudioModule()
{
  // Create an instance of jura::SamplerModule:
  juce::CriticalSection pluginLock;        // The AudioModule constructor expects pointers to...
  jura::MetaParameterManager metaManager;  // ...such objects
  SamplerModuleTest sampler(&pluginLock, &metaManager);


  SamplerEditorTest* editor = dynamic_cast<SamplerEditorTest*> (sampler.createEditor(0));

  //expectNotEquals((void*) editor, nullptr, "Failed to create jura::SamplerEditor");
  // ...that doesn't compile - let's do a workaround:
  bool ok = editor != nullptr; 
  //expectEquals(ok, true);    // this doesn't compile either
  //expectEquals(true, true);  // this neither - seems like we can't use bool for expectEquals?
  //int i_ok = (int)ok;
  expectEquals((int)ok, 1);       // OK - this finally compiles, but it's ugly!



  ok &= editor->testTreeViewNodeSelection(); expectEquals((int)ok, 1); 



  delete editor;  // clean up memory
}


jura::AudioModuleEditor* SamplerModuleTest::createEditor(int type)
{
  return new SamplerEditorTest(this);
}

bool SamplerEditorTest::testTreeViewNodeSelection()
{
  bool ok = true;

  // ToDo:
  // -Create an example sfz-string:
  //  -One group containign one region
  //  -No sample opcode
  //  -Have opcodes for cutoff, lfo1_freq, lfo1_cutoff
  // -Let the sampler set itself up from that sfz string
  // -Simulate clicking on the cutoff node in the TreeView
  // -Check, if the slider appears as it should
  // -Do the same for the other sliders for lfo1_freq, lfo1_cutoff

  juce::String sfzString = "<group><region>cutoff=1000 lfo1_freq=2.5 lfo1_cutoff=500";

  return ok;
}



/*

ToDo:
-Make subclasses SamplerModuleTest/SamplerEditorTest of SamplerModule/SamplerEditor here
-In SamplerModuleTest, override createEditor to return an object of class SamplerEditorTest
-In SamplerEditorTest have additional functions that test certain functionality. These tests need
 to access the (protected) widgets - that's hwy we need to create a subclass.
 
-Test selecting am opcode node in the TreeView and check, if the correct slider is displayed


*/