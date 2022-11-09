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
  jura::SamplerModule sampler(&pluginLock, &metaManager);


  jura::SamplerEditor* editor = dynamic_cast<jura::SamplerEditor*> (sampler.createEditor(0));

  //expectNotEquals((void*) editor, nullptr, "Failed to create jura::SamplerEditor");
  // ...that doesn't compile - let's do a workaround:
  bool ok = editor != nullptr; 
  //expectEquals(ok, true);    // this doesn't compile either
  //expectEquals(true, true);  // this neither - seems like we can't use bool for expectEquals?
  //int i_ok = (int)ok;
  expectEquals((int)ok, 1);       // OK - this finally compiles, but it's ugly!




  int dummy = 0;
}

// 