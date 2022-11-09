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




  int dummy = 0;
}

// 