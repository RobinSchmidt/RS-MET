#include "../JuceLibraryCode/JuceHeader.h"

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  //jura::FileManager::companyName = "Soundemote"; // test
  // If I remember correctly, this will determine the root folder for all the support files for all
  // products of given company.


  // Uncomment the relevant line of code for which AudioModule you want to get built as plugin. 
  // Normally, we want to build ToolChain but by just changing, which line is not commented, we can
  // choose at compile time any other AudioModule:

  jura::ToolChain *dummy = nullptr; return createPluginWithMidi(dummy, 10);

  //jura::PhaseScope *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhaseScope2 *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhaseScopeMultiColor *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhasorFilter *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Ladder       *dummy = nullptr; return createPluginWithMidi(dummy);
  //jura::Enveloper    *dummy = nullptr; return createPluginWithMidi(dummy);
  //jura::SamplerModule *dummy = nullptr; return createPluginWithMidi(dummy, 10);

  //// just for testing the resizing constraints:
  //jura::AudioPluginWithMidiIn *plugIn = new jura::AudioPluginWithMidiIn(10);
  ////plugIn->setEditorSizeLimits(400, 300, 800, 600);  // test constrained size
  //plugIn->setEditorSizeLimits(400, 300, 400, 300);    // test fixed size
  //jura::ToolChain *module = new jura::ToolChain(&plugIn->plugInLock, &plugIn->metaParaManager);
  //module->setSaveAndRecallMetaParameters(true);
  //plugIn->setAudioModuleToWrap(module);
  //return plugIn;
}
