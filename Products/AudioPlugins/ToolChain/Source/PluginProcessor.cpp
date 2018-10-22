#include "../JuceLibraryCode/JuceHeader.h"

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  //jura::FileManager::companyName = "Soundemote";

  //jura::PhaseScope *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhaseScope2 *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhaseScopeMultiColor *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhasorFilter *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Ladder       *dummy = nullptr; return createPluginWithMidi(dummy);
  //jura::Enveloper    *dummy = nullptr; return createPluginWithMidi(dummy);

  //jura::ToolChain *dummy = nullptr; return createPluginWithMidi(dummy, 10);

  // just for testing the resizing constrainer:
  jura::AudioPluginWithMidiIn *plugIn = new jura::AudioPluginWithMidiIn(10);

  // this doesn't work:
  //juce::ComponentBoundsConstrainer* bc = new juce::ComponentBoundsConstrainer;
  //bc->setMaximumSize(600, 800);
  //bc->setMinimumSize(300, 400);
  //plugIn->setEditorBoundsConstrainer(bc);

  jura::ToolChain *module = new jura::ToolChain(&plugIn->plugInLock, &plugIn->metaParaManager);
  module->setSaveAndRecallMetaParameters(true);
  plugIn->setAudioModuleToWrap(module);
  return plugIn;
}
