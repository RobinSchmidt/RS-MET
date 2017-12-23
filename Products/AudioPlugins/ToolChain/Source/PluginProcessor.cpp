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

  jura::ToolChain *dummy = nullptr; return createPluginWithMidi(dummy, 10);
}
