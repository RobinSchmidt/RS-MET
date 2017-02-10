#include "../JuceLibraryCode/JuceHeader.h"

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  //jura::PhaseScope *dummy = nullptr; return createPluginWithoutMidi(dummy);
  jura::PhaseScope2 *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhaseScopeMultiColor *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::PhasorFilter *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Ladder       *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Enveloper    *dummy = nullptr; return createPluginWithMidi(dummy);


  //jura::ModuleChainer *dummy = nullptr; return createPluginWithoutMidi(dummy);
}
