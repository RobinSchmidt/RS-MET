#include "../JuceLibraryCode/JuceHeader.h"

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  jura::PhasorFilter *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Ladder       *dummy = nullptr; return createPluginWithoutMidi(dummy);
  //jura::Enveloper    *dummy = nullptr; return createPluginWithMidi(dummy);
}
