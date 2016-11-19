#include "../JuceLibraryCode/JuceHeader.h"

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  //jura::Ladder    *dummy = nullptr; return createPluginWithoutMidi(dummy);
  jura::Enveloper *dummy = nullptr; return createPluginWithMidi(dummy);
}
