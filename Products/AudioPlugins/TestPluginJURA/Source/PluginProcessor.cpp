#include "../JuceLibraryCode/JuceHeader.h"

// This call is all that is needed to create the actual plugin. To make other plugins later, we
// will just have to create several of such calls that return different classes of objects and we 
// can switch between building one or another plugin by merely commenting out a line of code and
// uncommenting another. So, we can use this single project to build many plugins - that's very
// convenient and maintenance friendly. :-) ...OK - we may have to manually rename the .dll after 
// the build but that's better than having to maintain a dozen of projects...

//AudioProcessor* JUCE_CALLTYPE createPluginFilter() { return new jura::Ladder(); }

AudioProcessor* JUCE_CALLTYPE createPluginFilter() 
{ 
  return new jura::AudioPlugin(new jura::Ladder());

  //return new jura::Ladder();  // old - before getting rid of AudioProcessor inheritance
}

