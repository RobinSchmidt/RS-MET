#include "../JuceLibraryCode/JuceHeader.h"

// This call is all that is needed to create the actual plugin. To make other plugins later, we
// will just have to create several of such calls that return different classes of objects and we 
// can switch between building one or another plugin by merely commenting out a line of code and
// uncommenting another. So, we can use this single project to build many plugins - that's very
// convenient and maintenance friendly. :-) ...OK - we may have to manually rename the .dll after 
// the build but that's better than having to maintain a dozen of projects...

AudioProcessor* JUCE_CALLTYPE createPluginFilter() { return new jura::Ladder(); }

// hmm...but i don't actually want to derive my audiomodules from juce::AudioProcessor because that
// craetes an overhead for modules that are not supposed to built as plugin on their own. Maybe 
// it's neverthe possible to write a similar one-liner, like
// return wrapIntoPlugin(new jura::Ladder); or something, we'll see
