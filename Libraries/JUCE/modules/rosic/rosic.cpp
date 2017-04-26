#ifdef ROSIC_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using. */
#error "Incorrect use of JUCE cpp file"
#endif

#include "rosic.h"

/** ToDo: include the cpp files in the order in which they depend on each other */

// analysis:
#include "analysis/rosic_CyclicAutoCorrelator.cpp"
//#include "analysis/rosic_EnvelopeFollower.cpp"
//#include "analysis/rosic_FormantPreserver.cpp"
//#include "analysis/rosic_FormantRemover.cpp"
//#include "analysis/rosic_InstantaneousEnvelopeDetector.cpp"
//#include "analysis/rosic_LevelDetector.cpp"


