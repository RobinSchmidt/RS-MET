#ifdef ROSIC_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using. */
#error "Incorrect use of JUCE cpp file"
#endif

#include "rosic.h"

/** The cpp files are included in the order in which they depend on each other. ToDo: reorder them 
in the library accordingly, such that files in one library folder depend only on other files in 
folders that are considered "above" in the hierarchy. */

#include "basics/GlobalFunctions.cpp"
#include "basics/rosic_ChannelMatrix2x2.cpp"
#include "basics/rosic_Constants.cpp"                // empty
#include "basics/rosic_FunctionTemplates.cpp"        // empty
#include "basics/rosic_HelperFunctions.cpp"
#include "basics/rosic_Interpolator.cpp"
#include "basics/rosic_NumberManipulations.cpp"      // empty
#include "infrastructure/rosic_MutexLock.cpp"        // used by sample buffer - move to basis (and/or SampleBuffer elsewhere)
#include "basics/rosic_SampleBuffer.cpp"
#include "basics/rosic_SamplePlaybackParameters.cpp"

//#include "basics/rosic_TabulatedFunction.cpp" // needs ExpressionEvaluator


// analysis:
#include "analysis/rosic_CyclicAutoCorrelator.cpp"  // no dependencies
//#include "analysis/rosic_EnvelopeFollower.cpp"
//#include "analysis/rosic_FormantPreserver.cpp"
//#include "analysis/rosic_FormantRemover.cpp"
//#include "analysis/rosic_InstantaneousEnvelopeDetector.cpp"
//#include "analysis/rosic_LevelDetector.cpp"


