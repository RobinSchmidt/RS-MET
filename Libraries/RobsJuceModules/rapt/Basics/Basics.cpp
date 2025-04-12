#include "Basics.h"

//#if !defined PLOTTING_COMPILED
//
//#include "Plotting.cpp" // inlcudes GNUPlotter.cpp, should not be defined in RAPT namspace
//#define PLOTTING_COMPILED
//
//#endif


namespace RAPT
{

#include "DebugTools.cpp"
#include "SortAndSearch.cpp"

}

//=================================================================================================
/*

ToDo:

- Maybe make a subfolder "DevTools" into which we then move DebugTools and the Plotting stuff. 
  Maybe DebugTools can be re-organized a bit, too. Maybe have files: Printing.h, Plotting.h, 
  Logging.h, ErrorHandling.h. There, we could also have the heap allocation logging stuff which is
  currently in the research repo in \Projects\CppExperiments. It could also contain stuff for 
  instrumentation, benchmarking, etc. (although the latter may better be placed outside the 
  library). We could also have a logging allocator there.

*/