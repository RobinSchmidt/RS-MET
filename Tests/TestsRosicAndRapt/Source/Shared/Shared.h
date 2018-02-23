#ifndef RAPT_EXECUTABLES_SHARED_H
#define RAPT_EXECUTABLES_SHARED_H

#include "Plotting/Plotting.h"

//#include "RaptLibraryCode/RaptTypedefInstantiations.h"
#include "Utilities/TestInputCreation.h"
#include "Utilities/FileWriting.h"
#include "Utilities/PerformanceTestTools.h"

#include "Prototypes/LindenmayerSystem.h"

//#include "rosic/rosic.h"

// for testing the callback performance (this is actually in jura, but anyway):
#define JUCE_API
#include "jura_framework/control/jura_Callbacks.h"

inline bool detectMemoryLeaks()
{
#ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
#else
  return false;
#endif
}


#endif