#ifndef RAPT_EXECUTABLES_SHARED_H
#define RAPT_EXECUTABLES_SHARED_H

#include "Plotting/Plotting.h"

//#include "RaptLibraryCode/RaptTypedefInstantiations.h"
#include "Utilities/TestInputCreation.h"
#include "Utilities/FileWriting.h"
#include "Utilities/PerformanceTestTools.h"

//#include "rosic/rosic.h"

inline bool detectMemoryLeaks()
{
#ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
#else
  return false;
#endif
}


#endif