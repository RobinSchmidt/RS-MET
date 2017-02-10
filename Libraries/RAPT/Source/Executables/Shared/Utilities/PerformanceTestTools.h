#ifndef PERFORMANCETESTTOOLS_INCLUDED
#define PERFORMANCETESTTOOLS_INCLUDED

#include <iostream>

#ifdef _MSC_VER
#include <intrin.h> // Or #include <ia32intrin.h> etc.
#endif

void printPerformanceTestResult(const char *testName, double numCyclesPerOperation);

/** This class implements a CPU cycle counter which is useful for performance tests. It does not 
literally count but instead measure the time stamps before and after a sequence of commands. Works 
only with MSVC at the moment. */

class ProcessorCycleCounter  
{

public:

  /** Resets the counter to zero. */
  inline void init()
  {  
    initTime = ReadTSC();
  }

  /** Returns the number of CPU-cycles since the last call to init(). */
  inline __int64 getNumCyclesSinceInit()
  {
    return ReadTSC() - initTime;
  }

protected:

  __int64 initTime;    // remembers the time at which the call to init() occured 

  /** Function for reading the time-stamp counter - this is taken from Agner Fog's optimization 
  tutorials (check there for a new version, maybe a version fo GCC). */
  __int64 ReadTSC()  // Returns time stamp counter
  {
#ifdef _MSC_VER
    int     dummy[4];      // For unused returns
    __int64 clock;         // Time 
    __cpuid(dummy, 0);     // Serialize
    clock = __rdtsc();     // Read time
    __cpuid(dummy, 0);     // Serialize again
    return clock;
#else
    return 0;
#endif

    // maybe convert the ifdefs into ifs: if( isMicrosoftCompiler() )
  }

};







#endif 