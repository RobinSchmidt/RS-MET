#ifndef PERFORMANCETESTTOOLS_INCLUDED
#define PERFORMANCETESTTOOLS_INCLUDED

#include <iostream>

#ifdef _MSC_VER
#include <intrin.h> // Or #include <ia32intrin.h> etc.
#include <WTypesbase.h> 
#endif

void printPerformanceTestResult(const char *testName, double numCyclesPerOperation);

/** Prevents the compiler from optimizing away the variable x. This is useful for performance 
tests which could give completetly nonsensical results, if the compiler can detect that the results 
of the operations, you want to measure, are not used anywhere. Just pass your results to this 
function and the compiler can't optimze away your operations anymore. If the datatype of your 
result does not support to be passed to cout via <<, you can just pass a pointer. */
template<class T>
inline void preventOptimization(T x)
{
  volatile bool f = false;
  if(f)
    cout << x;
}

/** This class implements a CPU cycle counter which is useful for performance tests. It does not 
literally count but instead measure the time stamps before and after a sequence of commands. Works 
only with MSVC at the moment. 

hmm...maybe we should update this using QueryPerformanceCounter, etc. - see here:
https://en.wikipedia.org/wiki/Time_Stamp_Counter
https://msdn.microsoft.com/en-us/library/windows/desktop/dn553408(v=vs.85).aspx

*/

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

//=================================================================================================

/** A cycle counter based on 
https://msdn.microsoft.com/en-us/library/windows/desktop/dn553408(v=vs.85).aspx */

class ProcessorCycleCounter2
{

public:

  inline void init()
  {  
    QueryPerformanceCounter(&startTime);
  }

  inline LONGLONG getNumCyclesSinceInit()
  {
    QueryPerformanceCounter(&endTime);
    return endTime.QuadPart - startTime.QuadPart;
  }

protected:

  LARGE_INTEGER startTime, endTime;

};





#endif 