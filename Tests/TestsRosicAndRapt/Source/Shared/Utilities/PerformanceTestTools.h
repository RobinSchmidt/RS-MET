#ifndef PERFORMANCETESTTOOLS_INCLUDED
#define PERFORMANCETESTTOOLS_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <functional>

#ifdef _MSC_VER
#include <intrin.h> // Or #include <ia32intrin.h> etc.
#include <WTypesbase.h> 
#endif

void printPerformanceTestResult(const char *testName, double numCyclesPerOperation);

void printPerformanceTestResult(const std::string& testName, double numCyclesPerOperation);


/** Prints the memory occupation of the given object in bytes. */
template<class T>
void printMemoryOccupation(const char *name, T& object);

void print(const char *name, int x);

/** Prevents the compiler from optimizing away the variable x. This is useful for performance 
tests which could give completetly nonsensical results, if the compiler can detect that the results 
of the operations, you want to measure, are not used anywhere. Just pass your results to this 
function and the compiler can't optimze away your operations anymore. If the datatype of your 
result does not support to be passed to cout via <<, you can just pass a pointer. */
template<class T>
inline void dontOptimize(T x)
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



/** (...not yet finished - just a stub). 

Class to perform performance analysis. You can hand it a bunch of std::function objects, each 
of which is supposed to perform a certain task whose performance you want to measure and it gives
you back a string with a report of the test results. The std::function objects you need to pass
should take an integer parameter for (a measure of) the input size and return a double that is a 
measure of the performance (like the number of CPU cycles taken (per datapoint, audio-sample, 
whatever). For example, each of these objects could perform the same task (matrix-multiply, FFT, 
filtering, ...) but use a different algorithm/implementation internally.

It runs the tests multiple times (how many is something you can set up), gathers the data and may 
perform statistical analysis afterwards. In the simplest case, you might be interested in how 
minimum, average or maximum running times of the various functions differ. On a more 
sophisticated level, you might be interested in variances and statistical significance of 
differences in the averages. This class is meant to do the appropriate statistical tests 

todo:
-computations for the statiscal tests should be factored out
-it should be possible to make plots that show how the performance depends on input size - with
 error-bars and multiple graphs
 
*/ 

class PerformanceAnalyzer
{

public:


protected:


  std::vector<string> names;                     // names of the tests
  std::vector<std::function<double(int)>> tests; // tests themselves

};




#endif 