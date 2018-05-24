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
should take an integer parameter for (a measure of) the input size. For example, each of these 
objects could perform the same task (matrix-multiply, FFT, filtering, ...) but use a different 
algorithm/implementation internally.

It runs the tests multiple times (how many is something you can set up), gathers the data and may 
perform statistical analysis afterwards. In the simplest case, you might be interested in how 
minimum, average or maximum running times of the various functions differ. On a more 
sophisticated level, you might be interested in variances and statistical significance of 
differences in the averages. This class is meant to do the appropriate statistical tests 

todo:
-computations for the statiscal tests should be factored out
-it should be possible to make plots that show how the performance depends on input size - with
 error-bars and multiple graphs   */ 

class PerformanceAnalyzer
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Adds a new test to be performed. You have to pass a name that will be used to identify the 
  test in analysis/comparisons - for example something like "FFT, radix-2, decimation in time". 
  The test function itself should then run one such FFT with a length given by its integer 
  parameter. */
  void addTest(std::function<void(int)>* testFunc, const std::string& testName)
  {
    tests.push_back(testFunc);
    names.push_back(testName);
  }

  /** You can pass an array of inputs sizes to the algorithms to be tested. */
  void setTestInputSizes(const std::vector<int>& newInputSizes)
  {
    inputSizes = newInputSizes;
  }

  /** Sets the number of times that each test should be run for each input size. */
  void setNumRunsPerTest(int numberOfRuns) { numRuns = numberOfRuns; }
    // maybe we should allow to set separate numbers of runs for different input sizes (for large
    // input sizes, we may want to run the test less often because it takes too mauch time to run
    // multiple times)

  /** Sets up, how large a deviation form the median is allowed before a datapoint is considered
  an outlier with invalid value and thrown away. For example, passing 0.5 and 3.0 will mean that 
  values below 0.5 times the median and above 3.0 times the median will be cosidered invalid 
  outliers. */
  void setOutlierRemovalLevels(double smallOutlierRatio, double largeOutlierRatio)
  {
    smallOutlier = smallOutlierRatio;
    largeOutlier = largeOutlierRatio;
  }

  /** Clears the array of test functions. */
  void clearTests() { tests.clear(); names.clear(); }

  /** Clears the array of input sizes. */
  void clearInputSizes() { inputSizes.clear(); }

  /** Puts the object back into its initial state, i.e. clears all arrays of test-functors, 
  input sizes, etc. */
  void init();

  //-----------------------------------------------------------------------------------------------
  // \name Measurement:

  /** Runs all tests. Needs to be called after setting up the object and before retrieving the test 
  results. */
  void runTests();


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry:

  /** Creates a report string that summarizes the results. */
  std::string getReport();

  // void plotResults() 
  // may be better handled by another, higher level class

  /** Returns the number of valid datapoints that were obtained for the test with given index and
  input size index. Due to removal of outliers, this value may be different from the value passed 
  to setNumRunsPerTest. */
  int getNumValidDatapoints(int testIndex, int sizeIndex);

  /** Returns the raw test results as nested vector (i.e. 3D array) where the first index is for 
  the test-function (i.e. the algorithm among a family of algorithms) and the second is for the 
  input-size and the third index is the index of the actual datapoint, like
  result[i][j][k] is the result of the k-th run of algorithm i with input size (indexed by) j */
  std::vector<std::vector<std::vector<double>>> getRawData() { return rawData; }

  /** Returns the mean values of the test results. */
  std::vector<std::vector<double>> getMeans();

  /** Returns the variances of the test results. */
  std::vector<std::vector<double>> getVariances();


protected:

  //void runTest(

  /** Resizes the rawData 3D-array to what is required to hold all the results. */
  void initResultArray();
    // ...maybe intialize to all zeros - but that's actually not necessarry - it will be soon 
    // after filled with valid data anyway

  std::vector<std::string> names;               // names of the tests
  std::vector<std::function<void(int)>*> tests; // tests themselves
  std::vector<int> inputSizes;


  int numRuns = 30; // later use an array to allow smaller numRuns for large input sizes
  double smallOutlier = -std::numeric_limits<double>::infinity();
  double largeOutlier =  std::numeric_limits<double>::infinity();


  std::vector<std::vector<std::vector<double>>> rawData;      
    // raw results, 1st index: function, 2nd: input-size, 3rd: datapoint

  std::vector<std::vector<double>> means, variances;

  ProcessorCycleCounter cpuCounter;

};
// test with 2 functions (like Ooura FFT, rsFFT), 5 input sizes (128...2048), 30 runs

#endif 