#include "MiscPerformanceTests.h"
using namespace RAPT;

void callbackPerformance()
{
  // compare the performance and memory occupation of various types of function-pointers: 
  // raw function pointer, member function pointer, functor, std::function with lambda
  // std::function wrapping a member function call, ...can we actually assign a member function
  // callback to a std::function? ...try it

  int numCalls = 1000;  // number of function calls
  ProcessorCycleCounter counter;
  int n;

  counter.init(); 
  for(n = 0; n < numCalls; n++)
  {
    // call....
  }

  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("raw function pointer", cycles/numCalls);

}

