#include "MiscPerformanceTests.h"
using namespace RAPT;

void testFunction(double x)
{
  dontOptimize(x);
}

class TestFunctor
{
public:
  void operator()(double x)
  {
    dontOptimize(x);
  }
};

class TestCalleeClass
{
public:
  void callbackTarget(double x)
  {
    dontOptimize(x);
  }
};

void callbackPerformance()
{
  // compare the performance and memory occupation of various types of function-pointers: 
  // raw function pointer, member function pointer, functor, std::function with lambda
  // std::function wrapping a member function call, ...can we actually assign a member function
  // callback to a std::function? ...try it

  int numCalls = 10000;  // number of function calls
  ProcessorCycleCounter counter;
  int n;
  double x = 1.0;

  // raw function pointer:
  void (*functionPointer) (double x);
  functionPointer = testFunction;
  counter.init(); 
  for(n = 0; n < numCalls; n++) functionPointer(x);
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("raw function pointer", cycles/numCalls);
  print("num bytes", sizeof(functionPointer));

  // functor:
  TestFunctor functor;
  counter.init(); 
  for(n = 0; n < numCalls; n++) functor(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("functor", cycles/numCalls);
  print("num bytes", sizeof(functor));

  // member function pointer:
  TestCalleeClass calleeObject;
  SpecificMemberFunctionCallback1<TestCalleeClass, void, double>
     memberFuncPtr(&calleeObject, &TestCalleeClass::callbackTarget);
  counter.init(); 
  for(n = 0; n < numCalls; n++) memberFuncPtr.call(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("member function pointer", cycles/numCalls);
  print("num bytes", sizeof(functor));

  // std::function with lambda:
  std::function<void(double)> f;
  f = [] (double x) { dontOptimize(x); };
  counter.init(); 
  for(n = 0; n < numCalls; n++) f(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("std::function with lambda", cycles/numCalls);
  print("num bytes", sizeof(f));

  // std::function with function pointer:
  f = testFunction;
  counter.init(); 
  for(n = 0; n < numCalls; n++) f(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("std::function with pointer", cycles/numCalls);
  print("num bytes", sizeof(f));

  // std::function with functor:
  f = functor;
  counter.init(); 
  for(n = 0; n < numCalls; n++) f(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("std::function with functor", cycles/numCalls);
  print("num bytes", sizeof(f));

  // std::function with member function pointer:
  f = memberFuncPtr;
  counter.init(); 
  for(n = 0; n < numCalls; n++) f(x);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("std::function with member pointer", cycles/numCalls);
  print("num bytes", sizeof(f));


  // results:                  cycles  bytes
  // function pointer:         23       8
  // functor:                   4       1
  // member function pointer:   8       1
  // std::function, lambda:    11      64
  // std::function, pointer:   26      64
  // std::function, functor:   11      64
  // std::function, member:    11      64
}

