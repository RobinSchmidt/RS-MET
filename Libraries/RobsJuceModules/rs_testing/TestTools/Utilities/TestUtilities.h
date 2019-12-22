#ifndef RS_TESTUTILITIES_H
#define RS_TESTUTILITIES_H

// todo: merge file with other utility files

//old:
//#include "../Common/Prototypes.h"

// new:
//#include "../RaptLibraryCode/RaptInstantiations.h"
#include "rosic/rosic.h"
#include "rs_testing/rs_testing.h"



bool runUnitTest(bool (*test)(), const std::string& name);



//bool detectMemoryLeaks();  // currently works only in MSVC

/** This function should be called on program startup when automatic detection of memory leaks
should be turned on. */
inline void checkForMemoryLeaksOnExit()
{
#if defined _MSC_VER
  int tmpFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG); // gets the current flag
  tmpFlag |= _CRTDBG_LEAK_CHECK_DF;                  // turns on leak checking
  //tmpFlag &= ~_CRTDBG_CHECK_CRT_DF;                  // turns off CRT block checking bit
  _CrtSetDbgFlag(tmpFlag);                           // set flag to the new value;
#endif
}

// helper functions to create some vectors useful for testing purposes (maybe move them to
// somewhere else):
std::vector<double> rsLinearRangeVector(     int N, double min, double max);
std::vector<double> rsExponentialRangeVector(int N, double min, double max);
std::vector<double> rsRandomVector(          int N, double min, double max, int seed = 0);
std::vector<double> rsApplyFunction(const std::vector<double>& v, double p,
  double (*f) (double, double));

// conversions to std::string:
std::string toString(int n);

// replace with own prng:
inline double random(double min, double max)
{
  double tmp = (1.0/RAND_MAX) * rand();  // between 0...1
  return RAPT::rsLinToLin(tmp, 0.0, 1.0, min, max);
}

// returns x^2 = x*x, useful for testing application of a unary function using a function pointer
//double rsSquare(double x);

//
template<class T>
T square(T x)
{
  return x*x;
}

template<class T>
void rsFillWithComplexRandomValues(std::complex<T>* x, size_t N, T min, T max,
  unsigned long seed = 0)
{
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(min, max);
  prng.setSeed(seed);
  for(size_t n = 0; n < N; n++)
    x[n] = std::complex<T>(prng.getSample(), prng.getSample());
}

template<class T> // convenience function for std::vector
void rsFillWithComplexRandomValues(std::vector<std::complex<T>>& x, T min, T max,
  unsigned long seed = 0)
{
  rsFillWithComplexRandomValues(&x[0], x.size(), min, max, seed);
}
template<class T> // yet more convenient function
std::vector<std::complex<T>> rsComplexRandomVector(int N, T min, T max, unsigned long seed = 0)
{
  std::vector<std::complex<T>> x(N);
  rsFillWithComplexRandomValues(x, min, max, seed);
  return x;
}


template<class T>
T rsMaxComplexError(std::complex<T>* target, std::complex<T>* actual, size_t N)
{
  T maxErr = T(0);
  for(size_t n = 0; n < N; n++)
    maxErr = RAPT::rsMax(maxErr, abs(target[n]-actual[n]));
  return maxErr;
}

template<class T>
bool rsAlmostEqual(std::vector<std::complex<T>>& x, std::vector<std::complex<T>>& y, T tolerance)
{
  RAPT::rsAssert(x.size() == y.size());
  T maxErr = rsMaxComplexError(&x[0], &y[0], x.size());
  return maxErr <= tolerance;
}


/** Applies the inner function to the value x and then the outer function to the result of that
inner function and returns the final result. This is known as function composition in
mathematics. */
template<class T, class F1, class F2>
T applyComposedFunction(T x, F1 innerFunction, F2 outerFunction)
{
  return outerFunction(innerFunction(x));
}

/** Returns true, if the function f maps the given argument x to itself. */
template<class T, class F>
bool mapsToItself(T x, F f)
{
  return x == f(x);
}

/** Checks, if the 2nd function is the inverse function of the first for the given input argument
x. */
template<class T, class F1, class F2>
bool mapsBack(T x, F1 forwardFunction, F2 maybeInverseFunction)
{
  return x == applyComposedFunction(x, forwardFunction, maybeInverseFunction);
  //return x == maybeInverseFunction(forwardFunction(x));
}
// maybe rename to isFunctionLocallyInverse, isFunctionInverseAt

/** Checks, if the 2nd function is the inverse function of the first for a given range of input
arguments between minValue and maxValue with given. */
template<class T, class F1, class F2>
bool isInverseFunction(F1 forwardFunc, F2 maybeInverseFunc, T minValue, T maxValue, T increment)
{
  T value = minValue;
  while(value < maxValue) {
    if( !mapsBack(value, forwardFunc, maybeInverseFunc) )
      return false;
    value += increment;
  }
  return true;
}





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

// get rid of that:
inline void appendTestResultToReport(std::string &reportString, const std::string &nameOfTest,
  bool result)
{
  if( result == true )
    reportString += nameOfTest + ": OK \n";
  else
    reportString += nameOfTest + ": !!! FAILED !!!\n";
}

/** Comparison function that compares with a given error tolerance and also returns true when the
involved numbers are NaNs or infinities. */
bool areNumbersEqual(double x, double y, double relativeTolerance);

/** Convenience function to convert a string to a window-type.  options: rc,hn,hm,bm,bh */
RAPT::rsWindowFunction::WindowType stringToWindowType(const std::string& wt);



/** Experimental - goal: resemble numpy/scipy/matplotlib functionality, so we may easily port such 
code to C++. */

template<class T>
class rsNumPy
{
public:

  std::vector<T> linspace(T min, T max, int N) { return RAPT::rsRangeLinear(T(min), T(max), N); }

  std::vector<T> sin(const std::vector<T>& x)
  {
    std::vector<T> y(x.size());
    RAPT::rsArrayTools::applyFunction(&x[0], &y[0], (int)x.size(), &::sin);
    return y;
  }

  T pi = PI;




};
// move this into the research repo
// https://docs.scipy.org/doc/numpy-1.11.0/numpy-ref-1.11.0.pdf





#endif
