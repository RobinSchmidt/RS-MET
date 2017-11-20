#ifndef RS_TESTUTILITIES_H
#define RS_TESTUTILITIES_H

#include "../Common/Prototypes.h"

bool detectMemoryLeaks();  // currently works only in MSVC

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
rsVectorDbl rsLinearRangeVector(     int N, double min, double max);
rsVectorDbl rsExponentialRangeVector(int N, double min, double max);
rsVectorDbl rsRandomVector(          int N, double min, double max, int seed = 0);

// conversions to std::string:
std::string toString(int n);

// returns x^2 = x*x, useful for testing application of a unary function using a function pointer
//double rsSquare(double x);

// 
template<class T>
T square(T x)
{
  return x*x;
}


#endif
