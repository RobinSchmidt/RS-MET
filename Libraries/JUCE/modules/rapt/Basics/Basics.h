#ifndef RAPT_BASICS_H_INCLUDED
#define RAPT_BASICS_H_INCLUDED

namespace RAPT
{

#include "Constants.h"
#include "TypeDefinitions.h"
#include "MacroDefinitions.h"
//#include "MathBasics.h"  // add from RSLib
#include "SortAndSearch.h"

// move to some other file (BasicFunctions or something):
template<class T>
inline void rsSwap(T& x, T& y)
{
  T t = x;
  x = y;
  y = t;
}

/** This function should be used to indicate a runtime error. */
inline void rsError(const char *errorMessage = nullptr)
{
  //printf("%s", errorMessage);
  RS_DEBUG_BREAK;
  // \todo have some conditional compilation code based on the DEBUG macro (trigger a break),
  // maybe open an error message box, etc.
}

/** This function should be used for runtime assertions. */
inline void rsAssert(bool expression, const char *errorMessage = nullptr)
{
  if( expression == false )
    rsError(errorMessage);
}

inline void rsAssertFalse() { rsAssert(false); }


// abs, sign, floor, ceil
template<class T> inline T rsSqrt(T x) { return std::sqrt(x); }
template<class T> inline T rsExp( T x) { return std::exp( x); }
template<class T> inline T rsLog( T x) { return std::log( x); }
template<class T> inline T rsSin( T x) { return std::sin( x); }
template<class T> inline T rsCos( T x) { return std::cos( x); }
template<class T> inline T rsTan( T x) { return std::tan( x); }





}

#endif