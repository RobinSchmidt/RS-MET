#ifndef RAPT_BASICS_H_INCLUDED
#define RAPT_BASICS_H_INCLUDED

namespace RAPT
{

#include "Constants.h"
#include "TypeDefinitions.h"
#include "MacroDefinitions.h"
#include "DebugTools.h"
//#include "MathBasics.h"  // add from RSLib
#include "SortAndSearch.h"

// move to some other file (BasicFunctions or something):

/** Swaps x and y via copy assignment operator. */
template<class T>
inline void rsSwapNaive(T& x, T& y)
{
  T t = x;
  x = y;
  y = t;
}

/** Swaps x and y via std::move. */
template<class T>
inline void rsSwapViaMove(T& x, T& y)
{
  T t(std::move(x));
  x = std::move(y);
  y = std::move(t);
}

/** Swaps x and y. Uses move assignment, i.e. calls rsSwapViaMove. */
template<class T>
inline void rsSwap(T& x, T& y)
{
  rsSwapViaMove(x, y);
}


template<class T> inline T rsUnityValue(T /*value*/) { return T(1);  }
template<class T> inline T rsZeroValue( T /*value*/) { return T(0);  }
template<class T> inline T rsIdentity(  T value) { return value; }

/** Returns true, if x is not-a-number, false otherwise. */
template<class T> inline bool rsIsNaN(T x)
{
  return x != x; // NaN is the only value that returns false for this comparison
}

/** Returns true, if x is plus or minus infinity, false otherwise. */ 
template<class T> inline bool rsIsInfinite(T x)
{
  return x == std::numeric_limits<T>::infinity() || x == -std::numeric_limits<T>::infinity();
}

/** Returns true, if x is a finite number, i.e. not NaN and not +-infinity, false otherwise. */
template<class T> inline bool rsIsFiniteNumber(T x)
{
  if(rsIsNaN(x) || rsIsInfinite(x))
    return false;
  return true;
}

template<class T> inline bool rsIsFiniteNonNegativeNumber(T x)
{
  return rsIsFiniteNumber(x) && x >= T(0);
}


/** Returns true, if the passed array of length N contains only finite numbers, false otherwise, */
template<class T> inline bool rsIsFiniteNumbers(T* x, int N)
{
  for(int i = 0; i < N; i++)
    if(!rsIsFiniteNumber(x[i]))
      return false;
  return true;
}


// abs, sign, floor, ceil
template<class T> inline T rsSqrt(T x) { return std::sqrt(x); }
template<class T> inline T rsExp( T x) { return std::exp( x); }
template<class T> inline T rsLog( T x) { return std::log( x); }
template<class T> inline T rsSin( T x) { return std::sin( x); }
template<class T> inline T rsCos( T x) { return std::cos( x); }
template<class T> inline T rsTan( T x) { return std::tan( x); }

template<class T> inline int rsFloorInt(T x) { return (int) floor(x); }
template<class T> inline int rsCeilInt( T x) { return (int) ceil(x);  }

/** Returns the absolute value of the input argument. It is intended to replace the standard
"abs" and "fabs" c-functions where genericity is desired. */
template <class T>
T rsAbs(T x)
{
  if( x < rsZeroValue(x) )
    return -x;
  else
    return  x;
}
inline double  rsAbs(double  x) { return fabs(x); }
inline float   rsAbs(float   x) { return fabs(x); }
//inline rsInt8  rsAbs(rsInt8  x) { return  abs(x); }
//inline rsInt16 rsAbs(rsInt16 x) { return  abs(x); }
inline rsInt32 rsAbs(rsInt32 x) { return  abs(x); }
//inline rsInt64 rsAbs(rsInt64 x) { return  abs(x); } // doesn't work with MinGW gcc 4.7
template<class T> inline T rsAbs(std::complex<T> z) { return abs(z); }


/** Squared absolute value of a complex number. */
template<class T> 
T rsAbsSquared(const std::complex<T>& z)
{
  return z.real()*z.real() + z.imag()*z.imag(); // == conj(z) * z
}

/*
template <class T>
bool rsGreater(const T& left, const T& right)
{
  return left > right;
}
*/



}

#endif