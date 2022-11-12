#ifndef RAPT_BASICFUNCTIONS_H_INCLUDED
#define RAPT_BASICFUNCTIONS_H_INCLUDED

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
  //rsSwapNaive(x, y);  // temporary, for debug
  rsSwapViaMove(x, y);
}


template<class T> inline T rsUnityValue(T /*value*/) { return T(1);  }
template<class T> inline T rsZeroValue( T /*value*/) { return T(0);  }

template<class TVal, class TTgt> 
inline TTgt rsConstantValue(TVal value, TTgt targetTemplate) { return (TTgt) value; }

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
// rename to rsAllFiniteNumbers

template<class T> inline T rsSqrt( T x) { return std::sqrt( x); }
template<class T> inline T rsExp(  T x) { return std::exp(  x); }
template<class T> inline T rsLog(  T x) { return std::log(  x); }
template<class T> inline T rsSin(  T x) { return std::sin(  x); }
template<class T> inline T rsCos(  T x) { return std::cos(  x); }
template<class T> inline T rsTan(  T x) { return std::tan(  x); }
template<class T> inline T rsFloor(T x) { return std::floor(x); }
template<class T> inline T rsCeil( T x) { return std::ceil( x); }
template<class T> inline T rsRound(T x) { return std::round(x); }
//template<class T> inline T rsSinh( T x) { return std::sinh( x); }
//template<class T> inline T rsCosh( T x) { return std::cosh( x); }
//template<class T> inline T rsTanh( T x) { return std::tanh( x); }
// the hyperbolic functions are already defined elsewhere (and differently!)

template<class T> inline T rsAtan2(T y, T x) { return std::atan2(y, x); }

// todo: 
// -sort them alphabetically, maybe use a shorthand #define for the common prefix
//  template<class T> inline T
// -merge code with RealFunctions.h
// -maybe wrap into a class rsMathFunctions or rsRealFunctions, do the same with the integer and 
//  complex functions
// -add more, see:
//  https://en.cppreference.com/w/cpp/header/cmath
//  https://en.cppreference.com/w/cpp/numeric/special_functions (C++17)
// -wrap and use modf, where appropriate - maybe it optimizes the splitting into int/frac
// -maybe take the arguments by const reference


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
template<class T> 
T rsAbsSquared(const T& x)
{
return x*x;
}
*/


template<class T> bool rsGreater(const T& a, const T& b) { return a > b; }
template<class T> bool rsLess(const T& a, const T& b)    { return a < b; }

/** Returns true, iff "left" has greater absolute value than "right" */
template <class T>
bool rsGreaterAbs(const T& left, const T& right)
{
  return rsAbs(left) > rsAbs(right);
}

template <class T>
bool rsGreaterAbs(const std::complex<T>& left, const std::complex<T>& right)
{
  return rsAbsSquared(left) > rsAbsSquared(right);
}

template <class T>
bool rsGreaterAbs(const T& left, const std::complex<T>& right)
{
  return rsAbs(left) > rsAbs(right);
}

template <class T>
bool rsGreaterAbs(const std::complex<T>& left, const T& right)
{
  return rsAbs(left) > rsAbs(right);
}

template <class T>
bool rsLessAbs(const std::complex<T>& left, const std::complex<T>& right)
{
  return rsAbsSquared(left) < rsAbsSquared(right);
}

template <class T>
bool rsLessAbs(const std::complex<T>& left, const T& right)
{
  return rsAbsSquared(left) < right*right;
}

template <class T>
bool rsIsCloseTo(const std::complex<T>& a, const std::complex<T>& b, const T& tol)
{
  std::complex<T> d = a-b;  // difference between a and b
  T m2 = rsAbsSquared(d);   // magnitude squared of difference
  return m2 <= tol*tol;
}


template <class T>
bool rsLessOrEqual(const T& left, const T& right)
{
  return left <= right;
}

/** Returns the biggest of the two values x and y where "biggest" means: has largest absolute 
value. (...could also be called rsBigger, but "biggest" may generalized to more than two values 
later and bigger may suggest something else) */
template <class T>
T rsBiggest(const T& x, const T& y)
{
  if( rsGreaterAbs(x, y) )
    return x;
  else
    return y;
}

/*
template <class T>
bool rsGreater(const T& left, const T& right)
{
return left > right;
}
*/

/** The maximum of two objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2)
{
  if(in1 > in2)
    return in1;
  else
    return in2;
}

/** The maximum of three objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2, T in3)
{
  return rsMax(rsMax(in1, in2), in3);
}

/** The maximum of four objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2, T in3, T in4)
{
  return rsMax(rsMax(in1, in2), rsMax(in3, in4));
}

/** Like rsMax but based on the "<" operator and swapping arguments (rather than using ">"). Can be
useful when a type just defines "<" but not ">". */
template <class T>
inline T rsMaxViaLess(T in1, T in2)
{
  if(in2 < in1)
    return in1;
  else
    return in2;
}

/** The minimum of two objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2)
{
  if(in1 < in2)
    return in1;
  else
    return in2;
}

/** The minimum of three objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2, T in3)
{
  return rsMin(rsMin(in1, in2), in3);
}

/** The minimum of four objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2, T in3, T in4)
{
  return rsMin(rsMin(in1, in2), rsMin(in3, in4));
}

/** Checks, if x is even. */
template<class T>
inline bool rsIsEven(T x) { return x % 2 == 0; } // maybe use bit-mask

/** Checks, if x is odd. */
template<class T>
inline bool rsIsOdd(T x) { return x % 2 != 0; }


#endif