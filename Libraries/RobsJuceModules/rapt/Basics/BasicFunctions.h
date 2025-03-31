#ifndef RAPT_BASICFUNCTIONS_H_INCLUDED
#define RAPT_BASICFUNCTIONS_H_INCLUDED

// In Math/Functions/BasicMathFunctions.h/cpp, there is some stuff that overlaps with the stuff 
// here. Consolidate that code into one file!

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

/** Returns a unity value of the given type. The idea is to use this template function to create
a unity-value inside other template functions where it might be required that the unity-value is
somehow parametrized. As an example, rsPow uses a unity-value as initializer for a multiplicative
accumulator. When it is invoked with a matrix-type, an explicit instantiation of rsUnity for the
matrix-type will be used to create an identity matrix with the required size (which is the same
as the size of "value"). If no explicit instantiation exists for the given type, it will fall
back to the default implementation, which returns T(1). It's also useful for modular integers to
create 1 with the same modulus as some other number. */
template<class T> inline T rsUnityValue(T /*value*/) { return T(1); }

/** Returns a zero value of the given type. @see rsUnityValue */
template<class T> inline T rsZeroValue( T /*value*/) { return T(0); }

/** Turns a given integer constant into another target type T using a value from that target type
as prototype. It is used, for example, to convert an integer into a modular integer. in this case, 
the prototype value is used to copy the modulus from the prototype into the result. */
template<class T> inline T rsIntValue(int value, T targetTemplate) { return T(value); }


template<class T> inline T rsIdentity(  T value) { return value; }
// identity function


/** Returns true, iff x is equal to zero. */
template<class T> inline bool rsIsZero(T x)
{
  return x == rsZeroValue(x);
}

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


// We wrap the math functions from the standard library such that in the code, we can call these 
// instead. The reason is that we want the code to be generic - we may provide different explicit
// specializations for selfmade numeric data types (like simd vectors, multiprecision numbers, 
// etc.):
template<class T> inline T rsSqrt( T x) { return std::sqrt( x); }
template<class T> inline T rsExp(  T x) { return std::exp(  x); }
template<class T> inline T rsExp2( T x) { return std::exp2( x); }
template<class T> inline T rsLog(  T x) { return std::log(  x); }
template<class T> inline T rsLog2( T x) { return std::log2( x); }
template<class T> inline T rsSin(  T x) { return std::sin(  x); }
template<class T> inline T rsCos(  T x) { return std::cos(  x); }
template<class T> inline T rsTan(  T x) { return std::tan(  x); }
template<class T> inline T rsFloor(T x) { return std::floor(x); }
template<class T> inline T rsCeil( T x) { return std::ceil( x); }
template<class T> inline T rsRound(T x) { return std::round(x); }
//template<class T> inline T rsTrunc(T x) { return std::trunc(x); }
//template<class T> inline T rsSinh( T x) { return std::sinh( x); }
//template<class T> inline T rsCosh( T x) { return std::cosh( x); }
//template<class T> inline T rsTanh( T x) { return std::tanh( x); }
// The hyperbolic functions are already defined elsewhere (and differently!). That's kinda bad! May
// these selfmade definitions should use further qualifiers like "Fast" (when using exp to compute
// tanh, for example) or "Approx" (when approximations are used) or _mPi_Pi when only inputs 
// in -pi..+pi are allowed (such that we can do away any range reduction code) etc..
// ToDo: Try to declare them as constexpr.

template<class T> inline T rsAtan2(T y, T x) { return std::atan2(y, x); }
template<class T> inline T rsPow(  T a, T b) { return std::pow(  a, b); }

//template<class T> inline T rsPow(  T x, T y) { return std::pow(  x, y); }
// Defining this gives compilation error "ambiguous call ..." because of rsPow(T, int). This is 
// bad. Try to fix this! Maybe the version with integer exponent should be renamed to rsPowInt.
// ...OK - done (?)


// todo: 
// -sort them alphabetically, maybe use a shorthand #define for the common prefix
//  template<class T> inline T
// -merge code with RealFunctions.h
// -maybe wrap into a class rsMathFunctions or rsRealFunctions, do the same with the integer and 
//  complex functions
// -add more, see:
//  https://en.cppreference.com/w/cpp/header/cmath
//  https://en.cppreference.com/w/cpp/numeric/special_functions (C++17)
//  https://www.programiz.com/cpp-programming/library-function/cmath/trunc
// -wrap and use modf, where appropriate - maybe it optimizes the splitting into int/frac
// -maybe take the arguments by const reference


template<class T> inline int rsFloorInt(T x) { return (int) floor(x); }
template<class T> inline int rsCeilInt( T x) { return (int) ceil(x);  }

/** Returns the absolute value of the input argument. It is intended to replace the standard
"abs" and "fabs" c-functions where genericity is desired. */
template <class T>
T rsAbs(T x)
{
  // This (default, fallback, generic) implementation is suitable for real number types (float, 
  // double, int, rational, etc.) but not for complex types because in the complex case, the return
  // value is not of type T but rather of the underlying real number type.

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




template <class T>
T rsReal(const T& z)  // This is for when z is already a real number type such as float
{
  return z;
}

template <class T>
T rsReal(const std::complex<T>& z)
{
  return std::real(z);
}

template <class T>
T rsImag(const std::complex<T>& z)
{
  return std::imag(z);
}

template<class T>
inline void rsSetComplex(std::complex<T>* z, const T& newReal, const T& newImag)
{
  z->real(newReal);
  z->imag(newImag);
}

// ToDo:
//
// - Add function rsConj(). For real numbers, it should just be the identity. For complex numbers,
//   it should negate the imaginary part and keep the real part as is. Maybe have an in-place 
//   version of it, too - i.e. one that manipulates the input rather than returning a result. It 
//   may be more efficient to use that on complex types that have a complicated real type (such as
//   arbitrary precision floats, for example - they may use heap memory, so copies may be expensive
//   whereas conjugating in place may be a matter of inverting a flag)


// Maybe this should go into a file BitTwiddling.h where we collect various low level 
// bit-twiddling functions:
inline unsigned long rsBitReverse(unsigned long number, unsigned long numBits)
{
  unsigned long result = 0;
  for(unsigned long n=0; n<numBits; n++)
  {
    // Leftshift the previous result by one and accept the new LSB of the current number on the
    // right:
    result   = (result << 1) + (number & 1);

    // Rightshift the number to make the second bit from the right to the new LSB:
    number >>= 1;
  }
  return result;
}
// Needs documentation. This function is used in rsOrderBitReversed() which in turn is used in FFT 
// routines.


//--------------------------------------------------------------------------------------------------
// Definitions of a function rsMaxNorm for different types. The main intention of this function is 
// to find the maximum absolute value in an arbitrarily nested template container type. For 
// example, one could have a vector of complex values or a matrix of matrices of complex vectors or
// whatever. The behavior of rsMaxNorm should always be: return the maximum absolute element of the
// bottommost type which is one of the primitive types, i.e. built in C++ types like int or float.
// This rsMaxNorm function shall then be used as basis to implement an rsIsNegligible() function 
// that can be used to determine if a value is so close to zero that it can be considered zero. We
// need this mostly to deal with inexact floating point comparisons for equality (to zero).
//
// My first attempt was to use two template parameters: one for the argument type and another for
// the return type and then somehow structure the implementations such that it all works. It was 
// quite tricky to predict the interactions between C++ template instantiations, type deductions, 
// and function overload resolution rules to make the overload set of the rsMaxNorm function 
// behave exactly the way I want it to. I couldn't make it work this way. It turned out that the
// trick is to use return type deduction, i.e. to use auto for the return value rather than 
// manually declaring it via another template parameter (we already need one template parameter 
// for the argument type). Return type deduction requires at least C++14.
//
// I'm not sure if we should call it rsMaxNorm. Maybe think about the name some more. For 
// mathematical usage of the term, see:
//
// https://en.wikipedia.org/wiki/Norm_(mathematics)
// https://en.wikipedia.org/wiki/Norm_(mathematics)#Maximum_norm_(special_case_of:_infinity_norm,_uniform_norm,_or_supremum_norm)
// https://math.stackexchange.com/questions/285398/what-is-the-norm-of-a-complex-number 


// Base cases for rsMaxNorm for built in primitive types:

inline unsigned int rsMaxNorm(unsigned int x) { return x;           }
inline int          rsMaxNorm(int          x) { return std::abs(x); }
inline float        rsMaxNorm(float        x) { return std::abs(x); }
inline double       rsMaxNorm(double       x) { return std::abs(x); }
// ToDo: add all primitive types like int64_t etc.


/** Implements the maximum norm for std::complex<T> where T can be either double or float. The 
maximum of a complex number x + i*y is defined as max(|x|,|y|). When viewing the complex plane as
a 2D real vector space, this is also called the infinity norm for this vector space. That's because
it's the limit of the p-norm: L_p(x,y) = (|x|^p + |y|^p)^(1/p) as p approaches infinity. */
template<class T>
T rsMaxNorm(const std::complex<T>& z)
{
  return std::max(std::abs(z.real()), std::abs(z.imag()));
}

/** Implements the maximum norm for a std::vector of some type T. The max-norm of a vector is 
defined recursively as the maximum of the max-norms of the vector's elements. */
template<class T>
auto rsMaxNorm(const std::vector<T>& v)
{
  auto max = rsMaxNorm(T(0));
  for(auto& e : v)
    max = rsMax(max, rsMaxNorm(e));
  return max;

  // Maybe use std::max instead of rsMax and maybe use std::accumulate instead of a loop.
}

template<class T>
auto rsMaxNorm(const std::list<T>& v)
{
  auto max = rsMaxNorm(T(0));
  for(auto& e : v)
    max = rsMax(max, rsMaxNorm(e));
  return max;
}

// It's annyoing that we need to duplicate the code for any container type for which we want to
// support the rsMaxNorm operation. But if we want to implement it generically for all sorts of
// containers like below, we get an error related to rsMatrix not defining value_type. Apparently,
// the compiler tries to invoke this template for rsMatrix. Maybe it's because the specific 
// implementation actually takes an rsMatrixView rather than an rsMatrix.
//
//template<class TCont>
//auto rsMaxNorm(const TCont& v)
//{
//  using T = TCont::value_type;
//  auto max = rsMaxNorm(T(0));
//  for(auto& e : v)
//    max = rsMax(max, rsMaxNorm(e));  // Maybe try to use std::max
//  return max;
//}
//
// To make that work, I think rsMatrix (or maybe already rsMatrixView) needs to implement the 
// following STL compatibility features: a "using value_type = T;" definition, definition of the
// iterator type and begin() and end() functions. Maybe more. Another possibility could be to
// define an explicit specialization for rsMatrix itself such that the compiler selects that 
// instead of the generic container template. It could invoke the definition for rsMatrixView by 
// an upcast (cast to baseclass reference). Try that! It would be the less invasive solution and 
// therefore perhaps preferable over modifying rsMatrix(View). At the moment, it's fine as is 
// because I currently don't really need a max-norm function for any STL containers except 
// std::vector. The implementation for std::list is just there for testing purposes. So, for the 
// time being, it's fine. But maybe it's something to change later.

template<class T>
auto rsMaxNorm(const T* p, int N)
{
  auto max = rsMaxNorm(T(0));
  for(int i = 0; i < N; i++)
    max = rsMax(max, rsMaxNorm(p[i]));
  return max;

  // ToDo:
  //
  // - Try it with a type T that requires a prototype for correct initialization. Maybe something
  //   like rsMultiVector or rsModularInteger - although, for the latter, the notion of a 
  //   maximum-norm may be mathematically questionable and maybe for the former as well. But we may
  //   need to implement it, if we wnat to do linear algebra with them. But maybe we can directly
  //   implement rsIsNegligible or maybe we don't need it for rsModularInteger if rsIsZero is 
  //   correctly implemented. We'll see.
}




//-------------------------------------------------------------------------------------------------
// Under construction

// The rsIsNegligible() stuff is under construction. The intention is to provide a general 
// infrastructure to identify negligible values such as floating point numbers below a numerical 
// roundoff error threshold. This shall then be used to prune nonzero data that should actually be
// considered zero - for example, trailing polynomial coefficients after an addition of two 
// polynomials. Such things are especially needed when we deal with matrices of elements of a more
// complicated type (like polynomials or rational functions) and want to do linear algebra on them.
// If we don't prune/canonicalize the intermediate results after every arithmetic operation, the 
// expressions may blow up in e.g. Gaussian elimination. Matrices of rational functions occurr, for
// example, as transfer function matrices in state space filters and feedback delay networks. To 
// compute these, we need to run Gauss-Jordan matrix inversion on such matrices of transfer 
// functions. For such purposes, the negligibility/pruning/canonicalization infrastructure is 
// needed.

template<class TVal, class TTol> 
inline bool rsIsNegligible(TVal val, TTol tol)
{
  return rsLessOrEqual(rsAbs(val), tol);

  // This is the default implementation of the negligibility test. A value is considered negligible
  // if its absolute value is less than or equal to a given tolerance threshold. ...TBC...

  // ToDo: 
  //
  // - Figure out if we may need an rsAbs() function that takes two template parameters - one for 
  //   the argument type and one for the return type. The type of an absolute value may in general 
  //   differ from the argument type. For example, for complex arguments, the return type may be 
  //   the underlying real type. Maybe we need to call it like rsAbs<TVal, TTol>(val) here because 
  //   the return type (here TTol) cannot be inferred from the passed argument. So, maybe we need 
  //   to revisit the implementation(s) of rsAbs to make that work. Maybe we should call the 
  //   template parameters TVal, TAbs or TArg, TRes (for result). And maybe we should pass the 
  //   argument by const reference because at some point, the argument may be something like a 
  //   matrix (in which case we would return the max-abs value of all elements). We'll see...
  //
  // - Replace call to rsAbs by rsMaxNorm
}

template<class TVal> 
inline bool rsIsNegligible(TVal val, rsEmptyType tol)
{
  return rsIsZero(val);

  // This is an explicit partial specialization of the rsIsNegligible() function for when the 
  // tolerance type is an empty dummy type. In such cases, a value is considered to be negligible
  // only when it is exactly equal to zero.
}




#endif