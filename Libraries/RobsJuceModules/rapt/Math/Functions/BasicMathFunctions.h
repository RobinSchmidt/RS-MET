#ifndef RAPT_BASICMATHFUNCTIONS_H_INCLUDED
#define RAPT_BASICMATHFUNCTIONS_H_INCLUDED

// todo: but maybe the content could also be absorbed into RealFunctions.h/cpp and
// Basics/BasicFunctions.h (the min/max stuff could go into Basics)


/** Returns the absolute value of the input argument. It is intended to replace the standard
"abs" and "fabs" c-functions where genericity is desired. */
//template <class T>
//T rsAbs(T x);

/** Clips the value x into the range min...max such that for the returned value y, we have:
min <= y <= max  */
template<class T>
inline T rsClip(T x, T min = T(-1), T max = T(+1));

///** Adds the 2 (signed) 32 bit integer values and clips the result to the desired range.
//Overflow is avoided by using 64-bit addition internally. */
//inline rsInt32 rsClippedSum(rsInt32 a, rsInt32 b,
//  rsInt32 min = RS_MIN(rsInt32), rsInt32 max = RS_MAX(rsInt32));
//
///** Subtracts the 2 (signed) 32 bit integer values and clips the result to the desired range.
//Overflow is avoided by using 64-bit addition internally. */
//inline rsInt32 rsClippedDifference(rsInt32 a, rsInt32 b,
//  rsInt32 min = RS_MIN(rsInt32), rsInt32 max = RS_MAX(rsInt32));

/** Exclusive "or": Returns true iff either of the operands but not both are true. */
inline bool rsExOr(bool a, bool b);

/** The inverse of "rsLinToExp" */
template<class T>
inline T rsExpToLin(T in, T inMin, T inMax, T outMin, T outMax);

/** The inverse of "rsLinToExpWithOffset" */
template<class T>
inline T rsExpToLinWithOffset(T in, T inMin, T inMax, T outMin, T outMax, T offset = 0.0);

/** Checks, if x is close to some target-value within some tolerance. */
//inline bool rsIsCloseTo(double x, double targetValue, double tolerance);
  // actually, we should not need this anymore due to the templated version

/** Checks, if x is close to some target-value within some tolerance. */
template<class T>
//inline bool rsIsCloseTo(T x, T targetValue, double tolerance);
inline bool rsIsCloseTo(T x, T targetValue, T tolerance);

/** Checks, if x is a power of 2. */
inline bool rsIsPowerOfTwo(unsigned int x);

/** Returns true when min <= x <= max, false otherwise. */
template<class T>
inline bool rsIsInRange(T x, T min, T max);

/** Splits a double into integer and fractional part, the integer part is returned in the
return value, the fractional part in the second argument. */
inline int rsIntAndFracPart(double x, double &frac);

/** Computes an integer power of x by successively multiplying x with itself. */
inline double rsIntegerPower(double x, int exponent);

/** Limits the value of an object.   ...redundant with clip */
//template <class T>
//inline T rsLimit(T x, T lowerBound, T upperBound);

/** Converts a value between inMin and inMax into a value between outMin and outMax where the
mapping is linear for the input and the output. Example: y = rsLinToLin(x, 0.0, 1.0, -96.0, 24.0)
will map the input x assumed to lie inside 0.0...1.0 to the range between -96.0...24.0. This
function is useful to convert between normalized parameter representations between 0.0...1.0 and 
the clear-text parameters. */
template<class T>
inline T rsLinToLin(T in, T inMin, T inMax, T outMin, T outMax)
{ return outMin + (outMax-outMin) * (in-inMin) / (inMax-inMin); }


/** Converts a value between inMin and inMax into a value between outMin and outMax where the
mapping of the output is exponential. Example: y = linToExp(x, 0.0, 1.0, 20.0, 20000.0) will map
the input x assumed to lie inside 0.0...1.0 to the range between 20.0...20000.0 where equal
differences in the input lead to equal factors in the output. Make sure that the outMin value is
greater than zero! */
template<class T>
inline T rsLinToExp(T in, T inMin, T inMax, T outMin, T outMax);

/** Same as linToExp but adds an offset afterwards and compensates for that offset by scaling the
offsetted value so as to hit the outMax correctly. */
template<class T>
inline T rsLinToExpWithOffset(T in, T inMin, T inMax, T outMin, T outMax, T offset = 0.0);



/** Computes the median of 3 values. */
template<class T>
inline T rsMedian(T x1, T x2, T x3);

/** Returns x if x is even, else x+1. */
template <class T>
inline T rsNextEvenNumber(T x);

// rsNextOddNumber

/** Returns a power of two which is greater than or equal to the input argument x where x is 
assumed to be >= 1. */
template <class T>
inline T rsNextPowerOfTwo(T x);

/** Returns the base raised to the power given by the exponent. It uses an algorithm based on
repeated squaring which has a complexity of O(log(exponent)).
for details, see: Jörg Arndt - Matters Computational, Ch.28.5
\todo bring back the naive implementation, to be used for small powers (there, it may be more
efficient - well..will it? questionable!), in rsPow, switch between rsPowSmall, rsPowBig 
depending on the size of the exponent */
template <class T>
T rsPow(const T& base, int exponent);

/** Generates a random number that is uniformly distributed between min and max (inclusive). The
underlying integer pseudo random number generator is a linear congruential with period length of
2^32. It is based on Numerical Recipies in C (2nd edition), page 284. You may pass a seed to the
first call to initialize it - otherwise it will use 0 as seed. A negative number (as in the
default argument) will indicate to not initialize the state and just generate a random number
based on the last state (which is the case for a typical call). */
inline double rsRandomUniform(double min = 0.0, double max = 1.0, int seed = -1);

/** Computes the scale factor and additive offset that needs to be applied to a variable x, such
that an input range of x-values from inMin to inMax is mapped linearly to an output range of
transformed x'-values from outMin to outMax. You compute the transformed variable by
x' = scale*x + shift, using the scale/shift coefficients computed by this function. */
//inline void rsRangeConversionCoefficients(double inMin, double inMax,
//  double outMin, double outMax, double *scale, double *shift);
template<class T>
inline void rsRangeConversionCoefficients(T inMin, T inMax, T outMin, T outMax, 
  T *scale, T *shift);

/** Rounds the given x to the closest integer. */
template <class T>
inline int rsRoundToInt(T x) { return (int) ::round(x); }

/** Returns +1 for x > 0, -1 for x < 0 and 0 for x == 0. */
template <class T>
inline T rsSign(T x) { return T(T(0) < x) - (x < T(0)); }

/** Calculates sine and cosine of x - this may (depending on datatype, compiler and instruction 
set) be more efficient than calling sin(x) and cos(x) seperately. */
template<class T>
inline void rsSinCos(T x, T* sinResult, T* cosResult);
// move to RealFunctions

/** Swaps in1 and in2, if in1 > in2 */
template <class T>
inline void rsSortAscending(T &in1, T &in2);

/** Squares the input value. */
template <class T>
inline T rsSquare(T x);

/** Swaps two objects of class T. */
template <class T>
inline void rsSwap(T &in1, T &in2);
// merge with Basics

/** Returns a unity value of the given type. The idea is to use this template function to create
a unity-value inside other template functions where it might be required that the unity-value is
somehow parametrized. As an example, rsPow uses a unity-value as initializer for a multiplicative
accumulator. When it is invoked with a matrix-type, an explicit instantiation of rsUnity for the
matrix-type will be used to create an identity matrix with the required size (which is the same
as the size of "value"). If no explicit instantiation exists for the given type, it will fall
back to the default implementation, which returns T(1). */
template<class T>
inline T rsUnityValue(T value);
// Merge with Basics

/** Wraps the number to the interval 0...length. */
inline double rsWrapAround(double numberToWrap, double length);

/** Periodically wraps a number into an interval between min and max, for example, an arbitrary
angle x may be wrapped into the interval 0...2*PI via xWrapped = wrapToInterval(x, 0.0, 2*PI) or
into the interval -PI...PI via xWrapped = wrapToInterval(x, -PI, PI). The left limit is included 
and the right limit is excluded. */
inline double rsWrapToInterval(double x, double min, double max);
// rename to rsWrap

/** Just outputs the constant value 0.0 for all inputs - used as default function pointer when
client code selects an invalid function-index, for example in the waveform-renderers. */
inline double rsZeroFunction(double x);

/** Returns a zero value of the given type. @see rsUnityValue */
template<class T>
inline T rsZeroValue(T value);
// Merge with basics

// \todo - is it somehow possible to get rid of the inlining?

//=================================================================================================
// implementations:

template<class T>
inline T rsClip(T x, T min, T max)
{
  rsAssert(min <= max);
  if(x < min)
    return min;
  if(x > max)
    return max;
  return x;
}

template<class T>
inline T rsExpToLin(T in, T inMin, T inMax, T outMin, T outMax)
{
  T tmp = log(in / inMin) / log(inMax / inMin);
  return outMin + tmp * (outMax - outMin);
}

template<class T>
inline T rsExpToLinWithOffset(T in, T inMin, T inMax, T outMin, T outMax, T offset)
{
  T tmp = in * (inMax + offset) / inMax;
  tmp -= offset;
  return rsExpToLin(tmp, inMin, inMax, outMin, outMax);

  //T tmp = linToExp(in, inMin, inMax, outMin, outMax);
  //tmp += offset;
  //tmp *= outMax/(outMax+offset);
  //return tmp;
}

template<class T>
inline bool rsIsCloseTo(T x, T targetValue, T tolerance)
{
  rsAssert(tolerance >= T(0), "tolerance must be non-negative");
  if(rsAbs(x - targetValue) <= tolerance)
    return true;
  else
    return false;
}

template<class T>
inline bool rsIsCloseTo(std::complex<T> x, std::complex<T> targetValue, std::complex<T> tolerance)
{
  rsAssert(tolerance.imag() == T(0), "tolerance is assumed to be a real number");
  rsAssert(tolerance.real() >= T(0), "tolerance must be non-negative");
  if( abs(x - targetValue) <= tolerance.real() )
    return true;
  else
    return false;
}

template<class T>
inline bool rsIsInRange(T x, T min, T max)
{
  if(x >= min && x <= max)
    return true;
  else
    return false;
}


template<class T>
inline T rsLinToExp(T in, T inMin, T inMax, T outMin, T outMax)
{
  // map input to the range 0.0...1.0:
  T tmp = (in - inMin) / (inMax - inMin);

  // map the tmp-value exponentially to the range outMin...outMax:
  return outMin * std::exp(tmp * (log(outMax / outMin)));
}

template<class T>
inline T rsLinToExpWithOffset(T in, T inMin, T inMax, T outMin, T outMax, T offset)
{
  T tmp = rsLinToExp(in, inMin, inMax, outMin, outMax);
  tmp += offset;
  tmp *= outMax / (outMax + offset);
  return tmp;
}

template<class T>
inline T rsMedian(T x1, T x2, T x3)
{
  if(x1 >= x2 && x1 >= x3) return rsMax(x2, x3);  // x1 is greatest
  if(x2 >= x1 && x2 >= x3) return rsMax(x1, x3);  // x2 is greatest
  if(x3 >= x1 && x3 >= x2) return rsMax(x1, x2);  // x3 is greatest
  rsError("we should always take one of the branches above");
  return 0;
}
// todo: test with all permutations of 1,2,3 and 1,2,2
// we could also do: if(x3>x2) swap(x2,x3); if(x2>x1) swap(x1,x2); return x2; ...that would
// amount to doing a 3-value bubble-sort and returning the middle - we would trade 9 
// comparisons for 2 comparisons and 2 swaps - maybe benchmark, which is better

template <class T>
inline T rsNextPowerOfTwo(T x)
{
  T accu = 1;
  while(accu < x)
    accu *= 2;
  return accu;
}

// general fallback version:
template<class T>
inline void rsSinCos(T x, T* sinResult, T* cosResult)
{
  *sinResult = sin(x);
  *cosResult = cos(x);
}

// explicit specialization for double:
inline void rsSinCos(double x, double* sinResult, double* cosResult)
{
#if MSC_X86_ASM
  double s, c;        // do we need these intermediate variables?
  __asm fld x
  __asm fsincos
  __asm fstp c
  __asm fstp s
  *sinResult = s;
  *cosResult = c;
#else
  *sinResult = sin(x);
  *cosResult = cos(x);
#endif
}
// here's some info, how to (not) optimize this in X64 builds:
// http://www.gamedev.net/topic/598105-how-to-implement-sincos-on-vc-64bit/

// taken from rosic - todo: edit coding style and use them as explicit specializations 
// for efficiency:
/*
INLINE int rosic::roundToInt(double const x)
{
  int n;
#if defined(__unix__) || defined(__GNUC__)
  // 32-bit Linux, Gnu/AT&T syntax:
  __asm ("fldl %1 \n fistpl %0 " : "=m"(n) : "m"(x) : "memory");
#else
  // 32-bit Windows, Intel/MASM syntax:
  __asm fld qword ptr x;
  __asm fistp dword ptr n;
#endif
  return n;
}

INLINE int rosic::roundToIntSSE2(double const x)
{
  return _mm_cvtsd_si32(_mm_load_sd(&x));
}
*/

#endif
