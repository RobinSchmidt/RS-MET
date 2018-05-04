#ifndef RAPT_BASICFUNCTIONS_H_INCLUDED
#define RAPT_BASICFUNCTIONS_H_INCLUDED

/** Returns the absolute value of the input argument. It is intended to replace the standard
"abs" and "fabs" c-functions where genericity is desired. */
//template <class T>
//T rsAbs(T x);

/** Clips the value x into the range min...max such that for the returned value y, we have:
min <= y <= max  */
template<class T>
T rsClip(T x, T min = T(-1), T max = T(+1));

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
inline double rsExpToLin(double in, double inMin, double inMax, double outMin, double outMax);

/** The inverse of "rsLinToExpWithOffset" */
inline double rsExpToLinWithOffset(double in, double inMin, double inMax, double outMin,
  double outMax, double offset = 0.0);

/** Checks, if x is close to some target-value within some tolerance. */
inline bool rsIsCloseTo(double x, double targetValue, double tolerance);
  // actually, we should not need this anymore due to the templated version

/** Checks, if x is close to some target-value within some tolerance. */
template<class T>
inline bool rsIsCloseTo(T x, T targetValue, double tolerance);

/** Checks, if x is even. */
template<class T>
inline bool rsIsEven(T x) { return x % 2 == 0; } // maybe use bit-mask

/** Checks, if x is odd. */
template<class T>
inline bool rsIsOdd(T x) { return x % 2 != 0; }

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
inline double rsLinToLin(double in, double inMin, double inMax, double outMin, double outMax)
{ return outMin + (outMax-outMin) * (in-inMin) / (inMax-inMin); }

/** Converts a value between inMin and inMax into a value between outMin and outMax where the
mapping of the output is exponential. Example: y = linToExp(x, 0.0, 1.0, 20.0, 20000.0) will map
the input x assumed to lie inside 0.0...1.0 to the range between 20.0...20000.0 where equal
differences in the input lead to equal factors in the output. Make sure that the outMin value is
greater than zero! */
inline double rsLinToExp(double in, double inMin, double inMax, double outMin, double outMax);

/** Same as linToExp but adds an offset afterwards and compensates for that offset by scaling the
offsetted value so as to hit the outMax correctly. */
inline double rsLinToExpWithOffset(double in, double inMin, double inMax, double outMin,
  double outMax, double offset = 0.0);

/** The maximum of two objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2);

/** The maximum of three objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2, T in3);

/** The maximum of four objects on which the ">"-operator is defined. */
template <class T>
inline T rsMax(T in1, T in2, T in3, T in4);

/** The minimum of two objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2);

/** The minimum of three objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2, T in3);

/** The minimum of four objects on which the "<"-operator is defined. */
template <class T>
inline T rsMin(T in1, T in2, T in3, T in4);

/** Returns a power of two which is greater than or equal to the input argument. */
template <class T>
inline T rsNextPowerOfTwo(T x);

/** Returns the base raised to the power given by the exponent. It uses an algorithm based on
repeated squaring which has a complexity of O(log(exponent)).
for details, see: Jörg Arndt - Matters Computational, Ch.28.5
\todo bring back the naive implementation, to be used for small powers (there, it may be more
efficient), in rsPow, switch between rsPowSmall, rsPowBig depending on the size of the exponent
*/
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
int rsRoundToInt(T x);

/** Returns +1 for x > 0, -1 for x < 0 and 0 for x == 0. */
template <class T>
T rsSign(T x);

/** Swaps in1 and in2, if in1 > in2 */
template <class T>
inline void rsSortAscending(T &in1, T &in2);

/** Squares the input value. */
template <class T>
inline T rsSquare(T x);

/** Swaps two objects of class T. */
template <class T>
inline void rsSwap(T &in1, T &in2);

/** Returns a unity value of the given type. The idea is to use this template function to create
a unity-value inside other template functions where it might be required that the unity-value is
somehow parametrized. As an example, rsPow uses a unity-value as initializer for a multiplicative
accumulator. When it is invoked with a matrix-type, an explicit instantiation of rsUnity for the
matrix-type will be used to create an identity matrix with the required size (which is the same
as the size of "value"). If no explicit instantiation exists for the given type, it will fall
back to the default implementation, which returns T(1). */
template<class T>
inline T rsUnityValue(T value);

/** Wraps the number to the interval 0...length. */
inline double rsWrapAround(double numberToWrap, double length);

/** Periodically wraps a number into an interval between min and max, for example, an arbitrary
angle x may be wrapped into the interval 0...2*PI via xWrapped = wrapToInterval(x, 0.0, 2*PI) or
into the interval -PI...PI via xWrapped = wrapToInterval(x, -PI, PI). */
inline double rsWrapToInterval(double x, double min, double max);

/** Just outputs the constant value 0.0 for all inputs - used as default function pointer when
client code selects an invalid function-index, for example in the waveform-renderers. */
inline double rsZeroFunction(double x);

/** Returns a zero value of the given type. @see rsUnityValue */
template<class T>
inline T rsZeroValue(T value);

// \todo - is it somehow possible to get rid of the inlining?

#endif
