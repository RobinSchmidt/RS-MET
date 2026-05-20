#ifndef RAPT_FRACTION_H
#define RAPT_FRACTION_H

/** Class for representing fractions a.k.a. rational numbers, i.e. ratios of two integers. 
Numerator and denominator are kept as signed integers "num", "den". On construction and in 
arithmetic operations, fractions are always put into a canonical representation which is a 
reduced form where the minus sign (if any) is put into the numerator.  */

template<class T>  // T should be a signed int type
class rsFraction
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime


  /** Default constructor. Constructs the canonicial representation of zero which is given by the
  fraction 0/1. */
  rsFraction() 
  {
    // Implementation is empty because the initialization is done in the member declaration.
  }

  /** Constructor that initializes numerator and denominator to the given values and then 
  canonicalizes the result (i.e. reduces it to lowest terms and makes sure that the minus sign, if 
  any, is in the numerator). Can also be used to initialize a fraction from an integer in which 
  case the denominator defaults to 1. */
  rsFraction(T numerator, T denominator = T(1)) : num(numerator), den(denominator)
  { 
    canonicalize();
    // ToDo: 
    // -Maybe do a static assert to make sure that T is a signed integer type.
    // -Maybe provide a specialized constructor for converting integers to fractions. It should 
    //  avoid the potentially costly call to canonicalize() which, besides other thing, calls a gcd
    //  algorithm. The denominator is just 1 is such a case. No need to call a gcd algo and other 
    //  stuff.
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void set(T numerator, T denominator) { num = numerator; den = denominator; canonicalize(); }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getNumerator()   const { return num; }
  T getDenominator() const { return den; }

  double toDouble() const { return double(num) / double(den); }
  float  toFloat()  const { return (float) toDouble(); }

  explicit operator double() const { return toDouble(); }
  // It's not a good idea to allow implicit conversions to double but an explicit conversion 
  // operator can be useful. See class rsQuadraticField in the research repo for usage of that 
  // operator. 


  // These need tests:
  bool isZero()        const { return num == T(0); }
  bool isPositive()    const { return num >  T(0); }
  bool isNegative()    const { return num <  T(0); }
  bool isNonPositive() const { return num <= T(0); }
  bool isNonNegative() const { return num >= T(0); }
  bool isInteger()     const { return den == T(1); }

  // isInteger works because we always keep the representation canonical. Well, actually it seems 
  // they rely on a canonical representation (or at least on one or another aspect of it). For 
  // example, isPositive() assumes that the denominator is positive. A non-canonical fraction like
  // -3/-5 should mathematically also count as positive. Maybe we should insert assertions like
  // rsAssert(isCanonical()); everywhere.
  // ToDo:  isOne


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic operators

  // Unary plus and minus:
  rsFraction operator+() const { return rsFraction(+num, den); }
  rsFraction operator-() const { return rsFraction(-num, den); }
  // ToDo: Optimization: these may avoid the call to canonicalize() in the constructor. Maybe have a 
  // private factory function that bypasses the call to canonicalize.

  // +,-,*,/ where both arguments are fractions:
  rsFraction operator+(const rsFraction& b) const { return rsFraction(num*b.den + b.num*den, den * b.den); }
  rsFraction operator-(const rsFraction& b) const { return rsFraction(num*b.den - b.num*den, den * b.den); }
  rsFraction operator*(const rsFraction& b) const { return rsFraction(num * b.num, den * b.den); }
  rsFraction operator/(const rsFraction& b) const { return rsFraction(num * b.den, den * b.num); }

  // The same for integer right arguments:
  rsFraction operator+(const T& b) const { return rsFraction(num + b*den, den); }
  rsFraction operator-(const T& b) const { return rsFraction(num - b*den, den); }
  rsFraction operator*(const T& b) const { return rsFraction(num * b, den); }
  rsFraction operator/(const T& b) const { return rsFraction(num, den * b); }

  // Boilerplate for the +=, -=, *=, /= operators:
  rsFraction& operator+=(const rsFraction& b) { return *this = (*this) + b; }
  rsFraction& operator-=(const rsFraction& b) { return *this = (*this) - b; }
  rsFraction& operator*=(const rsFraction& b) { return *this = (*this) * b; }
  rsFraction& operator/=(const rsFraction& b) { return *this = (*this) / b; }
  rsFraction& operator+=(const T& b) { return *this = (*this) + b; }
  rsFraction& operator-=(const T& b) { return *this = (*this) - b; }
  rsFraction& operator*=(const T& b) { return *this = (*this) * b; }
  rsFraction& operator/=(const T& b) { return *this = (*this) / b; }




  //-----------------------------------------------------------------------------------------------
  // \name Comparison operators

  bool operator==(const rsFraction& b) const { return num == b.num && den == b.den; }
  bool operator!=(const rsFraction& b) const { return !(*this == b); }
  bool operator< (const rsFraction& b) const { return num * b.den <  b.num * den; }
  bool operator<=(const rsFraction& b) const { return num * b.den <= b.num * den; }
  bool operator> (const rsFraction& b) const { return num * b.den >  b.num * den; }
  bool operator>=(const rsFraction& b) const { return num * b.den >= b.num * den; }


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Reduces this number to lowest terms. */
  void reduce() { T gcd = rsGcd(num, den); num /= gcd; den /= gcd; }
  // ToDo: Verify and document that rsGcd does the right thing when the input is negative. 

  /** Brings this fraction into its canonical form by reducing it to lowest terms and ensuring that
  denominator is nonnegative. */
  void canonicalize() { reduce(); if(den < 0) { num = -num; den = -den; }  }
  // Actually, the denominator is supposed to be positive and not only "nonnegative". However,
  // this function here really does only ensure nonnegativity, so the documentation is actually
  // accurate. That the denominator is nonzero must be ensured elsewhere. But maybe we should
  // allow a denominator of zero to represent infinity and NaN (NaN would be when num _and_ den are
  // zero). Maybe it could be better to first do the potential sign-flip and then reduce()? The 
  // rationale being that the gcd algorithm may be tripped up by negative inputs. I'm not sure 
  // about that, though - it may work just fine. Also, what about when both, num and den are 
  // negative? Maybe that never happens? At least not in arithmetic operations? But it may happen
  // when the user passes such arguments to the constructor. Verify that we have unit tests for all
  // of these cases. 
  // Maybe rename to canonize(). That's much shorter and seems to be a legal English word, see:
  // https://www.merriam-webster.com/dictionary/canonize If doing so, do it also in opther classes 
  // with similar functionality, for example in rsSparseRationalFunction. Search the whole codebase
  // for such canonicalization (canonization) functions and name them all consistently.

  /** Numerator and denominator. They are always kept canonical, i.e. in reduced form and with 
  minus sign in numerator if the number is negative. */
  T num = 0, den = 1;

};

// Operators for integer left argument:
template<class T>
rsFraction<T> operator+(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator() + r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator-(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator() - r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator*(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator/(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator(), r.getNumerator()); }


template<class T>
inline bool rsIsBetterPivot(const rsFraction<T>& x, const rsFraction<T>& y)
{
  //return rsGreaterAbs(x, y); // Default implementation suitable for floating point types

  // A zero x is never a better pivot than any y:
  if(x.isZero())
    return false;

  // Any nonzero x is always better than a zero y:
  if(y.isZero())
    return true;

  // An x with smaller denominator is better than a y with larger denominator:
  if(x.getDenominator() < y.getDenominator())
    return true;

  // An x with larger denominator is worse that a y with smaller denominator:
  if(x.getDenominator() > y.getDenominator())
    return false;

  // When x and y have the same denominator, we compare the numerators. A smaller absolute value is
  // better (except when it's zero - but this has already been ruled out):
  return rsAbs(x.getNumerator()) < rsAbs(y.getNumerator());

  // The rationale behind this is that a good pivot should be "simple" in the sense that it is
  // unlikely to blow up the size of the denominators in subsequent arithmetic operations. It's a 
  // heuristic, though. If we will see a blow up of complexity depends on many more variables. But
  // it might be the most reasonable thing we can do when we take into account only two values. 
  // Doing any better would require us to analyze the whole matrix with all the interactions 
  // between all the elements. That may very well make pivoting so expensive that it totally 
  // dominates the cost of a matrix inversion, I guess (Verify! ...or at least justify)
}
// Needs tests. I'm not yet quite sure about the appropriateness of the applied criteria.


// Functions to compute floor, ceil, round, etc. So far, the implementations are rather naive and 
// can perhaps be optimized. They should also be documented:

template<class T>
rsFraction<T> rsTrunc(const rsFraction<T>& x)
{
  return x.getNumerator() / x.getDenominator();
  // Trunction is just integer division, i.e. floor-division.


  // I think, the code above implictly invokes the constructor call:
  // 
  //   return rsFraction(x.getNumerator() / x.getDenominator(), 1);
  //
  // Maybe it would be better to write it out explicitly for enhanced documentation value. 
  // Figure this out and perhaps change it!
}

template<class T>
rsFraction<T> rsFloor(const rsFraction<T>& x)
{
  if(x.isInteger())        // Integers stay as is
    return x;
  if(x.isNonNegative())    // Non-negative fractions are truncated
    return rsTrunc(x);
  return rsTrunc(x) - 1;   // Negative fractions need a -1 after truncation

  // ToDo: Optimize - maybe use:
  // return rsFraction<T>(x.getNumerator() / x.getDenominator() - 1, 1);
  // for the last case. This avoids the subtraction of fractions (which is expensive) and uses only
  // subtraction of integers.
}

template<class T>
rsFraction<T> rsCeil(const rsFraction<T>& x)
{
  if(x.isInteger())
    return x;
  if(x.isNonPositive())
    return rsTrunc(x);
  return rsTrunc(x) + 1;
}


// Some free functions that are relevant mainly in the context of matrices of fractions:

template<class T>
inline bool rsIsInvalidDivisor(const rsFraction<T>& p, const rsFraction<T>& tol)
{
  return p.isZero();
}

template<class T>
inline rsFraction<T> rsGetPivotingTolerance(const rsMatrixView<rsFraction<T>>& /*A*/)
{
  return rsFraction<T>(0, 1);
}

template<class T>
auto rsMaxNorm(const rsFraction<T>& q)
{
  return rsAbs(q);

  // This currently invokes the fallback implementation of rsAbs that checks if less then zero and
  // if so, returns -x otherwise returns x. Maybe we can provide a more efficient implementation 
  // of rsAbs specifically for rsFraction. We just need to take the abs of the numerator and leave 
  // the denominator as is. Maybe we can even bypass the canonicalization in the constructor. 
  // Perhaps by using a private factory function makeUnchecked(T newNum, T newDen) or something 
  // like that.
}


// ToDo:
// -Implement functions for truncation, floor, ceiling, rounding. Maybe as free functions rsTrunc, 
//  rsFloor, rsCeil, rsRound. I think, truncation can simply be done by returning num / den, i.e. 
//  doing integer division. See also Basics/BasicFunctions.h. There are the fallback versions for 
//  float, double, etc. We should implement explicit specializations for rsFraction here. There's
//  some prototype code in the unit test already which may be dragged over here.
// -Maybe implement rsPow - at least for integer exponents
// -Maybe use algorithms for the arithmetic operators that make overflow less likely (divide by gcd 
//  before computing products, use lcm in + and - instead of just computing products, etc.).
//  I think addition of a/b + c/d could be like (2 variants):
//    -add(a,b,c,d): g = gcd(b,d); A = a*(d/g); C = c*(b/g); D = b*(d/g); R = (A+C) / D;
//    -add(a,b,c,d): D = lcm(b,d); A = D/d; B = D/b; R = (A+C) / D;
//  Multiplications (a/b) * (c/d) could be like:
//    -mul(a,b,c,d): g = gcd(a,d); h = gcd(c,b); R = ((a/g)*(c/h)) / ((b/h)*(d/g);
//  Subtraction would be similar to addition (just with a - rather then a +) and division similar
//  to multiplication (roles of c and d need to be swapped, I think). Implement and verify those 
//  algorithms and test them for overflow behavior and performance compared to the naive algos! 
//  Keep the naive algos as prototypes for unit tests and reference. Test and debug with: 
//  (1/4)+(1/6), (2/3)*(3/2),(4/3)*(6/2),(9/5)*(5/6), ... - look at intermediate results. They 
//  should be smaller with the algorithms above. Try it also with larger numbers that actually do
//  produce overflow in the naive algos but work fine with the algos above. Document the conditions
//  under which overflow will occur. Maybe insert assertions that catch overflows. All of that
//  will uglify the code. Maybe move this simple and clean implementation into the prototypes 
//  section as rsFractionNaive and write a better (safer, faster, with higher range) but uglier 
//  production version.
// -Implement some functions like pow (with integer exponent)..maybe using the ^ operator - but 
//  care has to be taken to parenthesize expressions like (r^i) inside longer expressions due to 
//  C++ precendence rules
// -Maybe detect if overflow will happen and trigger an assert
// -In the Prototypes section, there's some stuff for converting between fractions and their
//  continued fraction representation - maybe drag that in. 
// -intAndFracPart via div and mod. r = n/d = i+f -> n = d*(i+f)
// -Maybe try to instantiate it for T = rsPolynomial. If that works at all, compare results to
//  rsRationalFunction...maybe that can even be rendered obsolete? ...but I don't think so, if only 
//  for efficiency reasons.
// -Try to implement reduce and canonicalize in a branchless way to admit T to be a SIMD type. But 
//  this requires a branchless implementation of rsGcd, or at least an implementation that runs the 
//  while loop until it has finished for all components.
// -Maybe represent +inf, -inf, -0, nan as 1/0, -1/0, 0/0, 0/(-1) and implement rsIsInf, 
//  rsIsNan for rsFraction ...but this requires even more branching in all arithmetic operators. 
//  There's a tradeoff to be made between efficiency and feature-set - but this class is not meant 
//  for realtime-dsp (that we do with floats). But actually, we only need a check if(den == 0) in 
//  all operation to detect an exception. This should be a highly predictable branch, so it may be
//  worth it. Maybe +-inf should be the result of operations that overflow in the numerator but not
//  in the numerator. +-0 should result when the denominator overflows but not the numerator, nan
//  results when both overflow or we subtract inf - inf, etc. Replicate IEEE floating point 
//  behavior. The code would become much more messy, though. Maybe it could be worthwhile to have 
//  different implementations for different purposes: fast (for production), safe (for mathematical
//  experiments in R&D), naive (for reference in unit tests), etc. Or maybe overflow shouldn't be 
//  interpreted as inf? Maybe +-inf should only be the result of dividing +-n by +-0 for any finite
//  n. We'll see...
// -Maybe integrate some stuff that deals with continued fractions - see the unit test
//  See:   https://www.youtube.com/watch?v=tBc_xcRzMxk  Continued Fraction Arithmetic
// -Document what happens when the user tries to create a fraction with zero as denominator.


// Notes:
// -Maybe it's sometimes convenient to keep fractions in unreduced form. It may be easier to spot 
//  patterns in sequences of unreduced rational numbers that come from some computation
//  -but this will be relevant only for research code, not production code
//  -maybe we could introduce a compile-time switch (maybe a boolean template parameter) that 
//   controls if we canonicalize or not
//  -enforced canonical representation is important for the == operator to work properly...maybe it 
//   should be implemented in a way that admits non-canonical representations? 
//   a/b == c/d  <->  a*d == b*c
//  -maybe we should have a sub- or baseclass rsUnreducedFraction or maybe rsFraction should not
//   enforce a canonical representation but another class rsRationalNumber should?
// -For a potentially better algorithm for computing the sum of two fractions, see the comment in 
//  implementation of the += operator here:
//  https://www.boost.org/doc/libs/1_77_0/boost/rational.hpp


#endif