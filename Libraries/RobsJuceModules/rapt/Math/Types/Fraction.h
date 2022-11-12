#ifndef RAPT_FRACTION_H
#define RAPT_FRACTION_H

/** Class for representing fractions a.k.a. rational numbers, i.e. ratios of two integers. 
Numerator and denominator are kept as signed integers "num", "den". On construction and in 
arithmetic operations, fractions are always put into a canonical representation which is a 
reduced form where the minus sign (if any) is put into the numerator. */

template<class T>  // T should be a signed int type
class rsFraction
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime


  rsFraction(T numerator = T(0), T denominator = T(1)) : num(numerator), den(denominator)
  { 
    canonicalize(); 
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void set(T numerator, T denominator) { num = numerator; den = denominator; canonicalize(); }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getNumerator()   const { return num; }
  T getDenominator() const { return den; }

  double toDouble() const { return double(num) / double(den); }

  //operator double() const { return double(num) / double(den); }
  // Perhaps it's not a good idea to allow implicit conversions to double. Client code should be
  // explicit.


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic operators

  // unary plus and minus:
  rsFraction operator+() const { return rsFraction(+num, den); }
  rsFraction operator-() const { return rsFraction(-num, den); }
  // optimization: these may avoid calling canonicalize()

  // +,-,*,/ where both arguments are fractions:
  rsFraction operator+(const rsFraction& b) const { return rsFraction(num*b.den + b.num*den, den * b.den); }
  rsFraction operator-(const rsFraction& b) const { return rsFraction(num*b.den - b.num*den, den * b.den); }
  rsFraction operator*(const rsFraction& b) const { return rsFraction(num * b.num, den * b.den); }
  rsFraction operator/(const rsFraction& b) const { return rsFraction(num * b.den, den * b.num); }

  // ..same for integer right arguments:
  rsFraction operator+(const T& b) const { return rsFraction(num + b*den, den); }
  rsFraction operator-(const T& b) const { return rsFraction(num - b*den, den); }
  rsFraction operator*(const T& b) const { return rsFraction(num * b, den); }
  rsFraction operator/(const T& b) const { return rsFraction(num, den * b); }

  // boilerplate for the +=, -=, *=, /= operators:
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

  /** Reduces to lowest terms and ensures that denominator is nonnegative. */
  void canonicalize() { reduce(); if(den < 0) { num = -num; den = -den; }  }



  T num, den;  // numerator and denominator (they are always kept canonical)

};

// operators for integer left argument:
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


// ToDo:
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
// -Implement functions to convert to double or float and/or operators for implicit conversion
//  (but maybe that's not a good idea - conversions should probably always be explicit)
// -In the Prototypes section, there's some stuff for converting between fractions and their
//  continued fraction representation - maybe drag that in. 
// -intAndFracPart via div and mod. r = n/d = i+f -> n = d*(i+f)
// -Maybe try to instantiate it for T = rsPolynomial. If that works at all, compare results to
//  rsRationalFunction...maybe that can even be rendered obsolete? ...but i don't think so, if only 
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
//  

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