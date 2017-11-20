#ifndef RS_BIGFLOAT_H
#define RS_BIGFLOAT_H

namespace RSLib
{

  /**

  This is a class for representing floating point numbers with arbitrary numeric precision in
  an arbitrary base up to 2^32. If arithmetical operators ar invoked with operands of different
  precision, the result willl have a precision equal to the greater of the two. arithmetic
  operations on numbers in different bases are not supported.

  \todo implement elementary functions (sqrt, sin, cos, tan, exp, log, asin, pow, floor, ceil,
  round, fmod, etc...), ideally as template functions
  -> use standard library functions as initial guess and improve iteratively (during debug, use a
     less precise initial guess - maybe a truncated taylor series)
  -> use newton iteration for sqrt
  -> use cordic for sin, cos, exp, log, asin, acos, etc. - nah, these require tables
  -> use identities for tan, sinh, pow, etc.

  \todo implement functions to calculate mathematical constants (pi, e, ...), maybe let the first
  instance that is instantiated with a particular precision settings calculate all these constants
  for later use in tanscendental functions

  \todo implement some special functions (bessel, elliptic, gamma, Si, ...)
  see: cephes and NR

  \todo implement functions returning constants that parallel the corresponding functions from
  std::numeric_limits (i.e. epsilon(), infinity(), min(), quiet_NaN, signaling_NaN(), ...)

  \todo: maybe have a subclass that carries information about the worst-case error with it

  */

  class RSLib_API rsBigFloat : public rsBigNumber
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Creates an arbitrary precsion number where the precision parameter determines
    the the number of digits for the fraction-part (aka mantissa) in the given base. */
    rsBigFloat(rsUint32 precision = 1, rsUint64 baseToUse = maxBase);

    /** Constructor. Creates an arbitrary precsion number from a string of decimal digits. */
    rsBigFloat(const char value[] /*, rsUint64 baseToUse = 10*/);
      // \todo generalize to allow digits in other bases to be passed in the string (2nd parameter
      // stringBase, defaults to 10

    /** Constructor. Creates an arbitrary precsion number from a C++ double precision number. */
    rsBigFloat(double value, rsUint32 precisionToUse = 2, rsUint64 baseToUse = maxBase);

    /** Constructor. Creates an arbitrary precision number from an rsBigInt. */
    rsBigFloat(const rsBigInt& bigInt);

    /** Copy constructor. */
    rsBigFloat(const rsBigFloat& other);


    /** \name Operators */

    rsBigFloat& operator=(const rsBigFloat& other);

    rsBigFloat operator-();

    bool operator==(const rsBigFloat& other) const;
    bool operator!=(const rsBigFloat& other) const;
    bool operator>=(const rsBigFloat& other) const;
    bool operator<=(const rsBigFloat& other) const;
    bool operator>(const rsBigFloat& other) const;
    bool operator<(const rsBigFloat& other) const;

    rsBigFloat operator+(const rsBigFloat& other);
    rsBigFloat operator-(const rsBigFloat& other);
    rsBigFloat operator*(const rsBigFloat& other);
    rsBigFloat operator/(const rsBigFloat& other);

    // hmm... in rsComplex, we have the binary operators defined with 2 parameters, like this:
    // rsComplex<RealType> operator*(const rsComplex<RealType> &z, const rsComplex<RealType> &w)

    rsBigFloat& operator+=(const rsBigFloat& other);
    rsBigFloat& operator-=(const rsBigFloat& other);
    rsBigFloat& operator*=(const rsBigFloat& other);
    rsBigFloat& operator/=(const rsBigFloat& other);


    /** \name Setup */

    /** Sets this number to zero. */
    void setZero();

    // setPrecision(rsUint32 newPrecision)
    //  should append zeros, when newPrecision > oldPrecision and round if
    //   newPrecision > oldPrecision

    inline void setExponent(rsInt32 newExponent)
    {
      exponent = newExponent;
      // \todo detect and handle over- and underflow
    }

    /** Sets the data directly - mainly for debugging. */
    void setData(rsUint32 *digits, rsUint32 numDigits, rsUint64 base, rsInt32 exponent,
                 bool negative, bool shouldNormalize = false);

    /** Sets the value from an array of digits using the current base and number of digits. If
    numNewDigits is greater than the number of digits of this number, rounding will be used. If it
    is less, zero padding will be used. */
    void setValue(rsUint32 *newDigits, rsUint32 numNewDigits, rsInt32 exponent, bool negative);


    /** \name Inquiry */


    /** Returns the exponent part of this number */
    inline rsInt32 getExponent() const
    {
      return exponent;
    }

    /** Returns true, iff the 1st digit of the fraction is nonzero. */
    inline bool isNormalized() const
    {
      return digits[0] != 0;
    }

    /** Returns true, iff this number is either normalized or zero.. */
    inline bool isNormalizedOrZero() const
    {
      return isNormalized() || isZero();
    }

    /** Returns true, iff the absolute value of this number is small compared to the absolute value
    of x.
    \todo: state explicitly, what "small" means - perhaps, it would make sense to define y being
    small compared to x, iff adding y to x can change x only in the last digit (except in cases
    where carry occurs), or maybe even better, let the number of digits in y that may change be a
    parameter that defaults to 1. */
    inline bool isSmallComparedTo(const rsBigFloat& x)
    {
      rsAssert(isNormalizedOrZero() && x.isNormalizedOrZero() );
        // in case of unnormalized  inputs, we can't rely on the exponents - we must actually
        // compare digits

      return isZero() || (x.exponent - this->exponent >= (rsInt32)x.numDigits - 1);
        // check, if this criterion makes sense, we probably should define "small" as within 1
        // unit in the last place
    }

    /** Returns the size of the unit in the last place ("ULP") of this number, taking the exponent
    into account. Note that the returned value is only precise within double-precision. */
    double getUnitInTheLastPlace()
    {
      return pow((double)base, (double)exponent - (double)numDigits);
    }

    /** Returns -1 if a < b, 0 if a == b, +1 if a > b. */
    static int compare(const rsBigFloat& a, const rsBigFloat& b);

    /** Returns the reciporcal of one unit in the last place "ULP" for the given number of digits
    and base - the higher it (the reciprocal) is, the higher is the precision of the number
    representation. */
    static rsBigInt unitInTheLastPlaceReciprocal(rsUint32 numDigits, rsUint32 base)
    {
      return rsPow(rsBigInt(base), numDigits);
    }

    /** Returns the precision of a number representation with given number of digits and base in
    the unit of bits. In general (i.e. if the base is not a power of two), this will be a
    non-integer number. */
    static double precisionInBits(rsUint32 numDigits, rsUint64 base)
    {
      return (double) numDigits * rsLog2((double) base);
    }

    /** Returns the number of required digits, such that the precision is at least equivalent to a
    binary representation with "desiredPrecisionInBits" bits.
    \todo: explain this in more depth, maybe in terms of "ULPs", also, specify in which ranges the
    desiredPrecision holds (the distribution of the absolute precision is different in different
    bases) */
    static rsUint32 numRequiredDigits(double desiredPrecisionInBits, rsUint64 base)
    {
       return rsCeilInt(desiredPrecisionInBits / rsLog2((double) base));
    }


    /** \name Conversion */

    /** Returns a string with the digits of this number in some base. */
    rsString toString(rsUint32 numDigits, rsUint32 stringBase = 10) const;

    /** Sets the value of this number from a decimal c-string into an rsBigFloat to the base 10,
    using as much digits as present in the string. The function also supports exponent-notation,
    like 1.23e2 for 123,  etc.
    \todo allow string to be given in a base other than 10, then use this other base for internal
    representaion as well (should be straightforward - see comments in implementation) */
    void fromString(const char numberString[] /*, rsUint64 baseToUse = 10*/);

    /** Returns this number with another precision setting. */
    rsBigFloat withPrecision(rsUint32 precisionToUse) const;
      // replace with a function setPrecision that directly operates on "this" instead of returning
      // another

    /** Returns this number in another base.
    \todo
    If the optional parameter "exact" is true (which is the default), the converted number will be
    as exact as possible - that is: the converted value will differ from this value by at most one
    half of a unit-in-the last place (of the new representation). If false is passed, we will use
    reduced precision in the internal calculations (determined by the minimum precision of the
    source- and target number) - background: if, for example, we want to convert number with 1000
    decimal digits to 53 binary digits (IEEE double) or vice versa, it would be overkill to do a
    conversion that is exact up to the last place. This is especially true, if a 1000-decimal-digit
    number is obtained from a 53-binary-digit number and then used as an initial approximation for
    some function evaluation with 1000 decimal digits precision (initial approximations based on
    built-in C++ IEEE-double-precision functions are regularly done in our iterative,
    high-precision function evaluation routines).  */
    rsBigFloat toBase(rsUint64 newBase, rsUint32 numDigitsToUse, bool exact = true) const;

    /** Returns an rsBigFloat variable that mimics the behavior of the built-in C++ type for double
    precision floating point values, i.e. it uses base 2 and has a precision of 53 (binary) digits
    for the fraction.
    \todo: also set up the minExponent/maxExponent to -1022, +1023 -> we need this then as member
    variable ...but maybe not  */
    rsBigFloat mimicBuiltInType(double value) const;

    /** Sets the value of this number from a C++ double precision floating point number.
    maybe rename to setFromDouble. */
    void fromDouble(double value, rsUint32 precisionToUse = 2, rsUint64 baseToUse = maxBase);

    /** Returns the value of this number as a C++ double precision floating point number. */
    double toDouble() const;


    /** \name Computation of mathematical constants */

    /** Computes and returns pi in the given precision and base. If you need use it repeatedly, you
    should call this function only once and store and re-use the result. */
    static rsBigFloat computePi(rsUint32 precisionToUse, rsUint64 baseToUse);







    static const rsInt32 minExponent = -1022; // use minimum value of rsInt32 later
    static const rsInt32 maxExponent =  1023; // use maximum value of rsInt32 later



    // 3 different versions of the addition routine and the 4th "add" function can switch between
    // them at compile-time (i.e. comment/uncomment calls to the other functions). the 3rd
    // routine handles all cases uniformly, but comes at the cost of always using a full
    // zero-extended addition whereas the 1st and 2nd are for 2 complementary cases depending
    // on whether or not zero extension is required for the operands (this occurs only if a
    // to-be-rightshifted operand is being subtracted from a positive other operand). i don't know,
    // which strategy (switching between 2 cases or treating all cases uniformly) is best, we need
    // some performance tests:
    static void subtractZeroExtended(     const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& d);
    static void addOrSubtractWithTailCopy(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& r);
    static void addOrSubtractAllCases(    const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& r);
    static void add(                      const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& s);


    static void multiply(                 const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& p);
    static void divide(                   const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& q);




  protected:

    /** Normalizes the fraction such that the 1st digit is nonzero and adjusts the exponent
    accordingly. */
    void normalize();

    void incrementExponent(rsUint32 amount);

    void decrementExponent(rsUint32 amount);

    void handleRoundingOverflow();




    /** \name Data */

    rsInt32   exponent;




  };




  /** Under construction... */
  rsBigFloat RSLib_API rsSqrt(const rsBigFloat& x);
  rsBigFloat RSLib_API rsAGM(const rsBigFloat& a, const rsBigFloat& b);



}

#endif
