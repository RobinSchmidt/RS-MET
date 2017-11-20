#ifndef RS_BIGINT_H
#define RS_BIGINT_H

namespace RSLib
{

  /**

  This is a class for representing (signed) integers with a larger range of values than the 
  built-in C++ integer types. It uses a digit representation in a positional number system with an 
  arbitrary base up to 2^32. That means, each digit can represent values from 0 to base-1 and is 
  internally represented as an unsigned 32-bit integer.

  Internally, an rsBigInt number will always use just as many digits as are strictly required to
  represent the number, i.e. there are never any leading zero digits.

  \todo rewrite the comparison operators like it is done in rsBigFloat - that seems neater, and 
  maybe they can even be factored out into the baseclass rsBigNumber (using a purely virtual
  compare-function there

  \todo have functions that parallel those in functions/IntegerFunctions.h - but that would mean 
  code-duplication - better: templatize the IntegerFunctions, so they can take an rsBigInt as well

  \todo write a class rsRational<T IntType> that represents rational numbers

  */

  class RSLib_API rsBigInt : public rsBigNumber
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. You may initialize the number by passing some 64-bit integer. */
    rsBigInt(rsInt64 initialValue = 0, rsUint64 baseToUse = maxBase);

    /** Constructor. Creates a BigInt from an array of (base 2^32) digits. If the array has leading
    zero digits, these will be scrapped off in out internal representation. */
    rsBigInt(rsUint32 *digits, rsUint32 numDigits, bool negative, rsUint64 baseToUse = maxBase);

    /** Copy constructor. */
    rsBigInt(const rsBigInt& other);


    /** \name Operators */

    rsBigInt& operator=(const rsBigInt& other);

    rsBigInt operator-() const;

    bool operator==(const rsBigInt& other) const;
    bool operator!=(const rsBigInt& other) const;
    bool operator>=(const rsBigInt& other) const;
    bool operator<=(const rsBigInt& other) const;
    bool operator>(const rsBigInt& other) const;
    bool operator<(const rsBigInt& other) const;

    rsBigInt operator+(const rsBigInt& other);
    rsBigInt operator-(const rsBigInt& other);
    rsBigInt operator*(const rsBigInt& other);
    rsBigInt operator/(const rsBigInt& other);
    rsBigInt operator%(const rsBigInt& other);

    rsBigInt& operator+=(const rsBigInt& other);
    rsBigInt& operator-=(const rsBigInt& other);
    rsBigInt& operator*=(const rsBigInt& other);
    rsBigInt& operator/=(const rsBigInt& other);
    rsBigInt& operator%=(const rsBigInt& other);

    rsBigInt& operator++();
    rsBigInt& operator--();


    /** \name Setup */

    /** Sets the value of this BigInt from a (signed) 64-bit integer number. */
    void setValue(rsInt64 newValue);


    /** \name Inquiry */

    /** Returns quotient and remainder of integer division. */
    void divMod(const rsBigInt& divisor, rsBigInt& quotient, rsBigInt& remainder) const;

    /** Returns a string representation of this number. For the meaning of the "base" parameter,
    @see fromString */
    rsString toString(rsUint64 base = 10) const;


    /** \name Misc */

    /** Creates a BigInt from a string of digits (and possibly a leading '-'). The digits may be
    given in any base up to 36 where we adopt the convention from hexadecimal notation for the 
    digit-values above 9: A=10, B=11, C=12, D=13, E=14, F=15 and generalize it to the rest of
    the alphabet: G=16, H=17, ..., Z=35. Lowercase letters will also work. If no base is given, 
    decimal digits will be assumed (base parameter defaults to 10). You may also pass a base to be
    used for the internal representation. */
    static rsBigInt fromString(const rsString& s, rsUint64 base = 10, 
                               rsUint64 baseToUse = maxBase);

  };

  /** Converts an integer represented as string of digits (possibly with a minus sign) in some
  base into the representation in another base. The bases should be between 2 and 36 (inclusive).
  It's a convenience function that comes as a little spin-off of the rsBigInt class. It's not the 
  most efficient way to to do it but anyway. */                 
  rsString RSLib_API rsConvertIntegerIntoBase(rsString value, rsUint64 oldBase, rsUint64 newBase);

  /** Explicit template instantiation, to be used by the rsPow template-function. */
  rsBigInt RSLib_API rsUnityValue(rsBigInt value);

}

#endif