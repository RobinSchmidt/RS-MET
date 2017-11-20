#ifndef RS_BIGNUMBER_H
#define RS_BIGNUMBER_H

namespace RSLib
{

  /**

  Baseclass for arbitrary-size integers and arbitrary-precision floating-point numbers, providing
  common data and member-functions as well as static functions for performing arithmetic on arrays
  of digits.

  Such digit arrays may either represent integers or fractions in an arbitrary base (up to 2^32).
  An integer x, represented as an array XYZ (XYZ shall denote the 3 array elements of a 3-digit
  value) in a given base B is given by: x = X*B^2 + Y*B^1 + Z*B^0. If the array represents a
  fraction, it is given by: x = X/B^1 + Y/B^2 + Z/B^3. These functions constitute the computational
  core for big number arithmetic, but may also be used their own right.

  References:
  (1): The Art of Computer Programming, Vol. 2, 2nd Ed, Ch. 4

  \todo reorder functions in the C++ file to the same order as in the header

  \todo provide optimized versions of add, subtract, ... for the case when the base is 2^32 - in
  this case, the modulo operations are often implicit, call these addMaxBase etc., rename the
  generic versions addGeneric, .etc. then, in the add functions, switch between generic and
  optimized version according to the base - maybe not, because when using very large numbers, we
  want to do NTT-multiplication anyway which requires a prime base, so the special case of 2^32
  is not used for very large numbers anyway

  \todo optimize division: get rid of the internal memory allocation by passing an (optional)
  paramater divide(..., rsUint32 *workArea = nullptr) which should point to an array of size
  3*(N+1), if it's not null. in the division algorithm, the pointer is checked and if it's
  non-null, it's used for the intermediate arrays, oherwise, internal memory allocation will occur
  do this also for divideFractions, require enough space for the function it self and the
  internally use divide - pass through part of the memory

  \todo implement number-theoretic-transform multiplication algorithm, use it, when it's more
  efficient (switch, based on the number of digits of the operands) - refer to "Matters
  Computational" by Jörg Arndt for research.

  \todo maybe ditch the default argument of 10 for the base - we'll see



  */

  class RSLib_API rsBigNumber
  {

  public:

    /** \name Construction/Destruction */

    /** Destructor.  */
    ~rsBigNumber();



    /** \name Operators */



    /** \name Setup */



    /** \name Inquiry */

    /** Returns the base in which this number is represented. */
    inline rsUint64 getBase() const
    {
      return base;
    }

    /** Returns the number of 32-bit words used to represent this number. */
    inline rsUint32 getNumDigits() const
    {
      return numDigits;
    }

    /** Copies the digits of this number into the passed buffer. The buffer should be large enough
    to hold the number of digits that this BigInt object has. */
    inline void copyDigitsToBuffer(rsUint32 *buffer) const
    {
      rsCopyBuffer(digits, buffer, numDigits);
    }

    /** Returns true, iff this number is zero. */
    inline bool isZero() const
    {
      return rsIsAllZeros(digits, numDigits);
    }

    /** Returns true, if this number is negative, false otherwise. */
    inline bool isNegative() const
    {
      return negative;
    }

    /** Returns true, if this number has the same sign as the other one, false otherwise. */
    inline bool hasSameSignAs(const rsBigNumber& other) const
    {
      return this->negative == other.negative;
    }


    /** \name Arithmetic */

    /** Assuming that a, b are arrays of digits to the given base, this function computes the
    sum a + b and stores it in s: s = a + b. It is assumed that a, b and s have N places - you
    should prepend leading zeros for small numbers, if necessary. The sum s may have N+1 digits due
    to a carry digit in the leftmost place. This carry digit will be the returned in the return
    value (it will be either 0 or 1). s may point to the same location as b or a. */
    static rsUint32 add(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *s, rsUint64 base = 10);

    /** Adds the single-digit number b to a, stores it in s and returns whether or not a carry
    occurred. a and s may point to the same location. */
    static rsUint32 add(rsUint32 *a, rsUint64 b, rsUint32 N, rsUint32 *s, rsUint64 base = 10);

    /** Subtracts a from b and stores the difference in d: d = a - b. a, b and d have N places and
    it is assumed that a >= b. The return value will be 1, if there was a borrow in the leftmost
    place, 0 otherwise. If it is 1, the result is invalid (occurs, when a < b). d may point to the
    same location as b or a.*/
    static rsUint32 subtract(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *d,
                             rsUint64 base = 10);

    /** Subtracts the single-digit number b from a. a and d may point to the same location.  */
    static rsUint32 subtract(rsUint32 *a, rsUint32 b, rsUint32 N, rsUint32 *d, rsUint64 base = 10);

    /** Multiplies a and b and stores the product in p: p = a * b. a and b have N digits and the
    product p has 2*N digits. If a and b represent fractions instead of integers, the 1st half of
    p will contain the digits at the same positions as the inputs and the 2nd half will contain
    digits at the less significant places. Note that p must point to a different location than a
    and b (the function can't be used 'in-place') */
    static void multiply(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *p, rsUint64 base = 10);

    /** Version of the multiplication function for when a and b have a different number of digits
    Na, Nb, respectively. The product p will then have has Na + Nb digits.  */
    static void multiply(rsUint32 *a, rsUint32 Na, rsUint32 *b, rsUint32 Nb, rsUint32 *p,
                         rsUint64 base = 10);

    /** Multiplies the N-digit number a by the single-digit number b and returns the result in p
    and the carry-digit for the leftmost place in the return value. a and p have N places and a
    and p may point to the same location. */
    static rsUint32 multiply(rsUint32 *a, rsUint32 N, rsUint64 b, rsUint32 *p, rsUint64 base = 10);

    /** Divides a by b and stores the quotient in q and the remainder in r. a, b, q, and r have N
    places. The function can be used in-place. a and b may point to the same locations as q and r
    in both possible orders) */
    static void divide(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *q, rsUint32 *r,
                       rsUint64 base = 10);

    /** Divides the N-digit number a by the single-digit number b and returns the result in q and
    the remainder in the return value. a and q have N places. a and q may point to the same
    location. */
    static rsUint32 divide(rsUint32 *a, rsUint32 N, rsUint64 b, rsUint32 *q, rsUint64 base = 10);

    /** Divides the fractions represented by the length-N arrays a and b and stores the result in
    the length 2*N array q. The (most significant) 1st half of q will contain digits before the
    point and the second half contains the digits after the point. Note, that the result may have
    been rounded due to the fact, that division may produce nonterminating digit sequences. This,
    in turn, may cause rounding overflow. The return value informs, if this was the case. If so,
    the newDigits sequence will contain all zeros and you should imagine that a leading 1 should be
    prepended - so the result actually represents unity.

    naaah! this is wrong - if the result is unity, we will have overflow into q[N-1] - we can never
    have overflow into q[-1]. the function should be void: overflow of the final result is
    impossible
    */
    static bool divideFractions(rsUint32 *a, rsUint32 *b, rsUint32 *q, rsUint32 N,
                                rsUint64 base = 10);

    /** Version of the fraction-division function for when a and b may have a different number of
    digits Na, Nb, respectively. The quotient q will then have Nq digits where Nq should
    be >= Na + Nb. The Nb leading digits in q represent the part before the point, the remaining
    Nq-Nb digits represent the part after the point. */
    static bool divideFractions(rsUint32 *a, rsUint32 Na, rsUint32 *b, rsUint32 Nb,
                                rsUint32 *q, rsUint32 Nq, rsUint64 base = 10);



    /** \name Conversion */

    /** Converts an unsigned integer into a digit-array-representation in the given base. */
    static void digitsFromInt(rsUint64 value, rsUint32 *digits, rsUint32 N, rsUint64 base = 10);

    /** Converts a digit-array-representation in the given base into an unsigned integer. */
    static rsUint64 digitsToInt(rsUint32 *digits, rsUint32 N, rsUint64 base = 10);

    /** Converts a digit-array-representation of an integer in the given (old) base into a
    representation in another (new) base. The resulting digit array as returned as rsArray (mainly
    because the number of digits in the new digit-array is not known beforehand the rsArray class
    conveniently encapsulates this number). */
    static rsArray<rsUint32> changeBaseForInteger(rsUint32 *digits, rsUint32 numDigits,
                                                  rsUint64 oldBase, rsUint64 newBase);

    /** Converts a digit-array-representation of a fraction in the given (old) base into a
    representation in another (new) base. This may require a rounding operation
    @see divideFractions how to deal with that. */
    static bool changeBaseForFraction(rsUint32 *oldDigits, rsUint32 oldNumDigits, rsUint64 oldBase,
                                      rsUint32 *newDigits, rsUint32 newNumDigits, rsUint64 newBase);



    /** \name Rounding */

    /** Rounds the digit array in-place and returns whether or not rounding overflow occurred. */
    static bool round(rsUint32 *digits, rsUint32 oldNumDigits, rsUint32 newNumDigits,
                      rsUint64 base, bool zeroOutTail = true);

    /** If some fraction x is represented by the 5-digit-array ABCDE and we want to round it to
    only 3 digits, we call the sequence ABC to be retained the "head" and the sequence DE to be cut
    off the "tail". This function takes the head- and tail sequences seperately and determines
    whether upward- or downward-rounding is desired. If upward rounding is desired, the head
    sequence will be incremented by 1, otherwise it will remain as is. In the former case,
    rounding-overflow may occur. The return value informs, if this was the case. Sometimes (notably
    in base-conversion), it's convenient to allow different bases for head and tail, that's why
    there are separate base parameters for each. */
    static bool roundHeadAccordingToTail(rsUint32* head, rsUint32 headLength, rsUint64 headBase,
                                         rsUint32 *tail, rsUint32 tailLength, rsUint64 tailBase);

    /** Indicates, based on the tail-sequence and the last digit of the head-sequence, whether the
    head should be rounded up (incremented by 1) or down (left as is). In ambiguous cases (where
    the tail is exactly equal to 1/2), the midwayRoundingRule() is applied. */
    static bool upwardRoundRequired(rsUint32 *tail, rsUint32 tailLength, rsUint64 tailBase,
                                    rsUint32 lastHeadDigit, rsUint64 headBase);

    /** If some fractional value x = 0.XYZ, where XYZ is an array of digits (of length 3 in this
    example), this function returns +1 if x > 1/2, -1 if x < 1/2 and 0 if x == 1/2, independently
    from the base. For example, in base 10, it would return 0, iff X=5,Y=0,Z=0, return -1 iff X<5
    and return +1 iff (X>=5 and (Y>0 or Z>0)). In base 8, the same would be true, if you replace
    5 with 4. This function is useful in determining a desired rounding direction: a return value
    of +1 would indicate to round upward, -1 to round downward and for 0, it would be undecided.
    The latter case occurs if some number would be exactly midway between 2 representable numbers
    in which special rules, like round-to-even or round-to-odd could be applied. */
    static int compareToOneHalf(rsUint32 *digits, rsUint32 numDigits, rsUint64 base);

    /** Given the last digit of the head-sequence, this function returns whether or not to round up
    in cases, where the tail sequence is exactly 1/2. It applies the following rule: if base/2 is
    even round to odd, if base/2 is odd round to even, according to the recommendation in (1)
    page 222. */
    static bool shouldRoundUpAtMidway(rsUint32 lastHeadDigit, rsUint64 base);


    /** \name Misc */

    /** Converts an ASCII character into a rsUint32 number (used by fromString). */
    static rsUint32 rsCharToUint32(char c);

    /** Converts a rsUint32 number character into an ASCII (used by toString). */
    static char rsUint32ToChar(rsUint32 x);

    /** Maximum value for the base that can be used. 4294967296 = 2^32, which is one more than the
    largest representable number for an rsUint32. */
    static const rsUint64 maxBase = 4294967296;

    /** Base, that is to be used, when number-theoretic-transform (NTT) multiplication should be
    used. 3489660929 = 13 * 2^28 + 1 is a prime number suitable for NTTs up to length 2^28. It is
    the largest prime of the form k * 2^n + 1 that is <= maxBase (is that true?). The base being a
    prime representable in this particular form is required to make it suitable for NTTs.
    For details, see Jörg Arndt - Matters Computational, Ch.26 (Number theoretic transforms).
    (this is not yet used - NTT-multiplication is not yet implemented). */
    static const rsUint64 nttBase = 3489660929;


    //

    /** Maximum number of digits, for which NTT-multiplication is possible when using nttBase as
    base. 134217728 = 2^28/2 = 2^27 is half of the maximum NTT-size. It's only half of it, because
    we must avoid cyclic convolution, so we must do the NTT with suitably zero-padded arrays.  */
    static const rsUint32 nttMaxNumDigits = 134217728;

    // 2147483647 = 2^31 - 1 is also a prime http://en.wikipedia.org/wiki/2147483647

    // from APFloat-doc: 3221225473 = 2^32 - 2^30 + 1 = 3 * 2^30 + 1 . this is smaller than
    // 3489660929 but allows for NTTs up to length 2^30


    // maybe have purely virtual functions add, subtract, ... that are called by the operators
    // maybe then it's possible to factor the operators out into the baseclass? ..but i don't
    // think so because they have to return a value of the subclass

  protected:

    /** Allocates memory for the desired number of digits. */
    void allocateMemory();

    /** Frees the memory for the digits. */
    void freeMemory();


    /** \name Data */

    rsUint32  *digits;
    rsUint32  numDigits;
    rsUint64  base;
    bool      negative;

  };

}

#endif
