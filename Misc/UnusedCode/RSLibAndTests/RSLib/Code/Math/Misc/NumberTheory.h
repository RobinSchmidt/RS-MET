#ifndef RS_NUMBERTHEORY_H
#define RS_NUMBERTHEORY_H

namespace RSLib
{

  /** Returns y = rsIntSqrt(x), the integer square root of x such that
  y*y <= x and (y+1)*(y+1) > x. */
  template <class T>
  T rsIntSqrt(T x);
    // move to IntegerFunctions

  /** Fills the array "primes" with all primes up to "upperLimit" (including "upperLimit", if it's
  a prime itself). It uses an optimized version of the the sieve of Erathostenes. */
  template<class T>
  void rsFindPrimesUpTo(rsArray<T> &primes, T upperLimit);

  /** This function fills the boolean array primeFlags with true/false, where a flag at index "i"
  will be set to true, if start+2*i is a prime, false otherwise. This means, each entry of the
  flag-array corresponds to an odd number start+2*i which implies that start should be itself an
  odd number (we ignore even numbers as these are not prime anyway, except for 2 which should be
  treated separately by the caller). The function requires some primes to be known beforehand and
  passed as "primeTable". This table should go (at least) up to a prime which is greater than or
  equal to rsIntSqrt(start+2*(numFlags-1)). */
  template<class T>
  void rsGetPrimeFlags(bool *primeFlags, T start, rsUint32 numFlags, T *primeTable);

  /** For each true entry with index "i" in the "primeFlags" array of length "numFlags", this
  function appends start+2*i to the "primeTable" starting at "writeStart" with appending new
  values. We use +2*i and not just +i, because it is assumed that the array contains the flags for
  odd numbers only, which also implies that "start" should be odd. The function returns the number
  of items written. The last parameter "tableLength" should be the total length of "primeTable" and
  is required only to avoid writing beyond the end of the table. */
  template<class T>
  rsUint32 rsAppendFlaggedPrimes(bool *primeFlags, T start, rsUint32 numFlags, T *primeTable,
                                 rsUint32 writeStart, rsUint32 tableLength);

  /** Fills the passed "primes" buffer with the 1st "numPrimes" primes. The algorithm needs only
  a constant amount of auxiliary memory determined by the "bufferSize" parameter (the algorithm
  will allocate a temporary buffer of values of type T of that length). Reasonable values are
  1024...32768. In a performance test, 32768 seemed to be the sweet spot (on that particular
  machine), so this value will be used as default value, when 0 is passed. It is not  directly
  used as default value for the parameter, so we may later apply some heuristic to choose a
  machine-dependent default value at runtime.
  \todo: introduce a parameter partiallyFilledUpTo with default argument 0 - to be used to
  conveniently expand a primeTable that is already partially filled. */
  template<class T>
  void rsFillPrimeTable(T *primes, rsUint32 numPrimes, rsUint32 bufferSize = 0);

  /** Given a number x, this function computes the prime-factors of x and their repsective
  exponents such that x = product-of factors[i]^exponents[i]. The optional primeTable parameter,
  if not null, should point to a table of primes up to at least the highest prime-factor in x.
  The largest number, that could possibly be a factor of x is sqrt(x), so if the primeTable
  contains all primes up to and including sqrt(x), you are on the safe side. If a null-pointer is
  passed, such a table will be created temporarily internally (which is expensive, but sometimes
  convenient for quick-and-dirty experimental code). The algorithm is based on simple
  trial-divsion, nothing fancy at all. */
  template<class T>
  void rsPrimeFactors(T x, rsArray<T>& factors, rsArray<T>& exponents,
                      rsArray<T> *primeTable = nullptr);

  /** Given 2 nonnegative integers x, y, (not both 0), this algorithm calculates a, b, g, such that
  a*x + b*y = g = gcd(x,y) by means of the extended GCD algorithm. */
  template<class T>
  void rsEGCD(T x, T y, T& a, T& b, T& g);

  /** Computes the modular power of the given base with respect to the given modulus.
  \todo - get rid of this function, implement a class for modular arithmetic instead and then use
  the general rsPow function. ..or maybe let the function take a T-type for the exponent as well. */
  template <class T>
  T rsModularPow(const T& base, const T& modulus, rsUint64 exponent);

  /** Computes the modular inverse of x in the given modulus m (if existent, otherwise returns
  zero). */
  template <class T>
  T rsModularInverse(const T& x, const T& m);

  /** Computes the modular for prime moduli p (using the powering algorithm). */
  template <class T>
  T rsPrimeModularInverse(const T& x, const T& p);

  /** Alternative algorithm to compute modular inverse for prime modulus p. */
  template <class T>
  T rsPrimeModularInverse2(const T& x, const T& p);


  /** Given a number of moduli, this function computes the weights in the weighted sum in the
  chinese remainder theorem. These are values which are typically precomputed for a set of moduli
  in order to be applied to multiple sets of remainders. */
  template <class T>
  T rsChineseRemainderWeights(T* moduli, T* weights, rsUint32 count);

  /** Given a set of remainders, a set of precomputed weights and the product-modulus (product of
  the moduli with respect to which the remainders are taken), this function returns the result of
  applying the chinese remainder theorem to these remainders. */
  template <class T>
  T rsApplyChineseRemainderTheorem(T* remainders, T* weights, T modulus, rsUint32 count);

  /** Given a set of remainders (i.e. numbers modulo some modulus) and the corresponding moduli
  (which are assumed to be pairwise coprime - i.e. no pair has a common divisor other than 1),
  this function returns the solution R to the set of equations: R % m[i] == r[i]. The solution
  is actually a whole equivalence of numbers which are congruent modulo M, where M is the product
  of the moduli. The function returns the smallest positive representative. It's a convenience
  function, using rsChineseRemainderWeights and rsApplyChineseRemainderTheorem, thereby merging the
  precomputations and the actual application into one call. If the theorem is to be applied to
  different sets of remainders with respect to the same set of moduli and efficiency is a concern,
  you should do the precomputaion (of the weights) only once. */
  template <class T>
  T rsChineseRemainderTheorem(T* remainders, T* moduli, rsUint32 count);


  // egcd, inverse and crt should be called with signed types as template-parameter


  // \todo: write a class rsPrimeFactorization that allows for multiplication and division:
  // a factorization is multiplied by another by adding exponents of corresponding prime-factors,
  // division is done by subtracting exponents (if a result of such a subtraction is < 0, this
  // indicates, that the number is not divisible by the divisor). question: is it possible in this
  // case, to obtain the correct integer-division result and remainder? but we may simply allow for
  // negative exponents as well - then, a number is integer, iff none of the expoenents is
  // negative. addition and subtraction would require to convert into an an actual number,
  // add/subtract and factor the result - perhaps, we should not provide operators for that

  // easy test for checking that p is nonprime: if p is prime, then 24 divides p^2-1 ->
  // if 24 doesn't divide p^2-1, p can't be prime. however, if p is nonprime, 24 might still
  // divide p^2-1 (an example is p=25,p=35,p=49) - do such nonprimes occur often?
  // or can we say something else about them, draw conclusions from the primality of the quotient,
  // maybe?
  // http://puzzles.nigelcoldwell.co.uk/fifteen.htm
  // can be used as preliminary check to sort out many non-primes before doing a more elaborate
  // primality check. i order to avoid overflow in the square, we could first check if either
  // p-1 or p+1 is divisible by 4 and if so, check if either p-1 or p+1 is divisible by 3, if
  // so, p^2-1 is divisible by 24 (p+1 could still overflow but only if p=maxInt)

}

#endif
