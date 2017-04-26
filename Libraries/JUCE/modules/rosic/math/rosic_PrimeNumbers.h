#ifndef rosic_PrimeNumbers_h
#define rosic_PrimeNumbers_h

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"
#include <math.h>

namespace rosic
{

  /**

  This is a class which can check a number for primality or obtain a prime 
  number which is closest to a given number. It is based on precomputed table of
  the first 10000 prime-numbers (the 10000th prime is 104729) - it is only 
  intended to work for numbers lower/equal this number (at the moment).

  */

  class PrimeNumbers
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    PrimeNumbers();
    ~PrimeNumbers();

    //---------------------------------------------------------------------------------------------
    // processing:


    /** Returns true if the input number is prime, returns false if it is composite or out of 
    range. */
    static bool isPrime(int someNumber);

    /** Returns the prime-number which is closest to and lower than (or euqal to) the 
    input-number. */
    static int findClosestLowerPrime(int someNumber);

    /** Returns the prime-number which is closest to (or euqal to) the input-number - if two 
    prime-numbers are equallly close, it will return the lower of them. */
    static int findClosestPrime(int someNumber);

    /** Finds the index of the cloesest prime inside our array of primes. */
    static int findClosestLowerPrimeIndex(int someNumber);

    /** Returns the (index+1)th prime number, that is: getPrime(0)==2, getPrime(1)==3, etc. The 
    index has to be <= 10000, otherwise an access violation will occur. */
    static int getPrime(int index) { return primeArray[index]; }

    //=============================================================================================

  protected:

    // this is the length of our array:
    static const intA primeArrayLength = 10000;

    // this is precomputed table with the first 10000 prime-numbers:
    static intA primeArray[10000];

  };

}  // end namespace rosic

#endif // rosic_PrimeNumbers_h
