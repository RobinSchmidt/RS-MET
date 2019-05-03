#ifndef PrimeNumbers_h
#define PrimeNumbers_h

#include "Definitions.h"
#include "MoreMath.h"
using namespace MoreMath;

/**

This is a class which can check a number for primality or obtain a prime 
number which is closest to a given number. It is based on precomputed table of
the first 10000 prime-numbers (the 10000th prime is 104729) - it is only 
intended to work for numbers lower/equal this number (at the moment).

*/

class PrimeNumbers
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 PrimeNumbers();
 ~PrimeNumbers();

 //---------------------------------------------------------------------------
 // processing:

 bool isPrime(int someNumber);
 /**< Returns true if the input number is prime, returns false if it is 
      composite or out of range. */

 int findClosestLowerPrime(int someNumber);
 /**< Returns the prime-number which is closest to and lower than (or euqal to)
      the input-number. */

 //===========================================================================

protected:

 // this is the length of our array:
 static const intA primeArrayLength = 10000;

 // this is precomputed table with the first 10000 prime-numbers:
 static intA primeArray[10000];

};

#endif // PrimeNumbers_h
