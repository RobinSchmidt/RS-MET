#include "RandomNumberGenerator.h"

//----------------------------------------------------------------------------
// construction/destruction:

RandomNumberGenerator::RandomNumberGenerator()
{
 a = 0;
 x = 0;
 c = 0;
 m = 32768;  // 32768=2^15
}

RandomNumberGenerator::~RandomNumberGenerator()
{

}

//----------------------------------------------------------------------------
// setup:

void RandomNumberGenerator::setFactor(unsigned long newFactor)
{
 unsigned long tmp = newFactor;

 // restrict the range of the result to 0...m:
 tmp %= m;

 // make sure that the condition a % 4 = 1 holds, maintaining the 
 // condition a <= m:
 unsigned long tmp2 = (unsigned long) tmp % 4;
 if( tmp2 == 0 )
  a = (unsigned long) ((tmp+1) % m);
 else if( tmp2 == 1 )
  a = (unsigned long) tmp;
 else if( tmp2 == 2 )
  a = (unsigned long) ((tmp+3) % m);
 else if( tmp2 == 3 )
  a = (unsigned long) ((tmp+2) % m);

 // make sure, that a assumes some minimum value:
 if( a < 4000 )
  a += 4000;
}

void RandomNumberGenerator::setState(unsigned long newState)
{
 x = newState % m;
}

void RandomNumberGenerator::setAdditiveConstant(unsigned long newConstant)
{
 unsigned long tmp = newConstant;

 // restrict the range of the result to 0...m:
 tmp %= m;

 // make sure that c is odd and c <= m:
 if( tmp%2 == 1 ) // tmp is odd and can directly be used as c
  c = (unsigned long) tmp;
 else             // tmp is even, we need to add 1 (mod m) to make it odd
  c = (unsigned long) ((tmp+1) % m);

 // make sure, that c assumes some minimum value:
 if( c < 3000 )
  c += 3000;
}

void RandomNumberGenerator::setModulus(unsigned long newModulus)
{
 m = newModulus;
}

unsigned long RandomNumberGenerator::getRandomNumber()
{
 x = (a*x + c) % m; // calculate output and update the state
 return x;
}

void RandomNumberGenerator::setFactorFromString(char *theString)
{
 // caculate an intermediate number as the product of the first 8 characters
 // (by interpreting them as 8-bit numbers), the result will be always in the 
 // range 0...255^8 < (2^64)-1
 unsigned long long tmp = 1; // can hold numbers in the range 0...(2^64)-1
 for(int i=0; i<8; i++)
  tmp *= (unsigned long long) theString[i]; // multiplicative accumulation

 // restrict the range of the result to 0...m:
 tmp %= m;

 // make sure that the condition a % 4 = 1 holds, maintaining the 
 // condition a <= m:
 unsigned long tmp2 = (unsigned long) tmp % 4;
 if( tmp2 == 0 )
  a = (unsigned long) ((tmp+1) % m);
 else if( tmp2 == 1 )
  a = (unsigned long) tmp;
 else if( tmp2 == 2 )
  a = (unsigned long) ((tmp+3) % m);
 else if( tmp2 == 3 )
  a = (unsigned long) ((tmp+2) % m);

 // make sure, that a assumes some minimum value:
 if( a < 4000 )
  a += 4000;
}

void RandomNumberGenerator::setStateFromStringAndNumber(char *theString, 
                                                     unsigned long theNumber)
{
 // caculate an intermediate number as the product of the first 8 characters
 // (by interpreting them as 8-bit numbers), the result will be always in the 
 // range 0...255^8 < (2^64)-1
 unsigned long long tmp = 1; // can hold numbers in the range 0...(2^64)-1
 for(int i=0; i<8; i++)
 {
  tmp *= (unsigned long long) theString[i]; // multiplicative accumulation
 }

 // restrict result to 0...m:
 tmp %= m;

 // add the number, and restrict the result again to 
 // derive the state:
 x = (unsigned long) (tmp+theNumber % m);
}

void RandomNumberGenerator::setAdditiveConstantFromString(char *theString)
{
 // caculate an intermediate number as the product of the first 8 characters
 // (by interpreting them as 8-bit numbers), the result will be always in the 
 // range 0...255^8 < (2^64)-1
 unsigned long long tmp = 1; // can hold numbers in the range 0...(2^64)-1
 for(int i=0; i<8; i++)
  tmp *= (unsigned long long) theString[i]; // multiplicative accumulation

 // restrict the range of the result to 0...m:
 tmp %= m;

 // make sure that c is odd and c <= m:
 if( tmp%2 == 1 ) // tmp is odd and can directly be used as c
  c = (unsigned long) tmp;
 else             // tmp is even, we need to add 1 (mod m) to make it odd
  c = (unsigned long) ((tmp+1) % m);

 // make sure, that c assumes some minimum value:
 if( c < 3000 )
  c += 3000;
}