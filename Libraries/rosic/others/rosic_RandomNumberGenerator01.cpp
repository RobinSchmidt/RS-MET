#include "rosic_RandomNumberGenerator01.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RandomNumberGenerator01::RandomNumberGenerator01()
{
  a = 0;
  x = 0;
  c = 0;
  m = 32768;  // 32768=2^15

  setFactor(36739);
  setAdditiveConstant(17389);
  setState(0);
}

RandomNumberGenerator01::~RandomNumberGenerator01()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void RandomNumberGenerator01::setFactor(unsigned long newFactor)
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

void RandomNumberGenerator01::setState(unsigned long newState)
{
  x = newState % m;
}

void RandomNumberGenerator01::setAdditiveConstant(unsigned long newConstant)
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

void RandomNumberGenerator01::setModulus(unsigned long newModulus)
{
  m = newModulus;
}

void RandomNumberGenerator01::setStateFromString(char *theString)
{
  if( theString != NULL )
  {
    // caculate an intermediate number based on the characters in the string and use it as seed:
    int length = strlen(theString);
    unsigned long long tmp = 1; // can hold numbers in the range 0...(2^64)-1
    for(int i=0; i<=length; i++)
      tmp += (unsigned long long) (i+1) * abs(theString[i]); // accumulation

    // restrict the result to the permitted range - this gives the state:
    x = (unsigned long) (tmp % m);
  }
}

//-------------------------------------------------------------------------------------------------
// random number generation:

unsigned long RandomNumberGenerator01::getRandomNumber()
{
  x = (a*x + c) % m; // calculate output and update the state
  return x;
}
