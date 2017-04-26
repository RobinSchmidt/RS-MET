#include "rosic_RandomNumberGenerator02.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RandomNumberGenerator02::RandomNumberGenerator02()
{
  x = 0;
  c = 0xC6A939C5;
  m = 32768;          // 32768=2^15
}

RandomNumberGenerator02::~RandomNumberGenerator02()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void RandomNumberGenerator02::setState(unsigned long newState)
{
  x = newState % m;
}

//-------------------------------------------------------------------------------------------------
// processing:

unsigned long RandomNumberGenerator02::getRandomNumber()
{
  x = (x + (x*x | c)) % m; // calculate output and update the state
  return x;
}




