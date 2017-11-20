using namespace RSLib;

unsigned long RSLib::rsBitReverse(unsigned long number, unsigned long numBits)
{
  unsigned long result = 0;
  for(unsigned long n=0; n<numBits; n++)
  {
    // leftshift the previous result by one and accept the new LSB of the current number on the
    // right:
    result   = (result << 1) + (number & 1);

    // rightshift the number to make the second bit from the right to the new LSB:
    number >>= 1;
  }
  return result;
}

void RSLib::rsSetExponentOfDouble(double *value, int exponent)
{
  rsUint64 *x = reinterpret_cast<rsUint64*> (value);
  rsUint64 biasedExponent = ((rsUint64) (exponent + 1023)) << 52;
  static const rsUint64 mask = 0x7ff0000000000000ULL; 
  *x &= ~mask;                    // zero out old exponent
  *x |= (mask & biasedExponent);  // set new exponent
}



/*
namespace RSLib
{
   
  unsigned long rsBitReverse(unsigned long number, unsigned long numBits)
  {
    unsigned long result = 0;
    for(unsigned long n=0; n<numBits; n++)
    {
      // leftshift the previous result by one and accept the new LSB of the current number on the
      // right:
      result   = (result << 1) + (number & 1);

      // rightshift the number to make the second bit from the right to the new LSB:
      number >>= 1;
    }
    return result;
  }
}
*/
