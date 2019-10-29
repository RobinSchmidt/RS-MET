#ifndef rosic_HelperFunctions_h
#define rosic_HelperFunctions_h

// get rid - redundant with rapt

namespace rosic
{

  /** Returns an integer that is represented by the bit-reversed input number, taking into account
  only 'numBits' bits in the reversal process. */
  INLINE int bitReverse(int number, int numBits);

  /** Returns an integer that is represented by the bit-reversed input number, taking into account
  only 'numBits' bits in the reversal process. */
  INLINE unsigned long bitReverse(unsigned long number, unsigned long numBits);

  /** Finds the (non-integer) position of the upward zero crossing in the buffer which is nearest
  to 'searchStart'. The fractional part of the zero crossing is determined by linearly 
  interpolating betwen the points below and above zero and solving fo the zero-crossing of the
  connecting line. */
  double findNearestUpwardZeroCrossing(float* buffer, int length, double searchStart);

  /** Function that writes an error message to the console and triggers a debug-break. */
  void error(const char *message);

  //===============================================================================================
  // implementation fo the inline functions:

  INLINE int bitReverse(int number, int numBits)
  {
    int result = 0;
    for(int n=0; n<numBits; n++)
    {
      // leftshift the previous result by one and accept the new LSB of the current number on the 
      // right:
      result   = (result << 1) + (number & 1);

      // rightshift the number to make the second bit from the right to the new LSB:
      number >>= 1;
    }
    return result;
  }

  INLINE unsigned long bitReverse(unsigned long number, unsigned long numBits)
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

} // end namespace rosic

#endif // #ifndef rosic_HelperFunctions_h