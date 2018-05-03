#ifndef RS_NUMBERMANIPULATIONS_H
#define RS_NUMBERMANIPULATIONS_H

namespace RSLib
{

  // maybe rename to BitManipulations, maybe wrap a class around it (with static functions)

  /** \todo: check the linux versions....they have been inserted here sloppily in order to make it 
  compile when porting. */

  /** Converts a double precision floating point number to an integer using the FPU's current
  rounding mode (which is 'to nearest even integer' by default).  */
  /*
  RS_INLINE int toInt(double x)
  {
    int a;
    ASM(fld x);
    ASM(fistp a);
    return (a);
  }
  */








  /** Returns an integer that is represented by the bit-reversed input number, taking into account
  only 'numBits' bits in the reversal process. */
  RSLib_API unsigned long rsBitReverse(unsigned long number, unsigned long numBits);

  /** Returns the exponent of a 64 bit IEEE 754 floating point number. */
  RS_INLINE int rsExtractExponentFromDouble(double x)
  {
    return (int) (((*((reinterpret_cast<rsUint64 *>(&x)))&0x7FFFFFFFFFFFFFFFULL)>>52)-1023);
  }

  /** Returns the exponent of a 32 bit IEEE 754 floating point number. */
  RS_INLINE int rsExtractExponentFromFloat(float x)
  {
    return (((*((reinterpret_cast<rsUint32 *>(&x)))&0x7FFFFFFF)>>23)-127);
  }

  /** Fills an array of rsUint32 values with the bit pattern of a 64-bit datatype. */
  template <class T>
  void rsGetBits64(T value, rsUint32 bits[64])
  {
    // \todo convert endianness (maybe, test on different machines)
    rsUint64 *x   = reinterpret_cast<rsUint64*> (&value);
    rsUint64 mask = 0x8000000000000000ULL;
    for(int i = 0; i < 64; i++)
      bits[i] = (rsUint8) (((*x<<i) & mask) != 0);
  }

  /** Returns zero, if the i-th bit of x is zero, otherwise returns a value where only the i-th bit 
  is set. */
  template <class T>
  RS_INLINE T rsGetMaskedBit(T x, T i)
  {
    return x & T(1) << i;
  }

  /** Returns the number of true (1) bits in x. */
  RS_INLINE rsUint64 rsGetNumTrueBits64(rsUint64 x)
  {
    x -=  (x>>1)      & 0x5555555555555555UL;
    x  = ((x>>2)      & 0x3333333333333333UL) + (x & 0x3333333333333333UL);
    x  = ((x>>4) + x) & 0x0f0f0f0f0f0f0f0fUL;
    x *= 0x0101010101010101UL;
    return x>>56;
  }

  /** Returns true, iff the i-th bit of x is set ot 1 (true). */
  template <class T>
  RS_INLINE bool rsIsBitTrue(T x, T i)
  {
    return rsGetMaskedBit(x, i) != 0;
  }

  /** Sets bit pattern of a 64-bit datatype from an ofarray rsUint32 values (assumed to be either 0 
  or 1). */
  template <class T>
  void rsSetBits64(T *value, rsUint32 bits[64])
  {
    // \todo convert endianness (maybe, test on different machines)
    rsUint64 *x = reinterpret_cast<rsUint64*> (value);
    rsUint64 mask = 0x8000000000000000ULL;
    *x = 0;
    for(int i = 63; i >= 0; i--)
    {
      *x = *x >> 1;
      if( bits[i] != 0 )
        *x |= mask;
    }
  }

  /** Sets the i-th bit of x to false. */
  template <class T>
  RS_INLINE void rsSetBitFalse(T& x, T i)
  {
    x &= ~(T(1) << i);
  }

  /** Sets the i-th bit of x to true. */
  template <class T>
  RS_INLINE void  rsSetBitTrue(T& x, T i) 
  {
    x |= T(1) << i;
  }

  /** Switches the i-th bit of x. */  
  template <class T>
  RS_INLINE void rsSwitchBit(T& x, T i)
  {
    x ^= T(1) << i;
  }




  // \todo write analoguous functions for 32-bit types, factor out common code and then maybe also 
  // provide such functions for 16 -bit and 8-bit types


  /** Sets up the exponent of a 64 bit IEEE 754 floating point number. The exponent parameter 
  should be the desired exponent without bias, in the range: -1022 <= exponent <= 1023. */
  void rsSetExponentOfDouble(double *value, int exponent);

  /** Converts a 16 bit unsigned integer's endianess from big to little of vice versa. */
  RS_INLINE void rsSwapEndians16(rsUint16 &value)
  {
    value = ((value >> 8) & 0x00FF) | ((value << 8) & 0xFF00);
  }

  /** Converts a 32 bit unsigned integer's endianess from big to little of vice versa. */
  RS_INLINE void rsSwapEndians32(rsUint32 &value)
  {
    value = ((value >> 24) & 0x000000FF) | ((value >> 8)  & 0x0000FF00) 
           | ((value << 8) & 0x00FF0000) | ((value << 24) & 0xFF000000);
  }

  /** Applies rsSwapEndians16 to a buffer of numbers. */
  RS_INLINE void rsSwapEndiansInBuffer16(rsUint16 *buffer, int length)
  {
    for(int i = 0; i < length; i++)
      rsSwapEndians16(buffer[i]);
  }

  /** Returns true, if x is denormalized number, false otherwise. */
  RS_INLINE bool rsIsDenormal(double x)
  {
    if( x != 0.0 && fabs(x) < RS_MIN(double) ) 
      return true;
    return false;
  }
  // i have found this for single precision floats here http://musicdsp.org/files/denormal.pdf,
  // maybe that's more efficient, but we need to convert it to the corresponding double version,
  // (and not kill, but just return biased_exponent == 0 && abs_mantissa != 0):
  // void test_and_kill_denormal (float &val)
  //  {
  //   const int x = *reinterpret_cast <const int *> (&val); // needs 32-bit int
  //   const int abs_mantissa = x & 0x007FFFFF;
  //   const int biased_exponent = x & 0x7F800000;
  //   if (biased_exponent == 0 && abs_mantissa != 0)
  //   {
  //     val = 0;
  //   }
  // }

  /** Returns true, if x is not-a-number, false otherwise. */
  RS_INLINE bool rsIsNaN(double x)
  {
    return !(x == x); 
    // comparison of a NaN to any number (including NaN itself) always returns false, so x == x 
    // will return false if (and only if) x is NaN
  }

  /** Assuming, that the FPU is in 'to nearest even integer' rounding mode (which is the default),
  this function rounds to the nearest integer using upward rounding when the argument is exactly
  halfway between two integers (instead of returning the nearest even integer in this case).
  Argument x must satify (INT_MIN/2)ñ1.0 < x < (INT_MAX/2)+1.0.  */
  RS_INLINE int rsRoundToInt(double x)
  {
#  if defined RS_COMPILER_MS_X86
    const float round_to_nearest = 0.5f;
    int i;
     __asm
     {
       fld x;
       fadd st, st (0);
       fadd round_to_nearest;
       fistp i;
       sar i, 1;
     }
     return (i);
    /*
#  elif defined __GNUC__  // too many memory references for 'add', too many memory references for 'sar'
    const float round_to_nearest = 0.5f;
    int i;
    asm("fld x");
    asm("fadd st, st(0)");
    asm("fadd round_to_nearest");
    asm("fistp i");
    asm("sar i, 1");
    return (i);
     */
#  else
    double xFloor = floor(x);
    double xFrac  = x-xFloor;
    if( xFrac >= 0.5 )
      return (int) xFloor + 1;
    else
      return (int) xFloor;
#  endif
  }

  /** Returns the largest integer that is smaller than or equal to the given number. */
  RS_INLINE int rsFloorInt(double x)
  {
#  if defined RS_COMPILER_MS_X86
    const float round_towards_m_i = -0.5f;
    int i;
    __asm
    {
      fld x;
      fadd st, st (0);
      fadd round_towards_m_i;
      fistp i;
      sar i, 1;
    }
    return (i);
    /*
#  elif defined __GNUC__  // under Linux, it says: memory input 2 is not directly addressable
    const float round_towards_m_i = -0.5f;
    int i;
    __asm__ __volatile__
    (
    "fldl     %1             \n\t"
    "fadd     %%st(0), %%st  \n\t"
    "fadds    %2             \n\t"
    "fistpl   %0             \n\t"
    "sarl     %0             \n\t"
    : "=m"(i)
    : "m"(x), "m"(round_towards_m_i)
    : "memory"
    );
    return (i);
#  elif defined __GNUC__  // too many memory references for 'add', too many memory references for 'sar'
    const float round_towards_m_i = -0.5f;
    int i;
    asm("fld x");
    asm("fadd st, st (0)");
    asm("fadd round_towards_m_i");
    asm("fistp i");
    asm("sar i, 1");
    return (i);
     */
#  else
     return (int) floor(x);
#  endif
  }

  /** Returns the smallest integer that is greater than or equal to the given number. */
  RS_INLINE int rsCeilInt(double x)
  {
    return -rsFloorInt(-x);
  }

  /** \todo comment this */
  RS_INLINE int rsTruncateToInt(double x)
  {
#  if defined RS_COMPILER_MS_X86
    const float round_towards_m_i = -0.5f;
    int i;
    __asm
    {
      fld x;
      fadd st, st (0);
      fabs;
      fadd round_towards_m_i;
      fistp i;
      sar i, 1;
    }
    if(x < 0)
      i = -i;
    return (i);
    /*
#  elif defined __GNUC__  // too many memory references for 'add', too many memory references for 'sar'
    const float round_towards_m_i = -0.5f;
    int i;
    asm("fld x");
    asm("fadd st, st (0)");
    asm("fabs");
    asm("fadd round_towards_m_i");
    asm("fistp i");
    asm("sar i, 1");
    if (x < 0)
      i = -i;
    return (i);
     */
#  else
    if( x >= 0.0 )
      return (int) floor(x);
    else
      return - (int) floor(fabs(x));
#  endif
  }

}

#endif
