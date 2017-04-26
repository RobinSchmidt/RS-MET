#ifndef rosic_NumberManipulations_h
#define rosic_NumberManipulations_h

// rosic includes:
#include "GlobalDefinitions.h"
#include <math.h>

namespace rosic
{

  /** \todo: check the linux versions....they have been inserted here sloppily in order to make it compile when
  porting. */

/** Converts a double precision floating point number to an integer using the FPU's current
rounding mode (which is 'to nearest even integer' by default).  */
INLINE int toInt(double x)
{
  int a;

#if MSC_X86_ASM
  __asm
  {
    fld x;
    fistp a;
  }
#else
  a = (int)x;
#endif

  return (a);
}

/** Assuming, that the FPU is in 'to nearest even integer' rounding mode (which is the default),
this function rounds to the nearest integer using upward rounding when the argument is exactly
halfway between two integers (instead of returning the nearest even integer in this case).
Argument x must satify (INT_MIN/2)ñ1.0 < x < (INT_MAX/2)+1.0.  */
INLINE int roundToInt(double x)
{
  const float round_to_nearest = 0.5f;
  int i;

#if MSC_X86_ASM
  __asm
  {
    fld x;
    fadd st, st (0);
    fadd round_to_nearest;
    fistp i;
    sar i, 1;
  }
#else
// preliminary:
  i = (int)floor(x);
  int r = (int) (x-i);
  if(r >= 0.5)
    i++;
#endif

  return (i);
}

INLINE int floorInt(double x)
{
  const float round_towards_m_i = -0.5f;
  int i;

#if MSC_X86_ASM
  __asm
  {
    fld x;
    fadd st, st (0);
    fadd round_towards_m_i;
    fistp i;
    sar i, 1;
  }
#else
  i = (int)floor(x);
#endif

  return (i);
}

INLINE int ceilInt(double x)
{
  const float round_towards_p_i = -0.5f;
  int i;

#if MSC_X86_ASM
  __asm
  {
    fld x;
    fadd st, st (0);
    fsubr round_towards_p_i;
    fistp i;
    sar i, 1;
  }
#else
  return (int)ceil(x);
#endif

  return (-i);
}

INLINE int truncateToInt(double x)
{
  const float round_towards_m_i = -0.5f;
  int i;

#if MSC_X86_ASM
  __asm
  {
    fld x;
    fadd st, st (0);
    fabs;
    fadd round_towards_m_i;
    fistp i;
    sar i, 1;
  }
#else
  i = 0;
  DEBUG_BREAK; // this isn't implemented on linux as of yet
#endif

  if(x < 0)
  {
    i = -i;
  }
  return (i);
}

INLINE int extractExponent(double x)
{
  //return (int) (((*((reinterpret_cast<unsigned long long *>(&x)))&0x7FFFFFFFFFFFFFFF)>>52)-1023);
  //int dummy = (int) (((*((reinterpret_cast<unsigned long long *>(&x)))&0x7FFFFFFFFFFFFFFFULL)>>52)-1023);
  return (int)(((*((reinterpret_cast<UINT64 *>(&x)))&0x7FFFFFFFFFFFFFFFULL)>>52)-1023);
  //EXPOFDBL(value) (((*((reinterpret_cast<unsigned __int64 *>(&value)))&0x7FFFFFFFFFFFFFFF)>>52)-1023)
}

INLINE int extractExponent(float x)
{
  return (((*((reinterpret_cast<UINT32 *>(&x)))&0x7FFFFFFF)>>23)-127);
  //EXPOFFLT(value) (((*((reinterpret_cast<unsigned __int32 *>(&value)))&0x7FFFFFFF)>>23)-127)
}

INLINE double rAbs(double x)
{
  unsigned long long intAbsValue = *((unsigned long long*) &x) & 0x7FFFFFFFFFFFFFFFULL;
  return *((double*)&intAbsValue);
}


} // end namespace rosic

#endif // #ifndef rosic_NumberManipulations_h
