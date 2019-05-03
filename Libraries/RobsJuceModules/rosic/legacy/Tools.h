// maybe the command-line option /FORCE:MULTIPLE has to be added
// to the linker options

#ifndef Tools_h
#define Tools_h

#include <math.h>
#include <float.h>
#include "Definitions.h"

namespace Tools
{
 INLINE int roundToInt(double value);
 // rounds a double-value to the nearest integer.


 //---------------------------------------------------------------------------
 //implementation:

 // code is taken form Agner Fog - Optimizing software in C++
 INLINE int Tools::roundToInt(double const x)
 {
  int n;

  #if defined(__unix__) || defined(__GNUC__)
   // 32-bit Linux, Gnu/AT&T syntax:
   __asm ("fldl %1 \n fistpl %0 " : "=m"(n) : "m"(x) : "memory" );
  #else
   // 32-bit Windows, Intel/MASM syntax:
   __asm fld qword ptr x;
   __asm fistp dword ptr n;
  #endif

  return n;
 }

}

#endif // #ifndef Tools_h