#ifndef romos_FunctionModules_h
#define romos_FunctionModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{

  /** Clips the input signal to a range between some specified minimum and maximum. */
  class ClipperModule : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_3(ClipperModule);
  };

  /** Computes the sine and cosine of 2*pi times the input value (multiplication of the argument by 
  2*pi results in a sine and cosing wave with a period of unity). */
  class SinCosModule : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_1(SinCosModule);
  };

  // round, ceil, floor, abs, max, min
  // rawToDecibels, DecibelsToRaw, pitchToFrequency (have also version with custom mapping 
  // (microtuning)), LinearToExponential, etc.

  // Function (evaluates expression entered on the gui), FunctionNto1, Function1ToN, FunctionNToM

} 

#endif 
