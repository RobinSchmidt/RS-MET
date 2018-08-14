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

/* TriSaw shaper module, takes a Phasor input (range -1..+1) and converts it into a TriSaw 
oscillator signal. Depending on the parameters, various waveshapes are possible, among them
triangle, saw-up, saw-down, square/pulse and (approximated) sine and everything in between.
With the asymmetry parameter, for example, you can morph from saw-down through triangle to 
saw-up
Parameters:
  none
Inputs:
  0: Input signal (supposed to be a sawtooth in -1..+1 for a proper TriSaw output)
  1: Asymmetry
  2: Attack Bend
  3: Attack Sigmodity
  4: Decay Bend
  5: Decay Sigmodity
Outputs:
  0: the trisaw signal  */
class TriSawModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_6(TriSawModule);
};


// round, ceil, floor, abs, max, min
// rawToDecibels, DecibelsToRaw, pitchToFrequency (have also version with custom mapping 
// (microtuning)), LinearToExponential, etc.

// Function (evaluates expression entered on the gui), FunctionNto1, Function1ToN, FunctionNToM

}

#endif 
