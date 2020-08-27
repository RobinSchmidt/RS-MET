#ifndef romos_FunctionModules_h
#define romos_FunctionModules_h


//-------------------------------------------------------------------------------------------------

/** Clips the input signal to a range between some specified minimum and maximum. */
class ClipperModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_3(ClipperModule);
};
class ClipperTypeInfo : public ModuleTypeInfo
{
public:
  ClipperTypeInfo() {
    shortName    = "Clip";
    fullName     = "Clipper";
    description  = "Clips input signal to the range between Min and Max";
    category     = "Functions";
    createModule =  []()->Module* { return new ClipperModule; };
    //hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

/** Computes the sine and cosine of 2*pi times the input value (multiplication of the argument by
2*pi results in a sine and cosing wave with a period of unity). */
class SinCosModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(SinCosModule);
};
class SinCosTypeInfo : public ModuleTypeInfo
{
public:
  SinCosTypeInfo() {
    shortName    = "SinCos";
    fullName     = "SinCos"; // maybe rename to SineAndCosine
    description  = "Sine and cosine of 2*pi times the input (1-periodic)";
    category     = "Functions";
    createModule =  []()->Module* { return new SinCosModule; };
    //hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

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
  3: Attack Sigmoidity
  4: Decay Bend
  5: Decay Sigmoidity
Outputs:
  0: the trisaw signal  */
class TriSawModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_6(TriSawModule);
};
class TriSawTypeInfo : public ModuleTypeInfo
{
public:
  TriSawTypeInfo() {
    shortName    = "TriSaw";
    fullName     = "TriSaw";
    description  = "Morph between triangle, saw, square and sine";
    category     = "Functions";
    createModule =  []()->Module* { return new TriSawModule; };
    //hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class SaturatorModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_3(SaturatorModule);
  double (*sigmoid)(double) = tanh;  // (normalized) prototype sigmoid function
};
class SaturatorTypeInfo : public ModuleTypeInfo
{
public:
  SaturatorTypeInfo() {
    shortName    = "Saturator";
    fullName     = "Saturator";
    description  = "Saturation with adjustable width and center";
    category     = "Functions";
    createModule =  []()->Module* { return new SaturatorModule; };
    hasHeader    = true;
  }
};






// round, ceil, floor, abs, max, min
// rawToDecibels, DecibelsToRaw, pitchToFrequency (have also version with custom mapping 
// (microtuning)), LinearToExponential, etc.

// Function (evaluates expression entered on the gui), FunctionNto1, Function1ToN, FunctionNToM
// or maybe the N-to-M handles all cases, the GUI should have 4 entry fields:
// Inputs:   x; a; b
// Defaults: x=0; a=1; b=0;
// Outputs:  y; z
// Formulas: y = a*x + b; z = a*x - b
// especially for the Formulas field, it would be nice to be multi-line
// the class should automatically assign the names for the pins according to the Inputs/Outputs


#endif 
