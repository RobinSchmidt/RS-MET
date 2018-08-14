#ifndef romos_ArithmeticModules_h
#define romos_ArithmeticModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "../Algorithms/romos_FilterDesign.h"
#include "../Framework/romos_WorkArea.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{

//-------------------------------------------------------------------------------------------------
// arithmetic modules:

/** Outputs a constant number. */
class ConstantModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_1(ConstantModule);
public:
  virtual void clearVoiceBuffer(int voiceIndex);
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void setModuleName(const rosic::rsString& newName);
  virtual unsigned int getNumOutputPins() const { return 1; }
  virtual double getValue() const { return value; }
protected:
  double value;
};
class ConstantModuleTypeInfo : public ModuleTypeInfo
{
public:
  ConstantModuleTypeInfo() {
    shortName    = "Const";
    fullName     = "Constant";
    description  = "Outputs a constant value";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new ConstantModule; };
  }
};

/** Copies the signal in the input slot into the output slot. */
class IdentityModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_1(IdentityModule);
};
class IdentityModuleTypeInfo : public ModuleTypeInfo
{
public:
  IdentityModuleTypeInfo() {
    shortName    = "id";
    fullName     = "Identity";
    description  = "Passes input to output as is";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new IdentityModule; };
  }
};

/** Negates the input signal (multiplies by -1). */
class UnaryMinusModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_1(UnaryMinusModule);
};
class UnaryMinusTypeInfo : public ModuleTypeInfo
{
public:
  UnaryMinusTypeInfo() {
    shortName    = "-";
    fullName     = "Unary minus";
    description  = "Negates the input";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new UnaryMinusModule; };
    hasHeader = false;
  }
};

/** Computes 1/x - when x is zero, it returns zero as well. */
class ReciprocalModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_1(ReciprocalModule);
};
class ReciprocalTypeInfo : public ModuleTypeInfo
{
public:
  ReciprocalTypeInfo() {
    shortName    = "1/x";
    fullName     = "Reciprocal";
    description  = "Computes the reciprocal of the input, with zero checking";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new ReciprocalModule; };
    hasHeader = false;
  }
};

/** Adds two input signals. */
class AdderModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_2(AdderModule);  // 2 because it has two inputs
};
class AdderModuleTypeInfo : public ModuleTypeInfo
{
public:
  AdderModuleTypeInfo() {
    shortName    = "+";                          // shown inside the block 
    fullName     = "Adder";                      // shown in the menu (required!)
    description  = "Adds the two input signals"; // used for the description field
    category     = "Arithmetic";                 // used for categorization in the menu
    createModule =  []()->Module* { return new AdderModule; };
    // addInputPinInfo( "1st", "First argument", "First argument a of c = a + b");
    // addInputPinInfo( "2nd", "Second argument", "Second argument b of c = a + b");
    // addOutputPinInfo("Res", "Result", "Result c of c = a + b");
    hasHeader = false; // draw no header at top of the block
  }
};

/** Subtracts the second input signal from the first. */
class SubtractorModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_2(SubtractorModule);
};
class SubtractorTypeInfo : public ModuleTypeInfo
{
public:
    SubtractorTypeInfo() {
    shortName    = "-";
    fullName     = "Subtractor";
    description  = "Subtracts the second input from the first";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new SubtractorModule; };
    hasHeader = false;
  }
};

/** Multiplies two input signals. */
class MultiplierModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_2(MultiplierModule);
};
class MultiplierTypeInfo : public ModuleTypeInfo
{
public:
  MultiplierTypeInfo() {
    shortName    = "*";
    fullName     = "Multiplier";
    description  = "Multiplies the two inputs";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new MultiplierModule; };
    hasHeader = false;
  }
};
// make a multiplier with N inputs


/** Divides the first input signal by the second - if the second input signal is zero, the output 
will also be zero. */
class DividerModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_2(DividerModule);
};
class DividerTypeInfo : public ModuleTypeInfo
{
public:
  DividerTypeInfo() {
    shortName    = "/";
    fullName     = "Divider";
    description  = "Divides the first input by the second (if the 2nd is 0, it will return 0)";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new DividerModule; };
    hasHeader = false;
  }
};

/** Adds 3 input signals. */
class Adder3Module : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_3(Adder3Module);
};
class Adder3ModuleTypeInfo : public ModuleTypeInfo
{
public:
  Adder3ModuleTypeInfo() {
    shortName    = "+";
    fullName     = "Adder3";
    description  = "Adds three input signals together";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new Adder3Module; };
    hasHeader = false;
  }
};

/** Adds 4 input signals. */
class Adder4Module : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_4(Adder4Module);
};
class Adder4ModuleTypeInfo : public ModuleTypeInfo
{
public:
  Adder4ModuleTypeInfo() {
    shortName    = "+";
    fullName     = "Adder4";
    description  = "Adds four input signals together";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new Adder4Module; };
    hasHeader = false;
  }
};


/** Adds 5 input signals. */
class Adder5Module : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_5(Adder5Module);
};
class Adder5ModuleTypeInfo : public ModuleTypeInfo
{
public:
  Adder5ModuleTypeInfo() {
    shortName    = "+";
    fullName     = "Adder5";
    description  = "Adds five input signals together";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new Adder5Module; };
    hasHeader = false;
  }
};


/** Sums an arbitrary number (up to WorkArea::maxNumPins) of input signals. */
class AdderNModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_N(AdderNModule);

  virtual void connectInputPinTo(int inputPinIndex, Module *sourceModule, int sourceOutputPinIndex);
  virtual void disconnectInputPin(int inputPinIndex);
   // factor these out into a class VariableNumInputsModule or something
};
class AdderNModuleTypeInfo : public ModuleTypeInfo
{
public:
  AdderNModuleTypeInfo() {
    shortName    = "+";
    fullName     = "AdderN";
    description  = "Sums an arbitrary number of input signals";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new AdderNModule; };
    hasHeader = false;
  }
};

//-----------------------------------------------------------------------------------------------
// power ^, scaler (by fixed number) - 
// ScaleAndShift: *3 scales by 3, +2 shifts by 2, *3+2 scales by 3 and the adds 2
// default: *1+0


//-----------------------------------------------------------------------------------------------
// relations:

// equal, unequal, greater, less, greaterOrEqual, lessOrEqual, isZero, approximatelyEqual, within


//-----------------------------------------------------------------------------------------------
// logic:

// logic operarions: and, or, not, xor, nand, ....







//---------------------------------------------------------------------------------------------------------------------------------------
// misc modules:

// \todo:
// SampleAndHold: Ins: Trigger; Outs: OutValue
// Clock - outputs impulses at regular intervals: Ins: TickInterval (seconds), Start, Stop; Outs: Tick
// Switch - switch between a number of values on each received trigger in a round robin manner; 
//          Ins: Reset, Direction, Increment (lets step more than one steps per trigger), Step-Values; Outs: StepValueOut
// StopWatch - outputs the elapsed time since the last trigger; Ins: Start, Stop; Outs: ElapsedTime
// Timer: outputs a trigger impulse a certain time after it recieved a trigger impulse; Ins: Start, Stop, Time; Outs: DelayedTrigger
// Ramp: generates a linar ramp between two values at a given period/cycle-length Ins: Start, End, Period, Phase, Reset
//   ->generally, let time-domain stuff define rates in terms of cycle-length and audio-stuff in terms of frequency
// Pan, Balance, Crossfade, BinaryNoise (switch randomly between 0..1)
// generalized mod: mod(a, b) = a - floor(a/b) * b
// ExponentialSegments - see Elements of Computer Music (EoCM), page 183 (similar to gen4, variable number of segments, maybe with loop)
// LineSegments - similar to the above, but without shape control
// ADSR as special case of the above
// decibelsToAmpModDepth: see EoCM, page 188
// SineBank: bunch of sine oscillators, each with frequency and start-phase, global phase reset input (variable number of freq,phase pairs)
// weighted sum: y[n] = sum_k ( wk[n] * xk[n] ) variable number of  (wk, xk) pairs
// SineBankViaFormula: inputs: reference frequency, reset, numSines - the partial relative frequencies, amplitudes, start phases are 
// specified via a formula
// ModalFilterBankViaFormula: similar to SineBankViaFormula but with additional decay parameter for each partial and an global 
// decayScale input (maybe also a stereo version with spread values for each of the above parameters)
// LocalVoiceKiller - stops processing of the voice only for the immediate parent container - useful to stop transient-synthesis in an
// LA-alike synth - needs processing-flags for that in the ModuleContainer class and perhaps a corresponding LocalVoiceStarter module
// unless this starting (i.e. setting the flag back to true) is also done by NoteOn module


}

#endif 
