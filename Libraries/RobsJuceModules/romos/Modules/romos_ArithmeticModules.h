#ifndef romos_ArithmeticModules_h
#define romos_ArithmeticModules_h

//namespace romos
//{

//-------------------------------------------------------------------------------------------------

class ConstantModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(ConstantModule);
public:
  virtual void clearVoiceBuffer(int voiceIndex);
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void setModuleName(const std::string& newName);
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
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class IdentityModule : public AtomicModule
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
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class UnaryMinusModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(UnaryMinusModule);
};
class UnaryMinusTypeInfo : public ModuleTypeInfo
{
public:
  UnaryMinusTypeInfo() {
    shortName    = "-";
    fullName     = "UnaryMinus";
    description  = "Negates the input (multiplies by -1)";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new UnaryMinusModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class ScalerModule : public ModuleWithParameters //  public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(ScalerModule);
public:
  virtual void parameterChanged(int index);
protected:
  double multiplier = 1;
};
class ScalerTypeInfo : public ModuleTypeInfo
{
public:
  ScalerTypeInfo() {
    shortName    = "*a";
    fullName     = "Scaler";
    description  = "Multiplies input by adjustable number";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new ScalerModule; };
    hasHeader    = false;
    hasEditor    = true;
  }
};
// todo: make offset/shift module and ScaleAndShift y = a*x+b

//-------------------------------------------------------------------------------------------------

class ReciprocalModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(ReciprocalModule);
};
class ReciprocalTypeInfo : public ModuleTypeInfo
{
public:
  ReciprocalTypeInfo() {
    shortName    = "1/x";
    fullName     = "Reciprocal";
    description  = "Computes the reciprocal 1/x of the input. Returns 0, if x=0. ";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new ReciprocalModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class AdderModule : public AtomicModule
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
    hasHeader    = false; // draw no header at top of the block
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class SubtractorModule : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class MultiplierModule : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};
// todo: make a multiplier with N inputs

//-------------------------------------------------------------------------------------------------

class DividerModule : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------
/*
// todo:
class PowerModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_2(PowerModule);
};
class PowerTypeInfo : public ModuleTypeInfo
{
public:
  PowerTypeInfo() {
    shortName    = "^";
    fullName     = "Power";
    description  = "Raises the first input to the power of the second";
    category     = "Arithmetic";
    createModule =  []()->Module* { return new PowerModule; };
    hasHeader = false;
  }
};
*/

//-------------------------------------------------------------------------------------------------

class Adder3Module : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class Adder4Module : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class Adder5Module : public AtomicModule
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
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class AdderNModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_N(AdderNModule);

  virtual void connectInputPinTo(int inputPinIndex, Module *sourceModule, int srcOutPinIndex);
  virtual void disconnectInputPin(int inputPinIndex);
   // factor these out into a class VariableNumInputsModule or something
};
class AdderNModuleTypeInfo : public ModuleTypeInfo
{
public:
  AdderNModuleTypeInfo() {
    shortName    = "+";
    fullName     = "AdderN";
    description  = "Sums an arbitrary number of input signals"; // up to WorkArea::maxNumPins
    category     = "Arithmetic";
    createModule =  []()->Module* { return new AdderNModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-----------------------------------------------------------------------------------------------

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
// Ins: Reset, Direction, Increment (lets step more than one steps per trigger), Step-Values; Outs: StepValueOut
// Switch: if in1 > in2 return in3 else in4
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
// LA-alike synth - needs processing-flags for that in the ContainerModule class and perhaps a corresponding LocalVoiceStarter module
// unless this starting (i.e. setting the flag back to true) is also done by NoteOn module

//COMPLEMENT        // y = 1-x
//SCALER,           // y = c*x
//OFFSET            // y = x+c
//MATRIX,      // dynamic Ins/Outs

//
//}

#endif 
