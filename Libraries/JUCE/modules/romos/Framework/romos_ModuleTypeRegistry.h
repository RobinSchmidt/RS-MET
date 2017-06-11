#ifndef romos_ModuleTypeRegistry_h
#define romos_ModuleTypeRegistry_h

namespace romos
{

/** This class is used to map numeric module identifiers for module types to and from 
indetifiers represented as strings. The former are mainly used in the DSP code where fast access 
matters whereas the latter are mainly used for persistence purposes, where we want to have human 
readability. It is implemented as a singleton class.

\todo maybe use this also as ModuleFactory  */

class ModuleTypeRegistry
{

public:

  /** Defines a list of numeric identifiers to be used to quickly identify the type of a module 
  inside the DSP code. */
  enum moduleIdentifiers  // rename to moduleTypeCodes
  {
    // should all be nouns, optionally with appended specifiers (adjectives)

    UNKNOWN_MODULE_TYPE = 0,

    // Infrastructural:
    CONTAINER,
    TOP_LEVEL_MODULE,
    AUDIO_INPUT,
    AUDIO_OUTPUT,
    EVENT_INPUT,
    EVENT_OUTPUT,
    PARAMETER,
    SYSTEM_SAMPLE_RATE,

    // Events:
    NOTE_GATE,  // maybe have PolyGate/MonoGate ...but maybe just have the same M/P switch here as on all modules
    NOTE_ON_TRIGGER,
    NOTE_OFF_TRIGGER,
    VOICE_KILLER,
    VOICE_COMBINER,
    NOTE_FREQUENCY,
    NOTE_VELOCITY,

    // NOTE_ON_TRIGGER, NOTE_OFF_TRIGGER

    // Arithmetic:
    ADDER,
    SUBTRACTOR,
    MULTIPLIER,
    DIVIDER,
    //POWER        "^"
    CONSTANT,
    IDENTITY,    // rename to IDENTITY_AUDIO
    UNARY_MINUS,
    RECIPROCAL,       // y = 1/x
    //COMPLEMENT        // y = 1-x
    //SCALER,           // y = c*x
    //OFFSET            // y = x+c
    //LINEAR_POLYNOMIAL            // y = a*x + b 
    //QUADRATIC_POLYNOMIAL         // y = a*x^2 + b*x + c
    //CUBIC_POLYNOMIAL
    //QUARTIC_POLYNOMIAL
    //QUINTIC_POLYNOMIAL
    ADDER_3,
    ADDER_4,
    ADDER_5,     // maybe go upt to Adder8 - gives opportunities to optimize via parentheses (x1+x2) + (x3+x4) + 
    ADDER_N,     // dynamic Ins - maybe use the implicit summing at the single input - but for this we may also use the Identity module
    PRODUCT,     // dynamic Ins 
    MATRIX,      // dynamic Ins/Outs

    // Functions:
    CLIPPER,
    SIN_COS,
    FORMULA,     // dynamic Ins

    // CROSS_FADER, XY_FADER (mixes 4 inputs ...mmm maybe VECTOR_MIXER would be a better name)

    // Delays:
    UNIT_DELAY,
    DELAY_ROUNDING,
    DELAY_LINEAR,                // add allpass, cubic, etc...
    MULTI_TAP_DELAY_ROUNDING,    // dynamic Ins/Outs (Ins are tap delay-times, Outs are the tap-outputs)
    MULTI_TAP_DELAY_LINEAR,      // dynamic Ins/Outs 

    // Filters:      
    FIRST_ORDER_LOWPASS,
    FIRST_ORDER_FILTER,
    BIQUAD,
    // make a dynamic I/O version of all filters - they take the coefficients and a variable number of Ins/Outs for
    // multichannel processing
    BIQUAD_DESIGNER,
    LADDER_FILTER,



    // Generators:
    WHITE_NOISE,
    PERIODIC_LINEAR_RAMP,
    BANDLIMITED_IMPULSE_TRAIN,
    BLIT_SAW_OSCILLATOR,
    DUAL_BLIT_SAW_OSCILLATOR,
    // PHASE_MODULATION_2OPS, PHASE_MODULATION_3OPS, PHASE_MODULATION_4OPS, 

    // Modulators:
    ENVELOPE_ADSR,


    NUM_NON_TEST_MODULE_TYPES,

    // test-modules (for development, debugging, testing):
    TEST_GAIN,
    TEST_SUM_DIFF,
    TEST_WRAPPED_SUM_DIFF,
    TEST_SUMMED_DIFFS,
    TEST_MOVING_AVERAGE,
    TEST_LEAKY_INTEGRATOR,
    TEST_FILTER1,
    TEST_BIQUAD,
    TEST_ADDED_CONSTANTS,
    TEST_PIN_SORTING,
    TEST_BLIP,
    TEST_POLY_BLIP_STEREO,
    TEST_NOISE_FLUTE,
    //EXAMPLE_MOOG_FILTER,
    //TEST_CONTAINERIZE, 
    //TEST_UNCONTAINERIZE,
    //TEST_MINIMIZE_INS1,

    NUM_MODULE_TYPES
  };

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Given a module's numeric type identifier, this function returns the corresponding type-name 
  string. */
  rosic::rsString getModuleTypeStringFromIdentifier(int identifier);

  // \todo getModuleTypeDescription(int identifier) - maybe store these in a "OneWayAssociator" 
  // template class

  /** Given a module's type-name as string, this function returns the corresponding numeric 
  identifier. */
  int getModuleIdentifierFromTypeString(rosic::rsString typeString);

  /** Returns true when the passed identifier indicates an input- or output module 
  (audio or event). */
  static bool isIdentifierInputOrOutput(int typeIdentifier);

  /** Returns true when the module type given by typeIdentifier has an editable name, false 
  otherwise. Some basic modules do not allow their name to be edited such as adders, multipliers, 
  unit delays, etc. */
  static bool isModuleNameEditable(int typeIdentifier);

  // todo: make a similar function doesModuleBlockHaveHeader() - get rid of the flag in Module 
  // class

  /** Returns true, if the modules with given typeCode have a GUI editor - some very simple modules 
  such as adders, etc. don't. */
  static bool hasModuleTypeEditor(int typeCode);

  //-----------------------------------------------------------------------------------------------
  // access to the sole instance:

  /** Returns the sole instance of the ModuleTypeRegistry class (it's a singleton class). It will 
  also create and initialize the instance before, if necessary (i.e. if its hasn't been done 
  already). */
  static ModuleTypeRegistry* getSoleInstance();

  /** Deletes the sole instance of the ModuleTypeRegistry class, if existent. We can use this to 
  manually clean up heap memory which would otherwise be cleaned up automatically on shutdown. 
  Manual cleanup is preferable in order to not distract the memory leak detection of MSVC during 
  development/debugging. */
  static void deleteSoleInstance();

protected:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ModuleTypeRegistry();

  /** Destructor. */
  ~ModuleTypeRegistry();

  rosic::KeyValueMap<int, rosic::rsString> identifierNameMap;

private:

  static ModuleTypeRegistry* soleInstance;

};

/** Shorthand function for 
romos::ModuleTypeRegistry::getSoleInstance()->getModuleIdentifierFromTypeString(typeString); */
int getTypeId(rosic::rsString typeString);

/*
Module-Types to do:

for sure:

regular:

Difference:
y[n] = x[n] - x[n-1]

Clip:
if( in < min )
  out = min;
else if( in > max )
  out = max;
else
  out = in;

CenterClip:
if( in < min || in > max )
  out = in;
else
  out = centerValue;  // could be 0.5*(min+max) but maybe also set from outside

Select:
if( in1 >= threshold )
  out = in2;
else
  out = in3;

  // nah, make a multi-selector with variable number of ins and switch on round(in1)

Min, Max, Floor, Ceil, Round, Abs, Angle(x, y), Length(x, y), ToRA(x, y) (radius/angle), 
ToXY(r, a), Equal, InEqual, Greater, GreaterOrEqual, Less, LessOrEqual, And, Or, Not, Nand, Nor, 
Xor

associated with per-voice values (maybe have a PerVoiceValue Module):

associated with global values (maybe have a GlobalValue Module):
SampleRate, Tempo, Time

detection:
-zero-crossing finder (outputs a trigger signal and an offset to indicate subsample location of the 
 crossing)
-maximum-finder: x[n-1] > x[n] && x[n-1] > x[n-2] - output a trigger and an offset (found by 
 parabolic interpolation)
 -min-finder likewise
-Schmitt-Trigger (see CSound-book, page 568)

maybe:
SendAudio, ReceiveAudio, SendEvent, ReceiveEvent

UpSample, DownSample
-when connecting modules with different oversampling factors, automatic conversion is done by zero 
 insertion or discarding

*/

}

#endif
