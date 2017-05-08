#include "romos_ModuleTypeRegistry.h"
using namespace romos;
 
ModuleTypeRegistry* ModuleTypeRegistry::soleInstance = NULL;

ModuleTypeRegistry* ModuleTypeRegistry::getSoleInstance()
{
  if( soleInstance == NULL )
    soleInstance = new ModuleTypeRegistry;
  return soleInstance;
}

void ModuleTypeRegistry::deleteSoleInstance()
{ 
  if( soleInstance != NULL )
  {
    delete soleInstance;
    soleInstance = NULL;
  }
}
 
rosic::String ModuleTypeRegistry::getModuleTypeStringFromIdentifier(int identifier)
{
  bool wasFound;
  rosic::String result = identifierNameMap.getValueForKey(identifier, wasFound);
  if( wasFound )
    return result;
  else
    return rosic::String("UnknownModuleType");
}

int ModuleTypeRegistry::getModuleIdentifierFromTypeString(rosic::String typeString)
{
  bool wasFound;
  int result = identifierNameMap.getKeyForValue(typeString, wasFound);
  if( wasFound )
    return result;
  else
    return UNKNOWN_MODULE_TYPE;
}

bool ModuleTypeRegistry::isIdentifierInputOrOutput(int typeIdentifier)
{
  return typeIdentifier == AUDIO_INPUT || typeIdentifier == AUDIO_OUTPUT || 
         typeIdentifier == EVENT_INPUT || typeIdentifier == EVENT_OUTPUT;
}

bool ModuleTypeRegistry::isModuleNameEditable(int typeIdentifier)
{
  return typeIdentifier != ADDER && typeIdentifier != SUBTRACTOR && typeIdentifier != MULTIPLIER && typeIdentifier != DIVIDER
    &&   typeIdentifier != IDENTITY && typeIdentifier != UNARY_MINUS && typeIdentifier != ADDER_N
    &&   typeIdentifier != PRODUCT && typeIdentifier != UNIT_DELAY;
    // include <,>,=,!=,Not,etc.
}

bool ModuleTypeRegistry::hasModuleTypeEditor(int typeCode)
{
  switch( typeCode )
  {
    // for each type that doesn't have an editor, include a case here - we don't need to do anything in the case statement because without 
    // breaks, each case falls through to the "return false" statement
  case CONSTANT: 
  case IDENTITY: 
  case UNARY_MINUS: 
  case RECIPROCAL: 
  case ADDER: 
  case SUBTRACTOR:
  case MULTIPLIER:
  case DIVIDER:
  case UNIT_DELAY:
  case ADDER_N: 
  case VOICE_COMBINER: 
  case PERIODIC_LINEAR_RAMP:
  case CLIPPER:
  case SIN_COS:
  case AUDIO_INPUT: 
  case AUDIO_OUTPUT: 
    return false;

  default: 
    return true;  // all others that are not listed in the case-marks should have an editor
  }
}

ModuleTypeRegistry::ModuleTypeRegistry()
{
  // Infrastructural:
  identifierNameMap.insertKeyValuePair(CONTAINER,                 rosic::String("Container"));
  identifierNameMap.insertKeyValuePair(TOP_LEVEL_MODULE,          rosic::String("TopLevelModule"));
  identifierNameMap.insertKeyValuePair(AUDIO_INPUT,               rosic::String("AudioInput"));
  identifierNameMap.insertKeyValuePair(AUDIO_OUTPUT,              rosic::String("AudioOutput"));
  identifierNameMap.insertKeyValuePair(EVENT_INPUT,               rosic::String("EventInput"));
  identifierNameMap.insertKeyValuePair(EVENT_OUTPUT,              rosic::String("EventOutput"));
  identifierNameMap.insertKeyValuePair(PARAMETER,                 rosic::String("Parameter"));
  identifierNameMap.insertKeyValuePair(SYSTEM_SAMPLE_RATE,        rosic::String("SystemSampleRate"));


  identifierNameMap.insertKeyValuePair(NOTE_GATE,                 rosic::String("NoteGate"));
  identifierNameMap.insertKeyValuePair(NOTE_ON_TRIGGER,           rosic::String("NoteOnTrigger"));
  identifierNameMap.insertKeyValuePair(NOTE_OFF_TRIGGER,          rosic::String("NoteOffTrigger"));
  identifierNameMap.insertKeyValuePair(VOICE_KILLER,              rosic::String("VoiceKiller"));
  identifierNameMap.insertKeyValuePair(VOICE_COMBINER,            rosic::String("VoiceCombiner"));
  identifierNameMap.insertKeyValuePair(NOTE_FREQUENCY,            rosic::String("NoteFrequency"));
  identifierNameMap.insertKeyValuePair(NOTE_VELOCITY,             rosic::String("NoteVelocity"));

  // Arithmetic:

  identifierNameMap.insertKeyValuePair(CONSTANT,                  rosic::String("Constant"));
  identifierNameMap.insertKeyValuePair(IDENTITY,                  rosic::String("Identity"));
  identifierNameMap.insertKeyValuePair(UNARY_MINUS,               rosic::String("UnaryMinus"));
  identifierNameMap.insertKeyValuePair(RECIPROCAL,                rosic::String("Reciprocal"));
  identifierNameMap.insertKeyValuePair(ADDER,                     rosic::String("Adder"));
  identifierNameMap.insertKeyValuePair(SUBTRACTOR,                rosic::String("Subtractor"));
  identifierNameMap.insertKeyValuePair(MULTIPLIER,                rosic::String("Multiplier"));
  identifierNameMap.insertKeyValuePair(DIVIDER,                   rosic::String("Divider"));
  identifierNameMap.insertKeyValuePair(ADDER_3,                   rosic::String("Adder3"));
  identifierNameMap.insertKeyValuePair(ADDER_4,                   rosic::String("Adder4"));
  identifierNameMap.insertKeyValuePair(ADDER_5,                   rosic::String("Adder5"));
  identifierNameMap.insertKeyValuePair(ADDER_N,                   rosic::String("AdderN"));
  identifierNameMap.insertKeyValuePair(PRODUCT,                   rosic::String("Product"));  // MultiplierN
  identifierNameMap.insertKeyValuePair(MATRIX,                    rosic::String("Matrix"));

  // Functions:
  identifierNameMap.insertKeyValuePair(CLIPPER,                   rosic::String("Clipper"));
  identifierNameMap.insertKeyValuePair(SIN_COS,                   rosic::String("SinCos"));
  identifierNameMap.insertKeyValuePair(FORMULA,                   rosic::String("Formula"));

  // Delays:
  identifierNameMap.insertKeyValuePair(UNIT_DELAY,                 rosic::String("UnitDelay"));
  identifierNameMap.insertKeyValuePair(DELAY_ROUNDING,             rosic::String("DelayRounding"));
  identifierNameMap.insertKeyValuePair(DELAY_LINEAR,               rosic::String("DelayLinear"));  
  identifierNameMap.insertKeyValuePair(MULTI_TAP_DELAY_ROUNDING,   rosic::String("MultiTapDelayRounding"));
  identifierNameMap.insertKeyValuePair(MULTI_TAP_DELAY_LINEAR,     rosic::String("MultiTapDelayLinear"));

  // Filters:
  identifierNameMap.insertKeyValuePair(FIRST_ORDER_LOWPASS,        rosic::String("FirstOrderLowpass"));
  identifierNameMap.insertKeyValuePair(FIRST_ORDER_FILTER,         rosic::String("FirstOrderFilter"));
  identifierNameMap.insertKeyValuePair(BIQUAD,                     rosic::String("Biquad"));
  identifierNameMap.insertKeyValuePair(BIQUAD_DESIGNER,            rosic::String("BiquadDesigner"));
  identifierNameMap.insertKeyValuePair(LADDER_FILTER,              rosic::String("LadderFilter"));


  // Generators:
  identifierNameMap.insertKeyValuePair(PERIODIC_LINEAR_RAMP,      rosic::String("PeriodicLinearRamp"));
  identifierNameMap.insertKeyValuePair(WHITE_NOISE,               rosic::String("WhiteNoise"));
  identifierNameMap.insertKeyValuePair(BANDLIMITED_IMPULSE_TRAIN, rosic::String("BandlimitedImpulseTrain"));
  identifierNameMap.insertKeyValuePair(BLIT_SAW_OSCILLATOR,       rosic::String("BlitSaw"));
  identifierNameMap.insertKeyValuePair(DUAL_BLIT_SAW_OSCILLATOR,  rosic::String("DualBlitSaw"));


  // Modulators:
  identifierNameMap.insertKeyValuePair(ENVELOPE_ADSR,             rosic::String("EnvelopeADSR"));

  // Testmodules:
  identifierNameMap.insertKeyValuePair(TEST_GAIN,             rosic::String("Gain"));
  identifierNameMap.insertKeyValuePair(TEST_SUM_DIFF,         rosic::String("SumDiff"));
  identifierNameMap.insertKeyValuePair(TEST_WRAPPED_SUM_DIFF, rosic::String("WrappedSumDiff"));
  identifierNameMap.insertKeyValuePair(TEST_SUMMED_DIFFS,     rosic::String("SummedDiffs"));
  identifierNameMap.insertKeyValuePair(TEST_MOVING_AVERAGE,   rosic::String("MovingAverage"));
  identifierNameMap.insertKeyValuePair(TEST_LEAKY_INTEGRATOR, rosic::String("LeakyIntegrator"));
  identifierNameMap.insertKeyValuePair(TEST_FILTER1,          rosic::String("TestFilter1"));
  identifierNameMap.insertKeyValuePair(TEST_BIQUAD,           rosic::String("BiquadMacro"));
  identifierNameMap.insertKeyValuePair(TEST_ADDED_CONSTANTS,  rosic::String("AddedConstants"));
  identifierNameMap.insertKeyValuePair(TEST_PIN_SORTING,      rosic::String("PinSorting"));
  identifierNameMap.insertKeyValuePair(TEST_BLIP,             rosic::String("TestBlip"));
  identifierNameMap.insertKeyValuePair(TEST_POLY_BLIP_STEREO, rosic::String("PolyBlipStereo"));
  identifierNameMap.insertKeyValuePair(TEST_NOISE_FLUTE,      rosic::String("NoiseFlute"));


  //identifierNameMap.insertKeyValuePair(EXAMPLE_MOOG_FILTER,       rosic::String("ExampleMoogFilter"));
  //identifierNameMap.insertKeyValuePair(TEST_CONTAINERIZE,         rosic::String("TestContainerize"));
  //identifierNameMap.insertKeyValuePair(TEST_UNCONTAINERIZE,       rosic::String("TestUncontainerize"));
  //identifierNameMap.insertKeyValuePair(TEST_MINIMIZE_INS1,        rosic::String("TestMinimizeIns1"));




  //identifierNameMap.insertKeyValuePair(FORMULA_ARRAY,             rosic::String("FormulaArray")); // wassat?
  //identifierNameMap.insertKeyValuePair(MULTI_IN_FORMULA,          rosic::String("MultiInFormula"));  
  //identifierNameMap.insertKeyValuePair(UNARY_FORMULA,             rosic::String("UnaryFormula"));
}

ModuleTypeRegistry::~ModuleTypeRegistry()
{

}

int romos::getTypeId(rosic::String typeString)
{
  return romos::ModuleTypeRegistry::getSoleInstance()->getModuleIdentifierFromTypeString(typeString);
}