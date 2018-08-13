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
 
rosic::rsString ModuleTypeRegistry::getModuleTypeStringFromIdentifier(int identifier)
{
  bool wasFound;
  rosic::rsString result = identifierNameMap.getValueForKey(identifier, wasFound);
  if( wasFound )
    return result;
  else
    return rosic::rsString("UnknownModuleType");
}

int ModuleTypeRegistry::getModuleIdentifierFromTypeString(rosic::rsString typeString)
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
  case PHASOR:
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
  identifierNameMap.insertKeyValuePair(CONTAINER,                 rosic::rsString("Container"));
  identifierNameMap.insertKeyValuePair(TOP_LEVEL_MODULE,          rosic::rsString("TopLevelModule"));
  identifierNameMap.insertKeyValuePair(AUDIO_INPUT,               rosic::rsString("AudioInput"));
  identifierNameMap.insertKeyValuePair(AUDIO_OUTPUT,              rosic::rsString("AudioOutput"));
  identifierNameMap.insertKeyValuePair(EVENT_INPUT,               rosic::rsString("EventInput"));
  identifierNameMap.insertKeyValuePair(EVENT_OUTPUT,              rosic::rsString("EventOutput"));
  identifierNameMap.insertKeyValuePair(PARAMETER,                 rosic::rsString("Parameter"));
  identifierNameMap.insertKeyValuePair(SYSTEM_SAMPLE_RATE,        rosic::rsString("SystemSampleRate"));


  identifierNameMap.insertKeyValuePair(NOTE_GATE,                 rosic::rsString("NoteGate"));
  identifierNameMap.insertKeyValuePair(NOTE_ON_TRIGGER,           rosic::rsString("NoteOnTrigger"));
  identifierNameMap.insertKeyValuePair(NOTE_OFF_TRIGGER,          rosic::rsString("NoteOffTrigger"));
  identifierNameMap.insertKeyValuePair(VOICE_KILLER,              rosic::rsString("VoiceKiller"));
  identifierNameMap.insertKeyValuePair(VOICE_COMBINER,            rosic::rsString("VoiceCombiner"));
  identifierNameMap.insertKeyValuePair(NOTE_FREQUENCY,            rosic::rsString("NoteFrequency"));
  identifierNameMap.insertKeyValuePair(NOTE_VELOCITY,             rosic::rsString("NoteVelocity"));

  // Arithmetic:

  identifierNameMap.insertKeyValuePair(CONSTANT,                  rosic::rsString("Constant"));
  identifierNameMap.insertKeyValuePair(IDENTITY,                  rosic::rsString("Identity"));
  identifierNameMap.insertKeyValuePair(UNARY_MINUS,               rosic::rsString("UnaryMinus"));
  identifierNameMap.insertKeyValuePair(RECIPROCAL,                rosic::rsString("Reciprocal"));
  identifierNameMap.insertKeyValuePair(ADDER,                     rosic::rsString("Adder"));
  identifierNameMap.insertKeyValuePair(SUBTRACTOR,                rosic::rsString("Subtractor"));
  identifierNameMap.insertKeyValuePair(MULTIPLIER,                rosic::rsString("Multiplier"));
  identifierNameMap.insertKeyValuePair(DIVIDER,                   rosic::rsString("Divider"));
  identifierNameMap.insertKeyValuePair(ADDER_3,                   rosic::rsString("Adder3"));
  identifierNameMap.insertKeyValuePair(ADDER_4,                   rosic::rsString("Adder4"));
  identifierNameMap.insertKeyValuePair(ADDER_5,                   rosic::rsString("Adder5"));
  identifierNameMap.insertKeyValuePair(ADDER_N,                   rosic::rsString("AdderN"));
  identifierNameMap.insertKeyValuePair(PRODUCT,                   rosic::rsString("Product"));  // MultiplierN
  identifierNameMap.insertKeyValuePair(MATRIX,                    rosic::rsString("Matrix"));

  // Functions:
  identifierNameMap.insertKeyValuePair(CLIPPER,                   rosic::rsString("Clipper"));
  identifierNameMap.insertKeyValuePair(SIN_COS,                   rosic::rsString("SinCos"));
  identifierNameMap.insertKeyValuePair(FORMULA,                   rosic::rsString("Formula"));

  // Delays:
  identifierNameMap.insertKeyValuePair(UNIT_DELAY,                 rosic::rsString("UnitDelay"));
  identifierNameMap.insertKeyValuePair(DELAY_ROUNDING,             rosic::rsString("DelayRounding"));
  identifierNameMap.insertKeyValuePair(DELAY_LINEAR,               rosic::rsString("DelayLinear"));  
  identifierNameMap.insertKeyValuePair(MULTI_TAP_DELAY_ROUNDING,   rosic::rsString("MultiTapDelayRounding"));
  identifierNameMap.insertKeyValuePair(MULTI_TAP_DELAY_LINEAR,     rosic::rsString("MultiTapDelayLinear"));

  // Filters:
  identifierNameMap.insertKeyValuePair(FIRST_ORDER_LOWPASS,        rosic::rsString("FirstOrderLowpass"));
  identifierNameMap.insertKeyValuePair(FIRST_ORDER_FILTER,         rosic::rsString("FirstOrderFilter"));
  identifierNameMap.insertKeyValuePair(BIQUAD,                     rosic::rsString("Biquad"));
  identifierNameMap.insertKeyValuePair(BIQUAD_DESIGNER,            rosic::rsString("BiquadDesigner"));
  identifierNameMap.insertKeyValuePair(LADDER_FILTER,              rosic::rsString("LadderFilter"));


  // Generators:
  identifierNameMap.insertKeyValuePair(PHASOR,                    rosic::rsString("Phasor"));
  identifierNameMap.insertKeyValuePair(WHITE_NOISE,               rosic::rsString("WhiteNoise"));
  identifierNameMap.insertKeyValuePair(BANDLIMITED_IMPULSE_TRAIN, rosic::rsString("BandlimitedImpulseTrain"));
  identifierNameMap.insertKeyValuePair(BLIT_SAW_OSCILLATOR,       rosic::rsString("BlitSaw"));
  identifierNameMap.insertKeyValuePair(DUAL_BLIT_SAW_OSCILLATOR,  rosic::rsString("DualBlitSaw"));


  // Modulators:
  identifierNameMap.insertKeyValuePair(ENVELOPE_ADSR,             rosic::rsString("EnvelopeADSR"));

  // Testmodules:
  identifierNameMap.insertKeyValuePair(TEST_GAIN,             rosic::rsString("Gain"));
  identifierNameMap.insertKeyValuePair(TEST_SUM_DIFF,         rosic::rsString("SumDiff"));
  identifierNameMap.insertKeyValuePair(TEST_WRAPPED_SUM_DIFF, rosic::rsString("WrappedSumDiff"));
  identifierNameMap.insertKeyValuePair(TEST_SUMMED_DIFFS,     rosic::rsString("SummedDiffs"));
  identifierNameMap.insertKeyValuePair(TEST_MOVING_AVERAGE,   rosic::rsString("MovingAverage"));
  identifierNameMap.insertKeyValuePair(TEST_LEAKY_INTEGRATOR, rosic::rsString("LeakyIntegrator"));
  identifierNameMap.insertKeyValuePair(TEST_FILTER1,          rosic::rsString("TestFilter1"));
  identifierNameMap.insertKeyValuePair(TEST_BIQUAD,           rosic::rsString("BiquadMacro"));
  identifierNameMap.insertKeyValuePair(TEST_ADDED_CONSTANTS,  rosic::rsString("AddedConstants"));
  identifierNameMap.insertKeyValuePair(TEST_PIN_SORTING,      rosic::rsString("PinSorting"));
  identifierNameMap.insertKeyValuePair(TEST_BLIP,             rosic::rsString("TestBlip"));
  identifierNameMap.insertKeyValuePair(TEST_POLY_BLIP_STEREO, rosic::rsString("PolyBlipStereo"));
  identifierNameMap.insertKeyValuePair(TEST_NOISE_FLUTE,      rosic::rsString("NoiseFlute"));


  //identifierNameMap.insertKeyValuePair(EXAMPLE_MOOG_FILTER,       rosic::rsString("ExampleMoogFilter"));
  //identifierNameMap.insertKeyValuePair(TEST_CONTAINERIZE,         rosic::rsString("TestContainerize"));
  //identifierNameMap.insertKeyValuePair(TEST_UNCONTAINERIZE,       rosic::rsString("TestUncontainerize"));
  //identifierNameMap.insertKeyValuePair(TEST_MINIMIZE_INS1,        rosic::rsString("TestMinimizeIns1"));




  //identifierNameMap.insertKeyValuePair(FORMULA_ARRAY,             rosic::rsString("FormulaArray")); // wassat?
  //identifierNameMap.insertKeyValuePair(MULTI_IN_FORMULA,          rosic::rsString("MultiInFormula"));  
  //identifierNameMap.insertKeyValuePair(UNARY_FORMULA,             rosic::rsString("UnaryFormula"));
}

ModuleTypeRegistry::~ModuleTypeRegistry()
{

}

int romos::getTypeId(rosic::rsString typeString)
{
  return romos::ModuleTypeRegistry::getSoleInstance()->getModuleIdentifierFromTypeString(typeString);
}