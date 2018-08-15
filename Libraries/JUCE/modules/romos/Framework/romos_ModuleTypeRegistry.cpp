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
  identifierNameMap.insertKeyValuePair(TRISAW,                    rosic::rsString("TriSaw"));
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

//=================================================================================================

void ModuleTypeInfo::addInputPinInfo(const char* shortName, const char* fullName,
  const char* description)
{
  inputShortNames.push_back(shortName);
  inputFullNames.push_back(fullName);
  inputDescriptions.push_back(description);
}

void ModuleTypeInfo::addOutputPinInfo(const char* shortName, const char* fullName,
  const char* description)
{
  outputShortNames.push_back(shortName);
  outputFullNames.push_back(fullName);
  outputDescriptions.push_back(description);
}

//-------------------------------------------------------------------------------------------------

ModuleFactoryNew romos::moduleFactory;  // definition of the global object - causes memleak?
MemLeakTest romos::memLeakTest;         // clean up - move to a memleak demo project

ModuleFactoryNew::ModuleFactoryNew()
{
  registerStandardModules();
}

ModuleFactoryNew::~ModuleFactoryNew()
{
  clearRegisteredTypes();
}

romos::Module* ModuleFactoryNew::createModule(int id, const std::string& name, int x, int y, 
  bool polyphonic) const
{
  rassert(id >= 0 && id < typeInfos.size());  // id out of range
  // todo: if the id is out of range, return some kind of "Error" dummy module


  romos::Module* m = typeInfos[id]->createModule();
  m->typeInfo = typeInfos[id];
  setupModule(m, name, x, y, polyphonic);
  return m;
}

romos::Module* ModuleFactoryNew::createModule(const std::string& fullTypeName, 
  const std::string& name, int x, int y, bool polyphonic) const
{
  return createModule(getModuleId(fullTypeName), name, x, y, polyphonic);
}

romos::TopLevelModule* ModuleFactoryNew::createTopLevelModule(const std::string& name, 
  int x, int y, bool polyphonic) const
{
  TopLevelModule* tlm = new TopLevelModule();
  setupModule(tlm, name, x, y, polyphonic);
  return tlm;
}

void ModuleFactoryNew::deleteModule(romos::Module* moduleToDelete)
{
  moduleToDelete->cleanUp();
  delete moduleToDelete;
}

int ModuleFactoryNew::getModuleId(const std::string& fullTypeName) const
{
  for(int i = 0; i < typeInfos.size(); i++)
    if(typeInfos[i]->fullName == fullTypeName)
      return i;
  return -1;
}

void ModuleFactoryNew::registerModuleType(ModuleTypeInfo* info)
{
  rassert(!doesTypeExist(info->fullName)); // type with that name was already registered...
  info->id = (int) typeInfos.size();       // ...full module names must be unique
  typeInfos.push_back(info);
}

void ModuleFactoryNew::registerStandardModules()
{
  // todo: remove the "Module" from the class names where it appears

  // Arithmetic:
  registerModuleType(new ConstantModuleTypeInfo); // constant value
  registerModuleType(new IdentityModuleTypeInfo); // a
  registerModuleType(new UnaryMinusTypeInfo);     // -a
  registerModuleType(new ReciprocalTypeInfo);     // 1/a
  registerModuleType(new AdderModuleTypeInfo);    // a+b
  registerModuleType(new SubtractorTypeInfo);     // a-b
  registerModuleType(new MultiplierTypeInfo);     // a*b
  registerModuleType(new DividerTypeInfo);        // a/b
  // todo: a^b, 
  registerModuleType(new Adder3ModuleTypeInfo);   // a+b+c
  registerModuleType(new Adder4ModuleTypeInfo);   // a+b+c+d
  registerModuleType(new Adder5ModuleTypeInfo);   // a+b+c+d+e
  registerModuleType(new AdderNModuleTypeInfo);   // a+b+c+d+e+...

  
  // Comparison: a<b a<=b, a>b, a>=b,
  // Logic: a&b, a|b

  // Functions:
  registerModuleType(new ClipperTypeInfo);    // hardClip(x, Min, Max)
  registerModuleType(new SinCosTypeInfo);     // sin(2*pi*x),cos(2*pi*x)
  registerModuleType(new TriSawTypeInfo);

  // Delays:
  registerModuleType(new UnitDelayTypeInfo);  // y = x[n-1]

  // Sources
  registerModuleType(new WhiteNoiseTypeInfo);
  registerModuleType(new PhasorTypeInfo);
  registerModuleType(new BandlimitedImpulseTrainTypeInfo);
  registerModuleType(new DualBlitSawTypeInfo);

  // Filters:
  registerModuleType(new FirstOrderLowpassTypeInfo);
  registerModuleType(new FirstOrderFilterTypeInfo);
  registerModuleType(new BiquadTypeInfo);
  registerModuleType(new BiquadDesignerTypeInfo);
  registerModuleType(new LadderFilterTypeInfo);

  // Infrastructure:
  registerModuleType(new AudioInputTypeInfo);
  registerModuleType(new AudioOutputTypeInfo);
  registerModuleType(new ParameterModuleTypeInfo);
  registerModuleType(new ContainerModuleTypeInfo);

  registerModuleType(new NoteGateTypeInfo);
  registerModuleType(new NoteOnTriggerTypeInfo);
  registerModuleType(new NoteOffTriggerTypeInfo);
  registerModuleType(new NoteFrequencyTypeInfo);
  registerModuleType(new NoteVelocityTypeInfo);
  registerModuleType(new VoiceCombinerTypeInfo);
  registerModuleType(new VoiceKillerTypeInfo);

  registerModuleType(new SystemSampleRateTypeInfo);
  registerModuleType(new SystemSamplePeriodTypeInfo);


  // Modulation:
  registerModuleType(new EnvelopeADSRTypeInfo);

  // before starting using this, compare if the new names match the old ones (full and short)
  // register also programatically built containers...but maybe do this in liberty
}

void ModuleFactoryNew::registerPreBuiltContainers()
{
  // todo:
  // TestModuleBuilder::createGain, createSumDiff, createWrappedSumDiff, createSummedDiffs, 
  // createMovingAverage, createLeakyIntegrator, createTestFilter1, createBiquadMacro, 
  // createAddedConstants, createPinSortTest, createBlip, createPolyBlipStereo, createNoiseFlute
}

void ModuleFactoryNew::clearRegisteredTypes()
{
  for(int i = 0; i < typeInfos.size(); i++)
    delete typeInfos[i];
  typeInfos.clear();
}

void ModuleFactoryNew::setupModule(romos::Module* module, const std::string& name, 
  int x, int y, bool polyphonic) const
{
  // copy code from the old ModuleFactory here
}