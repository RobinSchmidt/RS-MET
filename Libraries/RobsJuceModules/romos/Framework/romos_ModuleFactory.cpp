//#include "romos_ModuleTypeRegistry.h"
//using namespace romos;

void ModuleTypeInfo::addInputPinInfo(
  const char* shortName, const char* fullName, const char* description)
{
  inputShortNames.push_back(shortName);
  inputFullNames.push_back(fullName);
  inputDescriptions.push_back(description);
}
// rename parameters to shortPinName, fullPinName, pinDescription to fix
// "hides class member" warning

void ModuleTypeInfo::addOutputPinInfo(
  const char* shortName, const char* fullName, const char* description)
{
  outputShortNames.push_back(shortName);
  outputFullNames.push_back(fullName);
  outputDescriptions.push_back(description);
}

//-------------------------------------------------------------------------------------------------

ModuleFactory moduleFactory;  // definition of the global object - causes memleak?

ModuleFactory::ModuleFactory()
{
  registerStandardModules();
}

ModuleFactory::~ModuleFactory()
{
  clearRegisteredTypes();
}

romos::Module* ModuleFactory::createModule(int id, const std::string& name, int x, int y,
  bool polyphonic)
{
  ensureTypeInfoArrayAllocated();
  rassert(id >= 0 && id < (int)typeInfos->size());  // id out of range
  // todo: if the id is out of range, return some kind of "Error" dummy module, maybe just a
  // constant module outputting a 0, named "Error"? ..if that's possible

  romos::Module* m = (*typeInfos)[id]->createModule();
  m->typeInfo = (*typeInfos)[id];
  setupModule(m, name, x, y, polyphonic);
  return m;
}

romos::Module* ModuleFactory::createModule(const std::string& fullTypeName,
  const std::string& name, int x, int y, bool polyphonic)
{
  int id = getModuleId(fullTypeName);
  return createModule(id, name, x, y, polyphonic);
}

romos::TopLevelModule* ModuleFactory::createTopLevelModule(const std::string& name,
  int x, int y, bool polyphonic) const
{
  TopLevelModule* tlm = new TopLevelModule();
  setupModule(tlm, name, x, y, polyphonic);
  return tlm;
}

void ModuleFactory::deleteModule(romos::Module* moduleToDelete)
{
  moduleToDelete->cleanUp();
  delete moduleToDelete;
}

int ModuleFactory::getModuleId(const std::string& fullTypeName)
{
  ensureTypeInfoArrayAllocated();
  for(size_t i = 0; i < typeInfos->size(); i++)
    if((*typeInfos)[i]->fullName == fullTypeName)
      return (int)i;
  return -1;
}

ModuleTypeInfo* ModuleFactory::getModuleTypeInfo(const std::string& fullTypeName)
{
  ensureTypeInfoArrayAllocated();
  for(size_t i = 0; i < typeInfos->size(); i++)
    if((*typeInfos)[i]->fullName == fullTypeName)
      return (*typeInfos)[i];
  return nullptr;
}

ModuleTypeInfo* ModuleFactory::getModuleTypeInfo(size_t id)
{
  if(id < typeInfos->size())
    return (*typeInfos)[id];
  return nullptr;
}

void ModuleFactory::registerModuleType(ModuleTypeInfo* info)
{
  ensureTypeInfoArrayAllocated();
  rassert(!doesTypeExist(info->fullName)); // type with that name was already registered...
  info->id = (int) typeInfos->size();       // ...full module names must be unique
  typeInfos->push_back(info);
}

void ModuleFactory::removeModuleType(const std::string& fullTypeName)
{
  ensureTypeInfoArrayAllocated();
  for(size_t i = 0; i < typeInfos->size(); i++)
    if((*typeInfos)[i]->fullName == fullTypeName) {
      delete (*typeInfos)[i];
      RAPT::rsRemove(*typeInfos, i);
      break;
    }
  for(size_t i = 0; i < typeInfos->size(); i++) // update type ids
    (*typeInfos)[i]->id = (int)i;
}

void ModuleFactory::registerStandardModules()
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
  //registerModuleType(new ScalerTypeInfo);         // a*const ...not yet finished
  // todo: a^b,
  registerModuleType(new Adder3ModuleTypeInfo);   // a+b+c
  registerModuleType(new Adder4ModuleTypeInfo);   // a+b+c+d
  registerModuleType(new Adder5ModuleTypeInfo);   // a+b+c+d+e
  registerModuleType(new AdderNModuleTypeInfo);   // a+b+c+d+e+...


  // Comparison: a<b a<=b, a>b, a>=b,
  // Logic: a&b, a|b

  // Functions:
  registerModuleType(new ClipperTypeInfo);    // hardClip(x, Min, Max)
  registerModuleType(new SaturatorTypeInfo);  // a+b*tanh(c*x+d), a,b,c,d adjusted by width/center
  registerModuleType(new SinCosTypeInfo);     // sin(2*pi*x),cos(2*pi*x)
  registerModuleType(new TriSawTypeInfo);
  //registerModuleType(new FormulaModule_1_1TypeInfo);
  //registerModuleType(new FormulaModule_N_1TypeInfo);
  registerModuleType(new FormulaModule_N_MTypeInfo);


  // Delays:
  registerModuleType(new UnitDelayTypeInfo);  // y = x[n-1]

  // Sources
  registerModuleType(new WhiteNoiseTypeInfo);
  registerModuleType(new PhasorTypeInfo);
  registerModuleType(new SineOscillatorTypeInfo);
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

  //registerModuleType(new TopLevelTypeInfo); // nope - should not be user creatable


  // Modulation:
  registerModuleType(new EnvelopeADSRTypeInfo);

  // before starting using this, compare if the new names match the old ones (full and short)
  // register also programatically built containers...but maybe do this in liberty
}

void ModuleFactory::registerPreBuiltContainers()
{
  // todo:
  // TestModuleBuilder::createGain, createSumDiff, createWrappedSumDiff, createSummedDiffs,
  // createMovingAverage, createLeakyIntegrator, createTestFilter1, createBiquadMacro,
  // createAddedConstants, createPinSortTest, createBlip, createPolyBlipStereo, createNoiseFlute
}

void ModuleFactory::clearRegisteredTypes()
{
  ensureTypeInfoArrayAllocated();
  for(size_t i = 0; i < typeInfos->size(); i++)
    delete (*typeInfos)[i];
  typeInfos->clear();
  delete typeInfos;
  typeInfos = nullptr;
}

void ModuleFactory::ensureTypeInfoArrayAllocated()
{
  if(typeInfos == nullptr)
    typeInfos = new std::vector<ModuleTypeInfo*>;
}

void ModuleFactory::setupModule(romos::Module* module, const std::string& name,
  int x, int y, bool polyphonic) const
{
  // copied from the old ModuleFactory::createModule

  module->initialize(); // this should set up the number of pins needed, etc.
  module->setPositionXY(x, y);
  module->setPolyphonic(polyphonic);
  module->allocateMemory();

  module->setModuleName(name);
  // must be called after allocateMemory because the Constant fills its output arrays with the
  // corresponding value

  module->assignProcessingFunctions();
  module->resetStateForAllVoices();
}
