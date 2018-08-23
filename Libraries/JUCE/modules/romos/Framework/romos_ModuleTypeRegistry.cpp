#include "romos_ModuleTypeRegistry.h"
using namespace romos;

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

ModuleFactoryNew::ModuleFactoryNew()
{
  registerStandardModules();
}

ModuleFactoryNew::~ModuleFactoryNew()
{
  clearRegisteredTypes();
}

romos::Module* ModuleFactoryNew::createModule(int id, const std::string& name, int x, int y, 
  bool polyphonic)
{
  ensureTypeInfoArrayAllocated();
  rassert(id >= 0 && id < typeInfos->size());  // id out of range
  // todo: if the id is out of range, return some kind of "Error" dummy module, maybe just a 
  // constant module outputting a 0, named "Error"? ..if that's possible

  romos::Module* m = (*typeInfos)[id]->createModule();
  m->typeInfo = (*typeInfos)[id];
  setupModule(m, name, x, y, polyphonic);
  return m;
}

romos::Module* ModuleFactoryNew::createModule(const std::string& fullTypeName, 
  const std::string& name, int x, int y, bool polyphonic)
{
  int id = getModuleId(fullTypeName);
  return createModule(id, name, x, y, polyphonic);
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

int ModuleFactoryNew::getModuleId(const std::string& fullTypeName)
{
  ensureTypeInfoArrayAllocated();
  for(int i = 0; i < typeInfos->size(); i++)
    if((*typeInfos)[i]->fullName == fullTypeName)
      return i;
  return -1;
}

ModuleTypeInfo* ModuleFactoryNew::getModuleTypeInfo(const std::string& fullTypeName)
{
  ensureTypeInfoArrayAllocated();
  for(int i = 0; i < typeInfos->size(); i++)
    if((*typeInfos)[i]->fullName == fullTypeName)
      return (*typeInfos)[i];
  return nullptr;
}

ModuleTypeInfo* ModuleFactoryNew::getModuleTypeInfo(size_t id)
{
  if(id < typeInfos->size())
    return (*typeInfos)[id];
  return nullptr;
}

void ModuleFactoryNew::registerModuleType(ModuleTypeInfo* info)
{
  ensureTypeInfoArrayAllocated();
  rassert(!doesTypeExist(info->fullName)); // type with that name was already registered...
  info->id = (int) typeInfos->size();       // ...full module names must be unique
  typeInfos->push_back(info);
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
  registerModuleType(new SaturatorTypeInfo);  // a+b*tanh(c*x+d), a,b,c,d adjusted by width/center
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

  //registerModuleType(new TopLevelTypeInfo); // nope - should not be user creatable


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
  ensureTypeInfoArrayAllocated();
  for(int i = 0; i < typeInfos->size(); i++)
    delete (*typeInfos)[i];
  typeInfos->clear();
  delete typeInfos;
  typeInfos = nullptr;
}

void ModuleFactoryNew::ensureTypeInfoArrayAllocated()
{
  if(typeInfos == nullptr)
    typeInfos = new std::vector<ModuleTypeInfo*>;
}

void ModuleFactoryNew::setupModule(romos::Module* module, const std::string& name, 
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