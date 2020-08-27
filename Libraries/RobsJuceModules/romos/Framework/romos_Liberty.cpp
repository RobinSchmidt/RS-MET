//#include "romos_Liberty.h"
//using namespace romos;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Liberty::Liberty()
{
  //populateModuleTypeRegistry();

  // these are useful only for development and testing - not for the Liberty product:
  //moduleFactory.removeModuleType("Formula_1_1");
  //moduleFactory.removeModuleType("Formula_N_1");
  // not needed anymore - they are not registered by default anymore

  //topLevelModule = (TopLevelModule*) ModuleFactory::createModule(
  //  ModuleTypeRegistry::TOP_LEVEL_MODULE, "Instrument", 0, 0, false);

  // later use:
  // topLevelModule = moduleTypeRegistry.createModule("TopLevelModule");
  // ...naahh - the top-level module needs special treatment - it should not be available in the
  // menu - maybe
  topLevelModule = moduleFactory.createTopLevelModule("Instrument", 0, 0, false);
}

Liberty::~Liberty()
{
#ifdef RS_BUILD_OLD_MODULE_FACTORY
  //moduleFactory.deleteModule(topLevelModule);
  romos::ModuleTypeRegistry::deleteSoleInstance();
#endif

  moduleFactory.deleteModule(topLevelModule);
}

//-------------------------------------------------------------------------------------------------
// setup:

void Liberty::setSampleRate(double newSampleRate)
{
  romos::processingStatus.setSystemSampleRate(newSampleRate);
}

void Liberty::reset()
{
  topLevelModule->resetStateForAllVoices();
  voiceAllocator.reset();
}

//-------------------------------------------------------------------------------------------------
// inquiry:
    

//-------------------------------------------------------------------------------------------------
// event handling:
    
void Liberty::noteOn(int key, int velocity)
{
  voiceAllocator.noteOn(key, velocity);
}

void Liberty::noteOff(int key)
{
  voiceAllocator.noteOff(key);
}

/*
void Liberty::resetAllVoices()
{
  voiceAllocator.reset();
}
*/

/*
void Liberty::populateModuleTypeRegistry()
{
  moduleTypeRegistry.registerPreBuiltContainers();
}
*/
