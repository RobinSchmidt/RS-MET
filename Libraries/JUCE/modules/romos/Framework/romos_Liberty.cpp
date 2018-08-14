#include "romos_Liberty.h"
#include "romos_ModuleFactory.h"
using namespace romos;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Liberty::Liberty()
{
  topLevelModule = (TopLevelModule*) ModuleFactory::createModule(ModuleTypeRegistry::TOP_LEVEL_MODULE, "Instrument", 0, 0, false);
  populateModuleTypeRegistry();
}

Liberty::~Liberty()
{
  ModuleFactory::deleteModule(topLevelModule);
  romos::ModuleTypeRegistry::deleteSoleInstance();
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

void Liberty::populateModuleTypeRegistry()
{
  ModuleTypeRegistry2* mtr = &moduleTypeRegistry;
  //mtr->registerModuleType(AdderTypeInfo);
  //mtr->registerModuleType(SubtractorTypeInfo);
}
