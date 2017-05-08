renamed to .x so the Projucer doesn't include it for building into the JuceLibraryCode. has
eventually to go to jura anyway

#include "romos_Liberty.h"
#include "framework/romos_ModuleFactory.h"
using namespace romos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

Liberty::Liberty()
{
  topLevelModule = (TopLevelModule*) ModuleFactory::createModule(ModuleTypeRegistry::TOP_LEVEL_MODULE, "Instrument", 0, 0, false);
}

Liberty::~Liberty()
{
  ModuleFactory::deleteModule(topLevelModule);
  romos::ModuleTypeRegistry::deleteSoleInstance();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:
    

//-----------------------------------------------------------------------------------------------------------------------------------------   
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
