//#include "romos_TopLevelModule.h"
//using namespace romos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

TopLevelModule::TopLevelModule()
{
  //moduleTypeIdentifier = ModuleTypeRegistry::TOP_LEVEL_MODULE;
  typeInfo = moduleFactory.getModuleTypeInfo("Container"); 

  ContainerModule::addAudioInputModule( "AudioIn1",   2,  2);
  ContainerModule::addAudioInputModule( "AudioIn2",   2,  6);
  ContainerModule::addAudioOutputModule("AudioOut1", 40,  2);
  ContainerModule::addAudioOutputModule("AudioOut2", 40,  6);
    // yes, we indeed want to call the baseclass implementation here - our own overriden functions do nothing in order to prevent adding
    // I/O modules by the user / client-code

  initInputBuffersAndModules();
}

TopLevelModule::~TopLevelModule()
{
  cleanUpInputBuffersAndModules();

  // clean up the I/O modules (the baseclass cleanUp will have called out overriden deleteChildModule() method, which leaves the I/O
  // modules behind:
  while( childModules.size() > 0 )
  { 
    Module *moduleToDelete = childModules[childModules.size()-1];

    for(unsigned int i = 0; i < childModules.size(); i++)
      childModules[i]->disconnectInputPinsWithInputFrom(moduleToDelete); 

    rosic::removeElementByValue(childModules, moduleToDelete);
    //childModules.removeElementByValue(moduleToDelete);

    moduleFactory.deleteModule(moduleToDelete);
  }

  outFrameStride = 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void TopLevelModule::disconnectAudioOutputModules()
{
  Module *outModule;

  outModule = getChildModule(getNumChildModules()-1);
  outModule->disconnectAllInputPins();

  outModule = getChildModule(getNumChildModules()-2);
  outModule->disconnectAllInputPins();
}

void TopLevelModule::deleteChildModule(Module *moduleToDelete, bool updateHasDelayedConnectionFlag)
{
  if( !(moduleToDelete->isInputModule() || moduleToDelete->isOutputModule()) )
    ContainerModule::deleteChildModule(moduleToDelete, updateHasDelayedConnectionFlag);
}

void TopLevelModule::sortChildModuleArray()
{
  sortModuleArrayByCoordinates(childModules);
  updateHasDelayedConnectionsFlag();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:



//-----------------------------------------------------------------------------------------------------------------------------------------  
// processing:


//-----------------------------------------------------------------------------------------------------------------------------------------
// internal functions:

void TopLevelModule::initInputBuffersAndModules()
{
  allocateInputBuffers();
  clearInputBuffers();
  setupInputModules();
}

void TopLevelModule::cleanUpInputBuffersAndModules()
{
  deleteInputBuffers();
  setupInputModules();
}

void TopLevelModule::allocateInputBuffers()
{
  inL = new double[processingStatus.getBufferSize()];
  inR = new double[processingStatus.getBufferSize()];
}

void TopLevelModule::deleteInputBuffers()
{
  delete[] inL;
  delete[] inR;
  inL = NULL;
  inR = NULL;
}

void TopLevelModule::clearInputBuffers()
{
  memset(inL, 0, processingStatus.getBufferSize() * sizeof(double));
  memset(inR, 0, processingStatus.getBufferSize() * sizeof(double));
}

void TopLevelModule::setupInputModules()
{
  Module *inModule;

  inModule = getChildModule(0);
  inModule->inputPins[0].outputPointer   = inL;
  inModule->inputPins[0].outputFrameSize = 1;

  inModule = getChildModule(1);
  inModule->inputPins[0].outputPointer   = inR;
  inModule->inputPins[0].outputFrameSize = 1;
}


