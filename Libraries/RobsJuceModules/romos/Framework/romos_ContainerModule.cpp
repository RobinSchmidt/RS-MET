//#include "romos_ContainerModule.h"
//using namespace romos;

//-------------------------------------------------------------------------------------------------
// the processing functions:

// fallback function to process the module in frames when necessary:
void processContainerBlockFrameWiseMono(Module *moduleAsVoid, int blockSize)
{
  // algorithm:
  // -adjust the input-pins of the container such that they point to some temporary input area
  // -gather input frame 0 into temporary input
  // -compute output frame 0
  // -store output frame 0 in temporary output (it's going to be overwritten in the subsequent loop)
  // -for n = 1 ... blocksize-1
  //   -copy input frame n into temporary input
  //   -call containers' per-frame processing function
  //   -copy output frame 0 into output frame n
  //  end for
  // -copy the temporary output frame into output frame 0
  // -restore the original input pin data for the container

  romos::ContainerModule *container = ((romos::ContainerModule *) moduleAsVoid);

  int    pinIndex;
  int    inFrameSize    = container->numInputs;
  int    outFrameSize   = container->outFrameStride;
  double *outputPointer = container->getOutputPointer(0);

  // temporarily modify the input-pin data in the container's input modules:
  for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
    container->childModules[pinIndex]->inputPins[0].outputPointer
    = &(WorkArea::tmpInFrame[pinIndex]);

  // compute 0th frame and buffer it in temp-variable - we'll need it later to restore the 0th
  // frame which will repeatedly be overwritten in the subsequent loop:
  for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
    WorkArea::tmpInFrame[pinIndex] = *(container->inputPins[pinIndex].outputPointer);
  container->processFrame(container, 0);
  for( pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
    WorkArea::tmpOutFrame[pinIndex] = outputPointer[pinIndex];  // maybe try memcopy instead

  for(int frameIndex = 1; frameIndex < blockSize; frameIndex++)
  {
    // copy input into temporary location:
    for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
      WorkArea::tmpInFrame[pinIndex] = * (container->inputPins[pinIndex].outputPointer
                                          + frameIndex * container->inputPins[pinIndex].outputFrameSize);

    // let the container process a single frame:
    container->processFrame(container, 0);

    // copy output from 0th frame to where it should end up:
    for(pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
      outputPointer[frameIndex*outFrameSize + pinIndex] = outputPointer[pinIndex];
  }

  // we have used the 0th output frame throughout the loop, thereby repeatedly overwriting it so we
  // must now restore the 0th frame of the block (for which we have a temporary variable):
  for(pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
    outputPointer[pinIndex] = WorkArea::tmpOutFrame[pinIndex];   // maybe try memcopy instead

  // restore the original pin data in the container's input modules:
  for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
    container->childModules[pinIndex]->inputPins[0].outputPointer
    = container->inputPins[pinIndex].outputPointer;


  // can we avoid this copy/restore stuff by calling something like
  // container->processFrame(container, frameIndex);
  // instead of:
  // container->processFrame(container, 0);
  // in the loop?

  // maybe it is possible to split the modules inside a container into those that are participating
  // in a feedback loop and those who don't and use per-sample processing only for those who do?
  // like:
  // process() {
  //   processModulesPreFeedback();    // uses block-wise processing
  //   processModulesWithFeedback();   // uses frame-wise processing
  //   processModulesPostFeedback(); } // uses block-wise processing
}

// little helper:
void copyMatrix(double *source, double *destination, int numRows, int numColumns)
{
  for(int i = 0; i < numRows*numColumns; i++)
    destination[i] = source[i];
}

// fallback function to process the module in frames when necessary:
void processContainerBlockFrameWisePoly(Module *moduleAsVoid, int blockSize)
{
  romos::ContainerModule *container = ((romos::ContainerModule *) moduleAsVoid);
  const int *playingVoiceIndices    = voiceAllocator.getPlayingVoiceIndices();
  int numPlayingVoices              = voiceAllocator.getNumPlayingVoices();
  int bufferSize                    = processingStatus.getBufferSize();

  int    pinIndex, playIndex, voiceIndex;
  int    inFrameSize    = container->numInputs;
  int    outFrameSize   = container->outFrameStride;
  double *outputPointer = container->getOutputPointer(0);
  int outVoiceStride =
    outFrameSize * processingStatus.getBufferSize() * (int) container->polyphonic;

  // temporarily modify the input-pin data in the container's input modules such that they point
  // to some global memory area - we need to set the voice-strides temporarily to different values
  // too because in the global memory area we have a lower voice-stride (because there, we don't
  // need to multiply with the buffersize - there's only one buffered frame in the global area):
  for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
  {
    container->childModules[pinIndex]->inputPins[0].outputPointer     = &(WorkArea::tmpInFramePoly[pinIndex]);
    container->childModules[pinIndex]->inputPins[0].outputVoiceStride = WorkArea::maxNumPins;
  }

  // compute 0th frame for each voice and buffer it:
  for(playIndex = 0; playIndex < numPlayingVoices; playIndex++)
  {
    voiceIndex = playingVoiceIndices[playIndex];
    for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
    {
      WorkArea::tmpInFramePoly[voiceIndex * WorkArea::maxNumPins + pinIndex] =
        *(container->inputPins[pinIndex].outputPointer
          + voiceIndex * container->inputPins[pinIndex].outputVoiceStride);
    }
  }
  container->processFrame(container, 0);
  for(playIndex = 0; playIndex < numPlayingVoices; playIndex++)
  {
    voiceIndex    = playingVoiceIndices[playIndex];
    outputPointer = container->audioOutputs + voiceIndex * outVoiceStride;
    for(pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
      WorkArea::tmpOutFramePoly[voiceIndex * WorkArea::maxNumPins + pinIndex]
      = outputPointer[pinIndex];
  }

  // compute the other frames and put them into their target locations:
  for(int frameIndex = 1; frameIndex < blockSize; frameIndex++)
  {
    for(playIndex = 0; playIndex < numPlayingVoices; playIndex++)
    {
      voiceIndex = playingVoiceIndices[playIndex];
      for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
      {
        WorkArea::tmpInFramePoly[voiceIndex * WorkArea::maxNumPins + pinIndex] =
          *(container->inputPins[pinIndex].outputPointer
            + voiceIndex * container->inputPins[pinIndex].outputVoiceStride
            + frameIndex * container->inputPins[pinIndex].outputFrameSize);
      }
    }
    container->processFrame(container, 0);
    for(playIndex = 0; playIndex < numPlayingVoices; playIndex++)
    {
      voiceIndex    = playingVoiceIndices[playIndex];
      outputPointer = container->audioOutputs + voiceIndex * outVoiceStride;
      for(pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
        outputPointer[frameIndex*outFrameSize + pinIndex] = outputPointer[pinIndex];
    }
  }

  // restore the buffered 0th frame:
  for(playIndex = 0; playIndex < numPlayingVoices; playIndex++)
  {
    voiceIndex    = playingVoiceIndices[playIndex];
    outputPointer = container->audioOutputs + voiceIndex * outVoiceStride;
    for(pinIndex = 0; pinIndex < outFrameSize; pinIndex++)
      outputPointer[pinIndex] = WorkArea::tmpOutFramePoly[voiceIndex * WorkArea::maxNumPins + pinIndex];
  }

  // restore the original pin data in the container's input modules:
  for(pinIndex = 0; pinIndex < inFrameSize; pinIndex++)
  {
    container->childModules[pinIndex]->inputPins[0].outputPointer     = container->inputPins[pinIndex].outputPointer;
    container->childModules[pinIndex]->inputPins[0].outputFrameSize   = container->inputPins[pinIndex].outputFrameSize; // not needed?
    container->childModules[pinIndex]->inputPins[0].outputVoiceStride = container->inputPins[pinIndex].outputVoiceStride;
  }
}

// gcc complains
void processContainerBlockFrameWiseMixed(Module *moduleAsVoid, int blockSize)
{

}


// frame-wise processing when we have monophonic as well as polyphonic child-modules (needs dispatch-logic for the accumulation in the
// connections):
void processContainerMixedMonoPoly(Module *module, int voiceIndex)
{
  ContainerModule *container = (ContainerModule*) module;
  romos::Module   *cm;  // currently vistited child module (is also the target module of the currently visited connection in the loops)
  for(unsigned int i = 0; i < container->childModules.size(); i++)
  {
    cm = container->childModules[i];
    cm->processFrame(cm, voiceIndex);
    //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
    // \todo: do the poly/mono, mono/poly conversion
  }
}

// frame-wise processing when all child-modules are monophonic:
void processContainerAllMono(Module *module, int voiceIndex)
{
  ContainerModule *container = (ContainerModule*) module;
  romos::Module   *cm;
  for(unsigned int i = 0; i < container->childModules.size(); i++)
  {
    cm = container->childModules[i];
    cm->processFrame(cm, voiceIndex);
    //romos::retrieveModuleState(cm);
    //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
  }
}

// frame-wise processing when all child-modules are polyphonic:
void processContainerAllPoly(Module *module, int voiceIndex)
{
  ContainerModule *container = (ContainerModule*) module;
  romos::Module   *cm;

  for(unsigned int i = 0; i < container->childModules.size(); i++)
  {
    cm = container->childModules[i];
    cm->processFrame(cm, voiceIndex);
    //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
  }
}

// block-wise processing when we have monophonic as well as polyphonic child-modules (needs dispatch-logic for the accumulation in the
// connections):
void processContainerMixedMonoPolyBlock(Module *module, int voiceIndex, int blockSize)
{
  if( ((ContainerModule*) module)->hasDelayedConnections )
    processContainerBlockFrameWiseMixed(module, blockSize);
  else
  {
    ContainerModule *container = (ContainerModule*) module;
    romos::Module   *cm;
    for(unsigned int i = 0; i < container->childModules.size(); i++)
    {
      cm            = container->childModules[i];
      cm->processBlock(cm, voiceIndex, blockSize);
      //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
      // \todo: do the poly/mono, mono/poly conversion
    }
  }
}

void processContainerAllMonoBlock(Module *module, int voiceIndex, int blockSize)
{
  if( ((ContainerModule*) module)->hasDelayedConnections )
    processContainerBlockFrameWiseMono(module, blockSize);
  else
  {
    ContainerModule *container = (ContainerModule*) module;
    romos::Module   *cm;
    for(unsigned int i = 0; i < container->childModules.size(); i++)
    {
      cm = container->childModules[i];
      cm->processBlock(cm, voiceIndex, blockSize);  // wrong! pass outs as 2nd argument

      //romos::retrieveModuleState(cm);               // may be uncommented for debugging
      //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
    }
  }
}

void processContainerAllPolyBlock(Module *module, int voiceIndex, int blockSize)
{
  if( ((ContainerModule*) module)->hasDelayedConnections )
    processContainerBlockFrameWisePoly(module, blockSize);
  else
  {
    ContainerModule *container = (ContainerModule*) module;
    romos::Module   *cm;  // currently vistited child module (is also the target module of the currently visited connection in the loops)
    for(unsigned int i = 0; i < container->childModules.size(); i++)
    {
      cm = container->childModules[i];
      cm->processBlock(cm, voiceIndex, blockSize);
      //romos::writeModuleStateToConsole(cm, true); // may be uncommented for debugging
    }
  }
}


//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ContainerModule::ContainerModule(const std::string& name, int x, int y, bool polyphonic)
: Module(name, x, y, polyphonic)
{
  tmpOutFrame          = NULL;
  //moduleTypeIdentifier = ModuleTypeRegistry::CONTAINER; // not strictly necessary, factory will set this up also

  if( name == std::string() )
    setModuleName(std::string("Container"));

  assignProcessingFunctions();
  updateHasDelayedConnectionsFlag();
}

ContainerModule::~ContainerModule()
{
  //cleanUp();
}


void ContainerModule::initialize()
{
  //tmpOutFrame = NULL; no! this sets the pointer to NULL after it has possibly been allocated
  updateHasDelayedConnectionsFlag();
}

void ContainerModule::cleanUp()
{
  freeMemory();
  deleteAllChildModules();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void ContainerModule::setPolyphonic(bool shouldBePolyphonic)
{
  Module::setPolyphonic(shouldBePolyphonic);

  //std::vector<romos::Module*> inputModules = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_INPUT);
  std::vector<romos::Module*> inputModules = getChildModulesWithType("AudioInput");

  for(unsigned int i = 0; i < inputModules.size(); i++)
    inputModules[i]->setPolyphonic(shouldBePolyphonic);

  //std::vector<romos::Module*> outputModules = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_OUTPUT);
  std::vector<romos::Module*> outputModules = getChildModulesWithType("AudioOutput");

  for(unsigned int i = 0; i < outputModules.size(); i++)
    outputModules[i]->setPolyphonic(shouldBePolyphonic);

  assignProcessingFunctions();
}

void ContainerModule::setPolyphonicRecursively(bool shouldBePolyphonic)
{
  Module::setPolyphonic(shouldBePolyphonic);

  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules[i]->isContainerModule() )
      ((ContainerModule *) childModules[i])->setPolyphonicRecursively(shouldBePolyphonic);
    else
      childModules[i]->setPolyphonic(shouldBePolyphonic);
  }

  assignProcessingFunctions();
}

void ContainerModule::connectInputPinTo(int inputPinIndex, Module *sourceModule, int sourceOutputPinIndex)
{
  Module::connectInputPinTo(inputPinIndex, sourceModule, sourceOutputPinIndex);
  getAudioInputModule(inputPinIndex)->connectInputPinTo(0, sourceModule, sourceOutputPinIndex);
}

void ContainerModule::disconnectInputPin(int inputPinIndex)
{
  Module::disconnectInputPin(inputPinIndex);
  getAudioInputModule(inputPinIndex)->disconnectInputPin(0);
}

romos::Module* ContainerModule::addAudioInputModule(std::string name, int x, int y,
  bool sortModuleArrayAfterInsertion)
{
  if( name.empty() )
    name = std::string("In") + std::to_string(getNumInputPins()+1);

  //Module *newModule = ModuleFactory::createModule(ModuleTypeRegistry::AUDIO_INPUT, name, x, y, this->isPolyphonic());
  //newModule->typeInfo = moduleFactory.getModuleTypeInfo("AudioInput");
    // can be deleted when we create the newModule with the newFactory later (but until then, we
    // need to manually set the typeInfo pointer here

  Module *newModule = moduleFactory.createModule("AudioInput", name, x, y, this->isPolyphonic());
  newModule->parentModule = this;
  rosic::appendElement(childModules, (romos::Module*) newModule);

  if( sortModuleArrayAfterInsertion == true )
    sortChildModuleArray();

  rosic::appendElement(inputPins, AudioInputPinData());
  numInputs++;

  //allocateAudioInputs();
  updateInputPointersAndInFrameStrides();

  return newModule;
}

romos::Module* ContainerModule::addAudioOutputModule(std::string name, int x, int y,
  bool sortModuleArrayAfterInsertion)
{
  if( name.empty() )
    name = std::string("Out") + std::to_string(getNumOutputPins()+1);

  //Module *newModule = ModuleFactory::createModule(ModuleTypeRegistry::AUDIO_OUTPUT, name, x, y, this->isPolyphonic());
  //newModule->typeInfo = moduleFactory.getModuleTypeInfo("AudioOutput");
  // can be deleted when we create the newModule with the newFactory later (but until then, we
  // need to manually set the typeInfo pointer here

  Module *newModule = moduleFactory.createModule("AudioOutput", name, x, y, this->isPolyphonic());

  newModule->parentModule = this;
  rosic::appendElement(childModules, (romos::Module*) newModule);
  if( sortModuleArrayAfterInsertion == true )
    sortChildModuleArray();

  //outFrameStride++;  // old
  outFrameStride = 1;

  allocateAudioOutputs();

  return newModule;
}

romos::Module* ContainerModule::addChildModule(Module *moduleToAdd,
  bool sortChildModuleArrayAfterInsertion)
{
  rassert( !hasAsDirectlyEmbeddedModule(moduleToAdd) ); // adding a child multiple times?


  rassert( !moduleToAdd->isInputModule() );             // use addAudioInputModule instead
  rassert( !moduleToAdd->isOutputModule() );            // use addAudioOutputModule instead
    // \todo: check here whether the module is an I/O module and call the proper function instead of requiring the client code to do it

  rosic::appendElement(childModules, moduleToAdd);
  moduleToAdd->parentModule = this;

  if( sortChildModuleArrayAfterInsertion == true )
    sortChildModuleArray();
    // todo: figure out where to insert and insert it directly in the right place - get rid of the
    // boolean parameter

  return moduleToAdd;
}

/*
romos::Module* ContainerModule::addChildModule(int moduleIdentifier, rosic::rsString name,
  int x, int y, bool polyphonic, bool sortChildModulesAfterInsertion)
{
  if( moduleIdentifier == ModuleTypeRegistry::AUDIO_INPUT )
    return addAudioInputModule( name, x, y, true);
  else if( moduleIdentifier == ModuleTypeRegistry::AUDIO_OUTPUT )
    return addAudioOutputModule(name, x, y, true);
  else
  {
    romos::Module *moduleToAdd = ModuleFactory::createModule(moduleIdentifier, name, x, y,
      polyphonic);

    addChildModule(moduleToAdd, sortChildModulesAfterInsertion);
    return moduleToAdd;
  }
}
*/

// new version - not yet tested:
romos::Module* ContainerModule::addChildModule(const std::string& fullTypeName,
  const std::string& name, int x, int y, bool poly, bool sortChildModulesAfterInsertion)
{
  //rassert(false); return nullptr; // does not yet work

  if( fullTypeName == "AudioInput" )
    return addAudioInputModule( name, x, y, true);
  else if( fullTypeName == "AudioOutput" )
    return addAudioOutputModule(name, x, y, true);
  else
  {
    std::string nameToUse;
    if(name == "")
      nameToUse = moduleFactory.getShortTypeName(fullTypeName);
    else
      nameToUse = name;
    romos::Module *moduleToAdd = moduleFactory.createModule(fullTypeName, nameToUse, x, y, poly);
    addChildModule(moduleToAdd, sortChildModulesAfterInsertion);
    return moduleToAdd;
  }
}




void ContainerModule::deleteChildModule(Module *moduleToDelete, bool updateHasDelayedConnectionFlag)
{

  if( rosic::containsElement(childModules, moduleToDelete) )
  {
    moduleToDelete->disconnectAllOutputPins();
    int index = rosic::findElement(childModules, moduleToDelete);

    if( moduleToDelete->isInputModule() )
    {
      rosic::removeElementByIndex(childModules, index);
      rosic::removeElementByIndex(inputPins,    index); // removes corresponding pin-data object (should have same index)
      numInputs--;
    }
    else if( moduleToDelete->isOutputModule() )
    {
      int numNonOutputModules = getNumChildModules() - getNumOutputPins();
      int outputIndex         = index - numNonOutputModules;

      std::vector<Module*> targetModules = getConnectedTargetModules();
      for(unsigned int targetIndex = 0; targetIndex < targetModules.size(); targetIndex++)
      {
        Module* targetModule = targetModules[targetIndex];
        for(unsigned int pinIndex = 0; pinIndex < targetModule->getNumInputPins(); pinIndex++)
        {
          if( targetModule->inputPins[pinIndex].sourceModule == this )
          {
            if( (int)targetModule->inputPins[pinIndex].outputIndex == outputIndex )
              targetModule->disconnectInputPin(pinIndex);
            else if( (int) targetModule->inputPins[pinIndex].outputIndex > outputIndex )
            {
              targetModule->inputPins[pinIndex].outputIndex -= 1;
              // moves up subsequent connection by one position - the rest of the update will be done in the updater-callback that is
              // spawned from allocateAudioOutputs()
            }
          }
        }
      }

      //outFrameStride--; // why that - this variable should always be unity -  obsolete now?
      rosic::removeElementByIndex(childModules, index);
      allocateAudioOutputs();
    }
    else
      rosic::removeElementByIndex(childModules, index);

    moduleFactory.deleteModule(moduleToDelete);  // old
    //moduleFactory.deleteModule(moduleToDelete);  // new - activate later

    if( updateHasDelayedConnectionFlag == true )
      updateHasDelayedConnectionsFlag();
  }
  else
    triggerRuntimeError(
    "ContainerModule::deleteChildModule called with module that isn't a child of the container for which it was called");
}

void ContainerModule::deleteAllChildModules()
{
  for(int i = (int) childModules.size()-1; i>=0; i--)
    deleteChildModule(childModules[i], false);
  updateHasDelayedConnectionsFlag();

  //while( childModules.getNumElements() > 0 )
  //  deleteChildModule( childModules[childModules.getNumElements()-1] );
  // while loop becaomes infinite when the overriden TopLevelModule::deleteChildModule is called
}

void ContainerModule::deleteModules(std::vector<Module*> modulesToDelete)
{
  for(unsigned int i = 0; i < modulesToDelete.size(); i++)
    deleteChildModule(modulesToDelete[i], false);
  updateHasDelayedConnectionsFlag();
}

void ContainerModule::setPolyphonyForModules(std::vector<Module*> modules, bool shouldBePolyphonic)
{
  for(unsigned int i = 0; i < modules.size(); i++)
    modules[i]->setPolyphonic(shouldBePolyphonic);
}

ContainerModule* ContainerModule::containerizeModules(std::vector<Module*> modulesToContainerize)
{
  // algorithm:
  /*
  -create a new container and add it as child here
  for all modules to be containerized:

   -remove the module as child here (2)
   -add the module as child to the new container (3)
    ->nah move this to the end ..update thsi comment

   for all of the module's input pins:
    if the source of the connection is not among the to-be-containerized modules:
     -create an input module inside the container
     -redirect the source of the pin to the new input module
     -connect the input pin of the container that corresponds to the new input module to the original source-module
    end if
   end for
   for all of the module's output pins:
    if the output pin has targets that are not among the to-be-containerized modules:
     -create an output module inside the container
     for all modules that are target of the current output pin:
      if the target of the connection is not among the to-be-containerized modules:
       -connect the created output-module to the to-be-containerized module
       -connect the original target-module to the output pin of the container that corresponds to the created output module
      end if
     end for
   end for
  end for
  -sort our array of child modules
  */

  // create a new container and add it as child here:
  int xMin, yMin, xMax, yMax, x, y;
  getExtremeCoordinates(modulesToContainerize, xMin, yMin, xMax, yMax);
  getMidpointCoordinates(modulesToContainerize, x, y);

  //ContainerModule *container = new ContainerModule("Container", x, y, isPolyphonic()); // old
  ContainerModule* container =
    (ContainerModule*) moduleFactory.createModule("Container", "Container", x, y, isPolyphonic());

  addChildModule(container, false);

  // we don't want to containerize I/O modules:
  //removeModulesOfType(modulesToContainerize, ModuleTypeRegistry::AUDIO_INPUT);
  //removeModulesOfType(modulesToContainerize, ModuleTypeRegistry::AUDIO_OUTPUT);
  removeModulesOfType(modulesToContainerize, "AudioInput");
  removeModulesOfType(modulesToContainerize, "AudioOutput");


  // loop over the to-becontainerized modules:
  romos::Module *module;
  int numContainerInputs  = 0;
  int numContainerOutputs = 0;
  for(unsigned int i = 0; i < modulesToContainerize.size(); i++)
  {
    module = modulesToContainerize[i];

    // for all of the module's input pins:

    for(unsigned int pinIndex = 0; pinIndex < module->getNumInputPins(); pinIndex++)
    {

      // if the source of the connection is not among the to-be-containerized modules:
      romos::Module *sourceModule = module->inputPins[pinIndex].sourceModule;
      if( sourceModule != NULL && !rosic::containsElement(modulesToContainerize, sourceModule) )
      {
        // create an input module inside the container:
        numContainerInputs++;
        int sourceOutIndex = module->inputPins[pinIndex].outputIndex; // temporary storage - value needed later
        romos::Module *inputModule = container->addAudioInputModule(
          std::string("In") + std::to_string(numContainerInputs), 1, 2*numContainerInputs, false);

        // redirect the source of the pin to the new input module:
        module->connectInputPinTo(pinIndex, inputModule, 0);

        // connect the input pin of the container that corresponds to the new input module to the original source-module:
        container->connectInputPinTo(numContainerInputs-1, sourceModule, sourceOutIndex);
        //int dummy  = 0;
      }
    }

    // for all of the module's output pins:
    for(unsigned int pinIndex = 0; pinIndex < module->getNumOutputPins(); pinIndex++)
    {
      // obtain a vector of modules that are targets of the current pin and will not go themselves into the container:
      std::vector<romos::Module*> outsidePinTargets = module->getConnectedTargetModulesOfPin(pinIndex);
      rosic::removeAllOccurencesOfElements(outsidePinTargets, modulesToContainerize);

      if( outsidePinTargets.size() > 0 )
      {
        // create output module:
        numContainerOutputs++;
        romos::Module *outputModule =
          container->addAudioOutputModule(std::string("Out") + std::to_string(numContainerOutputs),
                                          xMax+20, 2*numContainerOutputs, false);

        // re-connect the pins of the outside target modules to the output pin of the container and re-connect the original source to the
        // output module in the container:
        for(unsigned int targetModuleIndex = 0; targetModuleIndex < outsidePinTargets.size(); targetModuleIndex++)
        {
          romos::Module* targetModule = outsidePinTargets[targetModuleIndex];
          for(unsigned int targetInputIndex = 0; targetInputIndex < targetModule->getNumInputPins(); targetInputIndex++)
          {
            if( targetModule->inputPins[targetInputIndex].sourceModule == module )
            {
              int tmpOutIndex = targetModule->inputPins[targetInputIndex].outputIndex;
              outputModule->connectInputPinTo(0, module, tmpOutIndex);
              targetModule->connectInputPinTo(targetInputIndex, container, numContainerOutputs-1);
              //int dummy = 0;
            }
          }
        }
      }
    }

    // remove the module as child here:
    rosic::removeElementByValue(childModules, module);

    // add the module as child to the new container:
    container->addChildModule(module, false);
  }

  container->sortChildModuleArray();
  sortChildModuleArray();
  return container;
}

void ContainerModule::unContainerizeModules(std::vector<Module*> modulesToUnContainerize)
{
  romos::ContainerModule *container;
  for(unsigned int i = 0; i < modulesToUnContainerize.size(); i++)
  {
    container = dynamic_cast<romos::ContainerModule*> (modulesToUnContainerize[i]);
    if( container != NULL )
      unContainerize(container);
  }
}

void ContainerModule::unContainerize(ContainerModule *container)
{
  // algorithm:
  // for each child-of-container
  //  -check input pins for connection to an input module
  //  -if one is found, re-connect it to the output of the outlying source that has a conenction into the container
  // end for
  // for each child-of-this
  //  -check input pins if they have a connection coming from the container
  //  -if one is found, re-connect the child to the inner source module of the container
  //   (unless the source is an input module, in this case do nothing)
  // end for
  // make the container's non I/O child modules direct children of this


  // handle connections from children of "this" into the container:
  for(unsigned int childIndex = 0; childIndex < container->getNumChildModules(); childIndex++)
  {
    Module *childModule = container->getChildModule(childIndex);
    for(unsigned int pinIndex = 0; pinIndex < childModule->getNumInputPins(); pinIndex++)
    {
      Module *innerSourceModule = childModule->inputPins[pinIndex].sourceModule;
      if( innerSourceModule != NULL && innerSourceModule->isInputModule() )
      {
        // child has connection to an input module - find out what is conneted to the container's corresponding input pin and re-connect
        // the child directly to this source:
        int innerSourceIndex = container->getIndexOfChildModule(innerSourceModule); // also the index of the conatiners input pin
        Module *outerSourceModule  = container->inputPins[innerSourceIndex].sourceModule;
        int outerSourceOutPinIndex = container->inputPins[innerSourceIndex].outputIndex;
        childModule->connectInputPinTo(pinIndex, outerSourceModule, outerSourceOutPinIndex);
      }
    }
  }

  // handle connections from the container into other children of "this":
  for(unsigned int childIndex = 0; childIndex < this->getNumChildModules(); childIndex++)
  {
    Module *childModule = this->getChildModule(childIndex);
    for(unsigned int pinIndex = 0; pinIndex < childModule->getNumInputPins(); pinIndex++)
    {
      if( childModule->inputPins[pinIndex].sourceModule == container )
      {
        // child receives input from the container - find out the inner source inside the container and re-connect the child directly to
        // this inner source:
        int containerOutputPinIndex = childModule->inputPins[pinIndex].outputIndex;
        Module *outputModule        = container->getAudioOutputModule(containerOutputPinIndex);
        Module *innerSourceModule   = outputModule->inputPins[0].sourceModule;
        if( !innerSourceModule->isInputModule() )
        {
          int innerSourceOutputPinIndex = outputModule->inputPins[0].outputIndex;
          childModule->connectInputPinTo(pinIndex, innerSourceModule, innerSourceOutputPinIndex);
        }
        else
          childModule->disconnectInputPin(pinIndex); // should be done automatically when the container is deleted -> check this
      }
    }
  }

  // drag over the container's non I/O child modules to become child-modules of "this":
  std::vector<Module*> tmpChildModules = container->childModules; // we need an temporary array
  for(unsigned int childIndex = 0; childIndex < tmpChildModules.size(); childIndex++)
  {
    if( tmpChildModules[childIndex]->isInputModule() || tmpChildModules[childIndex]->isOutputModule() )
      continue;
    rosic::removeElementByValue(container->childModules, tmpChildModules[childIndex]);
    rosic::appendElement(       this->childModules,      tmpChildModules[childIndex]);
    tmpChildModules[childIndex]->parentModule = this;
  }

  deleteChildModule(container);
  sortChildModuleArray();        // because we've got a bunch of new children appended at the end
}

void ContainerModule::addAudioConnection(romos::Module *sourceModule, int outputIndex, romos::Module *targetModule, int inputIndex)
{
  targetModule->connectInputPinTo(inputIndex, sourceModule, outputIndex);
  updateHasDelayedConnectionsFlag();
}

void ContainerModule::addAudioConnection(AudioConnection *connectionToAdd)
{
  addAudioConnection(connectionToAdd->getSourceModule(), connectionToAdd->getSourceOutputIndex(),
                     connectionToAdd->getTargetModule(), connectionToAdd->getTargetInputIndex());
}

bool ContainerModule::deleteAudioConnection(Module* /*sourceModule*/, int /*outputIndex*/, Module* targetModule, int inputIndex)
{
  targetModule->disconnectInputPin(inputIndex);
  updateHasDelayedConnectionsFlag();
  return true; // make this void
}

void ContainerModule::deleteAudioConnection(AudioConnection connectionToDelete)
{
  connectionToDelete.getTargetModule()->disconnectInputPin(connectionToDelete.getTargetInputIndex());
  updateHasDelayedConnectionsFlag();
}

void ContainerModule::deleteAudioConnections(std::vector<AudioConnection> connectionsToDelete)
{
  for(unsigned int i = 0; i < connectionsToDelete.size(); i++)
    deleteAudioConnection(connectionsToDelete[i]);
  updateHasDelayedConnectionsFlag();
}

void ContainerModule::minimizeNumberOfAudioInputs()
{
  int oldNumInputs = getNumInputPins();
  int newNumInputs = oldNumInputs;
  for(int i=0; i<newNumInputs; i++)
  {
    for(int j=i+1; j<newNumInputs; j++)
    {

      // If two input pins i and j carry the same signal, one of them (say, j) might be superfluous and all pins that receive a signal from
      // pin j can be recofigured to receive the same signal from pin i. Then, pin j can be deleted. But this works only if none of the
      // target-pins of these inputs receives input signals from both pins i and j, because in this case, a target pin receives the same
      // signal twice and reconfiguration and deletion of one of the inputs would leave only one of these two contributions due to
      // disallowed multiple connections between the same pins - between each pair of pins there's either none or one connection but not
      // many). We check the latter condition first as it is faster to evaluate:
      /*
      if( doOutputPinsHaveDisjointTargets(getAudioInputModule(i), 0, getAudioInputModule(j), 0) )
      {
        if( doInputPinsCarrySameSignal(this, i, this, j) )
        {
          // reconfigure all modules that are the target of input-module j to receive now their input from input-module i, and delete
          // input module j:

          // ...

          int dummy = 0;
        }
      }
      */


      // If all targets of pin i and j involve both pins as source (they all do the sum i+j), one of them (say, j) might be superfluous
      // because the summing could have already been done inside the pin. All incoming connections into pin j (outside the container) could
      // be re-connected to pin i (which sums the signal at pin j into pin i) and pin j could then be deleted. But this works only if pin i
      // and j have no common source-pin, because in this case, the signal of the common source pin would occur twice in the sum i+j (before
      // reconfiguration), but only once the new re-configured pin i.
      /*
      if( doOutputPinsHaveSameTargets(getAudioInputModule(i), 0, getAudioInputModule(j), 0) )
      {
        if( !doInputPinsHaveCommonSource(this, i, this, j) )
        {


          int dummy = 0;
        }
      }
      */

    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

romos::AudioInputModule* ContainerModule::getAudioInputModule(int index) const
{
  int inputId = moduleFactory.getModuleId("AudioInput");
  int numSkipped = 0;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    //if( childModules.at(i)->getTypeIdentifierOld() == romos::ModuleTypeRegistry::AUDIO_INPUT )
    if( childModules.at(i)->typeInfo->id == inputId )
    {
      if( index == numSkipped )
        return (AudioInputModule*) childModules.at(i);
      numSkipped++;
    }
  }
  return NULL;
}

romos::AudioOutputModule* ContainerModule::getAudioOutputModule(int index) const
{
  int numSkipped = 0;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    //if( childModules.at(i)->getTypeIdentifierOld() == romos::ModuleTypeRegistry::AUDIO_OUTPUT )
    if( childModules.at(i)->getTypeName() == "AudioOutput" ) // optimize, use id
    {
      if( index == numSkipped )
        return (AudioOutputModule*) childModules.at(i);
      numSkipped++;
    }
  }
  return NULL;
}

int ContainerModule::getInputPinIndexOf(AudioInputModule *inputModule) const
{
  int numSkipped = 0;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    //if( childModules.at(i)->getTypeIdentifierOld() == romos::ModuleTypeRegistry::AUDIO_INPUT )
    if( childModules.at(i)->getTypeName() == "AudioInput" )  // optimize
    {
      if( childModules.at(i) == inputModule )
        return numSkipped;
      numSkipped++;
    }
  }
  rassert(0);  // the passed module is none of our input modules
  return 0;
}

int ContainerModule::getOutputPinIndexOf(AudioOutputModule *outputModule) const
{
  int numSkipped = 0;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    //if( childModules.at(i)->getTypeIdentifierOld() == romos::ModuleTypeRegistry::AUDIO_OUTPUT )
    if( childModules.at(i)->getTypeName() == "AudioOutput" )  // optimize
    {
      if( childModules.at(i) == outputModule )
        return numSkipped;
      numSkipped++;
    }
  }
  rassert(0);  // the passed module is none of our output modules
  return 0;
}

romos::Module* ContainerModule::getChildModule(int index) const
{
  //rassert( index >= 0 && index < (int) childModules.size() );  // index out of range
  if( index >= 0 && index < (int) childModules.size() )
    return childModules.at(index);
  else
  {
    romos::triggerRuntimeError("Index out of range in ContainerModule::getChildModule");
    return NULL;
  }
}

bool ContainerModule::areAllChildModulesMonophonic() const
{
  for(unsigned int i = 0; i < childModules.size(); i++ )
  {
    if( childModules.at(i)->isPolyphonic() )
      return false;
  }
  return true;
}

bool ContainerModule::areAllChildModulesPolyphonic() const
{
  for(unsigned int i = 0; i < childModules.size(); i++ )
  {
    if( !childModules.at(i)->isPolyphonic() )
      return false;
  }
  return true;
}

int ContainerModule::getContainerNestingDepth() const
{
  int result = 0;
  int tmp    = 0;
  for(unsigned int i = 0; i<childModules.size(); i++)
  {
    //if( childModules.at(i)->getTypeIdentifierOld() == ModuleTypeRegistry::CONTAINER )
    if( childModules.at(i)->getTypeName() == "Container" )  // optimize
    {
      tmp = childModules.at(i)->getContainerNestingDepth();
      if( tmp > result )
        result = tmp;
    }
  }
  return result + 1;
}

int ContainerModule::getIndexOfChildModule(romos::Module *childToFindIndexFor) //const
{
  int index = rosic::findElement(childModules, childToFindIndexFor);
  //int index = childModules.findElement(childToFindIndexFor);
  rassert( index != -1 ); // trying to find out the index of a supposd-to-be-a-child that actually isn't
  return index;
}

std::vector<romos::Module*> ContainerModule::getNonInOutChildModules() const
{
  std::vector<romos::Module*> result;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( !childModules.at(i)->isInputModule() && !childModules.at(i)->isOutputModule() )
      rosic::appendElement(result, childModules.at(i));
      //result.appendElement(childModules.getElement(i));
  }
  return result;
}

/*
std::vector<romos::Module*> ContainerModule::getChildModulesWithTypeOld(int typeIdentifier) const
{
  // to be removed when the new type-info system fully works
  std::vector<romos::Module*> result;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules.at(i)->getTypeIdentifierOld() == typeIdentifier )
      rosic::appendElement(result, childModules.at(i));
      //result.appendElement(childModules.getElement(i));
  }
  return result;
}
*/

std::vector<romos::Module*> ContainerModule::getChildModulesWithTypeId(int typeId) const
{
  std::vector<romos::Module*> result;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules.at(i)->typeInfo->id == typeId )
      rosic::appendElement(result, childModules.at(i));
  }
  return result;
}

std::vector<romos::Module*> ContainerModule::getChildModulesWithType(const std::string& type) const
{
  int id = moduleFactory.getModuleId(type);
  return getChildModulesWithTypeId(id);
}

std::vector<romos::Module*> ContainerModule::getConnectedTargetModulesOf(const romos::Module* sourceModule) const
{
  std::vector<romos::Module*> result;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules.at(i)->hasIncomingConnectionFrom(sourceModule) )
      rosic::appendElement(result, childModules.at(i));
  }
  return result;
}

std::vector<romos::Module*> ContainerModule::getConnectedTargetModulesOf(const romos::Module* sourceModule, int outputPinIndex) const
{
  std::vector<romos::Module*> result;
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules.at(i)->hasIncomingConnectionFrom(sourceModule, outputPinIndex) )
      rosic::appendElement(result, childModules.at(i));
  }
  return result;
}

rosic::rsString ContainerModule::getPinName(int kind, int direction, int pinIndex) const
{
  if( kind == AUDIO )
  {
    if( direction == INCOMING )
      return getAudioInputModule(pinIndex)->getName();
    else if( direction == OUTGOING )
      return getAudioOutputModule(pinIndex)->getName();
  }
  return rosic::rsString("");
}

/*
void ContainerModule::getPinDataForModule(romos::Module *module, int &kind, int &direction, int &pinIndex) const
{
  kind = direction = pinIndex = -1;
   // something to do...
}
*/

/*
bool ContainerModule::containsConnectionsWithImplicitDelay() const
{
  rosic::Array<AudioConnection*> connections;
  for(int i=0; i<childModules.getNumElements(); i++)
  {
    connections = childModules.getElement(i)->getIncomingAudioConnections();
    for(int j=0; j<connections.getNumElements(); j++)
    {
      if( connections[j]->hasImplicitDelay() )
        return true;
    }
  }
  return false;
}
*/

bool ContainerModule::isPositionOccupied(int &x, int &y) const
{
  for(unsigned int i = 0; i < childModules.size(); i++)
  {
    if( childModules.at(i)->getPositionX() == x && childModules.at(i)->getPositionY() == y )
      return true;
  }
  return false;
}

void ContainerModule::getNonOccupiedPositionNear(int &x, int &y) const
{
  while(true)
  {
    if( isPositionOccupied(x, y) )
      x++;  // maybe we can do something more elaborate here
    else
      return;
  }
}

/*
//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void ContainerModule::moduleMoved(Module* moduleThatHasMoved)
{
  //moduleMutex.lock();
  sortChildModuleArray();
  //moduleMutex.unlock();
}

void ContainerModule::numAudioInputsChanged(Module* moduleThatHasChanged)
{
  // \todo: possibly remove some audio-connections
  int dummy = 0;
}

void ContainerModule::numAudioOutputsChanged(Module* moduleThatHasChanged)
{
  // \todo: possibly remove some audio-connections
  int dummy = 0;
}

void ContainerModule::numEventInputsChanged(Module* moduleThatHasChanged)
{
  // \todo: possibly remove some event-connections
}

void ContainerModule::numEventOutputsChanged(Module* moduleThatHasChanged)
{
  // \todo: possibly remove some event-connections
}
*/

void ContainerModule::childPolyphonyChanged(Module *childThatHasChangedPolyphony)
{
  for(unsigned int i = 0; i < childModules.size(); i++)
    childModules[i]->updateInputPointersAndInFrameStrides();


  // if the child that has changed its polyphony is connected to the outside, we must inform the outlying container too, so
  // it can update itself as well:
  if( childThatHasChangedPolyphony->isConnectedToAudioOutput() )
  {
    if( parentModule != NULL )
      parentModule->childPolyphonyChanged(this);
       // later, when the I/O modules are derived from PointerRedirectModule, pass the "childThatHasChangedPolyphony" instead of "this"
  }
  // obsolete? relevant only with Proxy I/O modules?

  assignProcessingFunctions();
}

void ContainerModule::resetVoiceState(int voiceIndex)
{
  for(unsigned int i = 0; i < childModules.size(); i++)
    childModules[i]->resetVoiceState(voiceIndex);

  //int dummy = 0;
}

void romos::ContainerModule::freeMemory()
{
  delete[] audioOutputs;
  delete[] tmpOutFrame;
  audioOutputs = NULL;
  tmpOutFrame  = NULL;

  //DEBUG_BREAK; // check code above


  /*
  // needs to be updated:
  unsigned int i;
  for(i = 0; i < getNumInputPins(); i++)
  {
    getAudioInputModule(i)->setAudioInputAddress(NULL);
    getAudioInputModule(i)->setAudioOutputAddress(NULL);
    getAudioInputModule(i)->numInputs  = numInputs;
    getAudioInputModule(i)->numOutputs = numInputs;
    // silly to call getAudioInputModule(i) repeatedly - clean this up...
  }
  for(i = 0; i < numOutputs; i++)
  {
    getAudioOutputModule(i)->setAudioInputAddress(NULL);
    getAudioOutputModule(i)->setAudioOutputAddress(NULL);
    getAudioOutputModule(i)->numInputs  = numOutputs;
    getAudioOutputModule(i)->numOutputs = numOutputs;
  }
  */
}

/*
void romos::ContainerModule::allocateAudioInputs()
{
  romos::Module::allocateAudioInputs();
  setupPointersInInputModules(); // check if still necessary...yep - seems so
}
*/

void romos::ContainerModule::updateInputPointersAndInFrameStrides()
{
  romos::Module::updateInputPointersAndInFrameStrides();
  updatePointersInInputModules();
}

void romos::ContainerModule::updatePointersInInputModules()
{
  //DEBUG_BREAK;  // check the below
  for(unsigned int i = 0; i < inputPins.size(); i++)
  {
    AudioInputModule *inModule = getAudioInputModule(i);

    //inModule->inputPins[0] = inputPins[i];

    if(inModule != nullptr)
      inModule->inputPins[0] = inputPins[i];
      // copies all the pin-data from this container's i-th input pin into the input-pin of the
      // i-th input module
    else
    {
      // hmm...what else? in what sort of situation may this happen?
    }



  }


  /*
  // old:
  for(unsigned int i = 0; i < numInputs; i++)
  {
    AudioInputModule *inModule  = getAudioInputModule(i);
    inModule->audioInputs[0]    = audioInputs[i];
    inModule->inFrameStrides[0] = inFrameStrides[i];
  }
  */
}





int romos::ContainerModule::getRequiredOutputBufferSize() const
{
  return getNumOutputPins() * getRequiredOutputBufferSizePerPin();
}

int romos::ContainerModule::getRequiredOutputBufferSizePerPin() const
{
  return getNumVoices() * processingStatus.getBufferSize();
}

void romos::ContainerModule::allocateAudioOutputs()
{
  // maybe take out commonalities with baseclass method and invoke it...

  delete[] audioOutputs;
  audioOutputs = NULL;

  delete[] tmpOutFrame;
  tmpOutFrame = NULL;

  int sizeToAllocate = getRequiredOutputBufferSize();
  if( sizeToAllocate > 0)
  {
    audioOutputs = new double[sizeToAllocate];
    memset(audioOutputs, 0, sizeToAllocate*sizeof(double));

    tmpOutFrame = new double[getNumVoices() * getNumOutputPins()];
    memset(tmpOutFrame, 0, getNumVoices() * getNumOutputPins() * sizeof(double));
  }



  if( parentModule != NULL )
    parentModule->outputsWereReAllocated(this);

  updatePointersInOutputModules();
}

void romos::ContainerModule::updatePointersInOutputModules()
{
  for(unsigned int i = 0; i < getNumOutputPins(); i++)
  {
    AudioOutputModule *outModule = getAudioOutputModule(i);
    outModule->audioOutputs      = audioOutputs + i * getRequiredOutputBufferSizePerPin();
    outModule->outFrameStride    = 1;
  }
}

void romos::ContainerModule::assignProcessingFunctions()
{
  if( isPolyphonic() && areAllChildModulesPolyphonic() )
  {
    processFrame = &processContainerAllPoly;
    processBlock = &processContainerAllPolyBlock;
  }
  else if( !isPolyphonic() && areAllChildModulesMonophonic() )
  {
    processFrame = &processContainerAllMono;
    processBlock = &processContainerAllMonoBlock;
  }
  else
  {
    processFrame = &processContainerMixedMonoPoly;
    processBlock = &processContainerMixedMonoPolyBlock;
  }
}

void romos::ContainerModule::outputsWereReAllocated(Module* /*moduleThatHasReAllocated*/)
{
  for(unsigned int i = 0; i < childModules.size(); i++)
    childModules[i]->updateInputPointersAndInFrameStrides();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// internal functions:

bool ContainerModule::hasAsDirectlyEmbeddedModule(romos::Module *moduleToSearchFor)
{
  return rosic::containsElement(childModules, moduleToSearchFor);
  //return childModules.hasElement(moduleToSearchFor);
}

bool ContainerModule::hasAsDescendant(Module *moduleToSearchFor)
{
  if( hasAsDirectlyEmbeddedModule(moduleToSearchFor) )
    return true;
  else
  {
    for(unsigned int i=0; i<childModules.size(); i++)
    {
      ContainerModule *mc = dynamic_cast<ContainerModule*> (childModules[i]);
      if( mc != NULL && mc->hasAsDescendant(moduleToSearchFor) )
        return true;
    }
  }
  return false;
}

void ContainerModule::updateHasDelayedConnectionsFlag()
{
  for(unsigned int childIndex = getNumInputPins(); childIndex < childModules.size() - getNumOutputPins(); childIndex++)
  {
    if( childModules.at(childIndex)->hasDelayedIncomingConnection() )
    {
      hasDelayedConnections = true;
      return;
    }
  }
  hasDelayedConnections = false;
}

void ContainerModule::sortChildModuleArray()
{
  // retrieve some stuff before sorting which is needed to re-connect the container's pins in case the order of I/O modules changes:
  //std::vector<Module*> oldInputModules  = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_INPUT);
  //std::vector<Module*> oldOutputModules = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_OUTPUT);
  std::vector<Module*> oldInputModules  = getChildModulesWithType("AudioInput");
  std::vector<Module*> oldOutputModules = getChildModulesWithType("AudioOutput");

  std::vector<std::vector<Module*> > oldTargetModuleArrays; // one array for each output pin
  for(unsigned int outIndex = 0; outIndex < getNumOutputPins(); outIndex++)
    rosic::appendElement(oldTargetModuleArrays, getConnectedTargetModulesOfPin(outIndex));

  // the actual sorting:
  sortModuleArrayByCoordinates(childModules);

  // re-connect input pins, if necessary:
  //std::vector<Module*> newInputModules  = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_INPUT);
  std::vector<Module*> newInputModules  = getChildModulesWithType("AudioInput");
  for(unsigned int newIndex = 0; newIndex < newInputModules.size(); newIndex++)
  {
    unsigned int oldIndex = rosic::findElement(oldInputModules, newInputModules[newIndex]);
    if(oldIndex != newIndex)
    {
      Module *sourceModule     = oldInputModules[oldIndex]->inputPins[0].sourceModule;
      int    sourceOutPinIndex = oldInputModules[oldIndex]->inputPins[0].outputIndex;
      connectInputPinTo(newIndex, sourceModule, sourceOutPinIndex);
    }
  }

  // re-connect output pins, if necessary:
  //std::vector<Module*> newOutputModules = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_OUTPUT);
  std::vector<Module*> newOutputModules = getChildModulesWithType("AudioOutput");
  for(unsigned int newIndex = 0; newIndex < newOutputModules.size(); newIndex++)
  {
    unsigned int oldIndex = rosic::findElement(oldOutputModules, newOutputModules[newIndex]);
    if(oldIndex != newIndex)
    {
      std::vector<Module*> oldTargetsOfPin = oldTargetModuleArrays[oldIndex];
      for(unsigned int targetIndex = 0; targetIndex < oldTargetsOfPin.size(); targetIndex++)
      {
        Module *targetModule = oldTargetsOfPin[targetIndex];
        for(unsigned int pinIndex = 0; pinIndex < targetModule->inputPins.size(); pinIndex++)
        {
          if( targetModule->inputPins[pinIndex].outputIndex == oldIndex )
            targetModule->connectInputPinTo(pinIndex, this, newIndex);
          //int dummy = 0;
        }
      }
    }
  }

  // hmmm...why this?:
  updatePointersInInputModules();
  updatePointersInOutputModules();

  updateHasDelayedConnectionsFlag();
}


void ContainerModule::mapApparentSourceToProcessingSource(Module * &sourceModule, int &sourceOutputPinIndex)
{
  if( sourceModule != this )
    triggerRuntimeError("sourceModule != this in ContainerModule::mapApparentSourceToProcessingSource");

  //std::vector<romos::Module*> outputModules = getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_OUTPUT);
  std::vector<romos::Module*> outputModules = getChildModulesWithType("AudioOutput");
  if( sourceOutputPinIndex >= 0 && sourceOutputPinIndex < (int) outputModules.size() )
  {
    sourceModule         = outputModules[sourceOutputPinIndex];
    sourceOutputPinIndex = 0;
    sourceModule->mapApparentSourceToProcessingSource(sourceModule, sourceOutputPinIndex);
  }
  else
  {
    triggerRuntimeError("pin-index invalid in ContainerModule::mapApparentSourceToProcessingSource");
    sourceModule         = NULL;  // maybe have a global default module to connect to such that realease code doesn't crash
    sourceOutputPinIndex = 0;
  }
}

void ContainerModule::mapProcessingSourceToSourceApparent(Module*& /*sourceModule*/, int& /*sourceOutputPinIndex*/)
{
  DEBUG_BREAK; // not yet implemented
}



void ContainerModule::sortModuleArrayByCoordinates(std::vector<romos::Module*> &modulesToSort)
{
  std::sort(modulesToSort.begin(), modulesToSort.end(), modulePointerLessByXY);

  //std::vector<romos::Module*> modulesVector;
  //modulesToSort.toVectorSTL(modulesVector);
  //std::sort(modulesVector.begin(), modulesVector.end(), modulePointerLessByXY);
  //modulesToSort.fromVectorSTL(modulesVector);
}

void ContainerModule::removeModulesOfType(std::vector<romos::Module*> &modules, const std::string& typeName)
{
  int i = (int) modules.size()-1;
  while( i >= 0)
  {
    if( modules[i]->getTypeName() == typeName ) // optimize: us id
      rosic::removeElementByIndex(modules, i);
    i--;
  }
}
/*
void ContainerModule::removeModulesOfType(std::vector<romos::Module*> &modules, int typeCodeToRemove)
{
  int i = (int) modules.size()-1;
  while( i >= 0)
  {
    if( modules[i]->getTypeIdentifierOld() == typeCodeToRemove )
      rosic::removeElementByIndex(modules, i);
    i--;
  }
}
*/

