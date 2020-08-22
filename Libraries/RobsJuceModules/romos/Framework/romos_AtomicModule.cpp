//#include "romos_AtomicModule.h"
//#include "romos_AudioConnection.h" 
//#include "romos_ContainerModule.h"
//using namespace romos;

//=========================================================================================================================================
// AtomicModule:

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

AtomicModule::AtomicModule(const std::string& name, int x, int y, bool polyphonic) 
: Module(name, x, y, polyphonic)
{

}

AtomicModule::~AtomicModule()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------    
// setup of pins:
/*
// this old code is obsolete now
void AtomicModule::initInputPins(int numberOfPins, const char*, ...)
{
  va_list ap;
  va_start(ap, numberOfPins);
  for(int i = 1; i <= numberOfPins; i++)
  {
    rosic::appendElement(audioInputNames, rosic::rsString(va_arg(ap, const char*)));
    rosic::appendElement(inputPins,       AudioInputPinData());
  }
  numInputs += numberOfPins;
  va_end(ap);
  // this function crashes - but only in a release build? am i using the var_arg mechanism wrong?
  // ...anyway - maybe best to get rid of that function entirely and use the one below:
}
*/
void AtomicModule::initInputPins(const std::vector<std::string>& pinNames)
{
  for(size_t i = 0; i < pinNames.size(); i++)
  {
    rosic::appendElement(audioInputNames, rosic::rsString(pinNames[i]));
    rosic::appendElement(inputPins,       AudioInputPinData());
  }
  numInputs = (int) inputPins.size();
}
/*
void AtomicModule::initOutputPins(int numberOfPins, const char*, ...)
{
  va_list ap;
  va_start(ap, numberOfPins);
  for(int i=1; i<=numberOfPins; i++)
    rosic::appendElement(audioOutputNames, rosic::rsString(va_arg(ap, const char*)));
  outFrameStride += numberOfPins;
  va_end(ap);
}
*/
void AtomicModule::initOutputPins(const std::vector<std::string>& pinNames)
{
  for(size_t i = 0; i < pinNames.size(); i++)
    rosic::appendElement(audioOutputNames, rosic::rsString(pinNames[i]));
  outFrameStride += (int) pinNames.size();
}


/*

// maybe reactivate later when we implement variable I/O:

void romos::AtomicModule::setPinName(int kind, int direction, int pinIndex, const rosic::rsString &newName)
{
  rosic::Array<rosic::rsString>* pinNames = properties->getPinNameArray(kind, direction);
  if( pinIndex >= 0 && pinIndex < pinNames->getNumElements() )
    pinNames->replaceElement(pinIndex, newName);
  else
    DEBUG_BREAK;
}

void romos::AtomicModule::setNumAudioInputs(int newNumber) 
{ 
  //moduleMutex.lock();
  numAudioInputs = newNumber; 
  sendNumAudioInputsChangedNotification(); 
  //moduleMutex.unlock();
}
 
void romos::AtomicModule::setNumAudioOutputs(int newNumber) 
{ 
  //moduleMutex.lock();
  numAudioOutputs = newNumber; 
  sendNumAudioOutputsChangedNotification(); 
  //moduleMutex.unlock();
}
*/


//-----------------------------------------------------------------------------------------------------------------------------------------    
// inquiry about pins:

rosic::rsString AtomicModule::getPinName(int kind, int direction, int pinIndex) const
{
  if( direction == romos::INCOMING )
    return audioInputNames.at(pinIndex);
    //return audioInputNames.getElement(pinIndex);
  else
    return audioOutputNames.at(pinIndex);
    //return audioOutputNames.getElement(pinIndex);  // get rid of the dispatch - provide separate functions instead
}


//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void AtomicModule::addAudioInput(const char* shortName, const char* longName, 
  const char* description)
{
  rosic::appendElement(audioInputNames,        rosic::rsString(shortName));
  //rosic::appendElement(audioInputLongNames,    rosic::rsString(longName));
  //rosic::appendElement(audioInputDescriptions, rosic::rsString(description));
  rosic::appendElement(inputPins, AudioInputPinData());
  numInputs++;
  updateInputPointersAndInFrameStrides();
}

void AtomicModule::addAudioOutput(const char* pinName)
{
  rosic::appendElement(audioOutputNames, rosic::rsString(pinName));



  outFrameStride++;
  allocateAudioOutputs();
}

void romos::AtomicModule::deleteAudioInput(int index)
{
  rosic::removeElementByIndex(audioInputNames, index);
  rosic::removeElementByIndex(inputPins,       index);
  numInputs--;
  updateInputPointersAndInFrameStrides(); // not needed - the remaining stay the same?
  //allocateAudioInputs();  


  /*
  // old:
  // delete outside connections to the to-be-deleted input:
  deleteIncomingAudioConnectionsToPin(index); 

  // update pin-indices in the remaining outside connections (decrement all that come after the to-be-deleted pin):
  for(unsigned int i = 0; i < incomingAudioConnections.size(); i++)
  {
    romos::AudioConnection *c = incomingAudioConnections[i];
    int oldIndex = c->getTargetInputIndex();
    if( oldIndex > index ) 
      c->setTargetInputIndex(oldIndex-1);
  }

  //audioInputNames.removeElementByIndex(index);
  rosic::removeElementByIndex(audioInputNames, index);
  numInputs--;
  allocateAudioInputs();
  */
}

void AtomicModule::deleteAudioOutput(int index)
{
  DEBUG_BREAK;

  /*
  // needs to be updated:
  deleteOutgoingAudioConnectionsFromPin(index);  
  std::vector<AudioConnection*> outgoingAudioConnections = getOutgoingAudioConnections();
  for(unsigned int i = 0; i < outgoingAudioConnections.size(); i++)
  {
    romos::AudioConnection *c = outgoingAudioConnections[i];
    int oldIndex = c->getSourceOutputIndex();
    if( oldIndex > index ) 
      c->setSourceOutputIndex(oldIndex-1);
  }
  //audioOutputNames.removeElementByIndex(index);
  rosic::removeElementByIndex(audioOutputNames, index);
  numOutputs--;
  allocateAudioOutputs();
  */
}







//=========================================================================================================================================
// ParameterMixIn:

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

bool ParameterMixIn::setParameter(const rosic::rsString &parameterName, const rosic::rsString &newValue, bool callInternalCallback)
{
  int index = findIndexOfParameterWithName(parameterName);
  if( index != -1 )
  {
    setParameter(index, newValue, callInternalCallback);
    return true;
  }
  else
    return false;
}

void ParameterMixIn::setParameter(int index, const rosic::rsString &newValue, bool callInternalCallback)
{
  rassert( index >= 0 && index < getNumParameters() );
  parameters[index].value = newValue;  
  if( callInternalCallback == true )
    parameterChanged(index);
}

void ParameterMixIn::addParameter(const rosic::rsString &parameterName, const rosic::rsString &defaultValue)
{
  Parameter p;
  p.name         = parameterName;
  p.defaultValue = defaultValue;
  p.value        = defaultValue;
  rosic::appendElement(parameters, p);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

rosic::rsString ParameterMixIn::getParameterName(int index) const
{
  rassert( index >= 0 && index < getNumParameters() );
  return parameters[index].name;
}
    
rosic::rsString ParameterMixIn::getParameterValue(int index) const
{
  rassert( index >= 0 && index < getNumParameters() );
  return parameters[index].value;
}

rosic::rsString ParameterMixIn::getParameterValue(const rosic::rsString &parameterName) const
{
  int index = findIndexOfParameterWithName(parameterName);
  if( index != -1 )
    return getParameterValue(index);
  else
    return 0.0;
}

rosic::rsString ParameterMixIn::getParameterDefaultValue(int index) const
{
  rassert( index >= 0 && index < getNumParameters() );
  return parameters[index].defaultValue;
}

int ParameterMixIn::findIndexOfParameterWithName(const rosic::rsString &nameToFind) const
{
  for(int i = 0; i < (int) parameters.size(); i++)
  {
    if( parameters[i].name == nameToFind ) 
      return i;
  }
  return -1;
}














//=========================================================================================================================================
// ModuleProxy:


//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ModuleProxy::ModuleProxy(const std::string& name, int x, int y, bool polyphonic) 
: AtomicModule(name, x, y, polyphonic)
{

}

ModuleProxy::~ModuleProxy()
{

}

void ModuleProxy::initialize()
{ 
  initInputPins({ "" });
  initOutputPins({ "" });
  hasHeaderFlag = false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void ModuleProxy::mapApparentSourceToProcessingSource(Module * &sourceModule, int &sourceOutputPinIndex)
{
  if( sourceModule != this )
    triggerRuntimeError("sourceModule != this in ModuleProxy::mapApparentSourceToProcessingSource");

  if( inputPins[0].sourceModule != NULL ) 
  {
    sourceModule         = inputPins[0].sourceModule;
    sourceOutputPinIndex = inputPins[0].outputIndex;
    sourceModule->mapApparentSourceToProcessingSource(sourceModule, sourceOutputPinIndex);
    int dummy = 0;
  }
  else
  {
    // leave pointer and index as is
  }
}
 
void ModuleProxy::mapProcessingSourceToSourceApparent(Module * &sourceModule, int &sourceOutputPinIndex)
{
  DEBUG_BREAK; // not yet implemented
}

void ModuleProxy::assignProcessingFunctions()
{
  processFrame = & ModuleProxy::processFrameDummy;    
  processBlock = & ModuleProxy::processBlockDummy;    
}
