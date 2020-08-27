//#include "romos_AudioConnection.h"
//#include "romos_ContainerModule.h"
//using namespace romos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:
    
AudioConnection::AudioConnection(romos::Module *sourceModule, int outputIndex, romos::Module *targetModule, int inputIndex)
{
  this->sourceModule = sourceModule;
  this->targetModule = targetModule;
  outIndex           = outputIndex;
  inIndex            = inputIndex;  
  //updateSourcePointer();
  //updateTargetPointer();
}

AudioConnection::~AudioConnection()
{
  //if( targetModule != NULL )
  //  targetModule->removeIncomingAudioConnection(this);
  // nah, the conneciton should not call the source/target modules because the destructor of the target module deletes its connections
  //
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

/*
void AudioConnection::updateSourcePointer() 
{ 
  if( sourceModule != NULL )
    sourcePointer = sourceModule->getAudioOutputAddress() + outIndex; 
  else
    sourcePointer = NULL;
}

void AudioConnection::updateTargetPointer() 
{ 
  if( targetModule != NULL )
    targetPointer = targetModule->getAudioInputAddress() + inIndex; 
  else
    targetPointer = NULL;
}
*/

void AudioConnection::resetToNull()
{
  sourceModule = nullptr;
  targetModule = nullptr;
  outIndex = 0;
  inIndex  = 0;  
}
    
void AudioConnection::setSourceModule(romos::Module *newSourceModule)
{
  sourceModule = newSourceModule;
  //updateSourcePointer();
}

void AudioConnection::setTargetModule(romos::Module *newTargetModule)
{
  targetModule = newTargetModule;
  //updateTargetPointer();
}
    
void AudioConnection::setSourceOutputIndex(int newIndex)
{
  outIndex = newIndex;
  //updateSourcePointer();
}

void AudioConnection::setTargetInputIndex(int newIndex)
{
  inIndex = newIndex;
  //updateTargetPointer();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

/*
romos::Module* AudioConnection::getApparentSourceModule() const 
{
  if( sourceModule->getTypeIdentifier() != ModuleTypeRegistry::AUDIO_OUTPUT )
    return sourceModule;
  else
  {
    rassert( sourceModule->getParentModule() != NULL );
    return (romos::Module*) sourceModule->getParentModule();
  }
}

romos::Module* AudioConnection::getApparentTargetModule() const
{
  if( targetModule->getTypeIdentifier() != ModuleTypeRegistry::AUDIO_INPUT )
    return targetModule;
  else
  {
    rassert( targetModule->getParentModule() != NULL );
    return (romos::Module*) targetModule->getParentModule();
  }
}

int AudioConnection::getApparentTargetInputIndex() const
{
  if( targetModule->getTypeIdentifier() != ModuleTypeRegistry::CONTAINER )
    return outIndex;
  else
  {
    rassert( dynamic_cast<romos::AudioInputModule*> (targetModule) != NULL );
    return ((romos::ContainerModule*) targetModule)->getInputPinIndexOf( (romos::AudioInputModule*) targetModule );
  } 
}

int AudioConnection::getApparentSourceOutputIndex() const
{
  if( sourceModule->getTypeIdentifier() != ModuleTypeRegistry::CONTAINER )
    return outIndex;
  else
  {
    rassert( dynamic_cast<romos::AudioOutputModule*> (sourceModule) != NULL );
    return ((romos::ContainerModule*) sourceModule)->getOutputPinIndexOf( (romos::AudioOutputModule*) sourceModule );
  } 
}
*/

bool AudioConnection::hasSameSourcePinAs(AudioConnection *otherConnection)
{
  if( sourceModule == otherConnection->sourceModule && outIndex == otherConnection->outIndex )
    return true;
  else
    return false;
}

bool AudioConnection::hasSameTargetPinAs(AudioConnection *otherConnection)
{
  if( targetModule == otherConnection->targetModule && inIndex == otherConnection->inIndex )
    return true;
  else
    return false;
}

bool AudioConnection::connectsSamePinsAs(AudioConnection *otherConnection)
{
  if( hasSameSourcePinAs(otherConnection) && hasSameTargetPinAs(otherConnection) )
    return true;
  else
    return false;
}

bool AudioConnection::hasImplicitDelay() const
{ 
  return (getTargetModule() == getSourceModule()) || modulePointerLessByXY(getTargetModule(), getSourceModule()); 
}

bool AudioConnection::areConnectionsAlike(AudioConnection *connection1, AudioConnection *connection2, int criterion)
{
  if( criterion == INPUT_PINS_EQUAL )
    return connection1->hasSameSourcePinAs(connection2);
  else if( criterion == OUTPUT_PINS_EQUAL )
    return connection1->hasSameTargetPinAs(connection2);
  else if( criterion == ALL_PINS_EQUAL )
    return connection1->connectsSamePinsAs(connection2);
  else
  {
    triggerRuntimeError("Unknown comparison criterion in AudioConnection::areConnectionsAlike");  
    return false;
  }
}

bool AudioConnection::areConnectionArraysAlike(std::vector<AudioConnection*> connections1,
                                               std::vector<AudioConnection*> connections2, 
                                               int criterion)
{
  if( connections1.size() != connections2.size() )
    return false;

  // For each connection in the first array, find a matching connection (one that is "alike") in the second array. If a match is 
  // found, remove the pair of connections from their respective arrays and investigate the next connection in the first array. If, at the 
  // end of the process, both arrays are empty (i.e. numRemainingConnections == 0 ), then we have found for each connection in array 1 a 
  // match in array 2 and vice versa. We need this nested loop because we can't assume any ordering of the connections in the two arrays. 

  int numRemainingConnections = (int) connections1.size();
  for(int i = 0; i < numRemainingConnections; i++)
  {
    for(int j = 0; j < numRemainingConnections; j++)
    {
      if( areConnectionsAlike(connections1[i], connections2[j], criterion) )
      {
        rosic::removeElementByIndex(connections1, i);
        rosic::removeElementByIndex(connections2, j);
        //connections1.removeElementByIndex(i);
        //connections2.removeElementByIndex(j);
        i--;
        numRemainingConnections--;
        break; // jump out of the inner loop
      }
    }
  }

  if( numRemainingConnections == 0 )
    return true;
  else
    return false;
}
