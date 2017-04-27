#include "rosic_Module.h"
using namespace rosic;

//=================================================================================================
// class ModulatableParameter:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulatableParameter::ModulatableParameter(const char* name, double initialValue)
{
  nominalValue = instantaneousValue = initialValue;
  if( name != NULL )
  {
    int length = (int) strlen(name);
    this->name = new char[length+1];
    for(int c=0; c<=length; c++) // the <= is valid here, because we have one more cell allocated
      this->name[c] = name[c];
  }
  else
    this->name = "unknown parameter";  // maybe we need strcpy
}


//=================================================================================================
// class ModulationRouter:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulationRouter::ModulationRouter()
{

}

ModulationRouter::~ModulationRouter()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void ModulationRouter::establishNewConnection(ModulationSource *sourceToConnect,
                                              ModulatableParameter *parameterToConnect,
                                              double connectionStrength)
{
  mutex.lock();
  ModulationConnection newConnection(sourceToConnect, parameterToConnect, connectionStrength);
  connections.appendElement(newConnection);
  mutex.unlock();
}

void ModulationRouter::removeConnection(int index)
{
  mutex.lock();
  if( index < 0 || index >= connections.getNumElements() )
  {
    mutex.unlock();
    return;
  }
  connections[index].parameter->initInstantaneousValue(); // to leave it 'clean'
  connections.removeElementByIndex(index);
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulationRouter::applyModulations(bool initWithNominalValue)
{
  mutex.lock();

  // this can't be dragged into the accumulation loop because that would cause parameters that
  // appear in more than one connection to be re-initialized in each iteration:
  if( initWithNominalValue == true )
  {
    for(int i=0; i<connections.getNumElements(); i++)
      connections[i].parameter->initInstantaneousValue();
  }

  // the accumulation loop:
  for(int i=0; i<connections.getNumElements(); i++)
  {
    connections[i].parameter->addModulationSignal(
      connections[i].source->getSample() * connections[i].strength);
  }

  mutex.unlock();
}


//=================================================================================================
// class Module:

char* Module::getModulatableParameterName(int index)
{
  if( index >= 0 && index < modulatableParameters.getNumElements() )
    return modulatableParameters[index].getName();
  else
  {
    DEBUG_BREAK;  // parameter-index out of range
    return NULL;
  }
}






