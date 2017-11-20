//#include "rojue_AutomatableModule.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AutomatableModule::AutomatableModule()
{

}

AutomatableModule::~AutomatableModule()
{
  removeAllObservedParameters(true);
}
    
//-------------------------------------------------------------------------------------------------
// MIDI controller stuff:

void AutomatableModule::assignMidiController(const String& nameOfParameter, int controllerNumber)
{
  observedParameters.getLock().enter();
  Parameter *p;
  p = getParameterByName(nameOfParameter);
  if( p != NULL )
  {
    AutomatableParameter* ap = dynamic_cast<AutomatableParameter*> (p);
    if( ap != NULL )
      ap->assignMidiController(controllerNumber);
  }
  observedParameters.getLock().exit();
}

void AutomatableModule::setMidiController(int controllerNumber, int controllerValue)
{
  // loop through all the observed parameters and pass the controller value to them - the 
  // parameters themselves will take care to respond only to controller-numbers which are assigned
  // to them:
  observedParameters.getLock().enter();
  Parameter            *p;
  AutomatableParameter *ap;
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    p  = observedParameters[i];
    ap = dynamic_cast<AutomatableParameter*> (p);
    if( ap != NULL )
      ap->setMidiController(controllerNumber, controllerValue);
    //observedParameters[i]->setMidiController(controllerNumber, controllerValue);
  }
  observedParameters.getLock().exit();
}

void AutomatableModule::revertToDefaultMapping()
{
  observedParameters.getLock().enter();
  Parameter            *p;
  AutomatableParameter *ap;
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    p  = observedParameters[i];
    ap = dynamic_cast<AutomatableParameter*> (p);
    if( ap != NULL )
      ap->revertToDefaults(false, false);
    //observedParameters[i]->revertToDefaults();
  }
  observedParameters.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// retrieve pointers to the observed parameters:

Parameter* AutomatableModule::getParameterByName(const String& nameOfParameter) const
{
  Parameter* result = NULL;

  observedParameters.getLock().enter();
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    if( observedParameters[i]->hasName(nameOfParameter) )
      result = observedParameters[i];
  }
  observedParameters.getLock().exit();

  return result;
}

Parameter* AutomatableModule::getParameterByIndex(int indexOfParameter) const
{
  Parameter* result = NULL;

  observedParameters.getLock().enter();
  if( indexOfParameter < (int) observedParameters.size() )
    result = observedParameters[indexOfParameter];
  observedParameters.getLock().exit();

  return result;
}

int AutomatableModule::getIndexOfParameter(Parameter* parameterToRetrieveIndexOf) const
{
  int parameterIndex = -1;
  observedParameters.getLock().enter();
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    if( parameterToRetrieveIndexOf == observedParameters[i] )
      parameterIndex = i;
  }  
  observedParameters.getLock().exit();
  return parameterIndex;
}

int AutomatableModule::getNumParameters() const
{
  observedParameters.getLock().enter();
  int result = observedParameters.size();
  observedParameters.getLock().exit();

  return result;
}

//-------------------------------------------------------------------------------------------------
// add/remove observed parameters:

void AutomatableModule::addObservedParameter(Parameter *parameterToAdd, 
  void (AutomatableModule::*handlerFunction)(double value)  )
{
  observedParameters.getLock().enter();
  //observedParameters.push_back(parameterToAdd);
  observedParameters.addIfNotAlreadyThere(parameterToAdd);
  parameterToAdd->registerParameterObserver(this);
  //handlerFunctions.push_back(handlerFunction);
  observedParameters.getLock().exit();
}

void AutomatableModule::removeObservedParameter(Parameter *parameterToRemove, bool deleteObject)
{
  observedParameters.getLock().enter();
  int i=0;
  while( i < (int) observedParameters.size() )
  {
    if( observedParameters[i] == parameterToRemove )
    {
      // de-register ourselves as listener:
      parameterToRemove->deRegisterParameterObserver(this);

      // remove the pointer from the array:
      observedParameters.removeFirstMatchingValue(parameterToRemove);
      //observedParameters.removeValue(parameterToRemove);
      //observedParameters.erase(observedParameters.begin() + i);

      // remove the handler function from the function array:
      //handlerFunctions.erase(handlerFunctions.begin() + i);

      // free the memory associated with the pointer, if so chosen:
      if( deleteObject == true )
        delete parameterToRemove;

      
      i--; // because array has shrunken
    }
    i++;
  }
  observedParameters.getLock().exit();
}

void AutomatableModule::removeAllObservedParameters(bool deleteObjects)
{
  observedParameters.getLock().enter();
  Parameter *removee; // this is the currently removed parameter
  while( observedParameters.size() > 0 )
  {
    // retrieve the pointer to the parameter to be removed from the array:
    removee = observedParameters[0];

    // de-register ourselves as listener:
    removee->deRegisterParameterObserver(this);

    // remove the parameter from the array:
    observedParameters.remove(0);
    //observedParameters.erase(observedParameters.begin());

    // remove the handler function from the function array:    
    //handlerFunctions.erase(handlerFunctions.begin());

    // optionally delete the parameter object associated with the pointer that was removed:
    if( deleteObjects == true )
      delete removee;   
  }
  observedParameters.getLock().exit();
}

void AutomatableModule::parameterChanged(Parameter *parameterThatHasChanged)
{
  observedParameters.getLock().enter();
  int index = getParameterIndex(parameterThatHasChanged);

  //if( handlerFunctions[index] != NULL )
  //  handlerFunctions[index](parameterThatHasChanged->getValue());

  observedParameters.getLock().exit();
}

void AutomatableModule::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  removeObservedParameter(parameterThatWillBeDeleted, false);
}

int AutomatableModule::getParameterIndex(Parameter *parameterToLookFor)
{
  int result = -1;
  observedParameters.getLock().enter();
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    if( observedParameters[i] == parameterToLookFor )
    {
      observedParameters.getLock().exit();
      return i;
    }
  }
  observedParameters.getLock().exit();
  return -1;
}