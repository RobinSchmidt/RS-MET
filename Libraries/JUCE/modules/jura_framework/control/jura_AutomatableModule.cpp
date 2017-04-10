//-------------------------------------------------------------------------------------------------
// construction/destruction:

AutomatableModule::AutomatableModule(CriticalSection *lockToUse)
{
  lock = lockToUse;
}

AutomatableModule::~AutomatableModule()
{
  ScopedLock scopedLock(*lock);
  removeAllObservedParameters(true);
}
    


//-------------------------------------------------------------------------------------------------
// retrieve pointers to the observed parameters:

Parameter* AutomatableModule::getParameterByName(const String& nameOfParameter) const
{
  ScopedLock scopedLock(*lock);
  Parameter* result = nullptr;
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    if( observedParameters[i]->hasName(nameOfParameter) )
      result = observedParameters[i];
  }
  jassert(result != nullptr);   // parameter with given name doesn't exist
  return result;
}

Parameter* AutomatableModule::getParameterByIndex(int indexOfParameter) const
{
  ScopedLock scopedLock(*lock);
  Parameter* result = nullptr;
  if( indexOfParameter < (int) observedParameters.size() )
    result = observedParameters[indexOfParameter];
  return result;
}

int AutomatableModule::getIndexOfParameter(Parameter* parameterToRetrieveIndexOf) const
{
  ScopedLock scopedLock(*lock);
  int parameterIndex = -1;
  for(int i=0; i < (int) observedParameters.size(); i++)
  {
    if( parameterToRetrieveIndexOf == observedParameters[i] )
      parameterIndex = i;
  }  
  return parameterIndex;
}

int AutomatableModule::getNumParameters() const
{
  ScopedLock scopedLock(*lock);
  return (int)observedParameters.size();
}

//-------------------------------------------------------------------------------------------------
// add/remove observed parameters:

void AutomatableModule::addObservedParameter(Parameter *parameterToAdd)
{
  ScopedLock scopedLock(*lock);
  parameterToAdd->setMutexToUse(lock);
  parameterToAdd->registerParameterObserver(this);
  appendIfNotAlreadyThere(observedParameters, parameterToAdd);
}

void AutomatableModule::removeObservedParameter(Parameter *parameterToRemove, bool deleteObject)
{
  ScopedLock scopedLock(*lock);
  int i=0;
  while( i < (int) observedParameters.size() ) 
  {
    if( observedParameters[i] == parameterToRemove ) 
    {
      parameterToRemove->deRegisterParameterObserver(this);
      parameterToRemove->setMutexToUse(nullptr);
      removeFirstOccurrence(observedParameters, parameterToRemove);
      if( deleteObject == true )
        delete parameterToRemove;
      i--; // because array has shrunken
    }
    i++;
  }
}

void AutomatableModule::removeAllObservedParameters(bool deleteObjects)
{
  ScopedLock scopedLock(*lock);
  Parameter *removee; // this is the currently removed parameter
  while( observedParameters.size() > 0 )
  {
    removee = observedParameters[0];
    removee->deRegisterParameterObserver(this);
    removee->setMutexToUse(nullptr);
    remove(observedParameters, 0);
    if( deleteObjects == true )
      delete removee;   
  }
}

void AutomatableModule::parameterChanged(Parameter *parameterThatHasChanged)
{

}

void AutomatableModule::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  ScopedLock scopedLock(*lock);
  removeObservedParameter(parameterThatWillBeDeleted, false);
}

int AutomatableModule::getParameterIndex(Parameter *parameterToLookFor)
{
  ScopedLock scopedLock(*lock);
  int result = -1;
  for(int i=0; i < (int) observedParameters.size(); i++) {
    if( observedParameters[i] == parameterToLookFor )
      return i; }
  return -1;
}
