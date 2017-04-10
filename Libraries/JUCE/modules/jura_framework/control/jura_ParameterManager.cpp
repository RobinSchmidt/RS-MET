//-------------------------------------------------------------------------------------------------
// construction/destruction:

ParameterManager::ParameterManager(CriticalSection *lockToUse)
{
  lock = lockToUse;
}

ParameterManager::~ParameterManager()
{
  ScopedLock scopedLock(*lock);
  removeAllObservedParameters(true);
}
    


//-------------------------------------------------------------------------------------------------
// retrieve pointers to the observed parameters:

Parameter* ParameterManager::getParameterByName(const String& nameOfParameter) const
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

Parameter* ParameterManager::getParameterByIndex(int indexOfParameter) const
{
  ScopedLock scopedLock(*lock);
  Parameter* result = nullptr;
  if( indexOfParameter < (int) observedParameters.size() )
    result = observedParameters[indexOfParameter];
  return result;
}

int ParameterManager::getIndexOfParameter(Parameter* parameterToRetrieveIndexOf) const
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

int ParameterManager::getNumParameters() const
{
  ScopedLock scopedLock(*lock);
  return (int)observedParameters.size();
}

//-------------------------------------------------------------------------------------------------
// add/remove observed parameters:

void ParameterManager::addObservedParameter(Parameter *parameterToAdd)
{
  ScopedLock scopedLock(*lock);
  parameterToAdd->setMutexToUse(lock);
  parameterToAdd->registerParameterObserver(this);
  appendIfNotAlreadyThere(observedParameters, parameterToAdd);
}

void ParameterManager::removeObservedParameter(Parameter *parameterToRemove, bool deleteObject)
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

void ParameterManager::removeAllObservedParameters(bool deleteObjects)
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

void ParameterManager::parameterChanged(Parameter *parameterThatHasChanged)
{

}

void ParameterManager::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  ScopedLock scopedLock(*lock);
  removeObservedParameter(parameterThatWillBeDeleted, false);
}

int ParameterManager::getParameterIndex(Parameter *parameterToLookFor)
{
  ScopedLock scopedLock(*lock);
  int result = -1;
  for(int i=0; i < (int) observedParameters.size(); i++) {
    if( observedParameters[i] == parameterToLookFor )
      return i; }
  return -1;
}
