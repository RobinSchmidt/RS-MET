//-------------------------------------------------------------------------------------------------
// construction/destruction:

ParameterManager::ParameterManager(CriticalSection *lockToUse)
{
  lock = lockToUse;
}

ParameterManager::~ParameterManager()
{
  ScopedPointerLock scopedLock(lock);
  removeAllObservedParameters(true);
}



//-------------------------------------------------------------------------------------------------
// retrieve pointers to the observed parameters:

Parameter* ParameterManager::getParameterByName(const String& nameOfParameter) const
{
  ScopedLock scopedLock(*lock);
  Parameter* result = nullptr;
  for(int i=0; i < (int) parameters.size(); i++)
  {
    if( parameters[i]->hasName(nameOfParameter) )
      result = parameters[i];
  }
  jassert(result != nullptr);   // parameter with given name doesn't exist
  return result;
}

bool ParameterManager::hasParameterWithName(const juce::String& name) const
{
  ScopedLock scopedLock(*lock);
  for(int i=0; i < (int)parameters.size(); i++)
    if( parameters[i]->hasName(name) )
      return true;
  return false;
}

Parameter* ParameterManager::getParameterByIndex(int indexOfParameter) const
{
  ScopedLock scopedLock(*lock);
  Parameter* result = nullptr;
  if( indexOfParameter < (int) parameters.size() )
    result = parameters[indexOfParameter];
  return result;
}

int ParameterManager::getIndexOfParameter(Parameter* parameterToRetrieveIndexOf) const
{
  ScopedLock scopedLock(*lock);
  int parameterIndex = -1;
  for(int i=0; i < (int) parameters.size(); i++)
  {
    if( parameterToRetrieveIndexOf == parameters[i] )
      parameterIndex = i;
  }
  return parameterIndex;
}

int ParameterManager::getNumParameters() const
{
  ScopedLock scopedLock(*lock);
  return (int)parameters.size();
}

//-------------------------------------------------------------------------------------------------
// add/remove observed parameters:

void ParameterManager::addObservedParameter(Parameter *p)
{
  ScopedLock scopedLock(*lock);
  jassert(!hasParameterWithName(p->getName())); // parameters must have unique names
  p->setMutexToUse(lock);
  p->registerParameterObserver(this);
  appendIfNotAlreadyThere(parameters, p);
}

void ParameterManager::removeObservedParameter(Parameter *parameterToRemove, bool deleteObject)
{
  ScopedLock scopedLock(*lock);
  int i=0;
  while( i < (int) parameters.size() )
  {
    if( parameters[i] == parameterToRemove )
    {
      parameterToRemove->deRegisterParameterObserver(this);
      parameterToRemove->setMutexToUse(nullptr);
      removeFirstOccurrence(parameters, parameterToRemove);
      if( deleteObject == true )
        delete parameterToRemove;
      i--; // because array has shrunken
    }
    i++;
  }
}

void ParameterManager::removeParameter(const juce::String& name, bool deleteObject)
{
  ScopedLock scopedLock(*lock);
  removeObservedParameter(getParameterByName(name), deleteObject);
}

void ParameterManager::removeAllObservedParameters(bool deleteObjects)
{
  ScopedPointerLock scopedLock(lock);
  Parameter *removee; // this is the currently removed parameter
  while( parameters.size() > 0 )
  {
    removee = parameters[0];
    removee->deRegisterParameterObserver(this);
    removee->setMutexToUse(nullptr);
    remove(parameters, 0);
    if( deleteObjects == true )
      delete removee;
  }
}

void ParameterManager::parameterChanged(Parameter *parameterThatHasChanged)
{

}

void ParameterManager::parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  ScopedLock scopedLock(*lock);
  removeObservedParameter(parameterThatWillBeDeleted, false);
}

int ParameterManager::getParameterIndex(Parameter *parameterToLookFor)
{
  ScopedLock scopedLock(*lock);
  //int result = -1;
  for(int i=0; i < (int) parameters.size(); i++) {
    if( parameters[i] == parameterToLookFor )
      return i; }
  return -1;
}
