MetaControlledParameter::MetaControlledParameter(const juce::String& name, double min, double max,
  double defaultValue, int scaling, double interval)
  : Parameter(name, min, max, defaultValue, scaling, interval)
{

}

void MetaControlledParameter::setFromMetaValue(double newMetaValue, bool sendNotification, 
  bool callCallbacks)
{
  setProportionalValue(newMetaValue, sendNotification, callCallbacks);

  // Later, we can introduce a nonlinear mapping function of the range 0..1 to itself here before
  // calling setProportionalValue. Maybe we can use an object some kind of ParameterMapper class
  // for this. We could do: 
  // if(parameterMapper == nullptr) 
  //   setProportionalValue(newMetaValue, ..);
  // else
  //   setProportionalValue(parameterMapper->map(newMetaValue), ..)
  // Subclases of this mapper class could realize Elan's rational function or we could have a 
  // breakpoint-based mapping, etc.
}

void MetaControlledParameter::setMetaParameterManager(MetaParameterManager *newManager)
{
  metaParaManager = newManager;
}

void MetaControlledParameter::attachToMetaParameter(int index)
{
  jassert(metaParaManager != nullptr); // we need a reference to a MetaParameterManager object
  if(metaParaManager == nullptr)
    return;
  bool success = metaParaManager->attachParameter(this, index);
  if(success)
    metaIndex = index;
  else
    index = -1;
}

void MetaControlledParameter::detachFromMetaParameter()
{
  jassert(metaParaManager != nullptr);
  if(metaParaManager == nullptr)
    return;
  metaParaManager->detachParameter(this);
  metaIndex = -1;
}

String MetaControlledParameter::getMetaParameterName()
{
  if(metaParaManager == nullptr || metaIndex == -1)
    return String::empty;
  return metaParaManager->getMetaParameterName(metaIndex);
}

//-------------------------------------------------------------------------------------------------

MetaParameter::MetaParameter()
{

}

void MetaParameter::attachParameter(MetaControlledParameter* p)
{
  if(contains(params, p))
    return; // already there, nothing to do, avoid recursive callbacks

  //// old:
  //p->setFromMetaValue(metaValue, true, true);
  //p->registerParameterObserver(this);
  //appendIfNotAlreadyThere(params, p);

  // new:
  if(size(params) == 0)
    metaValue = p->getProportionalValue();
  else
    p->setFromMetaValue(metaValue, false, false);
  p->registerParameterObserver(this);
  appendIfNotAlreadyThere(params, p);
  p->notifyObservers();          // notifies host that MetaParameter has (possibly) changed
  p->callValueChangeCallbacks(); // might be relevant in other contexts

  // Desired behavoir: when there are already other Parameters attached to this MetaParameter, set
  // the newly attached Parameter to the current value of the MetaParameter. If there are currently
  // none attached, let the MetaParameter take over the vealue from the attched Parameter.
  // ...but somehow the host must get notified
}

void MetaParameter::detachParameter(MetaControlledParameter* p)
{
  p->deRegisterParameterObserver(this);
  removeFirstOccurrence(params, p);
}

void MetaParameter::setMetaValue(double newValue)
{ 
  if(newValue == metaValue)
    return; // avoid superfluous updates (some DAWs continuously send constant values)

  jassert(newValue >= 0.0 && newValue <= 1.0); // must be a normalized value in the range 0..1
  metaValue = newValue;
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++)  
    params[i]->setFromMetaValue(metaValue, true, true);
  localAutomationSwitch = true;
}

void MetaParameter::parameterChanged(Parameter* p)
{
  metaValue = p->getProportionalValue();
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++) {
    if(params[i] != p)
      params[i]->setFromMetaValue(metaValue, true, true); }
  localAutomationSwitch = true;
}

void MetaParameter::parameterIsGoingToBeDeleted(Parameter* p)
{
  p->deRegisterParameterObserver(this);
  for(int i = 0; i < size(params); i++) {
    if(params[i] == p) {
      remove(params, i);
      return; }}
  // removeFirstOccurrence(params, p); cant be used because of ambiguous template parameter, so we 
  // have to do the search-loop ourselves
}

//-------------------------------------------------------------------------------------------------

void MetaParameterManager::addMetaParamater(MetaParameter* metaParameterToAdd)
{
  appendIfNotAlreadyThere(metaParams, metaParameterToAdd);
}

bool MetaParameterManager::attachParameter(MetaControlledParameter* param, int index)
{
  detachParameter(param);
  if(index >= 0 && index < size(metaParams)) {
    metaParams[index]->attachParameter(param);
    return true; }
  else {
    jassertfalse; // index out of range
    return false; }
}

void MetaParameterManager::detachParameter(MetaControlledParameter* param)
{
  for(int i = 0; i < size(metaParams); i++)
    metaParams[i]->detachParameter(param);
}

MetaParameter* MetaParameterManager::getMetaParameter(int index)
{
  if(index < 0 || index >= size(metaParams))
    return nullptr;
  return metaParams[index];
}

String MetaParameterManager::getMetaParameterName(int index)
{
  if(index < 0 || index >= size(metaParams))
    return String::empty;
  return metaParams[index]->getName();
}

void MetaParameterManager::resetAllToDefaults()
{
  for(int i = 0; i < size(metaParams); i++)
    metaParams[i]->resetToDefaultValue();
}

bool MetaParameterManager::setMetaValue(int index, double newValue)
{
  if(index < 0 || index >= size(metaParams))
    return false;
  metaParams[index]->setMetaValue(newValue);
  return true;
}

bool MetaParameterManager::setMetaName(int index, const String& newName)
{
  if(index < 0 || index >= size(metaParams))
    return false;
  metaParams[index]->setName(newName);
  return true;
}