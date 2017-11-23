
// experimental: null objects (as in https://sourcemaking.com/design_patterns/null_object) to be
// used by default:
//MetaParameterManager nullMetaParameterManager;
//NormalizedParameterMapper nullParameterMapper;


MetaControlledParameter::MetaControlledParameter(const juce::String& name, double min, double max,
  double defaultValue, int scaling, double interval)
  : rsSmoothableParameter(name, min, max, defaultValue, scaling, interval)
{
  // todo: initialize the metaParaManager and mapper member to Null Objects
  // like: 
  // metaParaManager = nullMetaParameterManager;
  // mapper = nullNormalizedParameterMapper; // identity mapper
}

void MetaControlledParameter::setFromMetaValue(double newMetaValue, bool sendNotification, 
  bool callCallbacks)
{
  if(mapper == nullptr)
    setProportionalValue(newMetaValue, sendNotification, callCallbacks);
  else
    setProportionalValue(mapper->map(newMetaValue), sendNotification, callCallbacks);

  //setProportionalValue(newMetaValue, sendNotification, callCallbacks);

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
  if(newManager != nullptr)
    metaParaManager = newManager;
  else
  {
    metaParaManager = nullptr;  // preliminary

    //jassertfalse;    
    // todo: if a nullptr is passed, set the metaParaManager to a Null Object, like so:
    //metaParaManager = &nullMetaParameterManager;
  }

}

void MetaControlledParameter::setParameterMapper(NormalizedParameterMapper *newMapper)
{
  mapper = newMapper;

  // maybe we should call something like 
  // setFromMetaValue(getMetaValue());
  // here, so a user could re-configure this object at runtime with a new mapper object and have 
  // that reconfiguration immediately reflected by the parameter values. getMetaValue would be 
  // implemented by (somehow) quering the metaParaManager (if not null) or use the proportional 
  // value of this Parameter (if metaParaManager null) 
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

bool MetaParameter::detachParameter(MetaControlledParameter* p)
{
  p->deRegisterParameterObserver(this);
  return removeFirstOccurrence(params, p);
}

void MetaParameter::setMetaValue(double newValue)
{ 
  jassert(newValue >= 0.0 && newValue <= 1.0); // must be a normalized value in the range 0..1
  metaValue = newValue;
  setLocalAutomationSwitch(false); // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++)  
    params[i]->setFromMetaValue(metaValue, true, true);
  setLocalAutomationSwitch(true);
}

void MetaParameter::parameterChanged(Parameter* p)
{
  metaValue = p->getProportionalValue();
  setLocalAutomationSwitch(false); // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++) {
    if(params[i] != p)
      params[i]->setFromMetaValue(metaValue, true, true); }
  setLocalAutomationSwitch(true);
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
    updateMetaName(index);
    return true; }
  else {
    jassertfalse; // index out of range
    return false; }
}

void MetaParameterManager::detachParameter(MetaControlledParameter* param)
{
  for(int i = 0; i < size(metaParams); i++) {
    if(metaParams[i]->detachParameter(param)) // checks, if actual detachment took place...
      updateMetaName(i); }                    // ...if so, name should be updated
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
  for(unsigned int i = 0; i < metaObservers.size(); ++i)
    metaObservers[i]->metaNameChanged(this, index);
  return true;
}

void MetaParameterManager::updateMetaName(int index)
{
  if(autoUpdateMetaNames && index >= 0 && index < size(metaParams))
  {
    auto sz = size(metaParams[index]->params);
    auto idxStr = String(index);
    if (sz > 0)
    {
      String name = idxStr+": ";

      for (auto i = 0; i < sz; i++)
        name += metaParams[index]->params[i]->getName()+", ";

      name = name.substring(0, name.length()-2);

      setMetaName(index, name);
    }
    else
    {
      setMetaName(index, "Meta "+idxStr);
    }      
  }
}