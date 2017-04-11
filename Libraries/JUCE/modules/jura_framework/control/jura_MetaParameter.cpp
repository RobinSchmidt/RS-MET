MetaControlledParameter::MetaControlledParameter(const juce::String& name, double min, double max,
  double defaultValue, int scaling, double interval)
  : Parameter(name, min, max, defaultValue, scaling, interval)
{
  proportionalValue = valueToProportion(value);
}

void MetaControlledParameter::setProportionalValue(double newProportionalValue,
  bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);
  proportionalValue = newProportionalValue;
  setValue(proportionToValue(proportionalValue), sendNotification, callCallbacks);
}

double MetaControlledParameter::valueToProportion(double value)
{
  if(minValue >= maxValue)
    return 0.0;
  switch( scaling )
  {
  case Parameter::EXPONENTIAL: 
  {
    if( minValue > 0.0 )
      return jlimit(0.0, 1.0, log(value/minValue) / (log(maxValue)-log(minValue)) );
    else
      return 0.0;
  }
  default: return (value - minValue) / (maxValue - minValue); // LINEAR(_BIPOLAR)
  }
}

double MetaControlledParameter::proportionToValue(double prop)
{
  switch( scaling )
  {
  case Parameter::LINEAR:         return minValue + (maxValue - minValue) * prop;
  case Parameter::EXPONENTIAL:    return minValue * exp(prop*(log(maxValue)-log(minValue)));
  case Parameter::LINEAR_BIPOLAR: return minValue + (maxValue - minValue) * prop;
  default: return 0.0;
  }
  // maybe make valueToProportion/proportionToValue static methods and pass in min, max, scaling
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

//-------------------------------------------------------------------------------------------------

MetaParameter::MetaParameter()
{

}

void MetaParameter::attachParameter(MetaControlledParameter* p)
{
  if(contains(params, p))
    return; // already there, nothing to do, avoid recursive callbacks
  p->setProportionalValue(metaValue, true, true);
  p->registerParameterObserver(this);
  appendIfNotAlreadyThere(params, p);
}

void MetaParameter::detachParameter(MetaControlledParameter* p)
{
  p->deRegisterParameterObserver(this);
  removeFirstOccurrence(params, p);
}

void MetaParameter::setMetaValue(double newValue)
{ 
  jassert(newValue >= 0.0 && newValue <= 1.0); // must be a normalized value in the range 0..1
  metaValue = newValue;
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++)  
    params[i]->setProportionalValue(metaValue, true, true);
  localAutomationSwitch = true;
}

void MetaParameter::parameterChanged(Parameter* p)
{
  MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*>(p);
  metaValue = mcp->getProportionalValue();
  localAutomationSwitch = false; // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++) {
    if(params[i] != mcp)
      params[i]->setProportionalValue(metaValue, true, true); }
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
