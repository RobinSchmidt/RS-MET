rsMetaParameterMapper::rsMetaParameterMapper()
{
  initToIdentity();
}

/*
size_t rsMetaParameterMapper::addNode(double x, double y)
{
  x = clip(x, 0, 1); 
  y = clip(y, 0, 1);
  return RAPT::rsNodeBasedFunction<double>::addNode(x, y);
}

bool rsMetaParameterMapper::removeNode(size_t index)
{
  if(index == 0 || index == nodes.size()-1) 
    return false; // first and last node can't be removed
  return RAPT::rsNodeBasedFunction<double>::removeNode(index);
}
*/

size_t rsMetaParameterMapper::moveNode(size_t index, double x, double y)
{
  x = RAPT::rsClip(x, 0.0, 1.0); 
  y = RAPT::rsClip(y, 0.0, 1.0);             // x and y must be in 0..1
  if(index == 0)               x = 0;    // first node's x value is fixed at 0
  if(index == nodes.size()-1)  x = 1;    // last node's x value is fixed at 1
  return RAPT::rsNodeBasedFunction<double>::moveNode(index, x, y);
}
// without this override, it would still possible to have first or last nodes x coordinates 
// different from 0 and 1 when the 1st or last node is moved through the position of another 
// node - then the other one becomes first or last but retains its position - it's something
// about the constrainNode function being called after re-indexing or something ...i need to
// figure that out..or not...it works like it is now, it's just unelegant

bool rsMetaParameterMapper::isNodeRemovable(size_t index) const
{
  if(index == 0 || index == nodes.size()-1) 
    return false; // first and last node can't be removed
  return true;
}

size_t rsMetaParameterMapper::constrainNode(size_t i)
{
  //return i; // preliminary
  nodes[i].x = RAPT::rsClip(nodes[i].x, 0.0, 1.0); 
  nodes[i].y = RAPT::rsClip(nodes[i].y, 0.0, 1.0);     // x and y must be in 0..1
  i = RAPT::rsNodeBasedFunction<double>::moveNodeToSortedIndex(i);
  if(i == 0)               
    nodes[i].x = 0; // first node's x value is fixed at 0
  if(i == nodes.size()-1)  
    nodes[i].x = 1; // last  node's x value is fixed at 1
  return RAPT::rsNodeBasedFunction<double>::moveNodeToSortedIndex(i);
}

bool rsMetaParameterMapper::isDefaultMap() const
{
  if(nodes.size() == 2  && nodes[1].shapeType == RAPT::rsFunctionNode<double>::LINEAR
    && nodes[0].getX() == 0.0 && nodes[0].getY() == 0.0 
    && nodes[1].getX() == 1.0 && nodes[1].getY() == 1.0)
    return true;
  return false;
}

void rsMetaParameterMapper::initToIdentity()
{
  nodes.clear(); 
  appendNode(0, 0); 
  appendNode(1, 1); 
}

void rsMetaParameterMapper::initToFlat(double v)
{
  nodes.clear(); 
  appendNode(0, v); 
  appendNode(1, v); 
}

XmlElement* rsMetaParameterMapper::getStateAsXml(const juce::String& tagName) const
{
  typedef RAPT::rsFunctionNode<double> FN;
  XmlElement* mapXml = new XmlElement(tagName);
  for(size_t i = 0; i < nodes.size(); i++) {
    XmlElement* nodeXml = new XmlElement("Node");
    nodeXml->setAttribute("X", nodes[i].getX());
    nodeXml->setAttribute("Y", nodes[i].getY());
    if(nodes[i].getShapeType() != FN::LINEAR) {
      nodeXml->setAttribute("Shape", shapeIndexToString(nodes[i].getShapeType()));
      nodeXml->setAttribute("ShapeParameter", nodes[i].getShapeParameter());
    }
    mapXml->addChildElement(nodeXml);
  }
  return mapXml;
}

void rsMetaParameterMapper::setStateFromXml(const XmlElement& mapXml)
{
  typedef RAPT::rsFunctionNode<double> FN;
  nodes.clear();
  forEachXmlChildElementWithTagName(mapXml, nodeXml, "Node") {
    double x  = nodeXml->getDoubleAttribute("X", 0.0);
    double y  = nodeXml->getDoubleAttribute("Y", 0.0);
    double sp = nodeXml->getDoubleAttribute("ShapeParameter", 0.0);
    int shape = stringToShapeIndex(nodeXml->getStringAttribute("Shape", "Linear"));
    appendNode(x, y, shape, sp);
    //addNode(x, y); // can't be used because when the array is empty, the constraint checker 
                     // doesn't work properly
  }
  jassert(nodes.size() >= 2); // xml corrupted? it should have at least 2 nodes
  if(nodes.size() < 2) 
    initToIdentity();
}

String rsMetaParameterMapper::shapeIndexToString(int index)
{
  typedef RAPT::rsFunctionNode<double> FN;
  switch(index)
  {
  case FN::LINEAR:      return "Linear";
  case FN::EXPONENTIAL: return "Exponential";
  case FN::RATIONAL:    return "Rational";
  default:              return "Linear";
  }
}

int rsMetaParameterMapper::stringToShapeIndex(const String& s)
{
  typedef RAPT::rsFunctionNode<double> FN;
  if(s == "Linear")       return FN::LINEAR;
  if(s == "Exponential")  return FN::EXPONENTIAL;
  if(s == "Rational")     return FN::RATIONAL;
  return FN::LINEAR;
}


// experimental: null object (as in https://sourcemaking.com/design_patterns/null_object) to be
// used by default:
//MetaParameterManager nullMetaParameterManager;


//=================================================================================================

MetaControlledParameter::MetaControlledParameter(const juce::String& name, double min, double max,
  double defaultValue, int scaling, double interval)
  : rsSmoothableParameter(name, min, max, defaultValue, scaling, interval)
{
  // todo: initialize the metaParaManager and mapper member to Null Objects
  // like: 
  // metaParaManager = nullMetaParameterManager;
  // mapper = nullNormalizedParameterMapper; // identity mapper

  //unmappedValue = normalizedValue;
}

void MetaControlledParameter::setFromMetaValue(double newMetaValue, bool sendNotification, 
  bool callCallbacks)
{
  setNormalizedValue(newMetaValue, sendNotification, callCallbacks);
}

void MetaControlledParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  newValue = restrictValueToParameterRange(newValue);
  double y = mapper->unmap(newValue);
  double x = metaMapper.unmap(y);
  //jassert(x >= 0 && x <= 1);
  x = RAPT::rsClip(x, 0.0, 1.0);
  setNormalizedValue(x, sendNotification, callCallbacks);

  // preliminary
  //rsSmoothableParameter::setValue(newValue, sendNotification, callCallbacks);

  // todo: find normalized value that corresponds to desired "value" and call setNormalizedValue
  // with that value
}

void MetaControlledParameter::setNormalizedValue(double newNormalizedValue, bool sendNotification,
  bool callCallbacks)
{
  if(normalizedValue == newNormalizedValue)
    return;
  if(!needsSmoothing())
  {
    normalizedValue = newNormalizedValue;
    value = applyBothMaps(normalizedValue); // maybe apply restrictions here
    if( callCallbacks == true )    callValueChangeCallbacks(value);
    if( sendNotification == true ) notifyObservers();
  }
  else
    setNormalizedTargetValue(newNormalizedValue, sendNotification, callCallbacks);

  /*
  unmappedValue = newValue;
  // hmm...maybe the inherited normalizedValue should always store the value pre-mapping and we
  // may have to copy over the function body of rsSmoothableParameter setNormalizedValue here and
  // just change
  // smoothingManager->addSmootherFor(this, normalizedValue, oldNormalizedValue); 
  // into:
  // smoothingManager->addSmootherFor(this, mappedNormalizedValue, oldMappedNormalizedValue); 
  // the we can rid of this weird (fabs(new-old)<tol) stuff there too - it doesn't belong there

  rsSmoothableParameter::setNormalizedValue(mapper.map(newValue), sendNotification, callCallbacks);
  */
}

void MetaControlledParameter::setNormalizedTargetValue(double newTargetValue, bool sendNotification,
  bool callCallbacks)
{
  double oldNormalizedValue = normalizedValue;

  double tol = 1.e-7;
  if(fabs(oldNormalizedValue-newTargetValue) < tol)
  {
    // code is same as in if(!needsSmoothing()) branch in setNormalizedValue 
    // -> get rid of duplication
    normalizedValue = newTargetValue;
    value = applyBothMaps(normalizedValue); // maybe apply restrictions here
    if( callCallbacks == true )    callValueChangeCallbacks(value);
    if( sendNotification == true ) notifyObservers();
    return;
  }
  // When this parameter has an attached meta, this function gets called twice. First, when 
  // setting the parameter from a slider and a second time from the meta. In the second call,
  // getValue will return a value that is already the new value, but only up to roundoff, so
  // the if(value == newValue) check doesn't trigger. This would effectively disable smoothing,
  // so we need this additional check here.
  // -may need some more thorough checking, especially with regard to the tolerance value and if 
  //  we should also use a relative tolerance..
  // -we may also set the value to newValue and invoke a callback and notification
  // maybe we should do this in MetaControlledParameter


  // copy/pasted/edited from rsSmoothableParameter:
  shouldSendNotification = sendNotification;
  //Parameter::setNormalizedValue(newTargetValue, false, false);
  normalizedValue = newTargetValue;
  value = applyBothMaps(normalizedValue); // maybe apply restrictions here
  if(sendNotification)
    notifyObserversPreSmoothing();
  smoothingManager->addSmootherFor(this, 
    metaMapper.map(normalizedValue), metaMapper.map(oldNormalizedValue)); 
}

void MetaControlledParameter::metaMapChanged()
{
  value = applyBothMaps(normalizedValue); 
  callValueChangeCallbacks(value);
  notifyObservers();
}

void MetaControlledParameter::initMetaMapToFlat() 
{ 
  //double v = getNormalizedValue();
  // this is not good - when there is currently a non-identity meta-map applied, the value may be
  // different from what it should be - we should retrieve the completely mapped value (after both 
  // maps) and undo the 2nd map (the one from normalized to physical units), like this:
  double v = mapper->unmap(getValue());
  // or v = metaMapper->map(getNormalizedValue());
  metaMapper.initToFlat(v);
}

void MetaControlledParameter::saveToXml(XmlElement* xml) const
{
  if(!metaMapper.isDefaultMap()) {
    XmlElement* mapXml = metaMapper.getStateAsXml(getName() + "ParameterMap");
    if(!hasAttachedMeta())
      mapXml->setAttribute("NormalizedValue", getNormalizedValue());
    xml->addChildElement(mapXml);
    // If there is a user-defined mapping function, we cannot just store the (final, doubly mapped) 
    // parameter value in the xml and rely on recovering the normalized value via the inverse 
    // mapping because user-mappings may be non-invertible. So we must also store the normalized 
    // value in the xml - but only if there's no meta attached because otherwise, the attached meta
    // is responsible to store the normalized value and we want to avoid redundant storage of 
    // values because the behavior may be confusing when manually editing xml files).
  }
  else
    rsSmoothableParameter::saveToXml(xml); // stores the (mapped) value directly as xml-attribute
}

void MetaControlledParameter::recallFromXml(const XmlElement& xml) 
{
  auto mapXml = xml.getChildByName(getName() + "ParameterMap");
  if(mapXml != nullptr) {
    metaMapper.setStateFromXml(*mapXml);
    double normVal = mapXml->getDoubleAttribute("NormalizedValue", getNormalizedDefaultValue());
    setNormalizedValue(normVal, true, true);
    // this may be later changed again due to re-attaching the meta as stored in xml ...verify this
  }
  else {
    metaMapper.initToIdentity();
    rsSmoothableParameter::recallFromXml(xml);
  }

  metaMapChanged();
  // This is necesarry to update the slider in case of an attached map but no attached meta and the
  // normalized (user-map)mapped value happens to be equal to th normalized default value
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

void MetaControlledParameter::attachToMetaParameter(int index, bool flatMap)
{
  jassert(metaParaManager != nullptr); // we need a reference to a MetaParameterManager object
  if(metaParaManager == nullptr)
    return;
  bool success = metaParaManager->attachParameter(this, index, flatMap);
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

MetaParameter* MetaControlledParameter::getAttachedMeta()
{
  if(metaParaManager == nullptr || metaIndex == -1)  // is the || -1 redundant? (checked also in metaParaManager->getMetaParameter?)
    return nullptr;
  return metaParaManager->getMetaParameter(metaIndex);
}

String MetaControlledParameter::getMetaParameterName()
{
  if(metaParaManager == nullptr || metaIndex == -1)
    return String();
  return metaParaManager->getMetaParameterName(metaIndex);
}


//-------------------------------------------------------------------------------------------------

MetaParameter::MetaParameter()
{
  notifyPreSmoothing(true);
  notifyPostSmoothing(false);
}

void MetaParameter::attachParameter(MetaControlledParameter* p, bool flatMap)
{
  if(contains(params, p))
    return; // already there, nothing to do, avoid recursive callbacks

  //// old:
  //p->setFromMetaValue(metaValue, true, true);
  //p->registerParameterObserver(this);
  //appendIfNotAlreadyThere(params, p);

  // new:
  if(size(params) == 0)
    metaValue = p->getNormalizedValue();
  else {
    if(flatMap)
      p->initMetaMapToFlat();
    else
      p->setFromMetaValue(metaValue, false, false);
  }
  p->registerParameterObserver(this);
  appendIfNotAlreadyThere(params, p);
  if(!flatMap) {
    p->notifyObservers();               // notifies host that MetaParameter has (possibly) changed
    p->callValueChangeCallbacks(p->getValue()); // might be relevant in other contexts
  }

  // Desired behavoir: when there are already other Parameters attached to this MetaParameter, set
  // the newly attached Parameter to the current value of the MetaParameter. If there are currently
  // none attached, let the MetaParameter take over the vealue from the attched Parameter.
  // ...but somehow the host must get notified

  // Currently the "Meta attach (flat)" option has confusing behavior...
}

bool MetaParameter::detachParameter(MetaControlledParameter* p)
{
  //p->initMetaMapToIdentity(); // reset map on detach..hmm...maybe not
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
  metaValue = p->getNormalizedValue();
  setLocalAutomationSwitch(false); // so we don't call ourselves recursively
  for(int i = 0; i < size(params); i++) {
    if(params[i] != p)
      params[i]->setFromMetaValue(metaValue, true, true); }
  setLocalAutomationSwitch(true);
}

void MetaParameter::parameterWillBeDeleted(Parameter* p)
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

bool MetaParameterManager::attachParameter(MetaControlledParameter* param, int index, bool flatMap)
{
  detachParameter(param);
  if(index >= 0 && index < size(metaParams)) {
    metaParams[index]->attachParameter(param, flatMap);
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
    return String();
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