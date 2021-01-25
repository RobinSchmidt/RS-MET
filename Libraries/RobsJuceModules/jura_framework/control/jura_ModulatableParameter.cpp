
std::vector<ModulationSource*> ModulationParticipant::dummySources;
std::vector<ModulationTarget*> ModulationParticipant::dummyTargets;
std::vector<ModulationConnection*> ModulationParticipant::dummyConnections;

const std::vector<ModulationSource*>& ModulationParticipant::getAvailableModulationSources()
{
  if(modManager)
    return modManager->getAvailableModulationSources();
  else
    return dummySources;
}

const std::vector<ModulationTarget*>& ModulationParticipant::getAvailableModulationTargets()
{
  if(modManager)
    return modManager->getAvailableModulationTargets();
  else
    return dummyTargets;
}

const std::vector<ModulationConnection*>& ModulationParticipant::getModulationConnections()
{
  if(modManager)
    return modManager->getModulationConnections();
  else
    return dummyConnections;
}

void ModulationParticipant::registerModulationSource(ModulationSource* source)
{
  if(modManager)
    modManager->registerModulationSource(source);
}

void ModulationParticipant::deRegisterModulationSource(ModulationSource* source)
{
  if(modManager)
    modManager->deRegisterModulationSource(source);
}

void ModulationParticipant::registerModulationTarget(ModulationTarget* target)
{
  if(modManager)
    modManager->registerModulationTarget(target);
}

void ModulationParticipant::deRegisterModulationTarget(ModulationTarget* target)
{
  if(modManager)
    modManager->deRegisterModulationTarget(target);
}

//-------------------------------------------------------------------------------------------------

ModulationSource::~ModulationSource() 
{
  ModulationParticipant::deRegisterModulationSource(this);
}

juce::String ModulationSource::getModulationSourceDisplayName() const
{ 
  if(displayName == "")
    return modSourceName;
  else
    return displayName;
}

bool ModulationSource::hasConnectedTargets() const
{
  if(modManager){
    const std::vector<ModulationConnection*>& 
      allConnections = modManager->getModulationConnections();
    for(int i = 0; i < size(allConnections); i++){
      if(this == allConnections[i]->source)
        return true; }}
  return false;
}

//-------------------------------------------------------------------------------------------------

ModulationTarget::~ModulationTarget() 
{
  ModulationParticipant::deRegisterModulationTarget(this);
}

void ModulationTarget::addModulationSource(ModulationSource* source)
{
  if(modManager)
    modManager->addConnection(source, this);
}

void ModulationTarget::removeModulationSource(ModulationSource* source)
{
  if(modManager)
    modManager->removeConnection(source, this);
}

bool ModulationTarget::isConnectedTo(ModulationSource* source) const
{
  if(modManager)
    return modManager->isConnected(source, this);
  return false;
}

ModulationConnection* ModulationTarget::getConnectionTo(ModulationSource* source)
{
  if(modManager)
    return modManager->getConnectionBetween(source, this);
  return nullptr;
}

std::vector<ModulationSource*> ModulationTarget::getConnectedSources() const
{
  std::vector<ModulationSource*> result;
  if(modManager)
  {
    //// with this code, the sources in the returned array are ordered in the same way as they were
    //// registered with the modManager:
    //const std::vector<ModulationSource*>& allSources = modManager->getAvailableModulationSources();
    //for(int i = 0; i < size(allSources); i++) {
    //  if(this->isConnectedTo(allSources[i]))
    //    result.push_back(allSources[i]); }

    // ...but what we want instead is the connected sources to appear in the order of the 
    // connections (also, this code is more efficient):
    const std::vector<ModulationConnection*>& connections = modManager->getModulationConnections();
    for(int i = 0; i < size(connections); i++){
      if(connections[i]->target == this)
        result.push_back(connections[i]->source); }
  }
  return result;
}

std::vector<ModulationSource*> ModulationTarget::getDisconnectedSources() const
{
  std::vector<ModulationSource*> result;
  if(modManager)
  {
    const std::vector<ModulationSource*>& allSources = modManager->getAvailableModulationSources();
    for(int i = 0; i < size(allSources); i++)
    {
      if(!this->isConnectedTo(allSources[i]))
        result.push_back(allSources[i]);
    }
  }
  return result;
}

std::vector<ModulationConnection*> ModulationTarget::getConnections() const
{
  std::vector<ModulationConnection*> result;
  if(modManager)
  {
    const std::vector<ModulationConnection*>& 
      allConnections = modManager->getModulationConnections();
    for(int i = 0; i < size(allConnections); i++)
    {
      if(this == allConnections[i]->target)
        result.push_back(allConnections[i]);
    }
  }
  return result;
}

bool ModulationTarget::hasConnectedSources() const
{
  if(modManager){
    const std::vector<ModulationConnection*>& 
      allConnections = modManager->getModulationConnections();
    for(int i = 0; i < size(allConnections); i++){
      if(this == allConnections[i]->target)
        return true; }}
  return false;
}

//-------------------------------------------------------------------------------------------------

ModulationConnection::ModulationConnection(ModulationSource* _source, ModulationTarget* _target, 
  MetaParameterManager* metaManager)
{
  source   = _source; 
  target   = _target;
  mode        = target->getDefaultModulationMode();
  //sourceValue = &(source->modValue);    // obsolete
  //targetValue = &(target->modulatedValue);  // obsolete

  double depthMin = target->getDefaultModulationDepthMin();
  double depthMax = target->getDefaultModulationDepthMax();
  double depthDef = RAPT::rsClip(0.0, depthMin, depthMax);      // default depth
  depth = RAPT::rsClip(target->getInitialModulationDepth(), depthMin, depthMax);

  juce::String name = source->getModulationSourceName(); // should we use the displayName here?
  depthParam = new 
    MetaControlledParameter(name, depthMin, depthMax, depthDef, Parameter::LINEAR, 0.0);
  depthParam->setValue(depth, false, false);
  depthParam->setValueChangeCallback<ModulationConnection>(
    this, &ModulationConnection::setDepthMember);
  depthParam->setMetaParameterManager(metaManager);

  // what, if the user changes min/max at runtime in a way such that the default value is not 
  // within the min/max range anymore? do we need additional consistency checks?
}

ModulationConnection::~ModulationConnection()
{
  delete depthParam;
}

XmlElement* ModulationConnection::getAsXml()
{
  XmlElement* xml = new XmlElement("Connection");

  xml->setAttribute("Source",   source->getModulationSourceName());
  xml->setAttribute("Target",   target->getModulationTargetName());
  xml->setAttribute("Depth",    depth);
  xml->setAttribute("DepthMin", depthParam->getMinValue()); 
  xml->setAttribute("DepthMax", depthParam->getMaxValue());
  if(depthParam->getMetaParameterIndex() != -1)
    xml->setAttribute("DepthMeta", depthParam->getMetaParameterIndex());
    // later, we may also have to store the mapping function here (if any)

  switch(mode)
  {
  case ABSOLUTE:       xml->setAttribute("Mode", "Absolute");       break;
  case RELATIVE:       xml->setAttribute("Mode", "Relative");       break;
  case EXPONENTIAL:    xml->setAttribute("Mode", "Exponential");    break;
  case MULTIPLICATIVE: xml->setAttribute("Mode", "Multiplicative"); break;
  }
  return xml;

  // todo (pseudocode):
  // if(!modMap.isIdentity())
  //   save modulation map

  // maybe move this function into ModulationManager as 
  // getConnectionXml(ModulationConnection *c);
}

//-------------------------------------------------------------------------------------------------

ModulationManager::ModulationManager(CriticalSection* lockToUse)
{
  modLock = lockToUse;
}

ModulationManager::~ModulationManager() 
{
  ScopedLock scopedLock(*modLock);
  removeAllConnections();
  deRegisterAllSources();  // to avoid dangling pointers in the sources
  deRegisterAllTargets();  // dito for the targets
}

void ModulationManager::applyModulations()
{
  ScopedLock scopedLock(*modLock); 
  applyModulationsNoLock();
}

void ModulationManager::applyModulationsNoLock()
{
  size_t i;

  // compute output signals of all modulators:
  for(i = 0; i < availableSources.size(); i++)
    availableSources[i]->updateModulationValue(); 
    // maybe we should loop only over an array of "usedSources" as optimization? similar to the way
    // we only loop over "affectedTargets"

  // initialize modulation target values with their unmodulated values:
  for(i = 0; i < affectedTargets.size(); i++)
    affectedTargets[i]->initModulatedValue();

  // apply all modulations:
  for(i = 0; i < modulationConnections.size(); i++)
    modulationConnections[i]->apply();

  // let the targets do whatever work they have to do with the modulated value (typically, 
  // call setter-callbacks):
  for(i = 0; i < affectedTargets.size(); i++)
    affectedTargets[i]->doModulationUpdate();
}

void ModulationManager::addConnection(ModulationSource* source, ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 

  //// old:
  //jassert(!isConnected(source, target)); // there is already a connection between source and target
  //modulationConnections.push_back(new ModulationConnection(source, target, metaManager));
  //appendIfNotAlreadyThere(affectedTargets, target);

  // new:
  addConnection(new ModulationConnection(source, target, metaManager));
}

void ModulationManager::addConnection(ModulationConnection* connection)
{
  ScopedLock scopedLock(*modLock); 
  jassert(!isConnected(connection->source, connection->target)); // connection already exists
  modulationConnections.push_back(connection);
  appendIfNotAlreadyThere(affectedTargets, connection->target);

  sendModulationChangeNotificationFor(connection->target); // new
}

void ModulationManager::removeConnection(ModulationSource* source, ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  jassert(isConnected(source, target)); // trying to remove no-existent connection
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->source == source && modulationConnections[i]->target == target)
      removeConnection(i, false);
  }
  updateAffectedTargetsArray();
  jassert(!isConnected(source, target)); // there must have been more than one connection between
                                         // given source and target - that should not happen
}

void ModulationManager::removeConnectionsWith(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->source == source)
    {
      removeConnection(i, false);
      i--; // array was shrunken
    }
  }
  updateAffectedTargetsArray();
}

void ModulationManager::removeConnectionsWith(ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->target == target)
    {
      removeConnection(i, false);
      i--; // array was shrunken
    }
  }
  updateAffectedTargetsArray();
}

void ModulationManager::removeAllConnections()
{
  ScopedLock scopedLock(*modLock); 
  while(size(modulationConnections) > 0)
    removeConnection(size(modulationConnections)-1, false);
  updateAffectedTargetsArray();
}

void ModulationManager::removeConnection(int i, bool updateAffectTargets)
{
  ScopedLock scopedLock(*modLock); 
  ModulationTarget* t = modulationConnections[i]->target;
  ObservableModulationTarget* omt = 
    dynamic_cast<ObservableModulationTarget*> (modulationConnections[i]->target);
  delete modulationConnections[i];
  remove(modulationConnections, i);
  if(omt)
    omt->sendModulationsChangedNotification();
  if(!t->hasConnectedSources()) {  // avoids the target getting stuck at modulated value when last 
    t->initModulatedValue();       // modulator was removed
    t->doModulationUpdate();  }
  if(updateAffectTargets)
    updateAffectedTargetsArray();
}

//void ModulationManager::resetAllTargetRangeLimits()
//{
//  ScopedLock scopedLock(*modLock); 
//  for(int i = 0; i < size(availableTargets); i++)
//  {
//    availableTargets[i]->setModulationRangeMin(-INF);
//    availableTargets[i]->setModulationRangeMax(+INF);
//  }
//}

void ModulationManager::registerModulationSource(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock);
  //sourceValues.resize(sourceValues.size()+1);
  appendIfNotAlreadyThere(availableSources, source);
  source->setModulationManager(this);
}

void ModulationManager::deRegisterModulationSource(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
  jassert(contains(availableSources, source)); // source was never registered
  removeFirstOccurrence(availableSources, source);
  removeConnectionsWith(source);
  source->setModulationManager(nullptr); 
  //sourceValues.resize(availableSources.size());
}

void ModulationManager::deRegisterAllSources()
{
  ScopedLock scopedLock(*modLock); 
  while(size(availableSources) > 0)
    deRegisterModulationSource(availableSources[size(availableSources)-1]);
}

void ModulationManager::registerModulationTarget(ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  appendIfNotAlreadyThere(availableTargets, target);
  target->setModulationManager(this);
}

void ModulationManager::deRegisterModulationTarget(ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  jassert(contains(availableTargets, target)); // target was never registered
  removeFirstOccurrence(availableTargets, target);
  removeFirstOccurrence(affectedTargets,  target);
  removeConnectionsWith(target);
  target->setModulationManager(nullptr);
}

void ModulationManager::deRegisterAllTargets()
{
  ScopedLock scopedLock(*modLock); 
  while(size(availableTargets) > 0)
    deRegisterModulationTarget(availableTargets[size(availableTargets)-1]);
}

bool ModulationManager::isConnected(const ModulationSource* source, 
  const ModulationTarget* target) const
{
  ScopedLock scopedLock(*modLock); 
  return getConnectionBetween(source, target) != nullptr;

  /*
  for(int i = 0; i < size(modulationConnections); i++)
    if(modulationConnections[i]->source == source && modulationConnections[i]->target == target)
      return true;
  return false;
  */
}

ModulationConnection* ModulationManager::getConnectionBetween(const ModulationSource* source, 
  const ModulationTarget* target) const
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
    if(modulationConnections[i]->source == source && modulationConnections[i]->target == target)
      return modulationConnections[i];
  return nullptr;
}

/*
int ModulationManager::numRegisteredSourcesOfType(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
  int result = 0;
  for(int i = 0; i < size(availableSources); i++)
  {
    if(typeid(*source) == typeid(*availableSources[i])) // gives warning on mac
      result++;
  }
  return result;
}
*/

ModulationSource* ModulationManager::getSourceByName(const juce::String& sourceName)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(availableSources); i++)
  {
    if(availableSources[i]->getModulationSourceName() == sourceName)
      return availableSources[i];
  }
  return nullptr;
}

ModulationTarget* ModulationManager::getTargetByName(const juce::String& targetName)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(availableTargets); i++)
  {
    if(availableTargets[i]->getModulationTargetName() == targetName)
      return availableTargets[i];
  }
  return nullptr;
}

bool ModulationManager::needsToStoreRangeLimits()
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(affectedTargets); i++)
  {
    if(  affectedTargets[i]->getModulationRangeMin() != -INF
      || affectedTargets[i]->getModulationRangeMax() !=  INF)
      return true;
  }
  return false;
}

void ModulationManager::setMetaParameterManager(MetaParameterManager* managerToUse)
{
  ScopedLock scopedLock(*modLock); 
  metaManager = managerToUse;
  // todo: set it to the null object, in case a nullptr is passed
}

void ModulationManager::setVoiceManager(rsVoiceManager* managerToUse)
{
  ScopedLock scopedLock(*modLock); 
  voiceManager = managerToUse;
}

void ModulationManager::setStateFromXml(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*modLock); 
  jassert(xmlState.hasTagName("Modulations"));  // not the right kind of xml element

  // recall connections:
  removeAllConnections();
  forEachXmlChildElementWithTagName(xmlState, conXml, "Connection")
  {
    ModulationSource* source = getSourceByName(conXml->getStringAttribute("Source"));
    ModulationTarget* target = getTargetByName(conXml->getStringAttribute("Target"));
    if(source != nullptr && target != nullptr)
    {
      ModulationConnection* c = new ModulationConnection(source, target, metaManager);

      c->depthParam->detachFromMetaParameter();
      int meta = conXml->getIntAttribute("DepthMeta", -1);
      if(meta >= 0)
        c->depthParam->attachToMetaParameter(meta);

      double min = conXml->getDoubleAttribute("DepthMin", -1.0);
      double max = conXml->getDoubleAttribute("DepthMax", +1.0);
      double val = conXml->getDoubleAttribute("Depth",     0.0);
      c->setDepthRangeAndValue(min, max, val);

      typedef ModulationConnection::modModes MM;
      juce::String modeStr = conXml->getStringAttribute("Mode");
      if(modeStr == "Absolute")       c->setMode(MM::ABSOLUTE);
      if(modeStr == "Relative")       c->setMode(MM::RELATIVE);
      if(modeStr == "Exponential")    c->setMode(MM::EXPONENTIAL);
      if(modeStr == "Multiplicative") c->setMode(MM::MULTIPLICATIVE);
      addConnection(c);
    }
    else
      jassertfalse; // source and/or target with given name doesn't exist - patch corrupted?
  }

  // recall range limits:
  //resetAllTargetRangeLimits();
  XmlElement* xmlLimits = xmlState.getChildByName("RangeLimits");
  if(xmlLimits != nullptr)
  {
    forEachXmlChildElement(*xmlLimits, targetLimitsXml)
    {
      ModulationTarget* target = getTargetByName(targetLimitsXml->getTagName());
      if(target != nullptr)
      {
        target->setModulationRangeMin(targetLimitsXml->getDoubleAttribute("Min", -INF));
        target->setModulationRangeMax(targetLimitsXml->getDoubleAttribute("Max", +INF));
      }
    }
  }
}

XmlElement* ModulationManager::getStateAsXml()
{
  ScopedLock scopedLock(*modLock); 

  if(size(modulationConnections) == 0)
    return nullptr; // "Modulations" child-element will be absent from preset file

  XmlElement* xmlState = new XmlElement("Modulations"); // maybe make the name a function parameter for consistency with the recall in Chainer

  for(int i = 0; i < size(modulationConnections); i++)
    xmlState->addChildElement(modulationConnections[i]->getAsXml());

  if(!needsToStoreRangeLimits())
    return xmlState;

  XmlElement* xmlLimits = new XmlElement("RangeLimits");
  for(int i = 0; i < size(affectedTargets); i++)
  {
    double min = affectedTargets[i]->getModulationRangeMin();
    double max = affectedTargets[i]->getModulationRangeMax();
    //if(min != -INF || max != INF)
    {
      XmlElement* xmlTargetLimits = new XmlElement(affectedTargets[i]->getModulationTargetName());
      //if(min != -INF) 
        xmlTargetLimits->setAttribute("Min", min);
      //if(max !=  INF) 
        xmlTargetLimits->setAttribute("Max", max);
      xmlLimits->addChildElement(xmlTargetLimits);
    }
  }
  xmlState->addChildElement(xmlLimits);
  return xmlState;
}

void ModulationManager::updateAffectedTargetsArray()
{
  ScopedLock scopedLock(*modLock); 
  affectedTargets.clear();
  for(int i = 0; i < size(modulationConnections); i++)
    appendIfNotAlreadyThere(affectedTargets, modulationConnections[i]->target);
  // not sure, if this is the best (most efficient) way to do it
}

void ModulationManager::sendModulationChangeNotificationFor(ModulationTarget* target)
{
  ObservableModulationTarget* omt = dynamic_cast<ObservableModulationTarget*> (target);
  if(omt != nullptr)
    omt->sendModulationsChangedNotification();
}

//-------------------------------------------------------------------------------------------------

void ModulatableParameter::setNormalizedValue(double newValue, bool sendNotification, 
  bool callCallbacks)
{
  MetaControlledParameter::setNormalizedValue(newValue, sendNotification, callCallbacks);
  unmodulatedValue = value;
  modulatedValue   = unmodulatedValue;
  // should we iterate over the array here? i guess not because this is supposed to be a 
  // monophonic parameter, so the 0th entry is the only one that exists

  //callValueChangeCallbacks(value); // might result in redundant call since baseclass may already 
  // call it? maybe we should use "false" for the callCallbacks parameter in the baseclass call and
  // call it here manually only when there is no smoothing and no modulation
}

void ModulatableParameter::setSmoothedValue(double newValue)
{
  unmodulatedValue = mapper->map(newValue);
  modulatedValue   = unmodulatedValue;
  callValueChangeCallbacks(unmodulatedValue);

  // later do:
  //if( hasConnectedSources() )
  //  callValueChangeCallbacks(unmodulatedValue); 
  // to avoid a superfluous callback call with the unmodulated value - but for this to make sense, 
  // we need a more efficient hasConnectedSources function - ModulationTarget must keep track of
  // the number of connected sources

  //callValueChangeCallbacks(mapper->map(newValue)); 
    // why - maybe it's wrong? maybe only needed when there's no mod-connection?

  // old:                                                
  //modulatedValue = unmodulatedValue = value = newValue;                                             
  //callValueChangeCallbacks(value); // maybe use modulatedValue
}

juce::String ModulatableParameter::getModulationTargetName()
{
  if(ownerModule == nullptr)
  {
    jassertfalse; // You need to set the owner via setOwnerAudioModule, so we can use that to 
                  // generate a unique name.
    return String();
  }
  else
    return ownerModule->getAudioModulePath() + getName();
}

//=================================================================================================


void ModulationManagerPoly::applyVoiceModulations(int i)
{    
  // todo loop over all connections to compute the modulated value and call the callback for 
  // voice i

  modulatedValues.resize(affectedTargets.size());
  // we should somehow assure that this doe not cause a re-allocation - the capacity should be 
  // large enough - the worst that can happen is min(numAvailableTargets, numConections)
  // but how do we assign the array-slots to the targets and how do we figure out when two 
  // connections go to the same target? ...maybe for this, we should also have the array sorted
  // by target
  // or the targets somehow need to store an array index that they want to use - but that must be
  // uodated when the number of targets changes...maybe its beset to just have an array of
  // size numAvailableTargets - that seems most practical

  // hmm.. i think, we need to make sure that connections that go to the same target are 
  // adjacent in the array of connections..we should also maintain a variable 
  // numDistinctActiveTargets - we should nevertheless compute all values first and then apply
  // all callbacks in one go after that in order to not introduce any dependency of the result
  // on the actual ordering of the connections

  for(i = 0; i < modulationConnections.size(); i++)
  {
    //modulationConnections[i]->apply();

    // todo:
    // -retrieve unmodulated value from source


    // it would be helpful, if the connections would be sorted by their target - then we could 
    // process one target at a time: loop over all connections involving that target, when done
    // call the callback and then move on to the next target
    // ...does it actually make any difference, in which order the connections are processed?
    // maybe only when we have some sort of feedback modulation

    // if we can't assume any ordering, we need a buffer of target-values whose size equals
    // affectedTargets.size

  }



  int dummy = 0;
}

void ModulationManagerPoly::applyModulationsNoLock()
{
  if(voiceManager != nullptr)
    for(int i = 0; i < voiceManager->getNumActiveVoices(); i++)
      applyVoiceModulations(voiceManager->getActiveVoiceIndex(i));

  ModulationManager::applyModulationsNoLock();
}
// It should perhaps first apply all the polyphonic modulations and then call the baseclass 
// method updateModulationValue should just copy the value from the newest voice into 
// modulationValue, if the source is polyphonic.



void ModulationManagerPoly::addConnection(ModulationConnection* connection)
{
  ModulationManager::addConnection(connection); // preliminary

  // ToDo: do not just use push_back but instead insert it at location thta ensures that 
  // connections with the same target are neighbors in the array, also increments 
  // numDistinctActiveTargets if the was not connection with that target yet, maybe resize the
  // modulatedValues buffer

}

void ModulationManagerPoly::removeConnection(int index, bool updateAffectTargets)
{
  ModulationManager::removeConnection(index, updateAffectTargets); // preliminary

  // ToDo: possibly decrement numDistinctActiveTargets, if this was the last connection to that
  // target, maybe resize the modulatedValues buffer
}

