
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

bool ModulationTarget::isConnectedTo(ModulationSource* source)
{
  if(modManager)
    return modManager->isConnected(source, this);
  return false;
}

std::vector<ModulationSource*> ModulationTarget::getConnectedSources()
{
  std::vector<ModulationSource*> result;
  if(modManager)
  {
    const std::vector<ModulationSource*>& allSources = modManager->getAvailableModulationSources();
    for(int i = 0; i < size(allSources); i++)
    {
      if(this->isConnectedTo(allSources[i]))
        result.push_back(allSources[i]);
    }
  }
  return result;
}

std::vector<ModulationSource*> ModulationTarget::getDisconnectedSources()
{
  std::vector<ModulationSource*> result;
  if(modManager)
  {
    const std::vector<ModulationSource*>& allSources = modManager->getAvailableModulationSources();
    for(int i = 0; i < size(allSources); i++)
    {
      if(!this->isConnectedTo(allSources[i])) // ! is only difference to getConnectedSources, maybe
        result.push_back(allSources[i]);      // we can refactor to avoid code duplication
    }
  }
  return result;
}

std::vector<ModulationConnection*> ModulationTarget::getConnections()
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

//-------------------------------------------------------------------------------------------------

ModulationConnection::ModulationConnection(ModulationSource* _source, ModulationTarget* _target)
{
  source   = _source; 
  target   = _target;
  amount   = 0.0;
  //relative = false;
  relative = true;
  sourceValue = &(source->modValue);
  targetValue = &(target->modulatedValue);

  juce::String name = source->getModulationSourceName();
  amountParam = new MetaControlledParameter(name, -1.0, 1.0, 0.0, Parameter::LINEAR, 0.0);
  amountParam->setValueChangeCallback<ModulationConnection>(this, 
    &ModulationConnection::setAmount);

  //amountParam->setMetaParameterManager(metaManager);
  // we may need a MetaParameterManager* parameter to the constructor ..or a pointer to the
  // ModulationManager...which may have to be a subclass of MetaParameterManager
}

ModulationConnection::~ModulationConnection()
{
  delete amountParam;
}

XmlElement* ModulationConnection::getAsXml()
{
  XmlElement* xml = new XmlElement("Connection");

  xml->setAttribute("Source", source->getModulationSourceName());
  xml->setAttribute("Target", target->getModulationTargetName());
  xml->setAttribute("Depth",  amount);
  juce::String modeString;
  if(relative)
    modeString = "Relative";
  else
    modeString = "Absolute";
    // we do it like this in order to be able to define more modes later - like Multiplicative
    // or whatever

  xml->setAttribute("Mode", modeString);
  return xml;

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
}

void ModulationManager::applyModulations()
{
  ScopedLock scopedLock(*modLock); 
  applyModulationsNoLock();
}

void ModulationManager::applyModulationsNoLock()
{
  int i;

  // compute output signals of all modulators:
  for(i = 0; i < size(availableSources); i++)
    availableSources[i]->updateModulationValue();

  // initialize modulation target values with their unmodulated values:
  for(i = 0; i < size(affectedTargets); i++)
    affectedTargets[i]->initModulatedValue();

  // apply all modulations:
  for(i = 0; i < size(modulationConnections); i++)
    modulationConnections[i]->apply();

  // let the targets do whatever work they have to do with the modulated value (typically, 
  // call setter-callbacks):
  for(i = 0; i < size(affectedTargets); i++)
    affectedTargets[i]->doModulationUpdate();
}

void ModulationManager::addConnection(ModulationSource* source, ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  jassert(!isConnected(source, target)); // there is already a connection between source and target
  modulationConnections.push_back(new ModulationConnection(source, target));
  appendIfNotAlreadyThere(affectedTargets, target);

  // we also need a function that removes a target from our affectedTargets array in case it has
  // no incoming connections
}

void ModulationManager::removeConnection(ModulationSource* source, ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  jassert(isConnected(source, target)); // trying to remove no-existent connection
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->source == source && modulationConnections[i]->target == target)
    {
      delete modulationConnections[i];
      remove(modulationConnections, i);
    }
  }
  jassert(!isConnected(source, target)); // there must have been more than one connection between
                                         // given source and target - that should not happen
}

void ModulationManager::removeAllConnections()
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
    delete modulationConnections[i];
  modulationConnections.clear();
}

void ModulationManager::removeConnectionsWith(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->source == source)
    {
      delete modulationConnections[i];
      remove(modulationConnections, i);
      i--; // array was shrunken
    }
  }
}

void ModulationManager::removeConnectionsWith(ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
  {
    if(modulationConnections[i]->target == target)
    {
      delete modulationConnections[i];
      remove(modulationConnections, i);
      i--; // array was shrunken
    }
  }
}

void ModulationManager::registerModulationSource(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
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

bool ModulationManager::isConnected(ModulationSource* source, ModulationTarget* target)
{
  ScopedLock scopedLock(*modLock); 
  for(int i = 0; i < size(modulationConnections); i++)
    if(modulationConnections[i]->source == source && modulationConnections[i]->target == target)
      return true;
  return false;
}

int ModulationManager::numRegisteredSourcesOfType(ModulationSource* source)
{
  ScopedLock scopedLock(*modLock); 
  int result = 0;
  for(int i = 0; i < size(availableSources); i++)
  {
    if(typeid(*source) == typeid(*availableSources[i]))
      result++;
  }
  return result;
}

void ModulationManager::setStateFromXml(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*modLock); 
  jassertfalse; // not yet implemented
}

XmlElement* ModulationManager::getStateAsXml()
{
  ScopedLock scopedLock(*modLock); 

  if(size(modulationConnections) == 0)
    return nullptr; // "Modulations" child-element will be absent from preset file

  XmlElement* xmlState = new XmlElement("Modulations");
  ModulationConnection* c;
  for(int i = 0; i < size(modulationConnections); i++)
  {
    c = modulationConnections[i];
    XmlElement* connectionXml = c->getAsXml();
    xmlState->addChildElement(connectionXml);
    // these 3 lines can be almagamated into one
  }
  return xmlState;
}
