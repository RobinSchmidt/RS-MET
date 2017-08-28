
std::vector<ModulationSource*> ModulationParticipant::dummySources;
std::vector<ModulationTarget*> ModulationParticipant::dummyTargets;
std::vector<ModulationConnection> ModulationParticipant::dummyConnections;

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

const std::vector<ModulationConnection>& ModulationParticipant::getModulationConnections()
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

bool ModulationTarget::isConnectedTo(ModulationSource* source)
{
  if(modManager)
    return modManager->isConnected(source, this);
  return false;
}

std::vector<ModulationSource*> ModulationTarget::getDisconnctedSources()
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

//-------------------------------------------------------------------------------------------------

ModulationConnection::ModulationConnection(ModulationSource* _source, ModulationTarget* _target)
{
  source   = _source; 
  target   = _target;
  amount   = 0.0;
  relative = false;
  sourceValue = &(source->modValue);
  targetValue = &(target->modulatedValue);

  amountParam = new MetaControlledParameter("Amount", -1.0, 1.0, 0.0, Parameter::LINEAR, 0.0);

  //amountParam->setValueChangeCallback(...this, &ModulationConnection::setAmount)
  //amountParam->setMetaParameterManager(metaManager);
  // we may need a MetaParameterManager* parameter to the constructor ..or a pointer to the
  // ModulationManager...which may have to be a subclass of MetaParameterManager
}

ModulationConnection::~ModulationConnection()
{
  delete amountParam;
}

//-------------------------------------------------------------------------------------------------

bool ModulationManager::isConnected(ModulationSource* source, ModulationTarget* target)
{
  for(int i = 0; i < size(modulationConnections); i++)
    if(modulationConnections[i].source == source && modulationConnections[i].target == target)
      return true;
  return false;
}