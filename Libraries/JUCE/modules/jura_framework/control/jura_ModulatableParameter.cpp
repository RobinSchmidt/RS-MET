
const std::vector<ModulationSource*>& ModulationParticipant::getAvailableModulationSources()
{
  jassert(modManager);
  return modManager->getAvailableModulationSources();
}

const std::vector<ModulationTarget*>& ModulationParticipant::getAvailableModulationTargets()
{
  jassert(modManager);
  return modManager->getAvailableModulationTargets();
}

void ModulationParticipant::registerModulationSource(ModulationSource* source)
{
  jassert(modManager);
  modManager->registerModulationSource(source);
}

void ModulationParticipant::deRegisterModulationSource(ModulationSource* source)
{
  jassert(modManager);
  modManager->deRegisterModulationSource(source);
}

void ModulationParticipant::registerModulationTarget(ModulationTarget* target)
{
  jassert(modManager);
  modManager->registerModulationTarget(target);
}

void ModulationParticipant::deRegisterModulationTarget(ModulationTarget* target)
{
  jassert(modManager);
  modManager->deRegisterModulationTarget(target);
}

//-------------------------------------------------------------------------------------------------

ModulationSource::~ModulationSource() 
{
  for(int i = 0; i < size(targets); i++)
    targets[i]->removeModulationSource(this);
  if(modManager != nullptr)
    modManager->deRegisterModulationSource(this);
}

