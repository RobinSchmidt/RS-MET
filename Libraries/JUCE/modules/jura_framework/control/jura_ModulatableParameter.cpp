
const std::vector<ModulationSource*>& ModulationParticipant::getAvailableModulationSources()
{
  jassert(modManager); // must be assigned before its used
  return modManager->getAvailableModulationSources();
}

const std::vector<ModulationTarget*>& ModulationParticipant::getAvailableModulationTargets()
{
  jassert(modManager); // must be assigned before its used
  return modManager->getAvailableModulationTargets();
}

//-------------------------------------------------------------------------------------------------

ModulationSource::~ModulationSource() 
{
  for(int i = 0; i < size(targets); i++)
    targets[i]->removeModulationSource(this);
  if(modManager != nullptr)
    modManager->deRegisterModulationSource(this);
}

