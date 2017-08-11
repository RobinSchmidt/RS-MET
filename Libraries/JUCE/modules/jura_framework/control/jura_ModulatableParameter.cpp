
ModulationSource::~ModulationSource() 
{
  for(int i = 0; i < size(targets); i++)
    targets[i]->removeModulationSource(this);
  if(modManager != nullptr)
    modManager->deRegisterModulationSource(this);
}

//-------------------------------------------------------------------------------------------------

const std::vector<ModulationSource*>* ModulationTarget::getAvailableSources() 
{ 
  if(modManager != nullptr)
    return modManager->getAvailableSources();
  return nullptr;
}
