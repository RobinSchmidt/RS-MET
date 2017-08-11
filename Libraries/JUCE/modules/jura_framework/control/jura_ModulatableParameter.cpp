
ModulationSource::~ModulationSource() 
{
  for(int i = 0; i < size(targets); i++)
    targets[i]->removeModulationSource(this);
  if(modManager != nullptr)
    modManager->deRegisterModulationSource(this);
}

//-------------------------------------------------------------------------------------------------
