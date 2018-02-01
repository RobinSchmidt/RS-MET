DualFilterAudioModule::DualFilterAudioModule(
  CriticalSection *lockToUse, rosic::rsDualFilter* coreToUse)
  : AudioModule(lockToUse), filterCore(coreToUse)
{

}

AudioModuleEditor* DualFilterAudioModule::createEditor()
{
  return new jura::DualFilterEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

DualFilterEditor::DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit)
  : AudioModuleEditor(lockToUse), filterModule(filterToEdit)
{

}