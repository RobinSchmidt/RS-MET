DualFilterAudioModule::DualFilterAudioModule(
  CriticalSection *lockToUse, rosic::rsDualFilter* coreToUse)
  : AudioModule(lockToUse), filterCore(coreToUse)
{

}

//=================================================================================================

DualFilterEditor::DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit)
  : AudioModuleEditor(lockToUse), filterModule(filterToEdit)
{

}