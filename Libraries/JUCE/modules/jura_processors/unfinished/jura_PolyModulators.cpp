PolyModulatorsAudioModule::PolyModulatorsAudioModule(CriticalSection *lockToUse) 
  : AudioModule(lockToUse)
{

}

//=================================================================================================

PolyModulatorsEditor::PolyModulatorsEditor(
  CriticalSection* lockToUse, PolyModulatorsAudioModule* modulatorsToEdit)
  : AudioModuleEditor(lockToUse), modulatorsModule(modulatorsToEdit)
{

}