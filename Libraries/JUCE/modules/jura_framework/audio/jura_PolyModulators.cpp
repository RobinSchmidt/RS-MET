PolyModulatorsAudioModule::PolyModulatorsAudioModule(CriticalSection *lockToUse) 
  : AudioModule(lockToUse)
{

}

AudioModuleEditor* PolyModulatorsAudioModule::createEditor()
{
  return new jura::PolyModulatorsEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

PolyModulatorsEditor::PolyModulatorsEditor(
  CriticalSection* lockToUse, PolyModulatorsAudioModule* modulatorsToEdit)
  : AudioModuleEditor(lockToUse), modulatorsModule(modulatorsToEdit)
{

}