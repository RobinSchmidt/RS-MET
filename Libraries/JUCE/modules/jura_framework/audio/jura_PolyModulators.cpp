PolyModulatorsAudioModule::PolyModulatorsAudioModule(CriticalSection *lockToUse) 
  : PolySlotAudioModule(lockToUse)
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
  //setPresetSectionPosition(INVISIBLE);
}