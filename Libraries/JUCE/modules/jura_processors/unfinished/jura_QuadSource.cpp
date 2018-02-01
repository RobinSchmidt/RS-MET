QuadSourceAudioModule::QuadSourceAudioModule(
  CriticalSection *lockToUse, rosic::rsQuadSource* coreToUse)
  : AudioModule(lockToUse), sourceCore(coreToUse)
{

}

AudioModuleEditor* QuadSourceAudioModule::createEditor()
{
  return new jura::QuadSourceEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

QuadSourceEditor::QuadSourceEditor(CriticalSection* lockToUse, QuadSourceAudioModule* sourceToEdit)
  : AudioModuleEditor(lockToUse), sourceModule(sourceToEdit)
{

}