QuadSourceAudioModule::QuadSourceAudioModule(
  CriticalSection *lockToUse/*, rosic::rsQuadSource* coreToUse*/)
  : AudioModule(lockToUse)/*, sourceCore(coreToUse)*/
{
  addChildAudioModule(topLeftModule     = new PolySlotAudioModule(lock));
  addChildAudioModule(topRightModule    = new PolySlotAudioModule(lock));
  addChildAudioModule(bottomLeftModule  = new PolySlotAudioModule(lock));
  addChildAudioModule(bottomRightModule = new PolySlotAudioModule(lock));
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