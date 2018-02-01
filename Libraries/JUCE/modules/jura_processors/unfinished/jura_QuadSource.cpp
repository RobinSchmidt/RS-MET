QuadSourceAudioModule::QuadSourceAudioModule(
  CriticalSection *lockToUse, rosic::rsQuadSource* coreToUse)
  : AudioModule(lockToUse), sourceCore(coreToUse)
{

}

//=================================================================================================

QuadSourceEditor::QuadSourceEditor(CriticalSection* lockToUse, QuadSourceAudioModule* sourceToEdit)
  : AudioModuleEditor(lockToUse), sourceModule(sourceToEdit)
{

}