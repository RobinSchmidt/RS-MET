QuadSourceAudioModule::QuadSourceAudioModule(
  CriticalSection *lockToUse, rosic::rsQuadSource* coreToUse)
  : AudioModule(lockToUse), sourceCore(coreToUse)
{

}