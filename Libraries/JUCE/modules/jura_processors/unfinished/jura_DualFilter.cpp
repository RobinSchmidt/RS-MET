DualFilterAudioModule::DualFilterAudioModule(
  CriticalSection *lockToUse, rosic::rsDualFilter* coreToUse)
  : AudioModule(lockToUse), filterCore(coreToUse)
{

}