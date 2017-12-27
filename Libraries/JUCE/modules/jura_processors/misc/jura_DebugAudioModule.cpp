DebugAudioModule::DebugAudioModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("DebugAudioModule");
  createParameters();
}

void DebugAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);
}

void DebugAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{

}

void DebugAudioModule::setSampleRate(double newSampleRate)
{
 
}

void DebugAudioModule::reset()
{

}