DebugAudioModule::DebugAudioModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("DebugAudioModule");
  createParameters();
}

void DebugAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  outParam = new MetaControlledParameter("Output" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  outParam->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setOutputValue);
  addObservedParameter(outParam);
}

void DebugAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int i = 0; i < numChannels; i++)
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[i][n] = outValue;
  // todo: loop over child modules and let them add their value
}

void DebugAudioModule::processStereoFrame(double *left, double *right)
{
  *left = *right = outValue;
}

void DebugAudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  AudioModule::setSampleRate(newSampleRate);
}

void DebugAudioModule::reset()
{
  ScopedLock scopedLock(*lock);
  AudioModule::reset();
}