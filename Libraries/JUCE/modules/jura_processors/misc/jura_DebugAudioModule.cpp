DebugAudioModule::DebugAudioModule(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("DebugAudioModule");
  createParameters();
}

void DebugAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  leftParam = new MetaControlledParameter("Left" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  leftParam->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setLeftValue);
  addObservedParameter(leftParam);

  rightParam = new MetaControlledParameter("Right" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  rightParam->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setRightValue);
  addObservedParameter(rightParam);
}

void DebugAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int i = 0; i < numChannels; i++)
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[i][n] = values[i];
  // todo: loop over child modules and let them add their value
}

void DebugAudioModule::processStereoFrame(double *left, double *right)
{
  *left  = values[0];
  *right = values[1];;
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

void DebugAudioModule::setMidiController(int controllerNumber, float controllerValue)
{
  double v = controllerValue / 127.0; //  0..1
  v = 2*v-1;                          // -1..+1
  if(controllerNumber == 74)
    getParameterByName("Left")->setValue(v, true, true);
  else if(controllerNumber == 71)
    getParameterByName("Right")->setValue(v, true, true);
}