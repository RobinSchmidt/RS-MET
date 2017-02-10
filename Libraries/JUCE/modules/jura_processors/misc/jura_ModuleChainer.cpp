
ModuleChainer::ModuleChainer(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "ModuleChainer";
  setActiveDirectory(getApplicationDirectory() + "/ModuleChainerPresets");
}

AudioModuleEditor* ModuleChainer::createEditor()
{
  return nullptr;
}

void ModuleChainer::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);

  //...
}

void ModuleChainer::setSampleRate(double newSampleRate)
{
  //...
}

void ModuleChainer::noteOn(int noteNumber, int velocity)
{
  //...
}

void ModuleChainer::noteOff(int noteNumber)
{
  //...
}

void ModuleChainer::reset()
{
  //...
}

