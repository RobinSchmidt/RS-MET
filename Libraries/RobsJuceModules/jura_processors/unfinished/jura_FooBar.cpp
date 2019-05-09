FooBarModule::FooBarModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("FooBar");
  setModuleName("FooBar1");
}

AudioModuleEditor* FooBarModule::createEditor(int type)
{
  return new jura::FooBarEditor(this);
}

void FooBarModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    fooBarCore->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void FooBarModule::processStereoFrame(double *left, double *right)
{
  fooBarCore->getSampleFrameStereo(left, right);
}

void FooBarModule::setSampleRate(double newSampleRate)
{ 
  fooBarCore->setSampleRate(newSampleRate); 
}

void FooBarModule::reset()
{ 
  fooBarCore->reset(); 
}

//=================================================================================================

FooBarEditor::FooBarEditor(FooBarModule* fooBarToEdit)
  : AudioModuleEditor(fooBarToEdit->lock), fooBarModule(fooBarToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...
}

void FooBarEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}