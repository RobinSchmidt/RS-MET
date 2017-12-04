MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "MultiComp";
  setActiveDirectory(getApplicationDirectory() + "/Presets/MultiComp");

  createParameters();
}

void MultiCompAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef RAPT::rsLissajousOscillator3D<double> LO;
  LO* lo = &oscCore;

  typedef ModulatableParameter Param;
  Param* p;


  p = new Param("Renormalize", 0.0, 2.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setRenormalizationAmount);

  p = new Param("Clip", 0.0, 2.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setClipAmplitude);
}

void MultiCompAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    multiCompCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void MultiCompAudioModule::setSampleRate(double newSampleRate)
{
  multiCompCore.setSampleRate(newSampleRate);
}

void MultiCompAudioModule::reset()
{
  multiCompCore.reset();
}
