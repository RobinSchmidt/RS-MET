EllipseOscillatorAudioModule::EllipseOscillatorAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("EllipseOscillator");
  createParameters();
}

void EllipseOscillatorAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsEllipseOscillator EO;
  EO* eo = &oscCore;

  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setDetune);

  p = new Param("Shift", 0.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setA);

  p = new Param("Scale", 0.0, 1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setC);

  p = new Param("Rotation", -180.0, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  //p->setValueChangeCallback<EO>(eo, &EO::setRotation);



  // Renormalize, Mix, FreqScaleY, FreqShiftY
}

void EllipseOscillatorAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    oscCore.getSamplePair(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void EllipseOscillatorAudioModule::processStereoFrame(double *left, double *right)
{
  oscCore.getSamplePair(left, right);
}

void EllipseOscillatorAudioModule::setSampleRate(double newSampleRate)
{
  oscCore.setSampleRate(newSampleRate);
}

void EllipseOscillatorAudioModule::reset()
{
  oscCore.reset();
}

void EllipseOscillatorAudioModule::noteOn(int noteNumber, int velocity)
{
  oscCore.setFrequency(pitchToFreq(noteNumber)); // preliminary - use tuning table
  oscCore.reset();
}