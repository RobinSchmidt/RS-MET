RotationOscillatorAudioModule::RotationOscillatorAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Oscillator3D");
  createParameters();
}

void RotationOscillatorAudioModule::createParameters()
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


  p = new Param("FreqScaleX", 0.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyScalerX);

  p = new Param("FreqScaleY", 0.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyScalerY);

  p = new Param("FreqScaleZ", 0.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyScalerZ);

  double maxFreqOffset = 20;
  p = new Param("FreqOffsetX", -maxFreqOffset, maxFreqOffset, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyOffsetX);

  p = new Param("FreqOffsetY", -maxFreqOffset, maxFreqOffset, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyOffsetY);

  p = new Param("FreqOffsetZ", -maxFreqOffset, maxFreqOffset, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setFrequencyOffsetZ);

  p = new Param("PhaseX", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setPhaseX);

  p = new Param("PhaseY", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setPhaseY);

  p = new Param("PhaseZ", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setPhaseZ);


  p = new Param("ShiftX", -1, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setShiftX);

  p = new Param("ShiftY", -1, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setShiftY);

  p = new Param("ShiftZ", -1, +1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setShiftZ);


  p = new Param("ScaleX", -2, +2.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setScaleX);

  p = new Param("ScaleY", -2, +2.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setScaleY);

  p = new Param("ScaleZ", -2, +2.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setScaleZ);


  p = new Param("RotX", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setRotationX);

  p = new Param("RotY", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setRotationY);

  p = new Param("RotZ", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setRotationZ);



  p = new Param("OutRotX", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setOutputRotationX);

  p = new Param("OutRotY", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setOutputRotationY);

  p = new Param("OutRotZ", -180, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<LO>(lo, &LO::setOutputRotationZ);
}

void RotationOscillatorAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  //double x, y, z;
  for(int n = 0; n < numSamples; n++)
    oscCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void RotationOscillatorAudioModule::processStereoFrame(double *left, double *right)
{
  oscCore.getSampleFrameStereo(left, right);
}

void RotationOscillatorAudioModule::setSampleRate(double newSampleRate)
{
  oscCore.setSampleRate(newSampleRate);
}

void RotationOscillatorAudioModule::reset()
{
  oscCore.reset();
}

void RotationOscillatorAudioModule::noteOn(int noteNumber, int velocity)
{
  oscCore.setFrequency(RAPT::rsPitchToFreq(noteNumber));
  oscCore.reset();
}