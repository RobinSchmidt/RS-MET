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

  p = new Param("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setAmplitude);

  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setDetune);

  //p = new Param("A", -1.0, 1.0, 0.0, Parameter::LINEAR);
  p = new Param("LowHigh", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setA);
  // determines, how the speed varies over the full cycle: 0 constant speed, 
  // +-1: skips half a cycle (i think)

  p = new Param("InOut", -10.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setC);
  // determines, how speed varies over a half-cycles?

  p = new Param("Phase", -180.0, +180.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setRotationDegrees);
  // adjusts horizontal vs vertical extent of ellipse - at 45°, its a circle


  // A is Upper/Lower Bias, B is Left/Right Bias, C is Upper/Lower Pinch



  // Renormalize, Mix, FreqScaleY, FreqShiftY

  // we really need some rescaling of LowHigh and InOut
}

void EllipseOscillatorAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    oscCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void EllipseOscillatorAudioModule::processStereoFrame(double *left, double *right)
{
  oscCore.getSampleFrameStereo(left, right);
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
  //oscCore.reset();
}

//=================================================================================================

TriSawOscModule::TriSawOscModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("TriSawOscillator");
  createParameters();
}


void TriSawOscModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef RAPT::rsTriSawOscillator<double> TSO;
  TSO* tso = &oscCore;

  typedef ModulatableParameter Param;
  Param* p;

  /*
  p = new Param("Amplitude", -1.0, +1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setAmplitude);

  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<EO>(eo, &EO::setDetune);
  */
}