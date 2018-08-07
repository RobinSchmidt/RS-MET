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

  p = new Param("Asymmetry", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.80)); // test
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setAsymmetry);

  p = new Param("AttackBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setBending1);

  p = new Param("AttackSigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setSigmoid1);

  p = new Param("DecayBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setMapper(new rsParameterMapperRationalBipolar(-1, +1, 0.9));
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setBending2);

  p = new Param("DecaySigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSO>(tso, &TSO::setSigmoid2);

  // add AttackDrive/DecayDrive

  // for the envelope, add floor/ceil, set it up in terms of AttackTime/DecayTime

  // asymmetry and bedning need finer resolution at the ends / less in the center
  // maybe try tanh mapper or rational s-curve
  // ok, the mappers are currently commented - maybe thei shape paremeters should themselves be
  // availabe as user parameters - at least during development to figure out best (most musical) 
  // settings for them, and later maybe hard-coded - maybe to figure out the best curve, use a 
  // sinusoidal LFO and modulate the parameter with it and choose a setting in which the LFO sounds
  // "most sinusoidal" ...whatever that means


}

void TriSawOscModule::processStereoFrame(double *left, double *right)
{
  *left = *right = oscCore.getSample();
}

void TriSawOscModule::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  oscCore.setPhaseIncrement(freq/sampleRate);
}

void TriSawOscModule::reset()
{
  oscCore.reset();
}

void TriSawOscModule::noteOn(int noteNumber, int velocity)
{
  freq = pitchToFreq(noteNumber);
  oscCore.setPhaseIncrement(freq/sampleRate);
}