
PhasorFilter::PhasorFilter(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "PhasorFilter";
  setActiveDirectory(getApplicationDirectory() + "/PhasorFilterPresets");

  // maybe later we should make a convenience class in RAPT that includes the statemapper already
  filterCore.setStateMapper(&stateMapper);

  createStaticParameters();
  inOld = 0;
}

void PhasorFilter::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "Frequency", 20.0, 20000.0, 0.0, 1000.0,  
    Parameter::EXPONENTIAL, 74);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorFilter>(&filterCore, &RAPTPhasorFilter::setFrequency);

  p = new AutomatableParameter(lock, "Decay", 0.0001, 0.1, 0.0, 0.01, 
    Parameter::EXPONENTIAL, 71);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorFilter>(&filterCore, &RAPTPhasorFilter::setDecayTime);

  // nonlinearity parameters:

  p = new AutomatableParameter(lock, "SatInput", 0.1, 100, 0.0, 1.0, 
    Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, 
    &RAPTPhasorMapper::setInputSaturation);

  p = new AutomatableParameter(lock, "Same", -2.0, 2.0, 0.0, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, &RAPTPhasorMapper::setSameSquare);

  p = new AutomatableParameter(lock, "Other", -2.0, 2.0, 0.0, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, &RAPTPhasorMapper::setOtherSquare);

  p = new AutomatableParameter(lock, "Cross", -2.0, 2.0, 0.0, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, &RAPTPhasorMapper::setCrossProduct);

  p = new AutomatableParameter(lock, "Offset", -2.0, 2.0, 0.0, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, &RAPTPhasorMapper::setOffset);

  p = new AutomatableParameter(lock, "SatPreRenorm", 0.1, 100, 0.0, 1.0, 
    Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, 
    &RAPTPhasorMapper::setPreNormalizeSaturation);

  p = new AutomatableParameter(lock, "SatPostRenorm", 0.1, 100, 0.0, 1.0, 
    Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorMapper>(&stateMapper, 
    &RAPTPhasorMapper::setPostNormalizeSaturation);
  // increasing this at high resonance creates a resonant ping when the input signal is turned off



  // more parameters to come..

  //// make sure that the parameters are initially in sync with the audio engine:
  //for(int i = 0; i < (int)observedParameters.size(); i++)
  //  observedParameters[i]->resetToDefaultValue(true, true);
}

void PhasorFilter::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  double in;
  double **b = inOutBuffer;   
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    // preliminary - mix inputs to mono, compute mono result and assign both outputs
    in = 0.5 * (b[0][n] + b[1][n]);

    //b[0][n] = b[1][n] = filterCore.getSample(in);

    b[0][n] = b[1][n] = filterCore.getSample(in - inOld) + in;  // with highpass
    inOld = in; 
  }
}

void PhasorFilter::setSampleRate(double newSampleRate)
{
  filterCore.setSampleRate(newSampleRate); 
}

void PhasorFilter::reset()
{
  filterCore.reset();
  inOld = 0;
}