
PhasorFilter::PhasorFilter(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhasorFilter";
  setActiveDirectory(getApplicationDirectory() + "/PhasorFilterPresets");

  createStaticParameters();
}

void PhasorFilter::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0,  
    Parameter::EXPONENTIAL, 74);
  addObservedParameter(p);
  p->setValueChangeCallback<RAPTPhasorFilter>(&filterCore, &RAPTPhasorFilter::setFrequency);

  // more parameters to come..

  //// make sure that the parameters are initially in sync with the audio engine:
  //for(int i = 0; i < (int)observedParameters.size(); i++)
  //  observedParameters[i]->resetToDefaultValue(true, true);
}

void PhasorFilter::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  double **b = inOutBuffer;   
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    // preliminary - mix inputs to mono, compute mono result and assign both outputs
    b[0][n] = b[1][n] = filterCore.getSample( 0.5 * (b[0][n] + b[1][n]) );
  }
}

void PhasorFilter::setSampleRate(double newSampleRate)
{
  filterCore.setSampleRate(newSampleRate); 
}

void PhasorFilter::reset()
{
  filterCore.reset();
}