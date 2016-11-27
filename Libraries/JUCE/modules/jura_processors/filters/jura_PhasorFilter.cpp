
PhasorFilter::PhasorFilter(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhasorFilter";
  setActiveDirectory(getApplicationDirectory() + "/PhasorFilterPresets");

  createStaticParameters();
}

void PhasorFilter::createStaticParameters()
{
  //ScopedLock scopedLock(*plugInLock);

  //juce::Array<double> defaultValues;
  //AutomatableParameter* p;

  //p = new AutomatableParameter(plugInLock, "Resonance", 0.0, 1.0, 0.0, 0.2,  
  //  Parameter::LINEAR, 71);
  //addObservedParameter(p);
  //p->setValueChangeCallback<Ladder>(this, &Ladder::setResonance);


  //// make sure that the parameters are initially in sync with the audio engine:
  //for(int i = 0; i < (int)observedParameters.size(); i++)
  //  observedParameters[i]->resetToDefaultValue(true, true);
}

AudioModuleEditor* PhasorFilter::createEditor()
{
  return nullptr;
  // todo: get rid of this override and provide a baseclass method that creates a generic editor
  // which just checks the parameters of the module and creates an appropriate widget for each.
  // According to the number of parameters, it may itself decide which size it wants to have.
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