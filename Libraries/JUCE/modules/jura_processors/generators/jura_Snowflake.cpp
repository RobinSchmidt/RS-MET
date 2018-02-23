Snowflake::Snowflake(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Snowflake");
  createParameters();
}

void Snowflake::createParameters()
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
}

void Snowflake::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    oscCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void Snowflake::processStereoFrame(double *left, double *right)
{
  oscCore.getSampleFrameStereo(left, right);
}

void Snowflake::setSampleRate(double newSampleRate)
{
  oscCore.setSampleRate(newSampleRate);
}

void Snowflake::reset()
{
  oscCore.reset();
}

void Snowflake::noteOn(int noteNumber, int velocity)
{
  oscCore.setFrequency(pitchToFreq(noteNumber)); // preliminary - use tuning table
  //oscCore.reset();
}