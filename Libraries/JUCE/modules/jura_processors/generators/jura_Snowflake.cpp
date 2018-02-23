Snowflake::Snowflake(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Snowflake");
  createParameters();
}

void Snowflake::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::Snowflake SF;
  SF* sf = &core;

  typedef Parameter Param;
  Param* p;

  p = new Param("Iterations", 0, 10, 0, Parameter::INTEGER);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setNumIterations);

  /*
  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setDetune);
  */
}

void Snowflake::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    core.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}
void Snowflake::processStereoFrame(double *left, double *right)
{
  core.getSampleFrameStereo(left, right);
}

void Snowflake::setSampleRate(double newSampleRate)
{
  core.setSampleRate(newSampleRate);
}

void Snowflake::reset()
{
  core.reset();
}

void Snowflake::noteOn(int noteNumber, int velocity)
{
  core.setFrequency(pitchToFreq(noteNumber)); // preliminary - use tuning table
  //oscCore.reset();
}