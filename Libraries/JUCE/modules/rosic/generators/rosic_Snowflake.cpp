void Snowflake::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateIncrement();
}

void Snowflake::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  updateIncrement();
}

void Snowflake::updateWaveTable()
{
  renderer.render(seed, numIterations, x, y);
  updateIncrement();
  reset();
}

void Snowflake::reset()
{
  pos = 0;
}

void Snowflake::updateIncrement()
{

}