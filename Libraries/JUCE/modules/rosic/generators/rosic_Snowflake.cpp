Snowflake::Snowflake()
{
  // init to order 4 Koch snowflake:
  seed = "F--F--F";
  clearRules();
  addRule('F', "F+F--F+F");
  numIterations = 4;
  updateWaveTable();
}

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

void Snowflake::setNumIterations(int newNumIterations) 
{ 
  numIterations = newNumIterations; 
  updateWaveTable(); 
}

void Snowflake::updateWaveTable()
{
  renderer.render(seed, numIterations, x, y);
  tableLength = x.size()-1;
  updateIncrement();
  reset();
}

void Snowflake::reset()
{
  pos = 0;
}

void Snowflake::updateIncrement()
{
  inc = tableLength * frequency / sampleRate; // avoid division
}