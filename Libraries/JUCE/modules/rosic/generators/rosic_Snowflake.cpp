Snowflake::Snowflake()
{
  // init to order 4 Koch snowflake:
  axiom = "F--F--F";
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

void Snowflake::setRotation(double newRotation)
{
  rotator.setAngle((PI/180) * newRotation);
}

void Snowflake::setNumIterations(int newNumIterations) 
{ 
  numIterations = newNumIterations; 
  updateWaveTable();
}

void Snowflake::setAngle(double newAngle) 
{ 
  renderer.setAngle(newAngle); 
  updateWaveTable();
  // todo: make the angle modulatable - don't re-render the table but only render the L-system 
  // output string, reduce it to the relevant characters and render one point at a time instead of
  // pre-rendering a wavetable
}

void Snowflake::updateWaveTable()
{
  renderer.render(axiom, numIterations, x, y);
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