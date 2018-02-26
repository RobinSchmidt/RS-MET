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
  incUpToDate = false;
  //updateIncrement();
}

void Snowflake::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  incUpToDate = false;
  //updateIncrement();
}

void Snowflake::setRotation(double newRotation)
{
  rotator.setAngle((PI/180) * newRotation);
}

void Snowflake::setNumIterations(int newNumIterations) 
{ 
  if(newNumIterations == numIterations)
    return;
  numIterations = newNumIterations; 
  tableUpToDate = false; 
  //updateWaveTable();
  // maybe don't update here, instead set a "dirty" flag, check it in getSample and update the
  // table there
}

void Snowflake::setAngle(double newAngle) 
{ 
  renderer.setAngle(newAngle); 
  tableUpToDate = false; 
  //updateWaveTable();
  // todo: make the angle modulatable - don't re-render the table but only render the L-system 
  // output string, reduce it to the relevant characters and render one point at a time instead of
  // pre-rendering a wavetable
}

void Snowflake::clearRules() 
{ 
  renderer.clearRules(); 
  tableUpToDate = false; 
}

void Snowflake::Snowflake::addRule(char input, const std::string& output) 
{ 
  renderer.addRule(input, output); 
  tableUpToDate = false; 
}

void Snowflake::setAxiom(const std::string& newAxiom) 
{ 
  axiom = newAxiom; 
  tableUpToDate = false; 
}

void Snowflake::updateWaveTable()
{
  renderer.render(axiom, numIterations, x, y);
  tableLength = (int)x.size()-1;
  tableUpToDate = true;
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
  incUpToDate = true;
}