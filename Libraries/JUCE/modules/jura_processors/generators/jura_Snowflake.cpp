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

  // init to Koch snowflake:
  seed  = "F--F--F";
  rules = "F = F+F--F+F";

  p = new Param("Iterations", 0, 10, 0, Parameter::INTEGER);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setNumIterations);
  p->setValue(4, true, true); // initial order is 4
  // maybe limit the order by a heuristic estimation of the resulting length of the string and/or
  // the number of points. it should be some kind of exponential growth a * exp(b*order) + c, maybe
  // the a,b,c parameters can be estimated by checking the lengths after the 1st 3 iterations

  p = new Param("Angle", 0, 360, 0, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setAngle);
  p->setValue(60, true, true); // Koch snowflake needs 60°


  /*
  p = new Param("Tune", -60.0, +60.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setDetune);
  */
}

AudioModuleEditor* Snowflake::createEditor()
{
  return new SnowflakeEditor(this);
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

void Snowflake::setSeed(const juce::String& newSeed)
{
  seed = newSeed;
  // translate to std::string and pass to core
}

bool Snowflake::setRules(const juce::String& newRules)
{

  return true;
}

bool Snowflake::validateRuleString(const juce::String& newRules)
{
  return true;
}

//=================================================================================================

SnowflakeEditor::SnowflakeEditor(jura::Snowflake *snowFlake) : AudioModuleEditor(snowFlake)
{
  ScopedLock scopedLock(*lock);
  snowflakeModule = snowFlake;
  createWidgets();
  setSize(300, 200);
}

void SnowflakeEditor::createWidgets()
{
  typedef RSlider Sld;
  Sld* s;
  Parameter* p;

  addWidget( sliderIterations = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("Iterations") );
  //s->setSliderName("Iterations");
  s->setDescription("Number of L-system iterations");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);

  addWidget( sliderAngle = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("Angle") );
  //s->setSliderName("Angle");
  s->setDescription("Turning angle of turtle graphics");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString2);
}

void SnowflakeEditor::resized()
{
  AudioModuleEditor::resized();
  int m  = 4; // margin
  int x  = m;
  int y  = getPresetSectionBottom() + m;
  int w  = getWidth() - 2*m;
  int h  = getHeight();
  int wh = 16;   // widget height
  int dy = wh-2; // delta-y between widgets

  sliderIterations->setBounds(x, y, w, wh); y += dy;
  sliderAngle->setBounds(x, y, w, wh); y += dy;
}