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

  typedef ModulatableParameter Param;
  Param* p;

  // init to Koch snowflake:
  axiom = "F--F--F--";
  rules = "F = F+F--F+F";

  p = new Param("Iterations", 0, 10, 0, Parameter::INTEGER, 1); // should not be modulatable
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setNumIterations);
  p->setValue(4, true, true); // initial order is 4
  // maybe limit the order by a heuristic estimation of the resulting length of the string and/or
  // the number of points. it should be some kind of exponential growth a * exp(b*order) + c, maybe
  // the a,b,c parameters can be estimated by checking the lengths after the 1st 3 iterations
  // ...or maybe it makes more sense to model it as a * b^order + c, a should be proportional (or
  // equal?) to the number of edges in the initiator and b the number of edges in the generator?
  // ...but that works only for simple one-rule systems that only involve F

  p = new Param("TurningAngle", 0, 360, 0, Parameter::LINEAR); // rename to TurnAngle
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setAngle);
  p->setValue(60, true, true); // Koch snowflake needs 60°

  p = new Param("Amplitude", -1, 1, 1, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setAmplitude);

  p = new Param("Rotation", -180, 180, 0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setRotation);

  p = new Param("CyclicReset", 0, 10, 1, Parameter::INTEGER, 1);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setResetAfterCycles);

  p = new Param("UseTable", 0, 1, 0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<SF>(sf, &SF::setUseTable);

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
  ScopedLock scopedLock(*lock);
  core.setSampleRate(newSampleRate);
}

void Snowflake::reset()
{
  ScopedLock scopedLock(*lock);
  core.reset();
}

void Snowflake::noteOn(int noteNumber, int velocity)
{
  core.setFrequency(pitchToFreq(noteNumber)); // preliminary - use tuning table
  core.reset();
}

void Snowflake::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  core.setNumIterations(0); // avoid excessive re-rendering of complicated tables
  String tmp = xmlState.getStringAttribute("Axiom");
  setAxiom(tmp);
  tmp = xmlState.getStringAttribute("Rules");
  setRules(tmp);

  AudioModuleWithMidiIn::setStateFromXml(xmlState, stateName, markAsClean); 
    // this calls core->setNumIterations which triggers re-rendering...but not always - not when 
    // the number of iterations does not change from one patch to the next, so we trigger 
    // re-rendering (it may result in updating the wavetable twice - maybe try to optimize):
  core.setNumIterations((int)getParameterByName("Iterations")->getValue());
}

XmlElement* Snowflake::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement* xml = AudioModuleWithMidiIn::getStateAsXml(stateName, markAsClean);
  xml->setAttribute("Axiom", axiom);
  xml->setAttribute("Rules", rules);
  return xml;
}

void Snowflake::setAxiom(const juce::String& newAxiom)
{
  ScopedLock scopedLock(*lock);
  if(axiom == newAxiom)
    return;
  axiom = newAxiom;
  core.setAxiom(axiom.toStdString());
  //core.updateWaveTable();
}

bool Snowflake::setRules(const juce::String& newRules)
{
  ScopedLock scopedLock(*lock);
  String tmp = newRules.removeCharacters(" \n"); // remove space and newline
  if(!validateRuleString(tmp))
    return false;

  // Check if new rule-string represents the same set of rules and exit early, if so. It may differ 
  // in terms of spaces, newlines, etc. but that doesn't require re-parsing/re-rendering). We do 
  // this, because this function seems to be called back excessively often.
  bool exitEarly = tmp == rules.removeCharacters(" \n");
  rules = newRules; // we store the rule-string with spaces and newlines here
  if(exitEarly)
    return true; 

  // parse tmp-string and add one rule at a time (maybe move the parsing into rosic::Snowflake):
  core.clearRules();
  bool done = false;
  while(done == false)
  {
    String rule = tmp.upToFirstOccurrenceOf(";", false, false);
    if(rule.length() < 2)
      break;
    jassert(rule[1] == '='); // malformed rule (should be ruled out by validation)
    char input = rule[0];
    String output = rule.substring(2);
    core.addRule(input, output.toStdString());
    tmp = tmp.fromFirstOccurrenceOf(";", false, false);
    done = tmp.length() == 0;
  }
  //core.updateWaveTable();
  return true;
}

bool Snowflake::validateRuleString(const juce::String& newRules)
{
  String tmp = newRules.removeCharacters(" \n"); // remove space and newline
  if(tmp.length() == 0) return true;             // empty string is valid
  if(tmp.length() == 1) return false;            // 1-element string is never valid
  bool done = false;
  while(done == false)
  {
    String rule = tmp.upToFirstOccurrenceOf(";", false, false);
    if(rule[1] != '=')
      return false;
    tmp = tmp.fromFirstOccurrenceOf(";", false, false);
    done = tmp.length() == 0;
  }
  return true;
}

//=================================================================================================

SnowflakeEditor::SnowflakeEditor(jura::Snowflake *snowFlake) : AudioModuleEditor(snowFlake)
{
  ScopedLock scopedLock(*lock);
  snowflakeModule = snowFlake;
  createWidgets();
  setSize(400, 200);
}

void SnowflakeEditor::createWidgets()
{
  ScopedLock scopedLock(*lock);
  typedef RTextField Lbl;
  //typedef RSlider Sld;
  typedef ModulatableSlider Sld;
  Sld* s;
  Lbl* l;
  Parameter* p;

  addWidget( axiomLabel = l = new Lbl("Axiom:") );
  l->setDescription("L-system axiom (aka initiator, seed)");
  l->setDescriptionField(infoField);

  addWidget( axiomEditor = new RTextEditor );
  axiomEditor->setText(snowflakeModule->getAxiom(), false);
  axiomEditor->addListener(this);
  axiomEditor->setDescription(axiomLabel->getDescription());
  axiomEditor->setDescriptionField(infoField);

  addWidget( rulesLabel = l = new Lbl("Rules:") );
  l->setDescription("L-system rules (aka generator(s), productions, semicolon separated)");
  l->setDescriptionField(infoField);

  addWidget( rulesEditor = new RTextEditor );
  rulesEditor->setMultiLine(true);
  rulesEditor->setText(snowflakeModule->getRules(), false);
  rulesEditor->setReturnKeyStartsNewLine(true);
  rulesEditor->addListener(this);
  rulesEditor->setDescription(rulesLabel->getDescription());
  rulesEditor->setDescriptionField(infoField);

  addWidget( sliderIterations = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("Iterations") );
  s->setDescription("Number of L-system iterations");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);

  addWidget( sliderAngle = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("TurningAngle") );
  s->setDescription("Turning angle of turtle graphics");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString2);


  addWidget( sliderAmplitude = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("Amplitude") );
  s->setDescription("Amplitude of output signal");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString2);

  addWidget( sliderRotation = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("Rotation") );
  s->setDescription("Rotation applied to xy coordinates");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString2);

  addWidget( sliderReset = s = new Sld );
  s->assignParameter( p = snowflakeModule->getParameterByName("CyclicReset") );
  s->setDescription("Number of cycles after which turtle resets (0 for never)");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);
}

void SnowflakeEditor::resized()
{
  AudioModuleEditor::resized();
  int m  = 4; // margin
  int x  = m;
  int y  = getPresetSectionBottom() + m;
  int w  = getWidth() - 2*m;
  int h  = getHeight();
  int wh = 16;    // widget height
  int dy = wh+2;  // delta-y between widgets

  int lw = 48; // label width
  axiomLabel->setBounds( x,    y, lw,   wh);
  axiomEditor->setBounds(x+lw, y, w-lw, wh);
  y += dy;

  rulesLabel->setBounds( x,    y, lw,   3*wh);
  rulesEditor->setBounds(x+lw, y, w-lw, 3*wh);
  y += 3*dy;

  // put initiator and generator plots here

  sliderIterations->setBounds(x, y, w, wh); y += dy;
  sliderAngle->setBounds(x, y, w, wh); y += dy;

  // put result (2D and 1D plots here)

  sliderAmplitude->setBounds(x, y, w, wh); y += dy;
  sliderRotation->setBounds(x, y, w, wh); y += dy;
  sliderReset->setBounds(x, y, w, wh); y += dy;
}

void SnowflakeEditor::rTextEditorTextChanged(RTextEditor& ed)
{
  if(&ed == axiomEditor)
    snowflakeModule->setAxiom(ed.getText());
  else if(&ed == rulesEditor)
  {
    bool rulesValid = snowflakeModule->setRules(ed.getText());
    //if(!rulesValid)
    //  ed->indicateInvalidText(true);
    //else
    //  ed->indicateInvalidText(false);
  }
}

void SnowflakeEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  axiomEditor->setText(snowflakeModule->getAxiom(), false);
  rulesEditor->setText(snowflakeModule->getRules(), false);
}