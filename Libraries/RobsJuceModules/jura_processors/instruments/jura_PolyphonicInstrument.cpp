//-------------------------------------------------------------------------------------------------
// construction/destruction:

PolyphonicInstrumentAudioModule::PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock,
  rosic::PolyphonicInstrument *instrumentToWrap) : AudioModuleWithMidiIn(newPlugInLock)
{
  jassert(instrumentToWrap != nullptr); // you must pass a valid rosic-object to the constructor
  underlyingRosicInstrument = instrumentToWrap;
  createParameters();
}

PolyphonicInstrumentAudioModule::PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock)
  : AudioModuleWithMidiIn(newPlugInLock)
{

}

void PolyphonicInstrumentAudioModule::setInstrumentToWrap(rosic::PolyphonicInstrument *instrumentToWrap)
{
  underlyingRosicInstrument = instrumentToWrap;
  createParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void PolyphonicInstrumentAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  AudioModuleWithMidiIn::parameterRangeChanged(parameterThatHasChanged);

  if( underlyingRosicInstrument == NULL )
    return;

  // find out the index in the vector of the parameter that has been changed:
  int parameterIndex = getIndexOfParameter(parameterThatHasChanged);

  // parameterIndex now contains the index in the array of the parameter that has changed now set
  // up the signal processing:

  // \todo replace this dispatching switch statement by the new callback-system

  double value = parameterThatHasChanged->getValue();
  switch( parameterIndex )
  {
  case   0: underlyingRosicInstrument->setMasterLevel(        value);     break;
  case   1: underlyingRosicInstrument->setVoiceLevelByKey(    value);     break;
  case   2: underlyingRosicInstrument->setVoiceLevelByVel(    value);     break;
  case   3: underlyingRosicInstrument->setMasterLevelByVoices(value);     break;
  case   4: underlyingRosicInstrument->setMidSideRatio(       value);     break;
  case   5: underlyingRosicInstrument->setGlideMode(   0.0 != value);     break;
  case   6: underlyingRosicInstrument->setGlideTime(          value);     break;
  case   7: underlyingRosicInstrument->setMasterTuneA4(       value);     break;
  case   8: underlyingRosicInstrument->setPitchWheelRange(    value);     break;
  default:
    {
      // do nothing
    }

  } // end of switch( parameterIndex )

}

void PolyphonicInstrumentAudioModule::setStateFromXml(const XmlElement& xmlState,
  const juce::String& stateName, bool markAsClean)
{
  underlyingRosicInstrument->allNotesOff();
  AudioModuleWithMidiIn::setStateFromXml(xmlState, stateName, markAsClean);
}

/*
XmlElement* PolyphonicInstrumentAudioModule::getStateAsXml(XmlElement* xmlElementToStartFrom)
{
  if( underlyingRosicInstrument == NULL )
    return NULL;

  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String(T("InstrumentState")));
  else
    xmlState = xmlElementToStartFrom;

  // store the setting of the inherited AudioModule object:
  xmlState = AudioModule::getStateAsXml(xmlState);

  return xmlState;
}

void PolyphonicInstrumentAudioModule::setStateFromXml(const XmlElement &xmlState)
{
  if( underlyingRosicInstrument == NULL )
    return;

  // we want to preserve the clean state but when we set automatable parameters it will be set to
  // dirty - so we remember if it was clean and restore it after setting the parameters:
  bool presetIsClean = !underlyingRosicInstrument->isPresetDirty();

  // restore the settings of the inherited AudioModule object:
  AudioModule::setStateFromXml(xmlState);

  if( presetIsClean )
    underlyingRosicInstrument->markPresetAsClean();
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

void PolyphonicInstrumentAudioModule::createParameters()
{
  // Create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particular parameter has changed.
  // Hmm - that comment is old - does the order still matter? I think, it shouldn't. I think, it 
  // mattered back in the day when they were mapped directly to host-automatable parameters. 
  // -> Figure out and update the comment!


  using Param = AutomatableParameter;
  Param* p;

  p = new Param(lock, "MasterLevel", -36.0, 12.0, 0.1, 0.0, Parameter::LINEAR, 7);
  addObservedParameter(p);

  p = new Param(lock, "VoiceLevelByKey", -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new Param(lock, "VoiceLevelByVel", 0.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new Param(lock, "MasterLevelByVoices", 0.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new Param(lock, "MidSideRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);

  p = new Param(lock, "GlideSwitch", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  p = new Param(lock, "GlideTime", 5.0, 2000.0, 0.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Param(lock, "MasterTuneA4", 220.0, 880.0, 0.01, 440.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Param(lock, "PitchWheelRange", 0.0, 24.0, 0.1, 12.0, Parameter::LINEAR);
  addObservedParameter(p);

  // Make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
  // ToDo: Verify, if this is still needed. It might also be a remnant from older days. Ah - I see:
  // we have a method  parameterChanged() change here which then calls the actual setter in the
  // embedded DSP-object. I think, we should instead do things like:
  // using Instrum = rosic::PolyphonicInstrument;
  // p->setValueChangeCallback<Instrum>(instrum, &Instrum::setMasterLevel);
  // ...etc.
  // such that the value change callbacks are called directly. See the AcidDevil code for how it's
  // now supposed to be done.
}

//=================================================================================================


PolyphonicInstrumentEditor::PolyphonicInstrumentEditor(CriticalSection *newPlugInLock,
  PolyphonicInstrumentAudioModule* newInstrumentToEdit) : AudioModuleEditor(newInstrumentToEdit)
{
  ScopedLock scopedLock(*lock);

  // Assign the pointer members:
  jassert( newInstrumentToEdit != nullptr );
  jassert( newInstrumentToEdit->underlyingRosicInstrument != nullptr );
  instrumentEngine = newInstrumentToEdit->underlyingRosicInstrument;

  // Assign the instruments tuning table to the inherited TuningFileManager:
  TuningFileManager::assignTuningTable(&(instrumentEngine->tuningTable));


  // Factor out into createWidgets:

  // Create and setup the widgets:
  addWidget( levelSlider = new RSlider("VolumeSlider") );
  levelSlider->assignParameter(moduleToEdit->getParameterByName("MasterLevel") );
  levelSlider->setSliderName(juce::String("Level"));
  levelSlider->setDescription(juce::String("Master output level"));
  levelSlider->setDescriptionField(infoField);
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( levelByKeySlider = new RSlider("LevelByKeySlider") );
  levelByKeySlider->assignParameter(moduleToEdit->getParameterByName("VoiceLevelByKey") );
  levelByKeySlider->setSliderName(juce::String("Key"));
  levelByKeySlider->setDescription(juce::String("Key dependence of level (per voice)"));
  levelByKeySlider->setDescriptionField(infoField);
  levelByKeySlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( levelByVelSlider = new RSlider("LevelByVelSlider") );
  levelByVelSlider->assignParameter(moduleToEdit->getParameterByName("VoiceLevelByVel") );
  levelByVelSlider->setSliderName(juce::String("Vel"));
  levelByVelSlider->setDescription(juce::String("Velocity dependence of level (per voice)"));
  levelByVelSlider->setDescriptionField(infoField);
  levelByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( midSideRatioSlider = new RSlider("MidSideSlider") );
  midSideRatioSlider->assignParameter(moduleToEdit->getParameterByName("MidSideRatio") );
  midSideRatioSlider->setSliderName(juce::String("M/S"));
  midSideRatioSlider->setDescription("Mid/side adjustment");
  midSideRatioSlider->setDescriptionField(infoField);
  midSideRatioSlider->setStringConversionFunction(&ratioToString0);

  addWidget( numVoicesSlider = new RSlider("NumVoicesSlider") );
  numVoicesSlider->setSliderName(juce::String("Voices"));
  numVoicesSlider->setDescription(juce::String("Maximum number of playing voices"));
  numVoicesSlider->addListener(this);
  numVoicesSlider->setDescriptionField(infoField);
  numVoicesSlider->setStringConversionFunction(&valueToString0);
  numVoicesSlider->setRange(0.0, 16.0, 1.0, 8.0);

  addWidget( compSlider = new RSlider("CompSlider") );
  compSlider->assignParameter(moduleToEdit->getParameterByName("MasterLevelByVoices") );
  compSlider->setSliderName(juce::String("Comp"));
  compSlider->setDescription(juce::String("Compensation for cumulative loudness of playing voices"));
  compSlider->setDescriptionField(infoField);
  compSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( tuningLabel = new RTextField( juce::String("Tuning")) );
  tuningLabel->setDescription(juce::String("Name of current tuning file (if any)"));
  tuningLabel->setDescriptionField(infoField);
  tuningLabel->setNoBackgroundAndOutline(true);

  addWidget( tuningFileNameLabel = new RTextField() );
  tuningFileNameLabel->setDescription(juce::String("Name of current tuning file (if any)"));
  tuningFileNameLabel->setNoBackgroundAndOutline(false);
  tuningFileNameLabel->setDescriptionField(infoField);

  addWidget( tuningLoadButton = new RButton(juce::String("Load")) );
  tuningLoadButton->addRButtonListener(this);
  tuningLoadButton->setDescription(juce::String("Load tuning from a file"));
  tuningLoadButton->setDescriptionField(infoField);
  tuningLoadButton->setClickingTogglesState(false);
  tuningLoadButton->setToggleState(false, false);

  //addWidget( tuningMinusButton = new RButton(RButton::MINUS) );
  addWidget( tuningMinusButton = new RButton(RButton::ARROW_LEFT) );
  tuningMinusButton->addRButtonListener(this);
  tuningMinusButton->setDescription(juce::String("Skip to previous tuning-file in current directory"));
  tuningMinusButton->setDescriptionField(infoField);
  tuningMinusButton->setClickingTogglesState(false);
  tuningMinusButton->setToggleState(false, true);

  //addWidget( tuningPlusButton = new RButton(RButton::PLUS) );
  addWidget( tuningPlusButton = new RButton(RButton::ARROW_RIGHT) );
  tuningPlusButton->addRButtonListener(this);
  tuningPlusButton->setDescription(juce::String("Skip to next tuning-file in current directory"));
  tuningPlusButton->setDescriptionField(infoField);
  tuningPlusButton->setClickingTogglesState(false);
  tuningPlusButton->setToggleState(false, true);

  addWidget( masterTuneSlider = new RSlider("MasterTuneSlider") );
  masterTuneSlider->assignParameter(moduleToEdit->getParameterByName("MasterTuneA4") );
  masterTuneSlider->setSliderName(juce::String("A4"));
  masterTuneSlider->setDescription(juce::String("Master tuning frequency for note A4"));
  masterTuneSlider->setDescriptionField(infoField);
  masterTuneSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( wheelRangeSlider = new RSlider("WheelRangeSlider") );
  wheelRangeSlider->assignParameter(moduleToEdit->getParameterByName("PitchWheelRange") );
  wheelRangeSlider->setSliderName(juce::String("Wheel"));
  wheelRangeSlider->setDescription(juce::String("Range for pitch wheel"));
  wheelRangeSlider->setDescriptionField(infoField);
  wheelRangeSlider->setStringConversionFunction(&semitonesToStringWithUnit1);

  addWidget( glideButton = new RButton(juce::String("Glide")) );
  glideButton->addRButtonListener(this);
  glideButton->setDescription(juce::String("Switch glide on/off"));
  glideButton->setDescriptionField(infoField);
  glideButton->setClickingTogglesState(true);
  glideButton->setToggleState(false, false);

  addWidget( glideTimeSlider = new RSlider("GlideTimeSlider") );
  glideTimeSlider->assignParameter(moduleToEdit->getParameterByName("GlideTime") );
  glideTimeSlider->setSliderName(juce::String("Glide Time"));
  glideTimeSlider->setDescription(juce::String("Adjust the glide time"));
  glideTimeSlider->setDescriptionField(infoField);
  glideTimeSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  updateWidgetsAccordingToState();

}

//-------------------------------------------------------------------------------------------------
// setup:

void PolyphonicInstrumentEditor::setPresetSectionColourScheme(const WidgetColourScheme& newColourScheme)
{
  //stateWidgetSet->setWidgetColourScheme(newColourScheme);
}

void PolyphonicInstrumentEditor::setTuningSectionColourScheme(const WidgetColourScheme& newColourScheme)
{
  /*
  tuningLabel->setColourScheme(newColourScheme);
  tuningFileNameLabel->setColourScheme(newColourScheme);
  tuningLoadButton->setColourScheme(newColourScheme);
  tuningPlusButton->setColourScheme(newColourScheme);
  tuningMinusButton->setColourScheme(newColourScheme);
  // \todo: maybe create also a class TuningWidgetset
  */
}

void PolyphonicInstrumentEditor::setInfoFieldTextColour(const Colour newColour)
{
  Colour tb = Colours::transparentBlack;

  //WidgetColourScheme tmpColourScheme(tb, tb, tb, newColour, tb, tb);
  WidgetColourScheme tmpColourScheme;  // preliminary

  infoField->setColourScheme(tmpColourScheme);
}

//-------------------------------------------------------------------------------------------------
// widget callbacks:

void PolyphonicInstrumentEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  ScopedLock scopedLock(*lock);
  if( instrumentEngine == NULL )
    return;

  if( buttonThatWasClicked == glideButton )
  {
    instrumentEngine->setGlideMode(glideButton->getToggleState());
    moduleToEdit->markStateAsDirty();
  }
  else if( buttonThatWasClicked == tuningLoadButton )
  {
    TuningFileManager::openLoadingDialog();
    TuningFileManager::loadTuningFromFile( TuningFileManager::getActiveFile() );
    tuningFileNameLabel->setText(juce::String(instrumentEngine->tuningTable.getName()));
    moduleToEdit->markStateAsDirty();
  }
  else if( buttonThatWasClicked == tuningPlusButton )
  {
    TuningFileManager::loadNextFile();
    tuningFileNameLabel->setText(juce::String(instrumentEngine->tuningTable.getName()));
    moduleToEdit->markStateAsDirty();
  }
  else if( buttonThatWasClicked == tuningMinusButton )
  {
    TuningFileManager::loadPreviousFile();
    tuningFileNameLabel->setText(juce::String(instrumentEngine->tuningTable.getName()));
    moduleToEdit->markStateAsDirty();
  }
  else
  {
    // it must have been an inherited button:
    AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
  }
}

void PolyphonicInstrumentEditor::rSliderValueChanged(RSlider* sliderThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( instrumentEngine == NULL )
    return;
  if( sliderThatHasChanged == numVoicesSlider )
  {
    instrumentEngine->setNumPlayableVoices( (int) numVoicesSlider->getValue() );
    numVoicesSlider->setValue((double) instrumentEngine->getNumPlayableVoices());
  }
  moduleToEdit->markStateAsDirty();
}

void PolyphonicInstrumentEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);
  if( instrumentEngine == NULL )
    return;

  AudioModuleEditor::updateWidgetsAccordingToState();

  // update tuning widgets:
  masterTuneSlider->setValue( instrumentEngine->getMasterTuneA4(),        false);
  tuningFileNameLabel->setText(juce::String(instrumentEngine->tuningTable.getName()));

  // update global widgets:
  levelSlider->setValue(       instrumentEngine->getMasterLevel(),         false);
  levelByKeySlider->setValue(  instrumentEngine->getVoiceLevelByKey(),     false);
  levelByVelSlider->setValue(  instrumentEngine->getVoiceLevelByVel(),     false);
  midSideRatioSlider->setValue(instrumentEngine->getMidSideRatio(),        false);
  numVoicesSlider->setValue(   instrumentEngine->getNumPlayableVoices(),   false);
  compSlider->setValue(        instrumentEngine->getMasterLevelByVoices(), false);
  wheelRangeSlider->setValue(  instrumentEngine->getPitchWheelRange(),     false);
  glideTimeSlider->setValue(   instrumentEngine->getGlideTime(),           false);
  glideButton->setToggleState( instrumentEngine->isInGlideMode(),          false);

  stateWidgetSet->stateFileNameLabel->setText(moduleToEdit->getStateNameWithStarIfDirty());
  // maybe this call is redundant because it is also called in
  // AudioModuleEditor::updateWidgetsAccordingToState();
}

void PolyphonicInstrumentEditor::setInstrumentToEdit(rosic::PolyphonicInstrument* newInstrumentToEdit)
{
  ScopedLock scopedLock(*lock);
  instrumentEngine = newInstrumentToEdit;
}

void PolyphonicInstrumentEditor::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();

  x = 0;
  y = getHeadlineBottom();
  w = getWidth()/4;

  stateWidgetSet->stateLabel->setBounds(x+4, y+4, 52, 20);

  stateWidgetSet->statePlusButton->setBounds(w-4-20, y+4, 20, 20);
  stateWidgetSet->stateMinusButton->setBounds(stateWidgetSet->statePlusButton->getX()-20, y+4, 20, 20);
  stateWidgetSet->stateLoadButton->setBounds(stateWidgetSet->stateMinusButton->getX()-40-4, y+4, 40, 20);
  stateWidgetSet->stateSaveButton->setBounds(stateWidgetSet->stateLoadButton->getX()-40-4, y+4, 40, 20);

  y = stateWidgetSet->stateLabel->getBottom();
  stateWidgetSet->stateFileNameLabel->setBounds(stateWidgetSet->stateLabel->getX(), y+4, w-8, 20);

  x  = w;
  w /= 2;
  y = getHeadlineBottom();
  levelSlider->setBounds(x+4, y+4, w-8, 20);
  y += 24;
  levelByVelSlider->setBounds(x+4, y+4, w-8, 20);
  y -= 24;
  x += w;
  numVoicesSlider->setBounds(x+4, y+4, w-8, 20);
  y += 24;
  compSlider->setBounds(x+4, y+4, w-8, 20);

  x = getWidth()/2;
  y = getHeadlineBottom();
  w = getWidth()/4;

  tuningLabel->setBounds(x+4, y+4, 80, 20);

  tuningPlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-20, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40-4, y+4, 40, 20);

  y = tuningLabel->getBottom();
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y+4, w-8, 20);

  y = getHeadlineBottom();
  x  = x+w;
  glideButton->setBounds(x+4, y+4, 56, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w-glideButton->getWidth()-12, 20);

  w /= 2;
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w-8, 20);
  x += w;
  wheelRangeSlider->setBounds(x+4, y+4, w-8, 20);
}
