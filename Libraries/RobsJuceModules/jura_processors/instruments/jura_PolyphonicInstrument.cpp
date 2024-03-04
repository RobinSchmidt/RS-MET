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

void PolyphonicInstrumentAudioModule::setStateFromXml(const XmlElement& xmlState,
  const juce::String& stateName, bool markAsClean)
{
  underlyingRosicInstrument->allNotesOff();
  AudioModuleWithMidiIn::setStateFromXml(xmlState, stateName, markAsClean);

  // Maybe an even more drastic resetMidiState (which perhaps needs to be written) would be 
  // appropriate. A call to allNotesOff will trigger releases but that may not be enough. Or will 
  // it? Perhaps it actually does a hrd reset? But then the function should be renamed.
}

//-------------------------------------------------------------------------------------------------
// internal functions:


void PolyphonicInstrumentAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  //using Param = AutomatableParameter;
  using Param = jura::Parameter;
  Param* p;

  using Inst = rosic::PolyphonicInstrument;
  Inst* inst = underlyingRosicInstrument;

  p = new Param("MasterLevel", -36.0, 12.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setMasterLevel);

  p = new Param("VoiceLevelByKey", -24.0, 24.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setVoiceLevelByKey);

  p = new Param("VoiceLevelByVel", 0.0, 12.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setVoiceLevelByVel);

  p = new Param("MasterLevelByVoices", 0.0, 100.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setMasterLevelByVoices);

  p = new Param("MidSideRatio", 0.0, 1.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setMidSideRatio);

  p = new Param("NumVoices", 1.0, 16.0, 16, Parameter::LINEAR, 1.0);
  //p->setCallbackLock(lock);  // might be a good idea to acquire the lock b4 changing NumVoices?
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setNumPlayableVoices);
  // A function p->setCallbackLock(lock) does not yet exist but might be good to add.

  p = new Param("GlideSwitch", 0.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setGlideMode);

  p = new Param("GlideTime", 5.0, 2000.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setGlideTime);

  p = new Param("MasterTuneA4", 220.0, 880.0, 440.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setMasterTuneA4);

  p = new Param("PitchWheelRange", 0.0, 24.0, 12.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<Inst>(inst, &Inst::setPitchWheelRange);
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

  // Create the GUI widgets (sliders, buttons, etc.):
  createWidgets();
}

//-------------------------------------------------------------------------------------------------
// setup:

// Check if we need these - if not, get rid:
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

 
  //// Obsolete:
  //if( buttonThatWasClicked == glideButton )
  //{
  //  instrumentEngine->setGlideMode(glideButton->getToggleState());
  //  moduleToEdit->markStateAsDirty();
  //}


  if( buttonThatWasClicked == tuningLoadButton )
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


  //// Obsolete:
  //if( sliderThatHasChanged == numVoicesSlider )
  //{
  //  instrumentEngine->setNumPlayableVoices( (int) numVoicesSlider->getValue() );
  //  numVoicesSlider->setValue((double) instrumentEngine->getNumPlayableVoices());
  //}


  moduleToEdit->markStateAsDirty();
}

/*
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
*/

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

void PolyphonicInstrumentEditor::createWidgets()
{
  addWidget( levelSlider = new RSlider("VolumeSlider") );
  levelSlider->assignParameter(moduleToEdit->getParameterByName("MasterLevel") );
  levelSlider->setSliderName("Level");
  levelSlider->setDescription("Master output level");
  levelSlider->setDescriptionField(infoField);
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( levelByKeySlider = new RSlider("LevelByKeySlider") );
  levelByKeySlider->assignParameter(moduleToEdit->getParameterByName("VoiceLevelByKey") );
  levelByKeySlider->setSliderName("Key");
  levelByKeySlider->setDescription("Key dependence of level (per voice)");
  levelByKeySlider->setDescriptionField(infoField);
  levelByKeySlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( levelByVelSlider = new RSlider("LevelByVelSlider") );
  levelByVelSlider->assignParameter(moduleToEdit->getParameterByName("VoiceLevelByVel") );
  levelByVelSlider->setSliderName("Vel");
  levelByVelSlider->setDescription("Velocity dependence of level (per voice)");
  levelByVelSlider->setDescriptionField(infoField);
  levelByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( midSideRatioSlider = new RSlider("MidSideSlider") );
  midSideRatioSlider->assignParameter(moduleToEdit->getParameterByName("MidSideRatio") );
  midSideRatioSlider->setSliderName("M/S");
  midSideRatioSlider->setDescription("Mid/side adjustment");
  midSideRatioSlider->setDescriptionField(infoField);
  midSideRatioSlider->setStringConversionFunction(&ratioToString0);


  // New:
  addWidget( numVoicesSlider = new RSlider("NumVoicesSlider") );
  numVoicesSlider->assignParameter(moduleToEdit->getParameterByName("NumVoices") );
  numVoicesSlider->setSliderName("Voices");
  numVoicesSlider->setDescription("Maximum number of playing voices");
  numVoicesSlider->setDescriptionField(infoField);
  numVoicesSlider->setStringConversionFunction(&valueToString0);

  // Old:
  //addWidget( numVoicesSlider = new RSlider("NumVoicesSlider") );
  //numVoicesSlider->setSliderName("Voices");
  //numVoicesSlider->setDescription("Maximum number of playing voices");
  //numVoicesSlider->addListener(this);
  //numVoicesSlider->setDescriptionField(infoField);
  //numVoicesSlider->setStringConversionFunction(&valueToString0);
  //numVoicesSlider->setRange(0.0, 16.0, 1.0, 8.0);


  addWidget( compSlider = new RSlider("CompSlider") );
  compSlider->assignParameter(moduleToEdit->getParameterByName("MasterLevelByVoices") );
  compSlider->setSliderName("Comp");
  compSlider->setDescription("Compensation for cumulative loudness of playing voices");
  compSlider->setDescriptionField(infoField);
  compSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( tuningLabel = new RTextField( "Tuning") );
  tuningLabel->setDescription("Name of current tuning file (if any)");
  tuningLabel->setDescriptionField(infoField);
  tuningLabel->setNoBackgroundAndOutline(true);

  addWidget( tuningFileNameLabel = new RTextField() );
  tuningFileNameLabel->setDescription("Name of current tuning file (if any)");
  tuningFileNameLabel->setNoBackgroundAndOutline(false);
  tuningFileNameLabel->setDescriptionField(infoField);

  addWidget( tuningLoadButton = new RButton("Load") );
  tuningLoadButton->addRButtonListener(this);
  tuningLoadButton->setDescription("Load tuning from a file");
  tuningLoadButton->setDescriptionField(infoField);
  tuningLoadButton->setClickingTogglesState(false);
  tuningLoadButton->setToggleState(false, false);

  addWidget( tuningMinusButton = new RButton(RButton::ARROW_LEFT) );
  tuningMinusButton->addRButtonListener(this);
  tuningMinusButton->setDescription("Skip to previous tuning-file in current directory");
  tuningMinusButton->setDescriptionField(infoField);
  tuningMinusButton->setClickingTogglesState(false);
  tuningMinusButton->setToggleState(false, true);

  addWidget( tuningPlusButton = new RButton(RButton::ARROW_RIGHT) );
  tuningPlusButton->addRButtonListener(this);
  tuningPlusButton->setDescription("Skip to next tuning-file in current directory");
  tuningPlusButton->setDescriptionField(infoField);
  tuningPlusButton->setClickingTogglesState(false);
  tuningPlusButton->setToggleState(false, true);

  addWidget( masterTuneSlider = new RSlider("MasterTuneSlider") );
  masterTuneSlider->assignParameter(moduleToEdit->getParameterByName("MasterTuneA4") );
  masterTuneSlider->setSliderName("A4");
  masterTuneSlider->setDescription("Master tuning frequency for note A4");
  masterTuneSlider->setDescriptionField(infoField);
  masterTuneSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( wheelRangeSlider = new RSlider("WheelRangeSlider") );
  wheelRangeSlider->assignParameter(moduleToEdit->getParameterByName("PitchWheelRange") );
  wheelRangeSlider->setSliderName("Wheel");
  wheelRangeSlider->setDescription("Range for pitch wheel");
  wheelRangeSlider->setDescriptionField(infoField);
  wheelRangeSlider->setStringConversionFunction(&semitonesToStringWithUnit1);



  // New:
  addWidget( glideButton = new RButton("Glide") );
  glideButton->assignParameter(moduleToEdit->getParameterByName("GlideSwitch"));
  glideButton->setDescription("Switch glide on/off");
  glideButton->setDescriptionField(infoField);

  // Old:
  //addWidget( glideButton = new RButton("Glide") );
  //glideButton->addRButtonListener(this);
  //glideButton->setDescription("Switch glide on/off");
  //glideButton->setDescriptionField(infoField);
  //glideButton->setClickingTogglesState(true);
  //glideButton->setToggleState(false, false);



  addWidget( glideTimeSlider = new RSlider("GlideTimeSlider") );
  glideTimeSlider->assignParameter(moduleToEdit->getParameterByName("GlideTime") );
  glideTimeSlider->setSliderName("Glide Time");
  glideTimeSlider->setDescription("Adjust the glide time");
  glideTimeSlider->setDescriptionField(infoField);
  glideTimeSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);

  updateWidgetsAccordingToState(); // ToDo: Check, if this is needed and if so, document why

  // ToDo: 
  // -The glideButton should be connected to the "Glide" parameter
  // -Maybe clean this up - see AciDevilModuleEditor::createWidgets()
  // -In the unit test, filter for the orphaned buttons - we expect 3: tuningLoad, tuningPlus, tuningMinus
}