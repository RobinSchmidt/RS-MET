#include "RPolyphonicInstrumentEditor.h"

RPolyphonicInstrumentEditor::RPolyphonicInstrumentEditor(rosic::PolyphonicInstrument* newInstrumentToEdit) 
: RPlugInEngineEditor(newInstrumentToEdit)
{
  // cast the inherited pointer to a pluginEngine to a pointer to a rosic::PolyphonicInstrument:
  //instrumentEngine = NULL;
  //instrumentEngine = dynamic_cast<rosic::PolyphonicInstrument*>(plugInEngine);
  plugInEngine     = newInstrumentToEdit;
  instrumentEngine = newInstrumentToEdit;

  // assign the instruments tuning table to the inherited TuningFileManager:
  TuningFileManager::assignTuningTable(&(instrumentEngine->tuningTable));

  //---------------------------------------------------------------------------
  // create and setup global widgets:

  addAndMakeVisible( levelSlider = new RSlider(T("VolumeSlider")) );
  levelSlider->addListener(this);
  levelSlider->setSliderName(String(T("Level")));
  levelSlider->setDescription(String(T("Master output level")));
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  levelSlider->setRange(-36.0, 12.0, 0.01, 0.0);
  levelSlider->setDescriptionField(infoField);
  levelSlider->assignParameter(plugInEngine->getAutomatableParameterByName("MasterLevel") );
  levelSlider->setBoxesRelativeWidths(0.5, 0.65);

  addAndMakeVisible( levelByVelSlider = new RSlider(T("LevelByVelSlider")) );
  levelByVelSlider->addListener(this);
  levelByVelSlider->setSliderName(String(T("Vel")));
  levelByVelSlider->setDescription(String(T("Velocity dependence of level (per voice)")));
  levelByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  levelByVelSlider->setRange(0.0, 12.0, 0.1, 0.0);
  levelByVelSlider->setDescriptionField(infoField);
  levelByVelSlider->assignParameter(plugInEngine->getAutomatableParameterByName("VoiceLevelByVel") );
  levelByVelSlider->setBoxesRelativeWidths(0.5, 0.65);

  addAndMakeVisible( numVoicesSlider = new RSlider(T("NumVoicesSlider")) );
  numVoicesSlider->addListener(this);
  numVoicesSlider->setSliderName(String(T("Voices")));
  numVoicesSlider->setDescription(String(T("Maximum number of playing voices")));
  numVoicesSlider->setStringConversionFunction(&valueToString0);
  numVoicesSlider->setRange(0.0, 16.0, 1.0, 8.0);
  numVoicesSlider->setDescriptionField(infoField);

  addAndMakeVisible( compSlider = new RSlider(T("CompSlider")) );
  compSlider->addListener(this);
  compSlider->setSliderName(String(T("Comp")));
  compSlider->setDescription(String(T("Compensation for cumulative loudness of playing voices")));
  compSlider->setStringConversionFunction(&percentToStringWithUnit1);
  compSlider->setRange(0.0, 100.0, 0.1, 0.0);
  compSlider->setDescriptionField(infoField);
  compSlider->assignParameter(plugInEngine->getAutomatableParameterByName("MasterLevelByVoices") );
  compSlider->setBoxesRelativeWidths(0.5, 0.65);

  addAndMakeVisible( tuningLabel = new RLabel(String(T("tuningLabel")), String(T("Tuning:"))) );
  tuningLabel->setEditable(false, false, false);
  tuningLabel->setDescription(String(T("Name of current tuning file (if any)")));
  //tuningLabel->setColour(Label::backgroundColourId, Colours::transparentWhite);
  //tuningLabel->setColour(Label::outlineColourId, Colours::transparentWhite);
  tuningLabel->setJustificationType(Justification::centredLeft);
  tuningLabel->setDescriptionField(infoField);

  addAndMakeVisible( tuningFileNameLabel = new RTextField(String(T("tuningFileNameLabel")), 
    String::empty) );
  tuningFileNameLabel->setEditable(false, false, false);
  tuningFileNameLabel->setDescription(String(T("Name of current tuning file (if any)")));
  tuningFileNameLabel->setDescriptionField(infoField);

  addAndMakeVisible( tuningLoadButton = new RButton(String("Load")) );
  tuningLoadButton->addButtonListener(this);
  tuningLoadButton->setDescription(String(T("Load tuning from a file")));
  tuningLoadButton->setDescriptionField(infoField);
  tuningLoadButton->setClickingTogglesState(false);
  tuningLoadButton->setToggleState(false, false);

  addAndMakeVisible( tuningMinusButton = new RButton(RButton::MINUS) );
  tuningMinusButton->addButtonListener(this);
  tuningMinusButton->setDescription(String(T("Skip to previous tuning-file in current directory")));
  tuningMinusButton->setDescriptionField(infoField);
  tuningMinusButton->setClickingTogglesState(false);
  tuningMinusButton->setToggleState(false, true);

  addAndMakeVisible( tuningPlusButton = new RButton(RButton::PLUS) );
  tuningPlusButton->addButtonListener(this);
  tuningPlusButton->setDescription(String(T("Skip to next tuning-file in current directory")));
  tuningPlusButton->setDescriptionField(infoField);
  tuningPlusButton->setClickingTogglesState(false);
  tuningPlusButton->setToggleState(false, true);

  addAndMakeVisible( masterTuneSlider = new RSlider(T("MasterTuneSlider")) );
  masterTuneSlider->addListener(this);
  masterTuneSlider->setSliderName(String(T("A4")));
  masterTuneSlider->setDescription(String(T("Master tuning frequency for note A4")));
  masterTuneSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  masterTuneSlider->setRange(430.0, 450.0, 0.1, 440.0);
  masterTuneSlider->setDescriptionField(infoField);
  masterTuneSlider->assignParameter(plugInEngine->getAutomatableParameterByName("MasterTuneA4") );
  masterTuneSlider->setBoxesRelativeWidths(0.25, 0.75);
  masterTuneSlider->setScaling(RSlider::LINEAR_BIPOLAR);

  addAndMakeVisible( wheelRangeSlider = new RSlider(T("WheelRangeSlider")) );
  wheelRangeSlider->addListener(this);
  wheelRangeSlider->setSliderName(String(T("Wheel")));
  wheelRangeSlider->setDescription(String(T("Range for pitch wheel")));
  wheelRangeSlider->setStringConversionFunction(&semitonesToStringWithUnit2);
  wheelRangeSlider->setRange(0.0, 24.0, 0.01, 12.0);
  wheelRangeSlider->setDescriptionField(infoField);
  wheelRangeSlider->assignParameter(plugInEngine->getAutomatableParameterByName("PitchWheelRange") );
  wheelRangeSlider->setBoxesRelativeWidths(0.5, 0.65);

  addAndMakeVisible( glideButton = new RButton(String("Glide")) );
  glideButton->addButtonListener(this);
  glideButton->setDescription(String(T("Switch glide on/off")));
  glideButton->setDescriptionField(infoField);
  glideButton->setClickingTogglesState(true);
  glideButton->setToggleState(false, false);

  addAndMakeVisible( glideTimeSlider = new RSlider(T("GlideTimeSlider")) );
  glideTimeSlider->addListener(this);
  glideTimeSlider->setSliderName(String(T("Glide Time")));
  glideTimeSlider->setDescription(String(T("Adjust the glide time")));
  glideTimeSlider->setStringConversionFunction(&millisecondsToStringWithUnit2);
  glideTimeSlider->setRange(10.0, 500.0, 0.01, 50.0);
  glideTimeSlider->assignParameter(plugInEngine->getAutomatableParameterByName("GlideTime") );
  glideTimeSlider->setDescriptionField(infoField);
  //glideTimeSlider->setBoxesRelativeWidths(0.25, 0.75);

  // set up the inherited RobsEditorBase-object:
  enablePresetFileManagement(true);

  // register as listener to change-messages from the plugIn (to reflect automation on the GUI):
  //plugInToEdit->addChangeListener(this);

  //setSize (944, 756);

  updateWidgetsAccordingToState();
}

RPolyphonicInstrumentEditor::~RPolyphonicInstrumentEditor()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// widget callbacks:

void RPolyphonicInstrumentEditor::buttonClicked(Button *buttonThatWasClicked)
{
  if( instrumentEngine == NULL )
    return;

  PresetRemembererEditor::buttonClicked(buttonThatWasClicked);

  if( buttonThatWasClicked == glideButton )
  {
    instrumentEngine->setGlideMode(glideButton->getToggleState());
    setPresetDirty();
  }
  else if( buttonThatWasClicked == tuningLoadButton )
  {
    TuningFileManager::openLoadingDialog();
    TuningFileManager::loadTuningFromFile( TuningFileManager::getCurrentFile() );
    tuningFileNameLabel->setText(String(instrumentEngine->tuningTable.getName()), false);
    setPresetDirty();
  }
  else if( buttonThatWasClicked == tuningPlusButton )
  {
    TuningFileManager::incrementCurrentFile();
    TuningFileManager::loadTuningFromFile( TuningFileManager::getCurrentFile() );
    tuningFileNameLabel->setText(String(instrumentEngine->tuningTable.getName()), false);
    setPresetDirty();
  }
  else if( buttonThatWasClicked == tuningMinusButton )
  {
    TuningFileManager::decrementCurrentFile();
    TuningFileManager::loadTuningFromFile( TuningFileManager::getCurrentFile() );
    tuningFileNameLabel->setText(String(instrumentEngine->tuningTable.getName()), false);
    setPresetDirty();
  }
  else
  {
    // it was none of the button of this class, so it must have been some inherited buttons, 
    // therefore we avoid the call to setDirty() by leaving now (but update the widgets and 
    // preset-field before):
    updateWidgetsAccordingToState();
    updatePresetField();
    return;
  }
}

void RPolyphonicInstrumentEditor::changeListenerCallback (void *objectThatHasChanged)
{
  if( instrumentEngine == NULL )
    return;

  /*
  if( objectThatHasChanged == plugIn )
  {
    // currently, the only events which are of interest are changes of the cutoff-frequency and
    // resonance of the filter - update the respective widgets:
    filterEditor->updateWidgetsAccordingToState();
  }
  */

  setPresetDirty();
}

void RPolyphonicInstrumentEditor::labelTextChanged(Label *labelThatHasChanged)
{
  if( instrumentEngine == NULL )
    return;

  setPresetDirty();
}

void RPolyphonicInstrumentEditor::rSliderValueChanged(RSlider* sliderThatHasChanged)
{
  if( instrumentEngine == NULL )
    return;

  if( sliderThatHasChanged == levelSlider )
    instrumentEngine->setMasterLevel( levelSlider->getValue() );
  else if( sliderThatHasChanged == levelByVelSlider )
    instrumentEngine->setVoiceLevelByVel( levelByVelSlider->getValue() );
  else if( sliderThatHasChanged == numVoicesSlider )
  {
    instrumentEngine->setNumPlayableVoices( (int) numVoicesSlider->getValue() );
    numVoicesSlider->setValue((double) instrumentEngine->getNumPlayableVoices());
  }
  else if( sliderThatHasChanged == compSlider )
    instrumentEngine->setMasterLevelByVoices( compSlider->getValue() );
  else if( sliderThatHasChanged == glideTimeSlider )
    instrumentEngine->setGlideTime( glideTimeSlider->getValue() );
  else if( sliderThatHasChanged == masterTuneSlider )
    instrumentEngine->setMasterTuneA4( masterTuneSlider->getValue() );
  else if( sliderThatHasChanged == wheelRangeSlider )
    instrumentEngine->setPitchWheelRange( wheelRangeSlider->getValue() );

  setPresetDirty();
}

//-------------------------------------------------------------------------------------------------
// state-management:

XmlElement* RPolyphonicInstrumentEditor::getStateAsXml(const String &stateName) const
{
  if( instrumentEngine == NULL )
    return NULL;
  else 
    return polyphonicInstrumentStateToXml(instrumentEngine);

  //return plugIn->getStateAsXml();
}

bool RPolyphonicInstrumentEditor::setStateFromXml(const XmlElement &xmlState)
{
  if( instrumentEngine == NULL )
    return false;
  else
  {
    polyphonicInstrumentStateFromXml(instrumentEngine, xmlState);
    return true;
  }

  //thePlugIn->setStateFromXml(xmlState);

  /*
  // update the preset fields of the sub-editors (their content is obsolete by now):
  oscSectionEditor->setPresetName(String(T("recalled from parent")), true);
  filterEditor->setPresetName(    String(T("recalled from parent")), true);
  filterEnvEditor->setPresetName( String(T("recalled from parent")), true);
  ampEnvEditor->setPresetName(    String(T("recalled from parent")), true);
  */

  updateWidgetsAccordingToState();
  return true;
}

void RPolyphonicInstrumentEditor::updateWidgetsAccordingToState()
{
  if( instrumentEngine == NULL )
    return;
    
  // update tuning widgets:
  masterTuneSlider->setValue( instrumentEngine->getMasterTuneA4(),        false, false);    
  tuningFileNameLabel->setText(String(instrumentEngine->tuningTable.getName()), false);

  // update global widgets:
  levelSlider->setValue(      instrumentEngine->getMasterLevel(),         false, false);
  levelByVelSlider->setValue( instrumentEngine->getVoiceLevelByVel(),     false, false);
  numVoicesSlider->setValue(  instrumentEngine->getNumPlayableVoices(),   false, false);
  compSlider->setValue(       instrumentEngine->getMasterLevelByVoices(), false, false);
  wheelRangeSlider->setValue( instrumentEngine->getPitchWheelRange(),     false, false);
  glideTimeSlider->setValue(  instrumentEngine->getGlideTime(),           false, false);
  glideButton->setToggleState(instrumentEngine->isInGlideMode(),          false);
}

void RPolyphonicInstrumentEditor::setInstrumentToEdit(rosic::PolyphonicInstrument* newInstrumentToEdit)
{
  instrumentEngine = newInstrumentToEdit;
}

void RPolyphonicInstrumentEditor::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  headline->setBounds(0, 0, getWidth(), 20);
  headline->setJustificationType(Justification::centred);
  //webLink->setBounds(4, getBottom()-20, getWidth()-8, 20);
  webLink->setBounds(getWidth()-112, 0, 112-4, 20);

  x = 0;
  y = headline->getBottom();
  w = getWidth()/4;

  presetLabel->setBounds(x+4, y+4, 52, 20);

  presetPlusButton->setBounds(w-4-20, y+4, 20, 20);
  presetMinusButton->setBounds(presetPlusButton->getX()-20, y+4, 20, 20);
  presetLoadButton->setBounds(presetMinusButton->getX()-40-4, y+4, 40, 20);
  presetSaveButton->setBounds(presetLoadButton->getX()-40-4, y+4, 40, 20);

  y = presetLabel->getBottom();
  presetFileNameLabel->setBounds(presetLabel->getX(), y+4, w-8, 20);

  x  = w;
  w /= 2;
  y = headline->getBottom();
  levelSlider->setBounds(x+4, y+4, w-8, 20);
  y += 24;
  levelByVelSlider->setBounds(x+4, y+4, w-8, 20);
  y -= 24;
  x += w;
  numVoicesSlider->setBounds(x+4, y+4, w-8, 20);
  y += 24;
  compSlider->setBounds(x+4, y+4, w-8, 20);

  x = getWidth()/2;
  y = headline->getBottom();
  w = getWidth()/4;

  tuningLabel->setBounds(x+4, y+4, 80, 20);

  tuningPlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-20, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40-4, y+4, 40, 20);
  //presetSaveButton->setBounds(presetLoadButton->getX()-40-4, y+4, 40, 20);

  y = tuningLabel->getBottom();
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y+4, w-8, 20);

  y = headline->getBottom();
  x  = x+w;
  glideButton->setBounds(x+4, y+4, 56, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w-glideButton->getWidth()-12, 20);

  w /= 2;
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w-8, 20);
  x += w;
  wheelRangeSlider->setBounds(x+4, y+4, w-8, 20);

  infoLabel->setBounds(0, getHeight()-20, 40, 20);
  infoField->setBounds(infoLabel->getRight(), getHeight()-20, getWidth()-infoLabel->getRight(),20);

  renderBackgroundImage();
}
