
// construction/destruction:

EchoLabDelayLineAudioModule::EchoLabDelayLineAudioModule(CriticalSection *newPlugInLock, rosic::EchoLabDelayLine *echoLabDelayLineToWrap)
: AudioModule(newPlugInLock)
{
  jassert( echoLabDelayLineToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedEchoLabDelayLine = echoLabDelayLineToWrap;
  moduleName = juce::String("DelayLine");
  setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabDelayLinePresets") );
  initializeAutomatableParameters();

  inputEqualizerModule = new EqualizerAudioModule(lock, &echoLabDelayLineToWrap->inputEqualizer);
  inputEqualizerModule->setModuleName(juce::String("InputFilter"));
  inputEqualizerModule->setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabPresets/InputFilterPresets") );
  addChildAudioModule(inputEqualizerModule);

  feedbackEqualizerModule = new EqualizerAudioModule(lock, &echoLabDelayLineToWrap->feedbackEqualizer);
  feedbackEqualizerModule->setModuleName(juce::String("FeedbackFilter"));
  feedbackEqualizerModule->setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabPresets/FeedbackFilterPresets") );
  addChildAudioModule(feedbackEqualizerModule);
}

// setup:

void EchoLabDelayLineAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedEchoLabDelayLine == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: wrappedEchoLabDelayLine->setDelayTime(        value         );   break;
  case 1: wrappedEchoLabDelayLine->setGlobalGainFactor( value         );   break;
  case 2: wrappedEchoLabDelayLine->setFeedbackInPercent(value         );   break;
  case 3: wrappedEchoLabDelayLine->setPan(              value         );   break;
  case 4: wrappedEchoLabDelayLine->setPingPongMode(     value >= 0.5  );   break;
  case 5: wrappedEchoLabDelayLine->setMute(             value >= 0.5  );   break;
  } 

  //sendChangeMessage();
  markStateAsDirty();
}

// audio processing:

/*
void EchoLabDelayLineAudioModule::getSampleFrameStereo(double* inOutL, double* inOutR)
{ 
  *inOutL = 0.0;
  *inOutR = 0.0;
  jassertfalse; // this function is not supposed to be used
}

void EchoLabDelayLineAudioModule::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages) 
{ 
  buffer.clear();
  jassertfalse; // this function is not supposed to be used
} 
*/

// others:

void EchoLabDelayLineAudioModule::initializeAutomatableParameters()
{
  std::vector<double> defaultValues;
  AutomatableParameter* p;

  // #00:
  p = new AutomatableParameter(lock, "DelayTime", 0.0001, 4.25, 0.0001, 1.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.push_back(0.125);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  addObservedParameter( new AutomatableParameter(lock, "Amplitude",  -1.0,  1.0,  0.001, 0.5, Parameter::LINEAR) );  // #01
  addObservedParameter( new AutomatableParameter(lock, "Feedback",  -99.0, 99.0,  0.1,   0.0, Parameter::LINEAR) );  // #02
  addObservedParameter( new AutomatableParameter(lock, "Pan",        -1.0,  1.0,  0.01,  0.0, Parameter::LINEAR) );  // #03
  addObservedParameter( new AutomatableParameter(lock, "PingPong",    0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #04
  addObservedParameter( new AutomatableParameter(lock, "Mute",        0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #05

  //addObservedParameter( new AutomatableParameter(plugInLock, "Solo",        0.0,  1.0,  1.0,   0.0, Parameter::BOOLEAN));  // #05
    // nah - "Solo" is not a parameter - it's a GUI feature - when slo is on, alway the currently selected delayline will play solo

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

/*
void EchoLabDelayLineAudioModule::setStateFromXml(const XmlElement& xmlState, 
const juce::String& stateName, bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);


  int dummy = 0;

}
*/

void EchoLabDelayLineAudioModule::reset()
{
  //wrappedEchoLabDelayLine->reset();  // maybe we should indeed implement this rest-function...
}

//=================================================================================================


EchoLabDelayLineModuleEditor::EchoLabDelayLineModuleEditor(CriticalSection *newPlugInLock, 
  EchoLabDelayLineAudioModule* newEchoLabDelayLineAudioModule) 
  : AudioModuleEditor(newEchoLabDelayLineAudioModule)
{
  setHeadlineStyle(NO_HEADLINE);
  setPresetSectionPosition(INVISIBLE);

  delayLineModuleToEdit = NULL;

  inputEqualizerEditor = new EqualizerModuleEditor(newPlugInLock, NULL);
  inputEqualizerEditor->setLayout(EqualizerModuleEditor::SLIDERS_ABOVE);
  //inputEqualizerEditor->setPresetSectionPosition(AudioModuleEditor::INVISIBLE);
  inputEqualizerEditor->setLinkPosition(AudioModuleEditor::INVISIBLE);
  //inputEqualizerEditor->setHeadlineStyle(SUB_HEADLINE);
  inputEqualizerEditor->setHeadlineText(juce::String("Input Filter"));
  inputEqualizerEditor->setUseShortSliderNames(true);
  inputEqualizerEditor->setUseSmallComboBox(true);
  addChildEditor(inputEqualizerEditor);
  inputEqualizerEditor->setDescriptionField(infoField, true);

  feedbackEqualizerEditor = new EqualizerModuleEditor(newPlugInLock, NULL);
  feedbackEqualizerEditor->setLayout(EqualizerModuleEditor::SLIDERS_ABOVE);
  //feedbackEqualizerEditor->setPresetSectionPosition(AudioModuleEditor::INVISIBLE);
  feedbackEqualizerEditor->setLinkPosition(AudioModuleEditor::INVISIBLE);
  feedbackEqualizerEditor->setUseShortSliderNames(true);
  feedbackEqualizerEditor->setUseSmallComboBox(true);
  feedbackEqualizerEditor->setHeadlineText(juce::String("Feedback Filter"));
  addChildEditor(feedbackEqualizerEditor);
  feedbackEqualizerEditor->setDescriptionField(infoField, true);

  addWidget( timeSlider = new RSlider("TimeSlider") );
  timeSlider->setSliderName(juce::String("Time"));
  timeSlider->setDescription( juce::String("Delaytime for the selected delayline (in seconds or beats)") );
  timeSlider->setRange(0.01, 4.25, 0.01, 0.25);
  timeSlider->setDescriptionField(infoField);
  timeSlider->setStringConversionFunction(secondsToStringWithUnit3);
  //timeSlider->addListener(this);

  addWidget( gainSlider = new RSlider("GainSlider") );
  gainSlider->setSliderName(juce::String("Gain"));
  gainSlider->setDescription( juce::String("Gain for the selected delayline (raw amplitude)") );
  gainSlider->setRange(-1.0, 1.0, 0.01, 0.0);
  gainSlider->setScaling(Parameter::LINEAR_BIPOLAR);
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&valueToString3);
  //gainSlider->addListener(this);

  addWidget( panSlider = new RSlider("PanSlider") );
  panSlider->setSliderName(juce::String("Pan"));
  panSlider->setDescription( juce::String("Panorama position for the selected delayline (-1...+1)") );
  panSlider->setRange(-1.0, 1.0, 0.01, 0.0);
  panSlider->setScaling(Parameter::LINEAR_BIPOLAR);
  panSlider->setDescriptionField(infoField);
  panSlider->setStringConversionFunction(&valueToString2);
  //panSlider->addListener(this);

  addWidget( feedbackSlider = new RSlider("FeedbackSlider") );
  feedbackSlider->setSliderName(juce::String("Feedback"));
  feedbackSlider->setDescription( juce::String("Feedback for the selected delayline in %") );
  feedbackSlider->setRange(-99.0, 99.0, 0.1, 0.0);
  feedbackSlider->setScaling(Parameter::LINEAR_BIPOLAR);
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&percentToStringWithUnit2);
  feedbackSlider->setLayout(RSlider::NAME_ABOVE);
  //feedbackSlider->addListener(this);

  addWidget( pingPongButton = new RButton(juce::String("Ping Pong")) );
  pingPongButton->setDescription(
    juce::String("Switch delayline into ping-pong mode (alternating pan-positions for successive echos)"));
  pingPongButton->setDescriptionField(infoField);
  pingPongButton->addRButtonListener(this);

  addWidget( muteButton = new RButton(juce::String("Mute")) );
  muteButton->setDescription(
    juce::String("Mute the delayline"));
  muteButton->setDescriptionField(infoField);
  muteButton->addRButtonListener(this);

  addWidget( soloButton = new RButton(juce::String("Solo")) );
  soloButton->setDescription(juce::String("Listen to the delayline in solo-mode (mute all others)"));
  soloButton->setDescriptionField(infoField);
  soloButton->addRButtonListener(this);

  addWidget( flushButton = new RButton(juce::String("Flush")) );
  flushButton->setDescription(juce::String("Flush/clear the content of the delaylines"));
  flushButton->setDescriptionField(infoField); 
  flushButton->setClickingTogglesState(false);
  flushButton->addRButtonListener(this);

  guiLayoutRectangles.add(middleRectangle);

  updateWidgetsAccordingToState();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EchoLabDelayLineModuleEditor::setDelayLineModuleToEdit(EchoLabDelayLineAudioModule* newEchoLabDelayLineModuleToEdit)
{
  ScopedLock scopedLock(*lock);

  removeWatchedAudioModule(delayLineModuleToEdit);

  delayLineModuleToEdit = newEchoLabDelayLineModuleToEdit;
  if( delayLineModuleToEdit != NULL )
  {
    addWatchedAudioModule(delayLineModuleToEdit);

    inputEqualizerEditor   ->setEqualizerModuleToEdit(delayLineModuleToEdit->getInputEqualizerModule());
    feedbackEqualizerEditor->setEqualizerModuleToEdit(delayLineModuleToEdit->getFeedbackEqualizerModule());

    timeSlider    ->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("DelayTime")));
    gainSlider    ->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("Amplitude")));
    panSlider     ->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("Pan")));
    feedbackSlider->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("Feedback")));
    pingPongButton->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("PingPong")));
    muteButton    ->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("Mute")));
    soloButton    ->assignParameter(delayLineModuleToEdit->getParameterByName(juce::String("Solo")));

    timeSlider->setSliderName(juce::String("Time"));
    gainSlider->setSliderName(juce::String("Gain"));
  }
  else
  {
    inputEqualizerEditor   ->setEqualizerModuleToEdit(NULL);
    feedbackEqualizerEditor->setEqualizerModuleToEdit(NULL);

    timeSlider    ->assignParameter(NULL);
    gainSlider    ->assignParameter(NULL);
    panSlider     ->assignParameter(NULL);
    feedbackSlider->assignParameter(NULL);
    pingPongButton->assignParameter(NULL);
    muteButton    ->assignParameter(NULL);
    soloButton    ->assignParameter(NULL);
  }

  //updateWidgetVisibility();
  updateWidgetsAccordingToState();
}

void EchoLabDelayLineModuleEditor::setHueOffsetForFilterEditors(float hueOffset)
{
  ScopedLock scopedLock(*lock);
  inputEqualizerEditor->setCentralHue(editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(0));
  feedbackEqualizerEditor->setCentralHue(editorColourScheme.getCentralHue() + editorColourScheme.getHueOffset(0));
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

/*
void EchoLabDelayLineModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
sendChangeMessage(); // triggers update of the plot
}
*/

void EchoLabDelayLineModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);

  AudioModuleEditor::updateWidgetsAccordingToState();
  updateWidgetVisibility();

  if( delayLineModuleToEdit == NULL )
    return;

  if( delayLineModuleToEdit->wrappedEchoLabDelayLine->isInSyncMode() )
    timeSlider->setStringConversionFunction(&beatsToStringWithUnit4);
  else
    timeSlider->setStringConversionFunction(&secondsToStringWithUnitTotal4);
}

void EchoLabDelayLineModuleEditor::audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted)
{
  ScopedLock scopedLock(*lock);
  if( moduleToBeDeleted == delayLineModuleToEdit && delayLineModuleToEdit != NULL )
    setDelayLineModuleToEdit(NULL);
}

void EchoLabDelayLineModuleEditor::paint(Graphics &g)
{
  //g.fillAll(Colours::lavender);

  AudioModuleEditor::paint(g);


  int x = middleRectangle.getX(); // + middleRectangle.getWidth()/2;
  int y = middleRectangle.getY();
  //drawBitmapFontText(g, x+4, y+4, juce::String("Delay-Setup"), &boldFont16px, editorColourScheme.headline);
  drawBitmapFontText(g, x+4, y+4, juce::String("Delay-Setup"), 
    &BitmapFontRoundedBoldA10D0::instance, editorColourScheme.headline);
}


void EchoLabDelayLineModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int middleWidth = 112;
  int x = 0;
  int y = 0;
  int w = (getWidth() - middleWidth) / 2;
  int h = getHeight();

  inputEqualizerEditor   ->setBounds(x,               y, w, h);
  feedbackEqualizerEditor->setBounds(x+w+middleWidth, y, w, h);

  guiLayoutRectangles.clear();
  middleRectangle.setBounds(getWidth()/2-middleWidth/2, y, middleWidth, h);
  guiLayoutRectangles.add(middleRectangle);

  x = middleRectangle.getX();
  y = middleRectangle.getY()+32;
  w = middleRectangle.getWidth();

  timeSlider->setBounds(x+4, y, w-8, 16);
  y += 20;
  gainSlider->setBounds(x+4, y, w-8, 16);
  y += 20;
  panSlider->setBounds(x+4, y, w-8, 16);
  y += 24;
  feedbackSlider->setBounds(x+4, y, w-8, 32);
  y += 48;

  pingPongButton->setBounds(x+16, y, w-32, 16);
  y += 32;
  muteButton->setBounds(x+4, y, w/2-8, 16);
  soloButton->setBounds(x+w/2+4, y, w/2-8, 16);
  y += 20;
  x = muteButton->getX()+muteButton->getWidth()/2;
  w = soloButton->getX()+soloButton->getWidth()/2 - x;
  flushButton->setBounds(x+4, y, w-8, 16);
}

void EchoLabDelayLineModuleEditor::updateWidgetVisibility()
{
  ScopedLock scopedLock(*lock);

  bool visible = true;
  if( delayLineModuleToEdit == NULL )
    visible = false;

  timeSlider    ->setVisible(visible);
  gainSlider    ->setVisible(visible);
  panSlider     ->setVisible(visible);
  feedbackSlider->setVisible(visible);
  pingPongButton->setVisible(visible);
  muteButton    ->setVisible(visible);
  soloButton    ->setVisible(visible);
  flushButton   ->setVisible(visible);
}

