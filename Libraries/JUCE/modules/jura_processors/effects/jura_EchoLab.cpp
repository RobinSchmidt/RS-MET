
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

//=================================================================================================

// construction/destruction:

EchoLabAudioModule::EchoLabAudioModule(CriticalSection *newPlugInLock, rosic::EchoLab *delayDesignerToWrap)
  : AudioModule(newPlugInLock)
{
  jassert(delayDesignerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedEchoLab = delayDesignerToWrap;
  moduleName = juce::String("EchoLab");
  setActiveDirectory(getApplicationDirectory() + juce::String("/EchoLabPresets") );

  //inputFilterModule = new EqualizerAudioModule(NULL);
  //inputFilterModule->setModuleName(juce::String(T("InputFilter")));

  //feedbackFilterModule = new EqualizerAudioModule(NULL);
  //feedbackFilterModule->setModuleName(juce::String(T("FeedbackFilter")));

  initializeAutomatableParameters();
}

// setup:

void EchoLabAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  //if( wrappedEchoLab == NULL )
  //  return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedEchoLab->setDryWet(  value); break;
  case   1: wrappedEchoLab->setWetLevel(value); break;
  } // end of switch( parameterIndex )
}

XmlElement* equalizerStateToXml(Equalizer* equalizer, XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the releveant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("Equalizer")); 
  else
    xmlState = xmlElementToStartFrom;

  xmlState->setAttribute("GlobalGain", equalizer->getGlobalGain());

  // create an XmlElement for each band and add it as child-XmlElement:
  for(int i=0; i<equalizer->getNumBands(); i++)
  {
    XmlElement* bandState = new XmlElement(juce::String("Band"));

    bandState->setAttribute(juce::String("Frequency"), equalizer->getBandFrequency(i));
    bandState->setAttribute(juce::String("Gain"),      equalizer->getBandGain(i));
    bandState->setAttribute(juce::String("Bandwidth"), equalizer->getBandBandwidth(i));

    juce::String modeString;
    int mode = equalizer->getBandMode(i);
    switch( mode )
    {
    case TwoPoleFilter::PEAK:       modeString = juce::String("Peak/Dip");           break;
    case TwoPoleFilter::LOW_SHELF:  modeString = juce::String("Low Shelving");       break;
    case TwoPoleFilter::HIGH_SHELF: modeString = juce::String("High Shelving");      break;
    case TwoPoleFilter::LOWPASS6:   modeString = juce::String("Lowpass 6 dB/oct");   break;
    case TwoPoleFilter::LOWPASS12:  modeString = juce::String("Lowpass 12 dB/oct");  break;
    case TwoPoleFilter::HIGHPASS6:  modeString = juce::String("Highpass 6 dB/oct");  break;
    case TwoPoleFilter::HIGHPASS12: modeString = juce::String("Highpass 12 dB/oct"); break;
    case TwoPoleFilter::BANDREJECT: modeString = juce::String("Notch 2*6 dB/oct");   break;
    }
    bandState->setAttribute("Mode", modeString);

    xmlState->addChildElement(bandState);
  } 

  return xmlState;
}

XmlElement* echoLabDelayLineStateToXml(EchoLabDelayLine* delayLine, XmlElement* xmlElementToStartFrom)
{
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("DelayLine")); 
  else
    xmlState = xmlElementToStartFrom;

  // todo: store value only if different from default values (save space)
  xmlState->setAttribute("DelayTime", delayLine->getDelayTime()              );
  xmlState->setAttribute("Amplitude", delayLine->getGlobalGainFactor()       );
  xmlState->setAttribute("Feedback",  delayLine->getFeedbackInPercent()      );
  xmlState->setAttribute("Pan",       delayLine->getPan()                    );
  xmlState->setAttribute("PingPong",  delayLine->isInPingPongMode()          );
  xmlState->setAttribute("Mute",      delayLine->isMuted()                   );

  // store the embedded equalizer's state as child element (unless it's neutral):
  if( delayLine->feedbackEqualizer.getNumBands(0) > 0 || delayLine->feedbackEqualizer.getGlobalGain() != 0.0 )
  {
    XmlElement* feedbackEqState = new XmlElement(juce::String("FeedbackFilter"));
    feedbackEqState             = equalizerStateToXml(&(delayLine->feedbackEqualizer.equalizers[0]), feedbackEqState);
    xmlState->addChildElement(feedbackEqState);
  }
  if( delayLine->inputEqualizer.getNumBands(0) > 0 || delayLine->inputEqualizer.getGlobalGain() != 0.0 )
  {
    XmlElement* inputEqState = new XmlElement(juce::String("InputFilter"));
    inputEqState             = equalizerStateToXml(&(delayLine->inputEqualizer.equalizers[0]), inputEqState);
    xmlState->addChildElement(inputEqState);
  }

  return xmlState;
}

XmlElement* echoLabStateToXml(EchoLab *echoLab, XmlElement* xmlElementToStartFrom)
{
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("EchoLabState")); 
  else
    xmlState = xmlElementToStartFrom;

  // create an XmlElement for each delayline and add it as child-XmlElement:
  echoLab->acquireLock();
  xmlState->setAttribute("SyncDelayTimes", echoLab->isDelayTimeSynced());
  xmlState->setAttribute("DryWet",         echoLab->getDryWet() );
  xmlState->setAttribute("WetLevel",       echoLab->getWetLevel() );
  xmlState->setAttribute("Solo",           echoLab->getSoloedDelayLineIndex() );
  for(int i=0; i<echoLab->getNumDelayLines(); i++)
  {
    rosic::EchoLabDelayLine* delayLine = echoLab->getDelayLine(i);
    XmlElement* delayLineState = echoLabDelayLineStateToXml(delayLine, NULL);
    xmlState->addChildElement(delayLineState);
  } 
  echoLab->releaseLock();

  return xmlState;
}

XmlElement* EchoLabAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedEchoLab != NULL )
  {
    wrappedEchoLab->acquireLock();
    xmlState = echoLabStateToXml(wrappedEchoLab, xmlState);
    wrappedEchoLab->releaseLock();
  }

  // mmm...actually, we may get away with the baseclass implementation, i think

  return xmlState;
}

void EchoLabAudioModule::setStateFromXml(const XmlElement& xmlState, 
  const juce::String& stateName, bool markAsClean)
{
  jassertfalse;

  // this function must be re-implemented
  // we have to do more here - we need to create and set up all the required child-modules
  // ahh - and quite possibly we need legacy-preset conversion (but maybe not)

  // perhaps we just need to look at echoLabStateFromXml, define the right static parameters 
  // inside EchoLsbDelayLineAudioModule, create the child-modules here and then rely on the 
  // recursion to recall all the parameters for the individual delaylines

  // we perhaps need to drag over the EchoLabStateFromXml function from the old codebase 
  // first

  ScopedLock scopedLock(*lock);

  removeAllDelayLines();

  // create the audiomodules for the delaylines:
  for(int i = 0; i < xmlState.getNumChildElements(); i++)
  {
    if( xmlState.getChildElement(i)->hasTagName(juce::String("DelayLine")) )
    {
      juce::XmlElement *delayLineState = xmlState.getChildElement(i);
      double delayTime = delayLineState->getDoubleAttribute(juce::String("DelayTime"), 1.0);
      double amplitude = delayLineState->getDoubleAttribute(juce::String("Amplitude"), 0.5);
      addDelayLine(delayTime, amplitude);
    }
  }

  // at this point, the delaytimes are still correct

  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  // should set up the internal states of the individual delaylines (which are now child-modules)
  // -> we must perhaps override setStateFromXml in EchoLabDelayLineAudioModule
  // ...this seems not yet to work....



  /*
  // this is the old implementation:
  ScopedLock scopedLock(*plugInLock);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  if( wrappedEchoLab != NULL )
  {
  wrappedEchoLab->acquireLock();
  echoLabStateFromXml(wrappedEchoLab, xmlState);
  wrappedEchoLab->releaseLock();
  }
  */
}

/*
XmlElement EchoLabAudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
ScopedLock scopedLock(*plugInLock);

DEBUG_BREAK; 
// \todo: implement this function - should take care of updating the sub-states for the two filters (the new EqualizerAudioModule
// admits 2 channels
// hmmmm....may it's better to just define a mono version and stereo version of EqualizerAudioModule
// ...or maybe we can adapt the state save/recall in the EqualizerAudioModule so as to admit for both versions
// one with and one without the "Channel" child element

int xmlPatchFormatIndex = xmlState.getIntAttribute(T("PatchFormat"), patchFormatIndex);

return AudioModule::convertXmlStateIfNecessary(xmlState);
}
*/

void EchoLabAudioModule::setSampleRate(double newSampleRate)   
{ 
  ScopedLock scopedLock(*lock);
  wrappedEchoLab->setSampleRate(newSampleRate); 
}

void EchoLabAudioModule::setBeatsPerMinute(double newBpm)   
{ 
  ScopedLock scopedLock(*lock);
  wrappedEchoLab->setTempoInBPM(newBpm); 
}

int EchoLabAudioModule::addDelayLine(double newDelayTime, double newGainFactor)
{
  ScopedLock scopedLock(*lock);

  int index = wrappedEchoLab->addDelayLine(newDelayTime, newGainFactor);

  if( index != -1 )
  {
    EchoLabDelayLineAudioModule *newDelayLineModule 
      = new EchoLabDelayLineAudioModule(lock, wrappedEchoLab->getDelayLine(index));
    newDelayLineModule->getParameterByName(
      juce::String("DelayTime"))->setValue(newDelayTime,  true, true);
    newDelayLineModule->getParameterByName(
      juce::String("Amplitude"))->setValue(newGainFactor, true, true);
    delayLineModules.add(newDelayLineModule);
    addChildAudioModule(newDelayLineModule);
  }

  return index;
}

bool EchoLabAudioModule::removeDelayLine(int index)
{
  ScopedLock scopedLock(*lock);

  bool success = false;

  //if( index < 0 || index >= delayLineModules.size() )  // huh? shouldn't it be index >= 0 && index < size?
  if( index >= 0 && index < delayLineModules.size() )
  {
    EchoLabDelayLineAudioModule *delayLineModule = delayLineModules[index];
    delayLineModules.remove(index);
    removeChildAudioModule(delayLineModule, true);

    success = wrappedEchoLab->removeDelayLine(index);
    jassert( success  == true   ); 
    // we successfully removed the EqualizerAudioModules but not the underlying rosic-delayline object - something must have gone wrong
  }
  else
    jassertfalse; // invalid index

  return success;
}

void EchoLabAudioModule::removeAllDelayLines()
{
  ScopedLock scopedLock(*lock);
  while( delayLineModules.size() > 0 )
    removeDelayLine(delayLineModules.size()-1);
}

// inquiry:

EchoLabDelayLineAudioModule* EchoLabAudioModule::getDelayLineModule(int index) const
{
  jassert( index >= 0  && index < delayLineModules.size() );
  return delayLineModules[index];
}

// audio processing:

void EchoLabAudioModule::getSampleFrameStereo(double* inOutL, double* inOutR)
{ 
  jassertfalse; // nor good idea to use this function - performance hog
  ScopedLock scopedLock(*lock);
  wrappedEchoLab->getSampleFrameStereo(inOutL, inOutR); 
}

void EchoLabAudioModule::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages) 
{ 
  ScopedLock scopedLock(*lock);

  if( wrappedEchoLab == NULL || buffer.getNumChannels() < 1 )
  {
    jassertfalse;
    return;
  }

  float *left, *right;  
  left  = buffer.getWritePointer(0, 0);
  if( buffer.getNumChannels() < 2 )
    right = buffer.getWritePointer(0, 0);
  else
    right = buffer.getWritePointer(1, 0);
  wrappedEchoLab->processBlock(left, right, buffer.getNumSamples());
} 

// others:

void EchoLabAudioModule::reset() 
{ 
  ScopedLock scopedLock(*lock);
  wrappedEchoLab->resetDelayLines(); 
}

void EchoLabAudioModule::initializeAutomatableParameters()
{
  ScopedLock scopedLock(*lock);

  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  juce::Array<double> defaultValues;

  // this pointer will be used to temporarily store the addresses of the created Parameter-objects:
  AutomatableParameter* p;

  // #00:
  p = new AutomatableParameter(lock, "DryWet", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR); 
  addObservedParameter(p);

  // #01:
  p = new AutomatableParameter(lock, "WetLevel", -36.0, 6.0, 0.01, 0.0, Parameter::LINEAR); 
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}


