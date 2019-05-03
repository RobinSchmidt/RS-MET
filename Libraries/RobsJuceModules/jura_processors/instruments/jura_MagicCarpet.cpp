//=================================================================================================
// class DelayPhaserAudioModule

DelayPhaserAudioModule::DelayPhaserAudioModule(CriticalSection *newPlugInLock,
  rosic::DelayPhaser *newDelayPhaserToWrap) : AudioModule(newPlugInLock)
{
  jassert( newDelayPhaserToWrap != NULL ); // you must pass a valid rosic-object
  wrappedDelayPhaser = newDelayPhaserToWrap;
  setModuleTypeName("DelayPhaser");

  initializeAutomatableParameters();

  phaser1Module = new PhaserAudioModule(lock, &wrappedDelayPhaser->phaser1);
  phaser1Module->setModuleName(juce::String("Phaser1"));
  addChildAudioModule(phaser1Module);

  delayModule = new PingPongEchoAudioModule(lock, &wrappedDelayPhaser->delay);
  delayModule->setModuleName(juce::String(("Delay")));
  addChildAudioModule(delayModule);

  phaser2Module = new PhaserAudioModule(lock, &wrappedDelayPhaser->phaser2);
  phaser2Module->setModuleName(juce::String(("Phaser2")));
  addChildAudioModule(phaser2Module);
}

void DelayPhaserAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedDelayPhaser == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  //case  0: wrappedDelayPhaser->setDryWetRatio(      value); break;
  case  1: wrappedDelayPhaser->setFeedback1(   0.01*value); break;
  case  2: wrappedDelayPhaser->setFeedback2(   0.01*value); break;
  case  3: wrappedDelayPhaser->setFeedback3(   0.01*value); break;
  }
}

void DelayPhaserAudioModule::initializeAutomatableParameters()
{
  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "FeedbackPhaser1Delay", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "FeedbackDelayPhaser2", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "FeedbackGlobal",       -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class MagicCarpetAudioModule

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MagicCarpetAudioModule::MagicCarpetAudioModule(CriticalSection *newPlugInLock)
: PolyphonicInstrumentAudioModule(newPlugInLock)
{
  wrappedMagicCarpet = new rosic::MagicCarpet;
  setInstrumentToWrap(wrappedMagicCarpet);
  setModuleTypeName("MagicCarpet");

  oscSectionModule = new VectorSamplePlayerAudioModule(lock,
    &wrappedMagicCarpet->voiceArray[0].oscSection);
  oscSectionModule->setModuleName(juce::String(("Oscillators")));
  addChildAudioModule(oscSectionModule);

  filterModule = new FourPoleFilterAudioModule(lock,
    &wrappedMagicCarpet->voiceArray[0].filter);
  filterModule->setModuleName(juce::String(("Filter")));
  addChildAudioModule(filterModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock,
    &wrappedMagicCarpet->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String(("FilterEnvelope")));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock,
    &wrappedMagicCarpet->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String(("AmpEnvelope")));
  addChildAudioModule(ampEnvModule);

  equalizerModule = new EqualizerAudioModule(lock,
    &wrappedMagicCarpet->equalizer);
  equalizerModule->setModuleName(juce::String(("Equalizer")));
  addChildAudioModule(equalizerModule);

  delayPhaserModule = new DelayPhaserAudioModule(lock,
    &wrappedMagicCarpet->delayPhaser);
  delayPhaserModule->setModuleName(juce::String(("DelayPhaser")));
  addChildAudioModule(delayPhaserModule);
}

MagicCarpetAudioModule::~MagicCarpetAudioModule()
{
  delete wrappedMagicCarpet;
}

AudioModuleEditor* MagicCarpetAudioModule::createEditor(int type)
{
  return new MagicCarpetModuleEditor(lock, this);
}

//=================================================================================================
// class MagicCarpetFilterEditor:

MagicCarpetFilterEditor::MagicCarpetFilterEditor(CriticalSection *newPlugInLock, FourPoleFilterAudioModule* newFourPoleFilterAudioModule)
  : FourPoleFilterModuleEditor(newPlugInLock, newFourPoleFilterAudioModule)
{
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);

  addWidget( frequencyByKeySlider = new RSlider("FrequencyByKeySlider") );
  //frequencyByKeySlider->assignParameter( moduleToEdit->getParameterByName("FrequencyByKey") );
  //frequencyByKeySlider->setDescription(juce::String("Key dependency of the characteristic frequency"));
  frequencyByKeySlider->setDescription(juce::String("Not yet implemented"));
  frequencyByKeySlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( frequencyByVelSlider = new RSlider(("FrequencyByVelSlider")));
  //frequencyByVelSlider->assignParameter( moduleToEdit->getParameterByName("FrequencyByVel") );
  //frequencyByVelSlider->setDescription(juce::String(("Velocity dependency of the characteristic frequency")));
  frequencyByVelSlider->setDescription(juce::String(("Not yet implemented")));
  frequencyByVelSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( gainByKeySlider = new RSlider (("GainByKeySlider")) );
  //gainByKeySlider->assignParameter( moduleToEdit->getParameterByName("GainByKey") );
  //gainByKeySlider->setDescription(juce::String(("Key dependency of the gain")));
  gainByKeySlider->setDescription(juce::String(("Not yet implemented")));
  gainByKeySlider->setStringConversionFunction(&decibelsToStringWithUnit1);


  addWidget( gainByVelSlider = new RSlider (("GainByVelSlider")) );
  //gainByVelSlider->assignParameter( moduleToEdit->getParameterByName("GainByVel") );
  //gainByVelSlider->setDescription(juce::String(("Velocity dependency of the gain")));
  gainByVelSlider->setDescription(juce::String(("Not yet implemented")));
  gainByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  bandwidthSlider->setSliderName(juce::String(("BW")));
}

void MagicCarpetFilterEditor::resized()
{

  AudioModuleEditor::resized();
  int x = getHeadlineRight();
  int y = 0;
  int w = getWidth()-x;
  //int h = 16;

  modeComboBox->setBounds(x+4, y+4, w-8, 16);

  x = 0;
  y = modeComboBox->getBottom();
  w = getWidth();

  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  x = frequencySlider->getX();
  y = frequencySlider->getBottom();
  w = frequencySlider->getWidth()/2;
  frequencyByKeySlider->setBounds(x,     y-2, w-4, 16);
  frequencyByVelSlider->setBounds(x+w+4, y-2, w-4, 16);

  y = frequencyByKeySlider->getBottom();
  w = 2*getWidth()/3;
  gainSlider->setBounds(0+4, y+4, w-8, 16);
  x = gainSlider->getX();
  y = gainSlider->getBottom();
  w = gainSlider->getWidth()/2;
  gainByKeySlider->setBounds(x,     y-2, w-4, 16);
  gainByVelSlider->setBounds(x+w+4, y-2, w-4, 16);

  x = gainSlider->getRight()+4;
  y = gainSlider->getBottom()-8;
  bandwidthSlider->setBounds(x+4, y, w-8, 16);
}

//=================================================================================================
// class PhaserModuleEditorCompact

PhaserModuleEditorCompact::PhaserModuleEditorCompact(CriticalSection *newPlugInLock, PhaserAudioModule* newPhaserAudioModule)
  : PhaserModuleEditor(newPlugInLock, newPhaserAudioModule)
{
  setHeadlineStyle(AudioModuleEditor::NO_HEADLINE);
  setPresetSectionPosition(AudioModuleEditor::INVISIBLE);

  addWidget( onOffButton = new RButton(juce::String(("OnOffButton"))) );
  //onOffButton->assignParameter( moduleToEdit->getParameterByName(("Activated")) );
  onOffButton->setDescription(juce::String(("Not yet implemented")));
  onOffButton->updateWidgetFromAssignedParameter(false); // shouldn't this happen automatically?

  addWidget( secondOrderButton = new RButton(juce::String(("SecondOrderButton"))) );
  secondOrderButton->assignParameter( moduleToEdit->getParameterByName(("FilterMode")) );
  secondOrderButton->setButtonText(juce::String(("2nd")));
  secondOrderButton->setDescription(juce::String(("Switch allpass chain to 2nd order allpasses")));
  secondOrderButton->addRButtonListener(this); // to en/disable q-slider  (still to do)
}

void PhaserModuleEditorCompact::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();
  onOffButton->setBounds( x+4,     y+4, 60,    16);
  dryWetSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);

  // something to do here....

  /*
  y += 20;
  cycleLengthSlider->setBounds( x+4,     y+4, w/2-8, 16);
  tempoSyncButton->setBounds(cycleLengthSlider->getRight()-40, y+4-14, 40, 16);
  y += 14;
  depthSlider->setBounds(       x+4,     y+4, w/2-8, 16);
  y += 14;
  startPhaseSlider->setBounds(  x+4,     y+4, w/2-8, 16);
  y += 14;
  stereoPhaseSlider->setBounds( x+4,     y+4, w/2-8, 16);

  y = dryWetSlider->getBottom();
  frequencySlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 14;
  feedbackSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 14;
  stagesSlider->setBounds(x+w/2+4, y+4, w/2-8-44, 16);
  y += 14;
  qSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  secondOrderButton->setBounds(qSlider->getRight()-32, y+4-14, 32, 16);
  */
}

//=================================================================================================
// class PingPongEchoModuleEditorCompact

PingPongEchoModuleEditorCompact::PingPongEchoModuleEditorCompact(CriticalSection *newPlugInLock,
  PingPongEchoAudioModule* newPingPongEchoAudioModule)
  : PingPongEchoModuleEditor(newPlugInLock, newPingPongEchoAudioModule)
{
  setHeadlineStyle(AudioModuleEditor::NO_HEADLINE);
  setPresetSectionPosition(AudioModuleEditor::INVISIBLE);

  delayTimeSlider->setSliderName(juce::String(("Time")));

  addWidget( onOffButton = new RButton(juce::String(("OnOffButton"))) );
  onOffButton->assignParameter( moduleToEdit->getParameterByName(("Activated")) );
  onOffButton->setButtonText(juce::String(("Delay")));
  onOffButton->setDescription(juce::String(("Switch delay/echo effect on/off")));
  onOffButton->updateWidgetFromAssignedParameter(false); // shouldn't this happen automatically?
}

void PingPongEchoModuleEditorCompact::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();
  onOffButton->setBounds( x+4,     y+4, 60,    16);
  dryWetSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20;

  delayTimeSlider->setBounds( x+4,     y+4, w/2-8, 16);
  tempoSyncButton->setBounds(delayTimeSlider->getRight()-40, y+4-14, 40, 16);
  y += 14;
  panSlider->setBounds(       x+4,     y+4, w/2-8, 16);
  y += 14;
  pingPongButton->setBounds(  x+4+32,     y+4, w/2-8-64, 16);

  y = dryWetSlider->getBottom();
  feedbackSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 14;
  highDampSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 14;
  lowDampSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
}

//=================================================================================================
// class DelayPhaserModuleEditor

DelayPhaserModuleEditor::DelayPhaserModuleEditor(CriticalSection *newPlugInLock,
  DelayPhaserAudioModule* newDelayPhaserAudioModule)
  : AudioModuleEditor(newDelayPhaserAudioModule)
{
  jassert(newDelayPhaserAudioModule != NULL ); // you must pass a valid module here

  addWidget( dryWetSlider = new RSlider (("DryWetRatioSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&ratioToString0);

  addWidget( feedbackLabel = new RTextField( juce::String(("Feedback"))) );
  feedbackLabel->setJustification(Justification::centred);
  feedbackLabel->setDescription(("Feedback around various parts of the phaser->delay->phaser chain"));
  feedbackLabel->setDescriptionField(infoField);

  addWidget( feedback1Slider = new RSlider (("Feedback1Slider")) );
  feedback1Slider->assignParameter( moduleToEdit->getParameterByName("FeedbackPhaser1Delay") );
  feedback1Slider->setSliderName(juce::String(("Ph1->Dl")));
  feedback1Slider->setDescription(juce::String(("Feedback around first phaser and delay")));
  feedback1Slider->setDescriptionField(infoField);
  feedback1Slider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( feedback2Slider = new RSlider (("Feedback2Slider")) );
  feedback2Slider->assignParameter( moduleToEdit->getParameterByName("FeedbackDelayPhaser2") );
  feedback2Slider->setSliderName(juce::String(("Dl->Ph2")));
  feedback2Slider->setDescription(juce::String(("Feedback around delay and second phaser")));
  feedback2Slider->setDescriptionField(infoField);
  feedback2Slider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( feedback3Slider = new RSlider (("Feedback3Slider")) );
  feedback3Slider->assignParameter( moduleToEdit->getParameterByName("FeedbackGlobal") );
  feedback3Slider->setSliderName(juce::String(("Global")));
  feedback3Slider->setDescription(juce::String(("Feedback around first phaser, delay and second phaser")));
  feedback3Slider->setDescriptionField(infoField);
  feedback3Slider->setStringConversionFunction(&percentToStringWithUnit1);

  phaser1Editor = new PhaserModuleEditorCompact(lock, newDelayPhaserAudioModule->phaser1Module);
  phaser1Editor->setDescription(juce::String(("Settings for the 1st phaser (before the delay)")));
  phaser1Editor->onOffButton->setButtonText(juce::String(("Phaser 1")));
  phaser1Editor->onOffButton->setDescription(juce::String(("Switch 1st phaser (before the delay) on/off")));
  addChildEditor( phaser1Editor );

  delayEditor = new PingPongEchoModuleEditorCompact(lock, newDelayPhaserAudioModule->delayModule);
  delayEditor->setDescription(juce::String(("Settings for the Delay")));
  addChildEditor( delayEditor );

  phaser2Editor = new PhaserModuleEditorCompact(lock, newDelayPhaserAudioModule->phaser2Module);
  phaser2Editor->setDescription(juce::String(("Settings for the 2nd phaser (after the delay)")));
  phaser2Editor->onOffButton->setButtonText(juce::String(("Phaser 2")));
  phaser2Editor->onOffButton->setDescription(juce::String(("Switch 2nd phaser (after the delay) on/off")));
  addChildEditor( phaser2Editor );

  updateWidgetsAccordingToState();
}

void DelayPhaserModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  //int h = getHeight();
  y = getPresetSectionBottom()+4;

  dryWetSlider->setBounds(140, 4, getWidth()-140-4, 16);


  int phaserHeight = 88;
  phaser1Editor->setBounds(x, y, w, phaserHeight);
  y = phaser1Editor->getBottom();

  int delayHeight = 72;
  delayEditor->setBounds(x, y, w, delayHeight);
  y = delayEditor->getBottom();

  phaser2Editor->setBounds(x, y, w, phaserHeight);
  y = phaser2Editor->getBottom();
  feedbackLabel->setBounds(  x    +4, y+4, w/2-8, 16);
  feedback1Slider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20;
  feedback3Slider->setBounds(x    +4, y+4, w/2-8, 16);
  feedback2Slider->setBounds(x+w/2+4, y+4, w/2-8, 16);
}

//=================================================================================================
// class MagicCarpetModuleEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MagicCarpetModuleEditor::MagicCarpetModuleEditor(CriticalSection *newPlugInLock, MagicCarpetAudioModule* newMagicCarpetAudioModule)
  : PolyphonicInstrumentEditor(newPlugInLock, newMagicCarpetAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);

  setHeadlinePosition(AudioModuleEditor::TOP_LEFT);
  setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

  // assign the pointer to the rosic::MagicCarpet object to be used as aduio engine:
  jassert(newMagicCarpetAudioModule != NULL ); // you must pass a valid module here
  magicCarpetAudioModule = newMagicCarpetAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  oscSectionEditor = new VectorSamplePlayerEditor(lock, magicCarpetAudioModule->oscSectionModule);
  oscSectionEditor->setDescription(juce::String(("Settings for the oscillator section")));
  addChildEditor( oscSectionEditor );

  filterEditor = new MagicCarpetFilterEditor(lock, magicCarpetAudioModule->filterModule);
  filterEditor->setDescription(juce::String(("Settings for the filter")));
  addChildEditor( filterEditor );

  filterEnvEditor = new BreakpointModulatorEditorCompact(lock, magicCarpetAudioModule->filterEnvModule);
  filterEnvEditor->setHeadlineText(juce::String(("FilterEnvelope")));
  filterEnvEditor->setPresetSectionPosition(AudioModuleEditor::INVISIBLE);
  filterEnvEditor->setPopUpEditorBounds(-600, -300, 632, 300);
  filterEnvEditor->setDescription(juce::String(("Envelope generator for the filter frequency")));
  //filterEnvEditor->setColourScheme(ColourSchemeComponent::PURPLE, true);
  addChildEditor( filterEnvEditor );

  ampEnvEditor = new BreakpointModulatorEditorCompact(lock, magicCarpetAudioModule->ampEnvModule);
  ampEnvEditor->setHeadlineText(juce::String(("AmpEnvelope")));
  ampEnvEditor->setPopUpEditorBounds(-600, -300, 632, 300);
  ampEnvEditor->setPresetSectionPosition(AudioModuleEditor::INVISIBLE);
  ampEnvEditor->setDescription(juce::String(("Envelope generator for the amplitude")));
  //ampEnvEditor->setColourScheme(ColourSchemeComponent::PURPLE, true);
  addChildEditor( ampEnvEditor );

  equalizerEditor = new EqualizerModuleEditor(lock, magicCarpetAudioModule->equalizerModule);
  equalizerEditor->setDescription(juce::String(("Settings for the master equalizer")));
  //equalizerEditor->setLayout(EqualizerModuleEditor::SLIDERS_ABOVE);
  equalizerEditor->setLayout(EqualizerModuleEditor::SLIDERS_BELOW);
  equalizerEditor->setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);
  equalizerEditor->setHeadlineText(juce::String(("Equalizer")));
  addChildEditor( equalizerEditor );

  delayPhaserEditor = new DelayPhaserModuleEditor(lock, magicCarpetAudioModule->delayPhaserModule);
  delayPhaserEditor->setDescription(juce::String(("Settings for Phaser/Delay/Phaser effect")));
  delayPhaserEditor->setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);
  delayPhaserEditor->setLinkPosition(AudioModuleEditor::INVISIBLE);
  addChildEditor( delayPhaserEditor );

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(800, 600);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void MagicCarpetModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  PolyphonicInstrumentEditor::rButtonClicked(buttonThatWasClicked);
}

void MagicCarpetModuleEditor::changeListenerCallback(ChangeBroadcaster* objectThatHasChanged)
{

}

void MagicCarpetModuleEditor::updateWidgetsAccordingToState()
{
  if( magicCarpetAudioModule == NULL )
    return;
  if( magicCarpetAudioModule->wrappedMagicCarpet == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a
  // dirtification of the preset-state:
  bool presetIsDirty = magicCarpetAudioModule->isStateDirty();
  const MessageManagerLock mmLock;
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  oscSectionEditor->updateWidgetsAccordingToState();
  filterEditor->updateWidgetsAccordingToState();
  filterEnvEditor->updateWidgetsAccordingToState();
  ampEnvEditor->updateWidgetsAccordingToState();
  equalizerEditor->updateWidgetsAccordingToState();
  delayPhaserEditor->updateWidgetsAccordingToState();

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is
  // not enough to make that sure because some of the have AutomatableParameters associated with
  // them which themselves may dirtify the preset:
  if( presetIsDirty )
    magicCarpetAudioModule->markStateAsDirty();
  else
    magicCarpetAudioModule->markStateAsClean();
}


void MagicCarpetModuleEditor::resized()
{

  AudioModuleEditor::resized();
  headlineX = 8;
  headlineY = 8;
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = 0;
  w = getWidth()/5;
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+8, w-8, 40);

  x = w;
  y = stateWidgetSet->getY();
  levelSlider->setBounds(       x+4,     y, w-8,   20);
  y += 18;
  levelByKeySlider->setBounds(  x+4,     y, w/2-8, 16);
  levelByVelSlider->setBounds(  x+w/2+4, y, w/2-8, 16);


  /*
  x = 0;
  y = getHeadlineBottom();
  w = getWidth()/5;

  stateWidgetSet->setBounds(x+4, y+4, w-8, 40);

  x = w;
  y = getHeadlineBottom()+4;
  levelSlider->setBounds(       x+4,     y, w-8,   20);
  y += 18;
  levelByKeySlider->setBounds(  x+4,     y, w/2-8, 16);
  levelByVelSlider->setBounds(  x+w/2+4, y, w/2-8, 16);

  x = 2*w;
  y = getHeadlineBottom()+4;
  midSideRatioSlider->setBounds(x+4, y, w-8, 14);
  y += 16;
  numVoicesSlider->setBounds(   x+4, y, w-8, 14);
  y += 16;
  compSlider->setBounds(        x+4, y, w-8, 14);

  x = 3*w;
  y = getHeadlineBottom();
  glideButton->setBounds(x+4, y+4, 40, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w-glideButton->getWidth()-12, 20);
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w/2-8, 16);
  x += w/2;
  wheelRangeSlider->setBounds(x+4, y+4, w/2-8, 16);

  x = 4*w;
  y = getHeadlineBottom();

  tuningLabel->setBounds(x+4, y+4, 92, 20);
  tuningPlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-18, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40-4, y+4, 40, 20);

  y = tuningLabel->getBottom()-2;
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y, w-8, 20);

  x = 0;
  y = stateWidgetSet->getBottom()+6;
  w = getWidth()/2;
  h = 304;
  */

  x = 0;
  y = getPresetSectionBottom()+4;

  int synthWidth = 660;
  oscSectionEditor->setBounds(x, y, synthWidth, 328);

  x = 0;
  //w = 2*synthWidth/5;
  w = 2*synthWidth/4;
  y = oscSectionEditor->getBottom()+8;
  filterEditor->setBounds(x, y, w, 94);
  x = filterEditor->getRight();
  w = synthWidth-x;
  filterEnvEditor->setBounds(x, y, w, 94);

  y = filterEditor->getBottom();
  w = filterEnvEditor->getWidth();
  h = filterEnvEditor->getHeight();
  ampEnvEditor->setBounds(x, y, w, 94);

  x = filterEditor->getX();
  w = ampEnvEditor->getX()-x;
  //w = samplePlayerTopRightEditor->getRight()-x;
  performanceRect.setBounds(x, y, w, h);

  x = synthWidth;
  w = getWidth()-x;
  y = oscSectionEditor->getY()-12;
  //y = stateWidgetSet->getY();
  h = ampEnvEditor->getBottom()-y;

  effectRect.setBounds(x+8, y, w-8, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(performanceRect);
  guiLayoutRectangles.add(effectRect);

  x = effectRect.getX();
  y = effectRect.getY();
  w = effectRect.getWidth();
  equalizerEditor->setBounds(x, y, w, 200);

  y = equalizerEditor->getBottom();
  h = effectRect.getBottom()-y;
  delayPhaserEditor->setBounds(x, y, w, h);
}
