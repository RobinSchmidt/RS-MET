
//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombStereoizerAudioModule::CombStereoizerAudioModule(CriticalSection *newPlugInLock,
  rosic::CombStereoizer *stereoizerToWrap) : AudioModule(newPlugInLock)
{
  jassert(stereoizerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedCombStereoizer = stereoizerToWrap;
  setModuleTypeName("CombStereoizer");
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* CombStereoizerAudioModule::getStateAsXml(const juce::String& stateName,
  bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedCombStereoizer == NULL )
    return xmlState;

  // add attributes for the non-automatable parameters (the automatable ones are already taken care
  // of by AudioModule::getStateAsXml()):
  xmlState->setAttribute("ChannelSwitch", wrappedCombStereoizer->isWetSignalChannelSwitched() );

  return xmlState;
}

void CombStereoizerAudioModule::setStateFromXml(const XmlElement& xmlState,
  const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedCombStereoizer == NULL )
    return;

  // restore the values of the non-automatable parameters (the automatable ones are already taken
  // care of by automatableModuleStateFromXml():
  wrappedCombStereoizer->setSwitchWetLeftForRight(
    xmlState.getBoolAttribute("ChannelSwitch", false));
}

//-------------------------------------------------------------------------------------------------
// automation:

void CombStereoizerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedCombStereoizer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedCombStereoizer->setDryWetRatio(                  (float) value);   break;
  case   1: wrappedCombStereoizer->delayLine.setDelayInMilliseconds((float) value);   break;
  case   2: wrappedCombStereoizer->wetFilter.setLowpassCutoff(      (float) value);   break;
  case   3: wrappedCombStereoizer->wetFilter.setHighpassCutoff(     (float) value);   break;
  case   4: wrappedCombStereoizer->wetFilter.setAllpassFrequency(   (float) value);   break;
  } // end of switch( parameterIndex )
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void CombStereoizerAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created
  // Parameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "DelayInMilliseconds", 1.0, 100.0, 0.0, 10.0,
    Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "WetLowpass", 20.0, 20000.0, 0.0, 20000.0,
    Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "WetHighpass", 20.0, 20000.0, 0.0, 20.0,
    Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "WetAllpass", 20.0, 20000.0, 0.0, 20000.0,
    Parameter::EXPONENTIAL);
  addObservedParameter(p);


  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

CombStereoizerModuleEditor::CombStereoizerModuleEditor(CriticalSection *newPlugInLock,
  CombStereoizerAudioModule* newCombStereoizerAudioModule)
  : AudioModuleEditor(newCombStereoizerAudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String("CombStereoizer") );

  // assign the pointer to the rosic::CombStereoizer object to be used as aduio engine:
  jassert(stereoizerAudioModule != NULL ); // you must pass a valid module here
  stereoizerAudioModule = newCombStereoizerAudioModule;

  // create the widgets and assign the automatable parameters to them:
  addWidget( dryWetSlider = new RSlider("DryWetSlider") );
  //dryWetSlider->addListener(this);
  //dryWetSlider->setRange(0.0, 100.0, 0.1, 50.0);
  dryWetSlider->assignParameter( stereoizerAudioModule->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String("Dry/Wet"));
  dryWetSlider->setDescription(juce::String("Ratio between dry and wet signal"));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(ratioToString0);
  //automatableSliders.addIfNotAlreadyThere(dryWetSlider);

  addWidget( delaySlider = new RSlider("DelaySlider"));
  //delaySlider->addListener(this);
  //delaySlider->setRange(1.0, 100.0, 0.01, 10.0);
  //delaySlider->setScaling(Parameter::EXPONENTIAL);
  delaySlider->assignParameter( stereoizerAudioModule->getParameterByName("DelayInMilliseconds") );
  delaySlider->setSliderName(juce::String("Delay"));
  delaySlider->setDescription(juce::String("Delay time for the wet signal"));
  delaySlider->setDescriptionField(infoField);
  delaySlider->setStringConversionFunction(&millisecondsToStringWithUnit2);
  //automatableSliders.addIfNotAlreadyThere(delaySlider);

  addWidget( wetLowpassSlider = new RSlider("WetLowpassSlider"));
  //wetLowpassSlider->addListener(this);
  //wetLowpassSlider->setRange(20.0, 20000.0, 0.001, 20000.0);
  //wetLowpassSlider->setScaling(Parameter::EXPONENTIAL);
  wetLowpassSlider->assignParameter( stereoizerAudioModule->getParameterByName("WetLowpass") );
  wetLowpassSlider->setSliderName(juce::String("Lowpass"));
  wetLowpassSlider->setDescription(juce::String("Lowpass cutoff for the wet signal"));
  wetLowpassSlider->setDescriptionField(infoField);
  wetLowpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(wetLowpassSlider);

  addWidget( wetHighpassSlider = new RSlider("WetHighpassSlider"));
  //wetHighpassSlider->addListener(this);
  //wetHighpassSlider->setRange(20.0, 20000.0, 0.001, 20.0);
  //wetHighpassSlider->setScaling(Parameter::EXPONENTIAL);
  wetHighpassSlider->assignParameter( stereoizerAudioModule->getParameterByName("WetHighpass") );
  wetHighpassSlider->setSliderName(juce::String("Highpass"));
  wetHighpassSlider->setDescription(juce::String("Highpass cutoff for the wet signal"));
  wetHighpassSlider->setDescriptionField(infoField);
  wetHighpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(wetHighpassSlider);

  addWidget( wetAllpassSlider = new RSlider("WetAllpassSlider"));
  //wetAllpassSlider->addListener(this);
  //wetAllpassSlider->setRange(20.0, 20000.0, 0.001, 20000.0);
  //wetAllpassSlider->setScaling(Parameter::EXPONENTIAL);
  wetAllpassSlider->assignParameter( stereoizerAudioModule->getParameterByName("WetAllpass") );
  wetAllpassSlider->setSliderName(juce::String("Allpass"));
  wetAllpassSlider->setDescription(juce::String("Allpass cutoff for the wet signal."));
  wetAllpassSlider->setDescriptionField(infoField);
  wetAllpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(wetAllpassSlider);

  addWidget( swapChannelsButton = new RButton(juce::String("Swap")));
  swapChannelsButton->addRButtonListener(this);
  swapChannelsButton->setDescription(juce::String("Swap channels of wet signal left for right"));
  swapChannelsButton->setDescriptionField(infoField);
  swapChannelsButton->setClickingTogglesState(true);

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void CombStereoizerModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( stereoizerAudioModule == NULL )
    return;
  if( stereoizerAudioModule->wrappedCombStereoizer == NULL )
    return;

  if( buttonThatWasClicked == swapChannelsButton )
    stereoizerAudioModule->wrappedCombStereoizer->setSwitchWetLeftForRight(
      swapChannelsButton->getToggleState());
  /*
  else
  {
  AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
  return;
  }
  */
  stereoizerAudioModule->markStateAsDirty();
}


void CombStereoizerModuleEditor::comboBoxChanged(ComboBox *comboBoxThatHasChanged)
{
  if( stereoizerAudioModule == NULL )
    return;
  if( stereoizerAudioModule->wrappedCombStereoizer == NULL )
    return;

  //....

  stereoizerAudioModule->markStateAsDirty();
}

void CombStereoizerModuleEditor::rSliderValueChanged(RSlider* sliderThatHasChanged)
{
  // \todo: put this inot the base-class
  if( stereoizerAudioModule == NULL )
    return;
  stereoizerAudioModule->markStateAsDirty();
}

void CombStereoizerModuleEditor::updateWidgetsAccordingToState()
{
  if( stereoizerAudioModule == NULL )
    return;
  if( stereoizerAudioModule->wrappedCombStereoizer == NULL )
    return;

  // update the global widgets and automatable sliders:
  AudioModuleEditor::updateWidgetsAccordingToState();

  // update the non-automatable widgets:
  swapChannelsButton->setToggleState(
    stereoizerAudioModule->wrappedCombStereoizer->isWetSignalChannelSwitched(), false);
}

void CombStereoizerModuleEditor::resized()
{
  //linkPosition          = RIGHT_TO_HEADLINE;
  //presetSectionPosition = BELOW_HEADLINE;
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth()/2;
  //int h = getHeight();

  delaySlider->setBounds(x+4, y+4, w-8, 16);

  y = delaySlider->getBottom();
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);

  y = dryWetSlider->getBottom();
  swapChannelsButton->setBounds(x+w/2+4, y+4, w/2-8, 16);

  x = w;
  y = getPresetSectionBottom();

  wetLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  y = wetLowpassSlider->getBottom();
  wetHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y = wetHighpassSlider->getBottom();
  wetAllpassSlider->setBounds(x+4, y+4, w-8, 16);
}


