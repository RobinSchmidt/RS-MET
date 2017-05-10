
//-------------------------------------------------------------------------------------------------
// construction/destruction:

PitchShifterAudioModule::PitchShifterAudioModule(CriticalSection *newPlugInLock, 
  rosic::PitchShifterGrainAdaptive *pitchShifterToWrap) : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*lock);
  //jassert(pitchShifterToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if(pitchShifterToWrap != nullptr)
    wrappedPitchShifter = pitchShifterToWrap;
  else
  {
    wrappedPitchShifter = new rosic::PitchShifterGrainAdaptive;
    wrappedPitchShifterIsOwned = true;
  }

  moduleName = juce::String("PitchShifter");
  setActiveDirectory(getApplicationDirectory() + juce::String("/PitchShifterPresets"));
  createStaticParameters();
}

PitchShifterAudioModule::~PitchShifterAudioModule()
{
  if(wrappedPitchShifterIsOwned)
    delete wrappedPitchShifter;
}

AudioModuleEditor* PitchShifterAudioModule::createEditor()
{
  return new PitchShifterModuleEditor(lock, this);
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void PitchShifterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  std::vector<double> defaultValues;
  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "DetuneCoarse", -48.0, 48.0, 0.1, 0.0, 
    Parameter::LINEAR_BIPOLAR, 74);
  defaultValues.clear();
  defaultValues.push_back(2.03910002); // f/f0 =  9/ 8 = 1.125
  defaultValues.push_back(3.15641287); // f/f0         = 1.2
  defaultValues.push_back(3.86313714); // f/f0 =  5/ 4 = 1.25
  //defaultValues.push_back(4.54213948); // f/f0 = 13/10 = 1.3
  defaultValues.push_back(4.98044999); // f/f0 =  4/ 3 = 1.333...
  defaultValues.push_back(7.01955001); // f/f0 =  3/ 2 = 1.5
  defaultValues.push_back(8.84358713); // f/f0 =  5/ 3 = 1.666...
  defaultValues.push_back(9.68825906); // f/f0 =  7/ 4 = 1.75
  defaultValues.push_back(12.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDetuneCoarse);
    // p->setValueChangeCallback(wrappedPitchShifter, &PitchShifterGrainAdaptive::setDetuneCoarse); does not work because the MS compiler 
    // seems not to be able to infer, that the template argument is "PitchShifterGrainAdaptive" - it ambiguates it with the "PitchShifter" 
    // baseclass - that's why we must pass the template argument here explicitly

  p = new AutomatableParameter(lock, "DetuneFine", -200.0, 200.0, 0.1, 0.0, 
    Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDetuneFine);

  p = new AutomatableParameter(lock, "GrainLengthInMilliseconds", 1.0, 2000.0, 0.001, 
    18.1818181818181818, Parameter::EXPONENTIAL);
     // a grain length of 18.18... will tune the amplitude modulation artifacts to +-55 Hz with 
     // respect to the carrier frequency
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInMilliseconds);

  p = new AutomatableParameter(lock, "GrainLengthInPitchCycles", 0.25, 64.0, 0.01, 4.0, 
    Parameter::EXPONENTIAL);
  defaultValues.clear();
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(8.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(16.0);
  defaultValues.push_back(24.0);
  defaultValues.push_back(32.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInPitchCycles);

  p = new AutomatableParameter(lock, "GrainLengthInBeats", 0.125, 4.0, 0.001, 0.5, 
    Parameter::EXPONENTIAL);
  defaultValues.clear();
  defaultValues.push_back(0.125);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.75);
  defaultValues.push_back(1.0);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthInBeats);

  p = new AutomatableParameter(lock, "GrainLengthUnit", 0.0, 2.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String("ms"));
  p->addStringValue(juce::String("cycles"));
  p->addStringValue(juce::String("beats"));
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setGrainLengthUnit);

  p = new AutomatableParameter(lock, "Feedback", -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setFeedback);

  p = new AutomatableParameter(lock, "DryWet", 0.0, 100.0, 0.1, 100.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setDryWet);

  p = new AutomatableParameter(lock, "AntiAlias", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setAntiAliasing);

  p = new AutomatableParameter(lock, "Reverse", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setReversePlayback);

  p = new AutomatableParameter(lock, "Invert", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setNegativePolarity);

  p = new AutomatableParameter(lock, "FormantPreserve", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<PitchShifterGrainAdaptive>(wrappedPitchShifter, 
    &PitchShifterGrainAdaptive::setFormantPreserve);

  // make sure that the parameters are initially in sync with the audio engine:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameters[i]->resetToDefaultValue(true, true);
}

XmlElement PitchShifterAudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);

  // retrieve the patch format of the xml-file to enable different interpretations of the patch for 
  // backwards compatibility:
  int xmlPatchFormat = xmlState.getIntAttribute("PatchFormat", 0);
  if( xmlPatchFormat == 0 ) // this is an old preset
  {
    // we had formerly only one "GrainLength" parameter which was either interpreted as being in ms 
    // or in pitch-cycles depending on a boolean flag "GrainLengthAdaption" - now we have 
    // different parameters for the different units and a string-parameter "GrainLengthUnit" to 
    // select which value is to be used
    XmlElement convertedState = xmlState;
    double d = xmlState.getDoubleAttribute("GrainLength", 8.0);
    bool   b = xmlState.getBoolAttribute("GrainLengthAdaption", false);
    if( b == true )
    {
      convertedState.setAttribute("GrainLengthInPitchCycles", d);
      convertedState.setAttribute("GrainLengthUnit", "cycles");
    }
    else
    {
      convertedState.setAttribute("GrainLengthInMilliseconds", d);
      convertedState.setAttribute("GrainLengthUnit", "ms");
    }
    return convertedState;
  }
  else
    return xmlState;
}

//=================================================================================================


PitchShifterModuleEditor::PitchShifterModuleEditor(CriticalSection *newPlugInLock, PitchShifterAudioModule* newPitchShifterAudioModule) 
  : AudioModuleEditor(newPitchShifterAudioModule)
{
  ScopedLock scopedLock(*lock);
  // maybe we should avoid this lock here and instead have a function that connects the widgets with the parameters where we acquire
  // the lock - but maybe not

  // set the plugIn-headline:
  setHeadlineText( juce::String("PitchShifter") );

  // assign the pointer to the rosic::PitchShifter object to be used as aduio engine:
  jassert(newPitchShifterAudioModule != NULL ); // you must pass a valid module here
  pitchShifterModuleToEdit = newPitchShifterAudioModule;

  // create the widgets and assign the automatable parameters to them:
  addWidget( coarseSlider = new RSlider("CoarseSlider") );
  coarseSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName("DetuneCoarse") );
  coarseSlider->setSliderName(juce::String("Coarse"));
  coarseSlider->setDescription(juce::String("Coarse pitch shifting factor in semitones"));
  coarseSlider->setDescriptionField(infoField);
  coarseSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( fineSlider = new RSlider("FineSlider") );
  fineSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName("DetuneFine") );
  fineSlider->setSliderName(juce::String("Fine"));
  fineSlider->setDescription(juce::String("Fine pitch shifting factor in cents"));
  fineSlider->setDescriptionField(infoField);
  fineSlider->setStringConversionFunction(&centsToStringWithUnit2);

  addWidget( grainLengthInMillisecondsSlider = new RSlider("GrainLengthSlider") );
  grainLengthInMillisecondsSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName("GrainLengthInMilliseconds") );
  grainLengthInMillisecondsSlider->setSliderName(juce::String("Grain Length"));
  grainLengthInMillisecondsSlider->setDescription(juce::String("Length of the grains in milliseconds"));
  grainLengthInMillisecondsSlider->setDescriptionField(infoField);
  grainLengthInMillisecondsSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthInCyclesSlider = new RSlider("CyclesPerGrainSlider"));
  grainLengthInCyclesSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName("GrainLengthInPitchCycles") );
  grainLengthInCyclesSlider->setSliderName(juce::String("Grain Length"));
  grainLengthInCyclesSlider->setDescription(juce::String("Length of the grains in pitch cylces"));
  grainLengthInCyclesSlider->setDescriptionField(infoField);
  grainLengthInCyclesSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthInBeatsSlider = new RSlider("GrainLengthInBeatsSlider"));
  grainLengthInBeatsSlider->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName("GrainLengthInBeats"));
  grainLengthInBeatsSlider->setSliderName(juce::String("Grain Length"));
  grainLengthInBeatsSlider->setDescription(juce::String("Length of the grains in beats"));
  grainLengthInBeatsSlider->setDescriptionField(infoField);
  grainLengthInBeatsSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( grainLengthUnitComboBox = new RComboBox(juce::String("GrainLengthUnitComboBox")));
  grainLengthUnitComboBox->assignParameter( 
    pitchShifterModuleToEdit->getParameterByName("GrainLengthUnit"));
  grainLengthUnitComboBox->setDescription("Choose the unit for the grain length");
  grainLengthUnitComboBox->setDescriptionField(infoField);
  grainLengthUnitComboBox->registerComboBoxObserver(this); // to update visibility of the sliders

  addWidget( feedbackSlider = new RSlider("FeedbackSlider"));
  feedbackSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName("Feedback"));
  feedbackSlider->setSliderName(juce::String("Feedback"));
  feedbackSlider->setDescription(juce::String("Feeds the pitch-shifted output back to the input"));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider("DryWet"));
  dryWetSlider->assignParameter( pitchShifterModuleToEdit->getParameterByName("DryWet"));
  dryWetSlider->setSliderName(juce::String("Dry/Wet"));
  dryWetSlider->setDescription(juce::String("Ratio between dry and wet signal (in % wet)"));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( antiAliasButton = new RButton(juce::String("Anti-Alias")));
  antiAliasButton->assignParameter( pitchShifterModuleToEdit->getParameterByName("AntiAlias"));
  antiAliasButton->setDescription(juce::String("Switch anti-alias filter (for up-shifting) on/off"));
  antiAliasButton->setDescriptionField(infoField);
  antiAliasButton->setClickingTogglesState(true);

  addWidget( reverseButton = new RButton(juce::String("Reverse")));
  reverseButton->assignParameter( pitchShifterModuleToEdit->getParameterByName("Reverse"));
  reverseButton->setDescription(juce::String("Reverse playback of the grains"));
  reverseButton->setDescriptionField(infoField);
  reverseButton->setClickingTogglesState(true);

  addWidget( invertButton = new RButton(juce::String("Invert")));
  invertButton->assignParameter( pitchShifterModuleToEdit->getParameterByName("Invert"));
  invertButton->setDescription(juce::String("Invert polarity of wet (shifted) signal"));
  invertButton->setDescriptionField(infoField);
  invertButton->setClickingTogglesState(true);

  /*
  addWidget( formantPreserveButton = new RButton(juce::String(T("Formant"))) );
  formantPreserveButton->assignParameter( 
  pitchShifterModuleToEdit->getParameterByName(T("FormantPreserve")) );
  formantPreserveButton->setDescription(juce::String(T("Preserve formants")));
  formantPreserveButton->setDescriptionField(infoField);
  formantPreserveButton->setClickingTogglesState(true);

  addWidget( monoButton = new RButton(juce::String(T("Mono"))) );
  monoButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Mono")) );
  //monoButton->addRButtonListener(this);
  monoButton->setDescription(juce::String(T("Save CPU for mono signals")));
  monoButton->setDescriptionField(infoField);
  monoButton->setClickingTogglesState(true);
  */

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(400, 180);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void PitchShifterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  updateWidgetVisibility();
}

void PitchShifterModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::updateWidgetsAccordingToState();
  updateWidgetVisibility();
}

void PitchShifterModuleEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  x  = 0;
  w /= 2;
  y = getPresetSectionBottom()+8;

  coarseSlider->setBounds(x+4, y, w-8, 16);

  y = coarseSlider->getBottom();  
  fineSlider->setBounds(x+4, y+4, w-8, 16);

  y = fineSlider->getBottom();  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);

  y = feedbackSlider->getBottom();  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);

  x = w;
  y = coarseSlider->getY();

  grainLengthInBeatsSlider->setBounds(x+4, y, w-64, 16);
  grainLengthInMillisecondsSlider->setBounds(grainLengthInBeatsSlider->getBounds());
  grainLengthInCyclesSlider->setBounds(grainLengthInBeatsSlider->getBounds());
  grainLengthUnitComboBox->setBounds(grainLengthInBeatsSlider->getRight()+4, y, 
    w-grainLengthInBeatsSlider->getWidth()-12, 16);

  y = grainLengthInMillisecondsSlider->getBottom()+8; 

  reverseButton->setBounds(x+4,    y, w/2-8, 16);
  invertButton->setBounds(x+w/2+4, y, w/2-8, 16);
  y += 24;
  //formantPreserveButton->setBounds(x+4, y+4, w/2-8, 16);
  x = reverseButton->getX() + reverseButton->getWidth()/2;
  w = invertButton->getX()  + invertButton->getWidth()/2   - x;
  antiAliasButton->setBounds(x, y+4, w, 16);

  //infoLabel->setBounds(0, getHeight()-20, 40, 20);
  infoField->setBounds(4, getHeight()-20, getWidth()-4,20);
  webLink->setBounds(getWidth()-112, getHeight()-20, 112-4, 20);
}

void PitchShifterModuleEditor::updateWidgetVisibility()
{
  ScopedLock scopedLock(*lock);
  if( pitchShifterModuleToEdit == NULL )
    return;
  if( pitchShifterModuleToEdit->wrappedPitchShifter == NULL )
    return;

  // update the visibility for the 3 grain-length sliders:
  grainLengthInMillisecondsSlider->setVisible(false);  
  grainLengthInCyclesSlider->setVisible(false);
  grainLengthInBeatsSlider->setVisible(false);
  switch( pitchShifterModuleToEdit->wrappedPitchShifter->getGrainLengthUnit() )
  {
  case PitchShifterGrainAdaptive::MILLISECONDS: grainLengthInMillisecondsSlider->setVisible(true);  break;
  case PitchShifterGrainAdaptive::PITCH_CYCLES: grainLengthInCyclesSlider->setVisible(true);        break;
  case PitchShifterGrainAdaptive::BEATS:        grainLengthInBeatsSlider->setVisible(true);         break;
  }
}
