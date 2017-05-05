//-------------------------------------------------------------------------------------------------
// construction/destruction:

AlgoVerbAudioModule::AlgoVerbAudioModule(CriticalSection *newPlugInLock, rosic::AlgoVerb *algoVerbToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(algoVerbToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedAlgoVerb = algoVerbToWrap;
  moduleName      = juce::String(T("AlgoVerb"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/AlgoVerbPresets")) );
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void AlgoVerbAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedAlgoVerb == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedAlgoVerb->setDryWetRatio(           (float) value );   break;
  case   1: wrappedAlgoVerb->setLateReverbLevel(       (float) value );   break;
  case   2: wrappedAlgoVerb->fdn.setReferenceDelayTime(        value );   break;
  case   3: wrappedAlgoVerb->fdn.setInjectionVector(     (int) value );   break;
  case   4: wrappedAlgoVerb->fdn.setFeedbackMatrix(      (int) value );   break;
  case   5: wrappedAlgoVerb->fdn.setOutputVector(        (int) value );   break;
  case   6: wrappedAlgoVerb->fdn.setAllpassMode(         value>=0.5  );   break;


  /*
  case   2: wrappedAlgoVerb->setGrainLengthInMilliseconds(value         );   break;
  case   3: wrappedAlgoVerb->setGrainLengthInPitchCycles( value         );   break;
  case   4: wrappedAlgoVerb->setGrainLengthInBeats(       value         );   break;
  case   5: wrappedAlgoVerb->setGrainLengthUnit(    (int) value         );   break;
  case   6: wrappedAlgoVerb->setFeedback(                 value         );   break;
  case   7: wrappedAlgoVerb->setDryWet(                   value         );   break;
  */
  } // end of switch( parameterIndex )
  markStateAsDirty();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void AlgoVerbAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  juce::Array<double> defaultValues;
  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, "LateLevel", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, "ReferenceDelayTime", 5.0, 100.0, 1.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "InjectionVector", 0.0, 1.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String(T("AllOnes")));
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "FeedbackMatrix", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Identity")));
  p->addStringValue(juce::String(T("MinusIdentity")));
  p->addStringValue(juce::String(T("SeriesConnection")));
  p->addStringValue(juce::String(T("SeriesWithFeedback")));
  p->addStringValue(juce::String(T("Hadamard")));
  p->addStringValue(juce::String(T("MinusHadamard")));
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "OutputVector", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String(T("AllOnes")));
  p->addStringValue(juce::String(T("Out01")));
  p->addStringValue(juce::String(T("Out02")));
  p->addStringValue(juce::String(T("Out03")));
  p->addStringValue(juce::String(T("Out04")));
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "AllpassMode", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);




  /*
  p = new Parameter("Diffusion", 0.0, 100.0, 1.0, 50.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new Parameter("LatePreDelay", 5.0, 100.0, 1.0, 0.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  */


  //p = new Parameter("FeedbackMatrix", 0.0, 5.0, 1.0, 1.0, Parameter::STRING);
  //addObservedParameter(p);



  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}

//=================================================================================================

// construction/destruction:

AlgoVerbModuleEditor::AlgoVerbModuleEditor(CriticalSection *newPlugInLock, AlgoVerbAudioModule* newAlgoVerbAudioModule) 
  : AudioModuleEditor(newPlugInLock, newAlgoVerbAudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String(T("AlgoVerb")) );

  // assign the pointer to the rosic::AlgoVerb object to be used as aduio engine:
  jassert(newAlgoVerbAudioModule != NULL ); // you must pass a valid module here
  algoVerbModuleToEdit = newAlgoVerbAudioModule;


  addWidget( globalLabel = new RTextField( juce::String(T("Global"))) );
  globalLabel->setJustification(Justification::centredLeft);
  globalLabel->setDescription(T("Global Parameters"));
  globalLabel->setDescriptionField(infoField);

  addWidget( dryWetSlider = new RSlider (T("DryWet")) );
  dryWetSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("DryWetRatio")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( pingButton = new RButton(juce::String(T("Ping"))) );
  pingButton->assignParameter( algoVerbModuleToEdit->getParameterByName(T("Ping")) );
  pingButton->setDescription(juce::String(T("Feed an impulse into the reverberator to audition the impulse response")));
  pingButton->setDescriptionField(infoField);
  pingButton->setClickingTogglesState(true);




  addWidget( earlyLabel = new RTextField( juce::String(T("Early Reflections"))) );
  earlyLabel->setJustification(Justification::centredLeft);
  earlyLabel->setDescription(T("Parameters for the early reflection module"));
  earlyLabel->setDescriptionField(infoField);


  addWidget( earlyLabel = new RTextField( juce::String(T("Early Reflections"))) );
  earlyLabel->setJustification(Justification::centredLeft);
  earlyLabel->setDescription(T("Parameters for the early reflection module"));
  earlyLabel->setDescriptionField(infoField);


  //......




  addWidget( lateLabel = new RTextField( juce::String(T("Late Reverb"))) );
  lateLabel->setJustification(Justification::centredLeft);
  lateLabel->setDescription(T("Parameters for the late reverberation module"));
  lateLabel->setDescriptionField(infoField);

  addWidget( lateLevelSlider = new RSlider (T("LateLevelSlider")) );
  lateLevelSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LateLevel")) );
  lateLevelSlider->setSliderName(juce::String(T("Level")));
  lateLevelSlider->setDescription(juce::String(T("Overall level of the late reverberation")));
  lateLevelSlider->setDescriptionField(infoField);
  lateLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( latePingButton = new RButton(juce::String(T("Ping"))) );
  //latePingButton->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LatePing")) );
  latePingButton->setDescription(juce::String(T("Feed an impulse into the late reverb module to audition the impulse response")));
  latePingButton->setDescriptionField(infoField);
  latePingButton->setClickingTogglesState(false);
  latePingButton->addRButtonListener(this);


  addWidget( decayTimeSlider = new RSlider (T("DecayTimeSlider")) );
  decayTimeSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("DecayTime")) );
  decayTimeSlider->setSliderName(juce::String(T("DecayTime")));
  decayTimeSlider->setDescription(juce::String(T("Time for the tail to decay to -60 dB")));
  decayTimeSlider->setDescriptionField(infoField);
  decayTimeSlider->setStringConversionFunction(&rojue::secondsToStringWithUnitTotal4);

  addWidget( lowDecayScaleSlider = new RSlider (T("LowDecayScaleSlider")) );
  lowDecayScaleSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LowDecayScale")) );
  lowDecayScaleSlider->setSliderName(juce::String(T("LowDecayScale")));
  lowDecayScaleSlider->setDescription(juce::String(T("Scale factor for the decay time at low frequencies")));
  lowDecayScaleSlider->setDescriptionField(infoField);
  lowDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( highDecayScaleSlider = new RSlider (T("HighDecayScaleSlider")) );
  highDecayScaleSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("HighDecayScale")) );
  highDecayScaleSlider->setSliderName(juce::String(T("HighDecayScale")));
  highDecayScaleSlider->setDescription(juce::String(T("Scale factor for the decay time at high frequencies")));
  highDecayScaleSlider->setDescriptionField(infoField);
  highDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( lowCrossFreqSlider = new RSlider (T("LowCrossFreqSlider")) );
  lowCrossFreqSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LowCrossFreq")) );
  lowCrossFreqSlider->setSliderName(juce::String(T("LowCrossFreq")));
  lowCrossFreqSlider->setDescription(juce::String(T("Crossover frequency between low and mid frequencies")));
  lowCrossFreqSlider->setDescriptionField(infoField);
  lowCrossFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( highCrossFreqSlider = new RSlider (T("HighCrossFreqSlider")) );
  highCrossFreqSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("HighCrossFreq")) );
  highCrossFreqSlider->setSliderName(juce::String(T("HighCrossFreq")));
  highCrossFreqSlider->setDescription(juce::String(T("Crossover frequency between high and mid frequencies")));
  highCrossFreqSlider->setDescriptionField(infoField);
  highCrossFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( referenceDelayTimeSlider = new RSlider (T("ReferenceDelayTimeSlider")) );
  referenceDelayTimeSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("ReferenceDelayTime")) );
  referenceDelayTimeSlider->setSliderName(juce::String(T("DelayTime")));
  referenceDelayTimeSlider->setDescription(juce::String(T("Arrival time of the first echo/reflection (excluding pre-delay)")));
  referenceDelayTimeSlider->setDescriptionField(infoField);
  referenceDelayTimeSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( latePreDelaySlider = new RSlider (T("LatePreDelaySlider")) );
  latePreDelaySlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LatePreDelay")) );
  latePreDelaySlider->setSliderName(juce::String(T("PreDelay")));
  latePreDelaySlider->setDescription(juce::String(T("Pre-delay for the late reverberation")));
  latePreDelaySlider->setDescriptionField(infoField);
  latePreDelaySlider->setStringConversionFunction(&rojue::secondsToStringWithUnitTotal4);

  /*
  addWidget( densitySlider = new RSlider (T("DensitySlider")) );
  densitySlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("Density")) );
  densitySlider->setSliderName(juce::String(T("Density")));
  densitySlider->setDescription(juce::String(T("Density/packaging of the reflections")));
  densitySlider->setDescriptionField(infoField);
  densitySlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( diffusionSlider = new RSlider (T("DiffusionSlider")) );
  diffusionSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("Diffusion")) );
  diffusionSlider->setSliderName(juce::String(T("Diffusion")));
  diffusionSlider->setDescription(juce::String(T("Diffusion of the reflections")));
  diffusionSlider->setDescriptionField(infoField);
  diffusionSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);
  */




  addWidget( injectionVectorComboBox = new RComboBox(juce::String(T("InjectionVectorComboBox"))) );
  injectionVectorComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName(T("InjectionVector")) );
  injectionVectorComboBox->setDescription(T("Choose the injection vector for the FDN"));
  injectionVectorComboBox->setDescriptionField(infoField);
  injectionVectorComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( feedbackMatrixComboBox = new RComboBox(juce::String(T("FeedbackMatrixComboBox"))) );
  feedbackMatrixComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName(T("FeedbackMatrix")) );
  feedbackMatrixComboBox->setDescription(T("Choose the feedback matrix for the FDN"));
  feedbackMatrixComboBox->setDescriptionField(infoField);
  feedbackMatrixComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( outputVectorComboBox = new RComboBox(juce::String(T("OutputVectorComboBox"))) );
  outputVectorComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName(T("OutputVector")) );
  outputVectorComboBox->setDescription(T("Choose the output vector for the FDN"));
  outputVectorComboBox->setDescriptionField(infoField);
  outputVectorComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( allpassModeButton = new RButton(juce::String(T("AllpassMode"))) );
  allpassModeButton->assignParameter( moduleToEdit->getParameterByName(T("AllpassMode")) );
  allpassModeButton->setDescription(juce::String(T("Switches delaylines into allpass mode")));
  allpassModeButton->setDescriptionField(infoField);
  allpassModeButton->setClickingTogglesState(true);


  // graphical RT60 editor....



  addAndMakeVisible( impulseResponsePlot = new WaveformDisplay() );  // todo: use addPlot
  impulseResponsePlot->setDescription(juce::String(T("Shows the impulse response")));
  impulseResponsePlot->setDescriptionField(infoField);

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void AlgoVerbModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( algoVerbModuleToEdit == NULL )
    return;
  if( algoVerbModuleToEdit->wrappedAlgoVerb == NULL )
    return;



  if( buttonThatWasClicked == latePingButton )
    algoVerbModuleToEdit->wrappedAlgoVerb->fdn.feedInImpulse();


}

void AlgoVerbModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{

}


void AlgoVerbModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();

  // updatePlots();
}

void AlgoVerbModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  globalRect.setBounds(x, y+4, w, 24);

  y += 24;
  h  = (getHeight()-y)/2;
  w /= 2;
  earlyRect.setBounds(x,   y, w, h);
  lateRect.setBounds(x+w, y, w, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(globalRect);
  guiLayoutRectangles.add(earlyRect); 
  guiLayoutRectangles.add(lateRect);




  x = globalRect.getX();
  y = globalRect.getY();
  globalLabel->setBounds(x+4, y+2, 80, 16);
  dryWetSlider->setBounds(globalLabel->getRight()+4, y+4, 120, 16);


  x = earlyRect.getX();
  y = earlyRect.getY();
  earlyLabel->setBounds(x+4, y+4, 120, 16);



  x = lateRect.getX();
  y = lateRect.getY();
  lateLabel->setBounds(x+4, y+4, 80, 16);

  latePingButton->setBounds(x+w-44, y+4, 40, 16);
  lateLevelSlider->setBounds(lateLabel->getRight()+4, y+4, latePingButton->getX()-lateLabel->getRight()-8, 16);
  w  = lateRect.getWidth()/2;
  y += 20;
  referenceDelayTimeSlider->setBounds(x+4, y+4, w-8, 16);
  latePreDelaySlider->setBounds(x+w+4, y+4, w-8, 16);
  y += 20;
  injectionVectorComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  feedbackMatrixComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  outputVectorComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  allpassModeButton->setBounds(x+4, y+4, w-8, 16);



  y = lateRect.getBottom();
  h = getHeight()-y-20;
  impulseResponsePlot->setBounds(0, y, getWidth(), h);


  /*
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
  */
}

void AlgoVerbModuleEditor::updatePlots()
{
  //...
}