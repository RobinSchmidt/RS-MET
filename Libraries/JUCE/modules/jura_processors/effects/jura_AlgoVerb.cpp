//-------------------------------------------------------------------------------------------------
// construction/destruction:

AlgoVerbAudioModule::AlgoVerbAudioModule(CriticalSection *newPlugInLock, rosic::AlgoVerb *algoVerbToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(algoVerbToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedAlgoVerb = algoVerbToWrap;
  setModuleTypeName("AlgoVerb");
  initializeAutomatableParameters();
}

AlgoVerbAudioModule::AlgoVerbAudioModule(CriticalSection *newPlugInLock) : AudioModule(newPlugInLock)
{
  wrappedAlgoVerb = new rosic::AlgoVerb;
  wrappedAlgoVerbIsOwned = true;

  // get rid of duplication - make init function
  setModuleTypeName("AlgoVerb");
  initializeAutomatableParameters();
}

AlgoVerbAudioModule::~AlgoVerbAudioModule()
{
  if(wrappedAlgoVerbIsOwned)
    delete wrappedAlgoVerb;
}

AudioModuleEditor* AlgoVerbAudioModule::createEditor(int type)
{
  return new jura::AlgoVerbModuleEditor(lock, this); // get rid of passing the lock
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

  p = new AutomatableParameter(lock, "DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "LateLevel", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "ReferenceDelayTime", 5.0, 100.0, 1.0, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "InjectionVector", 0.0, 1.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String("AllOnes"));
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "FeedbackMatrix", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String("Identity"));
  p->addStringValue(juce::String("MinusIdentity"));
  p->addStringValue(juce::String("SeriesConnection"));
  p->addStringValue(juce::String("SeriesWithFeedback"));
  p->addStringValue(juce::String("Hadamard"));
  p->addStringValue(juce::String("MinusHadamard"));
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "OutputVector", 0.0, 5.0, 1.0, 4.0, Parameter::STRING);
  p->addStringValue(juce::String("AllOnes"));
  p->addStringValue(juce::String("Out01"));
  p->addStringValue(juce::String("Out02"));
  p->addStringValue(juce::String("Out03"));
  p->addStringValue(juce::String("Out04"));
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "AllpassMode", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
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
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

AlgoVerbModuleEditor::AlgoVerbModuleEditor(CriticalSection *newPlugInLock, 
  AlgoVerbAudioModule* newAlgoVerbAudioModule) 
  : AudioModuleEditor(newAlgoVerbAudioModule)
{
  // assign the pointer to the rosic::AlgoVerb object to be used as aduio engine:
  jassert(newAlgoVerbAudioModule != NULL ); // you must pass a valid module here
  algoVerbModuleToEdit = newAlgoVerbAudioModule;


  addWidget( globalLabel = new RTextField( juce::String("Global")) );
  globalLabel->setJustification(Justification::centredLeft);
  globalLabel->setDescription("Global Parameters");
  globalLabel->setDescriptionField(infoField);

  addWidget( dryWetSlider = new RSlider("DryWet") );
  dryWetSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String("Dry/Wet"));
  dryWetSlider->setDescription(juce::String("Ratio between dry and wet signal"));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&ratioToString0);

  addWidget( pingButton = new RButton(juce::String("Ping")) );
  //pingButton->assignParameter( algoVerbModuleToEdit->getParameterByName("Ping") ); // there's no "Ping" Parameter
  pingButton->setDescription(juce::String("Feed an impulse into the reverberator to audition the impulse response"));
  pingButton->setDescriptionField(infoField);
  pingButton->setClickingTogglesState(true);

  addWidget( earlyLabel = new RTextField( juce::String("Early Reflections")) );
  earlyLabel->setJustification(Justification::centredLeft);
  earlyLabel->setDescription("Parameters for the early reflection module");
  earlyLabel->setDescriptionField(infoField);

  addWidget( earlyLabel = new RTextField( juce::String("Early Reflections")) );
  earlyLabel->setJustification(Justification::centredLeft);
  earlyLabel->setDescription("Parameters for the early reflection module");
  earlyLabel->setDescriptionField(infoField);

  //......

  addWidget( lateLabel = new RTextField( juce::String("Late Reverb")) );
  lateLabel->setJustification(Justification::centredLeft);
  lateLabel->setDescription("Parameters for the late reverberation module");
  lateLabel->setDescriptionField(infoField);

  addWidget( lateLevelSlider = new RSlider("LateLevelSlider") );
  lateLevelSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("LateLevel") );
  lateLevelSlider->setSliderName(juce::String("Level"));
  lateLevelSlider->setDescription(juce::String("Overall level of the late reverberation"));
  lateLevelSlider->setDescriptionField(infoField);
  lateLevelSlider->setStringConversionFunction(decibelsToStringWithUnit1);

  addWidget( latePingButton = new RButton(juce::String("Ping")) );
  //latePingButton->assignParameter( algoVerbModuleToEdit->getParameterByName(T("LatePing")) );
  latePingButton->setDescription(juce::String("Feed an impulse into the late reverb module to audition the impulse response"));
  latePingButton->setDescriptionField(infoField);
  latePingButton->setClickingTogglesState(false);
  latePingButton->addRButtonListener(this);


  addWidget( decayTimeSlider = new RSlider("DecayTimeSlider") );
  decayTimeSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("DecayTime") );
  decayTimeSlider->setSliderName(juce::String("DecayTime"));
  decayTimeSlider->setDescription(juce::String("Time for the tail to decay to -60 dB"));
  decayTimeSlider->setDescriptionField(infoField);
  decayTimeSlider->setStringConversionFunction(secondsToStringWithUnitTotal4);

  addWidget( lowDecayScaleSlider = new RSlider("LowDecayScaleSlider") );
  lowDecayScaleSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("LowDecayScale") );
  lowDecayScaleSlider->setSliderName(juce::String("LowDecayScale"));
  lowDecayScaleSlider->setDescription(juce::String("Scale factor for the decay time at low frequencies"));
  lowDecayScaleSlider->setDescriptionField(infoField);
  lowDecayScaleSlider->setStringConversionFunction(valueToString2);

  addWidget( highDecayScaleSlider = new RSlider("HighDecayScaleSlider") );
  highDecayScaleSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("HighDecayScale") );
  highDecayScaleSlider->setSliderName(juce::String("HighDecayScale"));
  highDecayScaleSlider->setDescription(juce::String("Scale factor for the decay time at high frequencies"));
  highDecayScaleSlider->setDescriptionField(infoField);
  highDecayScaleSlider->setStringConversionFunction(valueToString2);

  addWidget( lowCrossFreqSlider = new RSlider("LowCrossFreqSlider") );
  lowCrossFreqSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("LowCrossFreq") );
  lowCrossFreqSlider->setSliderName(juce::String("LowCrossFreq"));
  lowCrossFreqSlider->setDescription(juce::String("Crossover frequency between low and mid frequencies"));
  lowCrossFreqSlider->setDescriptionField(infoField);
  lowCrossFreqSlider->setStringConversionFunction(hertzToStringWithUnitTotal5);

  addWidget( highCrossFreqSlider = new RSlider("HighCrossFreqSlider") );
  highCrossFreqSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("HighCrossFreq") );
  highCrossFreqSlider->setSliderName(juce::String("HighCrossFreq"));
  highCrossFreqSlider->setDescription(juce::String("Crossover frequency between high and mid frequencies"));
  highCrossFreqSlider->setDescriptionField(infoField);
  highCrossFreqSlider->setStringConversionFunction(hertzToStringWithUnitTotal5);

  addWidget( referenceDelayTimeSlider = new RSlider("ReferenceDelayTimeSlider") );
  referenceDelayTimeSlider->assignParameter( algoVerbModuleToEdit->getParameterByName("ReferenceDelayTime") );
  referenceDelayTimeSlider->setSliderName(juce::String("DelayTime"));
  referenceDelayTimeSlider->setDescription(juce::String("Arrival time of the first echo/reflection (excluding pre-delay)"));
  referenceDelayTimeSlider->setDescriptionField(infoField);
  referenceDelayTimeSlider->setStringConversionFunction(millisecondsToStringWithUnit2);

  addWidget( latePreDelaySlider = new RSlider("LatePreDelaySlider") );
  latePreDelaySlider->assignParameter( algoVerbModuleToEdit->getParameterByName("LatePreDelay") );
  latePreDelaySlider->setSliderName(juce::String("PreDelay"));
  latePreDelaySlider->setDescription(juce::String("Pre-delay for the late reverberation"));
  latePreDelaySlider->setDescriptionField(infoField);
  latePreDelaySlider->setStringConversionFunction(secondsToStringWithUnitTotal4);


  /*
  addWidget( densitySlider = new RSlider (T("DensitySlider")) );
  densitySlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("Density")) );
  densitySlider->setSliderName(juce::String(T("Density")));
  densitySlider->setDescription(juce::String(T("Density/packaging of the reflections")));
  densitySlider->setDescriptionField(infoField);
  densitySlider->setStringConversionFunction(percentToStringWithUnit0);

  addWidget( diffusionSlider = new RSlider (T("DiffusionSlider")) );
  diffusionSlider->assignParameter( algoVerbModuleToEdit->getParameterByName(T("Diffusion")) );
  diffusionSlider->setSliderName(juce::String(T("Diffusion")));
  diffusionSlider->setDescription(juce::String(T("Diffusion of the reflections")));
  diffusionSlider->setDescriptionField(infoField);
  diffusionSlider->setStringConversionFunction(percentToStringWithUnit0);
  */

  addWidget( injectionVectorComboBox = new RComboBox(juce::String("InjectionVectorComboBox")) );
  injectionVectorComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName("InjectionVector") );
  injectionVectorComboBox->setDescription("Choose the injection vector for the FDN");
  injectionVectorComboBox->setDescriptionField(infoField);
  injectionVectorComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( feedbackMatrixComboBox = new RComboBox(juce::String("FeedbackMatrixComboBox")) );
  feedbackMatrixComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName("FeedbackMatrix") );
  feedbackMatrixComboBox->setDescription("Choose the feedback matrix for the FDN");
  feedbackMatrixComboBox->setDescriptionField(infoField);
  feedbackMatrixComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( outputVectorComboBox = new RComboBox(juce::String("OutputVectorComboBox")) );
  outputVectorComboBox->assignParameter( algoVerbModuleToEdit->getParameterByName("OutputVector") );
  outputVectorComboBox->setDescription("Choose the output vector for the FDN");
  outputVectorComboBox->setDescriptionField(infoField);
  outputVectorComboBox->registerComboBoxObserver(this); // to update the plot

  addWidget( allpassModeButton = new RButton(juce::String("AllpassMode")) );
  allpassModeButton->assignParameter( moduleToEdit->getParameterByName("AllpassMode") );
  allpassModeButton->setDescription(juce::String("Switches delaylines into allpass mode"));
  allpassModeButton->setDescriptionField(infoField);
  allpassModeButton->setClickingTogglesState(true);

  // graphical RT60 editor....

  addAndMakeVisible( impulseResponsePlot = new WaveformDisplay() );  // todo: use addPlot
  impulseResponsePlot->setDescription(juce::String("Shows the impulse response"));
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