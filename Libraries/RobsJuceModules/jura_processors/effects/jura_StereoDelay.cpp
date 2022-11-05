//-------------------------------------------------------------------------------------------------
// construction/destruction:

StereoDelayAudioModule::StereoDelayAudioModule(CriticalSection *newPlugInLock,
  rosic::StereoDelay *stereoDelayToWrap)
: AudioModule(newPlugInLock)
{
  //jassert(stereoDelayToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if( stereoDelayToWrap != nullptr)
    wrappedStereoDelay = stereoDelayToWrap;
  else
  {
    wrappedStereoDelay = new rosic::StereoDelay;
    wrappedStereoDelayIsOwned = true;
  }
  setModuleTypeName("StereoDelay");
  initializeAutomatableParameters();
}

StereoDelayAudioModule::~StereoDelayAudioModule()
{
  if(wrappedStereoDelayIsOwned)
    delete wrappedStereoDelay;
}

AudioModuleEditor* StereoDelayAudioModule::createEditor(int type)
{

  return new StereoDelayModuleEditor(lock, this);
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* StereoDelayAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedStereoDelay == NULL )
    return xmlState;

  // add attributes for the non-automatable parameters (the automatable ones are already taken care
  // of by AudioModule::getStateAsXml()):
  xmlState->setAttribute(("DelayInBeatsL"), wrappedStereoDelay->getDelayInBeats(StereoDelay::LEFT)  );
  xmlState->setAttribute(("DelayInBeatsR"), wrappedStereoDelay->getDelayInBeats(StereoDelay::RIGHT) );

  return xmlState;
}

void StereoDelayAudioModule::setStateFromXml(const XmlElement& xmlState,
                                             const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedStereoDelay == NULL )
    return;

  // restore the values of the non-automatable parameters (the automatable ones are already taken
  // care of by automatableModuleStateFromXml():
  wrappedStereoDelay->setDelayInBeats(
    xmlState.getDoubleAttribute(("DelayInBeatsL"),   1.0), StereoDelay::LEFT  );
  wrappedStereoDelay->setDelayInBeats(
    xmlState.getDoubleAttribute(("DelayInBeatsR"),   1.0), StereoDelay::RIGHT );
}

//-------------------------------------------------------------------------------------------------
// automation:

void StereoDelayAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedStereoDelay == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedStereoDelay->setDelayScale(                value, StereoDelay::LEFT);                       break;
  case   1: wrappedStereoDelay->setDelayScale(                value, StereoDelay::RIGHT);                      break;
  case   2: wrappedStereoDelay->setInjection(                 value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case   3: wrappedStereoDelay->setInjection(                 value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case   4: wrappedStereoDelay->setInjection(                 value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case   5: wrappedStereoDelay->setInjection(                 value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case   6: wrappedStereoDelay->setDiffusorTimeInMilliseconds(value, StereoDelay::LEFT);                       break;
  case   7: wrappedStereoDelay->setDiffusorTimeInMilliseconds(value, StereoDelay::RIGHT);                      break;
  case   8: wrappedStereoDelay->setDiffusorAmount(            value, StereoDelay::LEFT);                       break;
  case   9: wrappedStereoDelay->setDiffusorAmount(            value, StereoDelay::RIGHT);                      break;
  case  10: wrappedStereoDelay->setLowpassCutoff(             value, StereoDelay::LEFT);                       break;
  case  11: wrappedStereoDelay->setLowpassCutoff(             value, StereoDelay::RIGHT);                      break;
  case  12: wrappedStereoDelay->setHighpassCutoff(            value, StereoDelay::LEFT);                       break;
  case  13: wrappedStereoDelay->setHighpassCutoff(            value, StereoDelay::RIGHT);                      break;
  case  14: wrappedStereoDelay->setCutoffScale(               value);                                          break;
  case  15: wrappedStereoDelay->setFeedback(                  value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case  16: wrappedStereoDelay->setFeedback(                  value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case  17: wrappedStereoDelay->setFeedback(                  value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case  18: wrappedStereoDelay->setFeedback(                  value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case  19: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::LEFT,  StereoDelay::LEFT);   break;
  case  20: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::LEFT,  StereoDelay::RIGHT);  break;
  case  21: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::RIGHT, StereoDelay::LEFT);   break;
  case  22: wrappedStereoDelay->setOutputMix(                 value, StereoDelay::RIGHT, StereoDelay::RIGHT);  break;
  case  23: wrappedStereoDelay->setWetDelayInBeats(           value, StereoDelay::LEFT);                       break;
  case  24: wrappedStereoDelay->setWetDelayInBeats(           value, StereoDelay::RIGHT);                      break;
  case  25: wrappedStereoDelay->setDryWet(                    value);                                          break;
  } // end of switch( parameterIndex )
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void StereoDelayAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created
  // Parameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "DelayScaleL", 0.5, 2.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "DelayScaleR", 0.5, 2.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #002:
  p = new AutomatableParameter(lock, "LeftInToLeftDelay", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #003:
  p = new AutomatableParameter(lock, "RightInToLeftDelay", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #004:
  p = new AutomatableParameter(lock, "LeftInToRightDelay", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #005:
  p = new AutomatableParameter(lock, "RightInToRightDelay", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #006:
  p = new AutomatableParameter(lock, "DiffusorTimeL", 1.0, 50.0, 0.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #007:
  p = new AutomatableParameter(lock, "DiffusorTimeR", 1.0, 50.0, 0.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #008:
  p = new AutomatableParameter(lock, "DiffusorAmountL", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #009:
  p = new AutomatableParameter(lock, "DiffusorAmountR", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #010:
  p = new AutomatableParameter(lock, "LowpassCutoffL", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #011:
  p = new AutomatableParameter(lock, "LowpassCutoffR", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #012:
  p = new AutomatableParameter(lock, "HighpassCutoffL", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #013:
  p = new AutomatableParameter(lock, "HighpassCutoffR", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #014:
  p = new AutomatableParameter(lock, "CutoffScale", 0.0625, 16.0, 0.0, 1.0, Parameter::EXPONENTIAL, 74);
  addObservedParameter(p);

  // #015:
  p = new AutomatableParameter(lock, "FeedbackLeftToLeft", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #016:
  p = new AutomatableParameter(lock, "FeedbackRightToRight", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #017:
  p = new AutomatableParameter(lock, "FeedbackLeftToRight", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #018:
  p = new AutomatableParameter(lock, "FeedbackRightToLeft", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #019:
  p = new AutomatableParameter(lock, "LeftDelayToLeftOut", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #020:
  p = new AutomatableParameter(lock, "LeftDelayToRightOut", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #021:
  p = new AutomatableParameter(lock, "RightDelayToLeftOut", -100.0, 100.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #022:
  p = new AutomatableParameter(lock, "RightDelayToRightOut", -100.0, 100.0, 0.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #023:
  p = new AutomatableParameter(lock, "WetDelayInBeatsL", 0.0, 4.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #025:
  p = new AutomatableParameter(lock, "WetDelayInBeatsR", 0.0, 4.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #025:
  p = new AutomatableParameter(lock, "DryWet", 0.0, 100.0, 0.0, 50.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

// construction/destruction:

StereoDelayModuleEditor::StereoDelayModuleEditor(CriticalSection *newPlugInLock, StereoDelayAudioModule* newStereoDelayAudioModule)
  : AudioModuleEditor(newStereoDelayAudioModule)
{
  // assign the pointer to the rosic::StereoDelay object to be used as aduio engine:
  jassert(newStereoDelayAudioModule != NULL ); // you must pass a valid module here
  stereoDelayAudioModule = newStereoDelayAudioModule;

  // global widgets:

  addWidget( dryWetSlider = new RSlider(("DryWetSlider")) );
  dryWetSlider->assignParameter( stereoDelayAudioModule->getParameterByName("DryWet") );
  dryWetSlider->setSliderName(juce::String(("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(("Ratio between dry and wet signal (in % wet)")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(dryWetSlider);

  addWidget( cutoffScaleSlider = new RSlider (("CutoffScaleSlider")) );
  cutoffScaleSlider->assignParameter( stereoDelayAudioModule->getParameterByName("CutoffScale") );
  cutoffScaleSlider->setSliderName(juce::String(("Cutoff Scale")));
  cutoffScaleSlider->setDescription(juce::String(("Scales the cutoff frequencies of the filters by a factor")));
  cutoffScaleSlider->setDescriptionField(infoField);
  cutoffScaleSlider->setStringConversionFunction(&valueToString4);
  //automatableSliders.addIfNotAlreadyThere(cutoffScaleSlider);

  // left delayline widgets:

  addWidget( delayLineLabelL = new RTextField( juce::String(("Left Delay"))) );
  delayLineLabelL->setDescription(("These are the parameters for the left delayline."));
  delayLineLabelL->setDescriptionField(infoField);

  addWidget( delayComboBoxL = new RSyncIntervalComboBox(juce::String(("DelayComboBoxL"))) );
  delayComboBoxL->registerComboBoxObserver(this);
  delayComboBoxL->setDescription(("Left delay in beats"));
  delayComboBoxL->setDescriptionField(infoField);

  addWidget( delayScaleSliderL = new RSlider (("DelayScaleSliderL")) );
  //delayScaleSliderL->addListener(this);
  //delayScaleSliderL->setRange(0.5, 2.0, 0.0001, 1.0);
  //delayScaleSliderL->setScaling(Parameter::EXPONENTIAL);
  delayScaleSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("DelayScaleL") );
  delayScaleSliderL->setSliderName(juce::String(("Delay Scale")));
  delayScaleSliderL->setDescription(juce::String(("Scales the left delay time by a factor")));
  delayScaleSliderL->setDescriptionField(infoField);
  delayScaleSliderL->setStringConversionFunction(&valueToString4);
  //automatableSliders.addIfNotAlreadyThere(delayScaleSliderL);

  addWidget( inputLabelL = new RTextField( juce::String(("Input:"))) );
  inputLabelL->setDescription(("Injection parameters for the left delayline."));
  inputLabelL->setDescriptionField(infoField);

  addWidget( inputSliderL2L = new RSlider (("InputSliderL2L")) );
  //inputSliderL2L->addListener(this);
  //inputSliderL2L->setRange(-100.0, 100.0, 0.1, 100.0);
  inputSliderL2L->assignParameter( stereoDelayAudioModule->getParameterByName("LeftInToLeftDelay") );
  inputSliderL2L->setSliderName(juce::String(("from left")));
  inputSliderL2L->setDescription(juce::String(("Amount by which the left input goes to the left delayline")));
  inputSliderL2L->setDescriptionField(infoField);
  inputSliderL2L->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(inputSliderL2L);

  addWidget( inputSliderR2L = new RSlider (("InputSliderR2L")) );
  //inputSliderR2L->addListener(this);
  //inputSliderR2L->setRange(-100.0, 100.0, 0.1, 0.0);
  inputSliderR2L->assignParameter( stereoDelayAudioModule->getParameterByName("RightInToLeftDelay") );
  inputSliderR2L->setSliderName(juce::String(("from right")));
  inputSliderR2L->setDescription(juce::String(("Amount by which the right input goes to the left delayline")));
  inputSliderR2L->setDescriptionField(infoField);
  inputSliderR2L->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(inputSliderR2L);

  addWidget( diffusorLabelL = new RTextField( juce::String(("Diffusor:"))) );
  diffusorLabelL->setDescription(("Diffusion parameters for the left delayline."));
  diffusorLabelL->setDescriptionField(infoField);

  addWidget( diffusorTimeSliderL = new RSlider (("DiffusorTimeSliderL")) );
  //diffusorTimeSliderL->addListener(this);
  //diffusorTimeSliderL->setRange(1.0, 50.0, 0.01, 10.0);
  //diffusorTimeSliderL->setScaling(Parameter::EXPONENTIAL);
  diffusorTimeSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("DiffusorTimeL") );
  diffusorTimeSliderL->setSliderName(juce::String(("Time")));
  diffusorTimeSliderL->setDescription(juce::String(("Delay time for the left allpass diffusor")));
  diffusorTimeSliderL->setDescriptionField(infoField);
  diffusorTimeSliderL->setStringConversionFunction(&millisecondsToStringWithUnit2);
  //automatableSliders.addIfNotAlreadyThere(diffusorTimeSliderL);

  addWidget( diffusorAmountSliderL = new RSlider (("DiffusorAmountSliderL")) );
  //diffusorAmountSliderL->addListener(this);
  //diffusorAmountSliderL->setRange(-100.0, 100.0, 0.1, 0.0);
  diffusorAmountSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("DiffusorAmountL") );
  diffusorAmountSliderL->setSliderName(juce::String(("Amount")));
  diffusorAmountSliderL->setDescription(juce::String(("Amount for the left allpass diffusor")));
  diffusorAmountSliderL->setDescriptionField(infoField);
  diffusorAmountSliderL->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(diffusorAmountSliderL);

  addWidget( filterLabelL = new RTextField( juce::String(("Filter:"))) );
  filterLabelL->setDescription(("Filter parameters for the left delayline."));
  filterLabelL->setDescriptionField(infoField);

  addWidget( lowpassSliderL = new RSlider (("LowpassSliderL")) );
  //lowpassSliderL->addListener(this);
  //lowpassSliderL->setRange(20.0, 20000.0, 0.001, 20000.0);
  //lowpassSliderL->setScaling(Parameter::EXPONENTIAL);
  lowpassSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("LowpassCutoffL") );
  lowpassSliderL->setSliderName(juce::String(("Lowpass")));
  lowpassSliderL->setDescription(juce::String(("Lowpass cutoff for the left delayline")));
  lowpassSliderL->setDescriptionField(infoField);
  lowpassSliderL->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(lowpassSliderL);

  addWidget( highpassSliderL = new RSlider (("HighpassSliderL")) );
  //highpassSliderL->addListener(this);
  //highpassSliderL->setRange(20.0, 20000.0, 0.001, 20.0);
  //highpassSliderL->setScaling(Parameter::EXPONENTIAL);
  highpassSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("HighpassCutoffL") );
  highpassSliderL->setSliderName(juce::String(("Highpass")));
  highpassSliderL->setDescription(juce::String(("Highpass cutoff for the left delayline")));
  highpassSliderL->setDescriptionField(infoField);
  highpassSliderL->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(highpassSliderL);

  addWidget( feedbackLabelL = new RTextField( juce::String(("Feedback:"))) );
  feedbackLabelL->setDescription(("Feedback parameters for the left delayline."));
  feedbackLabelL->setDescriptionField(infoField);

  addWidget( feedbackSliderL = new RSlider (("FeedbackSliderL")) );
  //feedbackSliderL->addListener(this);
  //feedbackSliderL->setRange(-100.0, 100.0, 0.1, 0.0);
  feedbackSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("FeedbackLeftToLeft") );
  feedbackSliderL->setSliderName(juce::String(("from left")));
  feedbackSliderL->setDescription(juce::String(("Amount of self-feedback for the left delayline")));
  feedbackSliderL->setDescriptionField(infoField);
  feedbackSliderL->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(feedbackSliderL);

  addWidget( crossFeedbackSliderL = new RSlider (("FeedbackSliderL")) );
  //crossFeedbackSliderL->addListener(this);
  //crossFeedbackSliderL->setRange(-100.0, 100.0, 0.1, 0.0);
  crossFeedbackSliderL->assignParameter( stereoDelayAudioModule->getParameterByName("FeedbackRightToLeft") );
  crossFeedbackSliderL->setSliderName(juce::String(("from right")));
  crossFeedbackSliderL->setDescription(juce::String(("Amount of cross-feedback from right to left")));
  crossFeedbackSliderL->setDescriptionField(infoField);
  crossFeedbackSliderL->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(crossFeedbackSliderL);

  addWidget( outputLabelL = new RTextField( juce::String(("Output:"))) );
  outputLabelL->setDescription(("Output parameters for the left delayline."));
  outputLabelL->setDescriptionField(infoField);

  addWidget( outputSliderL2L = new RSlider (("OutputSliderL2L")) );
  //outputSliderL2L->addListener(this);
  //outputSliderL2L->setRange(-100.0, 100.0, 0.1, 100.0);
  outputSliderL2L->assignParameter( stereoDelayAudioModule->getParameterByName("LeftDelayToLeftOut") );
  outputSliderL2L->setSliderName(juce::String(("to left")));
  outputSliderL2L->setDescription(juce::String(("Amount by which the left wet signal goes to the left output")));
  outputSliderL2L->setDescriptionField(infoField);
  outputSliderL2L->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(outputSliderL2L);

  addWidget( outputSliderL2R = new RSlider (("OutputSliderL2R")) );
  //outputSliderL2R->addListener(this);
  //outputSliderL2R->setRange(-100.0, 100.0, 0.1, 0.0);
  outputSliderL2R->assignParameter( stereoDelayAudioModule->getParameterByName("LeftDelayToRightOut") );
  outputSliderL2R->setSliderName(juce::String(("to right")));
  outputSliderL2R->setDescription(juce::String(("Amount by which the left wet signal goes to the right output")));
  outputSliderL2R->setDescriptionField(infoField);
  outputSliderL2R->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(outputSliderL2R);

  addWidget( outDelaySliderL = new RSlider (("OutputDelaySliderL")) );
  //outDelaySliderL->addListener(this);
  //outDelaySliderL->setRange(0.0, 4.0, 0.0, 0.0);
  outDelaySliderL->assignParameter( stereoDelayAudioModule->getParameterByName("WetDelayInBeatsL") );
  outDelaySliderL->setSliderName(juce::String(("Delay")));
  outDelaySliderL->setDescription(juce::String(("Delay for the left channels wet signal (aka 'Pre-Delay')")));
  outDelaySliderL->setDescriptionField(infoField);
  outDelaySliderL->setStringConversionFunction(&beatsToStringWithUnit4);
  //automatableSliders.addIfNotAlreadyThere(outDelaySliderL);

  // right delayline widgets:

  addWidget( delayLineLabelR = new RTextField( juce::String(("Right Delay"))) );
  delayLineLabelR->setDescription(("These are the parameters for the right delayline"));
  delayLineLabelR->setDescriptionField(infoField);

  addWidget( delayComboBoxR = new RSyncIntervalComboBox(juce::String(("DelayComboBoxR"))) );
  delayComboBoxR->registerComboBoxObserver(this);
  delayComboBoxR->setDescription(("Right delay in beats"));
  delayComboBoxR->setDescriptionField(infoField);

  addWidget( delayScaleSliderR = new RSlider (("DelayScaleSliderR")) );
  //delayScaleSliderR->addListener(this);
  //delayScaleSliderR->setRange(0.5, 2.0, 0.0001, 1.0);
  //delayScaleSliderR->setScaling(Parameter::EXPONENTIAL);
  delayScaleSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("DelayScaleR") );
  delayScaleSliderR->setSliderName(juce::String(("Delay Scale")));
  delayScaleSliderR->setDescription(juce::String(("Scales the right delay time by a factor")));
  delayScaleSliderR->setDescriptionField(infoField);
  delayScaleSliderR->setStringConversionFunction(&valueToString4);
  //automatableSliders.addIfNotAlreadyThere(delayScaleSliderR);

  addWidget( inputLabelR = new RTextField( juce::String(("Input:"))) );
  inputLabelR->setDescription(("Injection parameters for the right delayline."));
  inputLabelR->setDescriptionField(infoField);

  addWidget( inputSliderL2R = new RSlider (("InputSliderL2R")) );
  //inputSliderL2R->addListener(this);
  //inputSliderL2R->setRange(-100.0, 100.0, 0.1, 0.0);
  inputSliderL2R->assignParameter( stereoDelayAudioModule->getParameterByName("LeftInToRightDelay") );
  inputSliderL2R->setSliderName(juce::String(("from left")));
  inputSliderL2R->setDescription(juce::String(("Amount by which the left input goes to the right delayline")));
  inputSliderL2R->setDescriptionField(infoField);
  inputSliderL2R->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(inputSliderL2R);

  addWidget( inputSliderR2R = new RSlider (("InputSliderR2R")) );
  //inputSliderR2R->addListener(this);
  //inputSliderR2R->setRange(-100.0, 100.0, 0.1, 100.0);
  inputSliderR2R->assignParameter( stereoDelayAudioModule->getParameterByName("RightInToRightDelay") );
  inputSliderR2R->setSliderName(juce::String(("from right")));
  inputSliderR2R->setDescription(juce::String(("Amount by which the right input goes to the right delayline")));
  inputSliderR2R->setDescriptionField(infoField);
  inputSliderR2R->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(inputSliderR2R);

  addWidget( diffusorLabelR = new RTextField( juce::String(("Diffusor:"))) );
  diffusorLabelR->setDescription(("Diffusion parameters for the right delayline."));
  diffusorLabelR->setDescriptionField(infoField);

  addWidget( diffusorTimeSliderR = new RSlider (("DiffusorTimeSliderR")) );
  //diffusorTimeSliderR->addListener(this);
  //diffusorTimeSliderR->setRange(1.0, 50.0, 0.01, 10.0);
  //diffusorTimeSliderR->setScaling(Parameter::EXPONENTIAL);
  diffusorTimeSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("DiffusorTimeR") );
  diffusorTimeSliderR->setSliderName(juce::String(("Diffusor Time")));
  diffusorTimeSliderR->setDescription(juce::String(("Delay time for the right allpass diffusor")));
  diffusorTimeSliderR->setDescriptionField(infoField);
  diffusorTimeSliderR->setStringConversionFunction(&millisecondsToStringWithUnit2);
  //automatableSliders.addIfNotAlreadyThere(diffusorTimeSliderR);

  addWidget( diffusorAmountSliderR = new RSlider (("DiffusorAmountSliderR")) );
  //diffusorAmountSliderR->addListener(this);
  //diffusorAmountSliderR->setRange(-100.0, 100.0, 0.1, 0.0);
  diffusorAmountSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("DiffusorAmountR") );
  diffusorAmountSliderR->setSliderName(juce::String(("Diffusion")));
  diffusorAmountSliderR->setDescription(juce::String(("Amount for the right allpass diffusor")));
  diffusorAmountSliderR->setDescriptionField(infoField);
  diffusorAmountSliderR->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(diffusorAmountSliderR);

  addWidget( filterLabelR = new RTextField( juce::String(("Filter:"))) );
  filterLabelR->setDescription(("Filter parameters for the right delayline."));
  filterLabelR->setDescriptionField(infoField);

  addWidget( lowpassSliderR = new RSlider (("LowpassSliderR")) );
  //lowpassSliderR->addListener(this);
  //lowpassSliderR->setRange(20.0, 20000.0, 0.001, 20000.0);
  //lowpassSliderR->setScaling(Parameter::EXPONENTIAL);
  lowpassSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("LowpassCutoffR") );
  lowpassSliderR->setSliderName(juce::String(("Lowpass")));
  lowpassSliderR->setDescription(juce::String(("Lowpass cutoff for the right delayline")));
  lowpassSliderR->setDescriptionField(infoField);
  lowpassSliderR->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(lowpassSliderR);

  addWidget( highpassSliderR = new RSlider (("HighpassSliderR")) );
  //highpassSliderR->addListener(this);
  //highpassSliderR->setRange(20.0, 20000.0, 0.001, 20.0);
  //highpassSliderR->setScaling(Parameter::EXPONENTIAL);
  highpassSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("HighpassCutoffR") );
  highpassSliderR->setSliderName(juce::String(("Highpass")));
  highpassSliderR->setDescription(juce::String(("Highpass cutoff for the right delayline")));
  highpassSliderR->setDescriptionField(infoField);
  highpassSliderR->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  //automatableSliders.addIfNotAlreadyThere(highpassSliderR);

  addWidget( feedbackLabelR = new RTextField( juce::String(("Feedback:"))) );
  feedbackLabelR->setDescription(("Feedback parameters for the right delayline."));
  feedbackLabelR->setDescriptionField(infoField);

  addWidget( feedbackSliderR = new RSlider (("FeedbackSliderR")) );
  //feedbackSliderR->addListener(this);
  //feedbackSliderR->setRange(-100.0, 100.0, 0.1, 0.0);
  feedbackSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("FeedbackRightToRight") );
  feedbackSliderR->setSliderName(juce::String(("from right")));
  feedbackSliderR->setDescription(juce::String(("Amount of self-feedback for the right delayline")));
  feedbackSliderR->setDescriptionField(infoField);
  feedbackSliderR->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(feedbackSliderR);

  addWidget( crossFeedbackSliderR = new RSlider (("FeedbackSliderR")) );
  //crossFeedbackSliderR->addListener(this);
  //crossFeedbackSliderR->setRange(-100.0, 100.0, 0.1, 0.0);
  crossFeedbackSliderR->assignParameter( stereoDelayAudioModule->getParameterByName("FeedbackLeftToRight") );
  crossFeedbackSliderR->setSliderName(juce::String(("from left")));
  crossFeedbackSliderR->setDescription(juce::String(("Amount of cross-feedback from left to right")));
  crossFeedbackSliderR->setDescriptionField(infoField);
  crossFeedbackSliderR->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(crossFeedbackSliderR);

  addWidget( outputLabelR = new RTextField( juce::String(("Output:"))) );
  outputLabelR->setDescription(("Output parameters for the right delayline."));
  outputLabelR->setDescriptionField(infoField);

  addWidget( outputSliderR2L = new RSlider (("OutputSliderR2L")) );
  //outputSliderR2L->addListener(this);
  //outputSliderR2L->setRange(-100.0, 100.0, 0.1, 0.0);
  outputSliderR2L->assignParameter( stereoDelayAudioModule->getParameterByName("RightDelayToLeftOut") );
  outputSliderR2L->setSliderName(juce::String(("to left")));
  outputSliderR2L->setDescription(juce::String(("Amount by which the right wet signal goes to the left output")));
  outputSliderR2L->setDescriptionField(infoField);
  outputSliderR2L->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(outputSliderR2L);

  addWidget( outputSliderR2R = new RSlider (("OutputSliderR2R")) );
  //outputSliderR2R->addListener(this);
  //outputSliderR2R->setRange(-100.0, 100.0, 0.1, 100.0);
  outputSliderR2R->assignParameter( stereoDelayAudioModule->getParameterByName("RightDelayToRightOut") );
  outputSliderR2R->setSliderName(juce::String(("to right")));
  outputSliderR2R->setDescription(juce::String(("Amount by which the right wet signal goes to the right output")));
  outputSliderR2R->setDescriptionField(infoField);
  outputSliderR2R->setStringConversionFunction(&percentToStringWithUnit1);
  //automatableSliders.addIfNotAlreadyThere(outputSliderR2R);

  addWidget( outDelaySliderR = new RSlider (("OutputDelaySliderR")) );
  //outDelaySliderR->addListener(this);
  //outDelaySliderR->setRange(0.0, 4.0, 0.0, 0.0);
  outDelaySliderR->assignParameter( stereoDelayAudioModule->getParameterByName("WetDelayInBeatsR") );
  outDelaySliderR->setSliderName(juce::String(("Delay")));
  outDelaySliderR->setDescription(juce::String(("Delay for the right channels wet signal (aka 'Pre-Delay')")));
  outDelaySliderR->setDescriptionField(infoField);
  outDelaySliderR->setStringConversionFunction(&beatsToStringWithUnit4);
  //automatableSliders.addIfNotAlreadyThere(outDelaySliderR);

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(360, 392);
}

/*
StereoDelayModuleEditor::~StereoDelayModuleEditor()
{
deleteAllChildren();
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

/*
void StereoDelayModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
if( stereoDelayAudioModule == NULL )
return;

if( buttonThatWasClicked == reverseButton )
stereoDelayAudioModule->setReversePlayback( reverseButton->getToggleState() );
else if( buttonThatWasClicked == invertButton )
stereoDelayAudioModule->setNegativePolarity( invertButton->getToggleState() );
else if( buttonThatWasClicked == formantPreserveButton )
stereoDelayAudioModule->setFormantPreserve( formantPreserveButton->getToggleState() );
else if( buttonThatWasClicked == antiAliasButton )
stereoDelayAudioModule->setAntiAliasing( antiAliasButton->getToggleState() );
else
{
AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
return;
}
setPresetDirty();
}
*/

void StereoDelayModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  if( stereoDelayAudioModule == NULL )
    return;
  if( stereoDelayAudioModule->wrappedStereoDelay == NULL )
    return;

  if( rComboBoxThatHasChanged == delayComboBoxL )
    stereoDelayAudioModule->wrappedStereoDelay->setDelayInBeats( delayComboBoxL->getValue(), rosic::StereoDelay::LEFT );
  else if( rComboBoxThatHasChanged == delayComboBoxR )
    stereoDelayAudioModule->wrappedStereoDelay->setDelayInBeats( delayComboBoxR->getValue(), rosic::StereoDelay::RIGHT );
  stereoDelayAudioModule->markStateAsDirty();
}

void StereoDelayModuleEditor::rSliderValueChanged(RSlider* sliderThatHasChanged)
{
  // \todo: put this inot the base-class
  if( stereoDelayAudioModule == NULL )
    return;
  stereoDelayAudioModule->markStateAsDirty();
}

void StereoDelayModuleEditor::updateWidgetsAccordingToState()
{
  if( stereoDelayAudioModule == NULL )
    return;
  if( stereoDelayAudioModule->wrappedStereoDelay == NULL )
    return;

  // update the global widgets and automatable sliders:
  AudioModuleEditor::updateWidgetsAccordingToState();

  // update the non-automatable widgets:
  delayComboBoxL->setValue( stereoDelayAudioModule->wrappedStereoDelay
    ->getDelayInBeats(rosic::StereoDelay::LEFT),  false );
  delayComboBoxR->setValue( stereoDelayAudioModule->wrappedStereoDelay
    ->getDelayInBeats(rosic::StereoDelay::RIGHT), false );
}

void StereoDelayModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
  fillRectWithBilinearGradient(g, leftRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  fillRectWithBilinearGradient(g, rightRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  g.setColour(editorColourScheme.outline);
  g.drawRect(leftRectangle);
  g.drawRect(rightRectangle);
}

void StereoDelayModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth()/2;
  //int h = getHeight();

  y  = getHeadlineBottom();
  dryWetSlider->setBounds(x+4, y+8, w-8, 20);
  cutoffScaleSlider->setBounds(x+w+4, y+8, w-8, 20);

  y = cutoffScaleSlider->getBottom()+8;

  // the left rectangle:

  leftRectangle.setBounds(x+4, y, w-8, 332);
  x = leftRectangle.getX();
  w = leftRectangle.getWidth();

  delayLineLabelL->setBounds(x+4, y+4, 88, 20);
  delayComboBoxL->setBounds(delayLineLabelL->getRight(), y+4, w-88-8, 20);
  y = delayLineLabelL->getBottom();
  delayScaleSliderL->setBounds(x+4, y+4, w-8, 16);
  y = delayScaleSliderL->getBottom()+4;

  inputLabelL->setBounds(x+4, y, 88, 20);
  y = inputLabelL->getBottom();
  inputSliderL2L->setBounds(x+4, y, w-8, 16);
  y = inputSliderL2L->getBottom()-2;
  inputSliderR2L->setBounds(x+4, y, w-8, 16);
  y = inputSliderR2L->getBottom()+4;

  diffusorLabelL->setBounds(x+4, y, 88, 20);
  y = diffusorLabelL->getBottom();
  diffusorTimeSliderL->setBounds(x+4, y, w-8, 16);
  y = diffusorTimeSliderL->getBottom()-2;
  diffusorAmountSliderL->setBounds(x+4, y, w-8, 16);
  y = diffusorAmountSliderL->getBottom()+4;

  filterLabelL->setBounds(x+4, y, 88, 20);
  y = filterLabelL->getBottom();
  lowpassSliderL->setBounds(x+4, y, w-8, 16);
  y = lowpassSliderL->getBottom()-2;
  highpassSliderL->setBounds(x+4, y, w-8, 16);
  y = highpassSliderL->getBottom()+4;

  feedbackLabelL->setBounds(x+4, y, 88, 20);
  y = feedbackLabelL->getBottom();
  feedbackSliderL->setBounds(x+4, y, w-8, 16);
  y = feedbackSliderL->getBottom()-2;
  crossFeedbackSliderL->setBounds(x+4, y, w-8, 16);
  y = crossFeedbackSliderL->getBottom()+4;

  outputLabelL->setBounds(x+4, y, 88, 20);
  y = outputLabelL->getBottom();
  outputSliderL2L->setBounds(x+4, y, w-8, 16);
  y = outputSliderL2L->getBottom()-2;
  outputSliderL2R->setBounds(x+4, y, w-8, 16);
  y = outputSliderL2R->getBottom()-2;
  outDelaySliderL->setBounds(x+4, y, w-8, 16);

  // the right rectangle:

  rightRectangle.setBounds(leftRectangle.getRight()+8, leftRectangle.getY(),
    leftRectangle.getWidth(), leftRectangle.getHeight());
  x = rightRectangle.getX();
  y = rightRectangle.getY();

  delayLineLabelR->setBounds(x+4, y+4, 88, 20);
  delayComboBoxR->setBounds(delayLineLabelR->getRight(), y+4, w-88-8, 20);
  y = delayLineLabelR->getBottom();
  delayScaleSliderR->setBounds(x+4, y+4, w-8, 16);
  y = delayScaleSliderR->getBottom()+4;

  inputLabelR->setBounds(x+4, y, 88, 20);
  y = inputLabelR->getBottom();
  inputSliderL2R->setBounds(x+4, y, w-8, 16);
  y = inputSliderL2R->getBottom()-2;
  inputSliderR2R->setBounds(x+4, y, w-8, 16);
  y = inputSliderR2R->getBottom()+4;

  diffusorLabelR->setBounds(x+4, y, 88, 20);
  y = diffusorLabelR->getBottom();
  diffusorTimeSliderR->setBounds(x+4, y, w-8, 16);
  y = diffusorTimeSliderR->getBottom()-2;
  diffusorAmountSliderR->setBounds(x+4, y, w-8, 16);
  y = diffusorAmountSliderR->getBottom()+4;

  filterLabelR->setBounds(x+4, y, 88, 20);
  y = filterLabelR->getBottom();
  lowpassSliderR->setBounds(x+4, y, w-8, 16);
  y = lowpassSliderR->getBottom()-2;
  highpassSliderR->setBounds(x+4, y, w-8, 16);
  y = highpassSliderR->getBottom()+4;

  feedbackLabelR->setBounds(x+4, y, 88, 20);
  y = feedbackLabelR->getBottom();
  feedbackSliderR->setBounds(x+4, y, w-8, 16);
  y = feedbackSliderR->getBottom()-2;
  crossFeedbackSliderR->setBounds(x+4, y, w-8, 16);
  y = crossFeedbackSliderR->getBottom()+4;

  outputLabelR->setBounds(x+4, y, 88, 20);
  y = outputLabelR->getBottom();
  outputSliderR2L->setBounds(x+4, y, w-8, 16);
  y = outputSliderR2L->getBottom()-2;
  outputSliderR2R->setBounds(x+4, y, w-8, 16);
  y = outputSliderR2R->getBottom()-2;
  outDelaySliderR->setBounds(x+4, y, w-8, 16);
}



/*

Ideas:
-Improve the diffusors, perhaps as follows:
 -Start with the general 1st order allpass difference equation:
    y[n] = c*x[n] + x[n-1] - c * y[n-1]   (see DAFX, pg 39)
 -Replace the unit delay with a delay by k samples:
    y[n] = c*x[n] + x[n-k] - c * y[n-k]
 -Maybe allow for non-integer k, using delayline interpolation and use the warped allpass 
  interpolation algorithm or Thiran interpolation for that
    

*/