
//-------------------------------------------------------------------------------------------------
// construction/destruction:

LowFrequencyOscillatorAudioModule::LowFrequencyOscillatorAudioModule(CriticalSection *newPlugInLock, 
  rosic::LowFrequencyOscillator *newLowFrequencyOscillatorToWrap)
  : AudioModule(newPlugInLock)
{
  jassert( newLowFrequencyOscillatorToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedLowFrequencyOscillator = newLowFrequencyOscillatorToWrap;
  moduleName = juce::String("LowFrequencyOscillator");
  setActiveDirectory(getApplicationDirectory() 
    + juce::File::getSeparatorString() + juce::String("Presets/LowFrequencyOscillator"));

  /*
  audioFileManager.setPermissibleWildcardPatterns(juce::String(T("*.wav;*.flac;*.xml")));
  audioFileManager.setActiveDirectory(getApplicationDirectory() 
    + File::separatorString + juce::String(T("Samples")) 
    + File::separatorString + juce::String(T("SingleCycle"))
    + File::separatorString + juce::String(T("Patterns")) );  
    */

  waveTableModule = new WaveTableAudioModule(lock, wrappedLowFrequencyOscillator->waveTable);
  addChildAudioModule(waveTableModule);

  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* LowFrequencyOscillatorAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  /*
  if( wrappedLowFrequencyOscillator != NULL )
  {
    juce::String path = juce::String(wrappedLowFrequencyOscillator->getWaveformName());
    xmlState->setAttribute(T("WaveformName"), path);
  }
  */
  return xmlState;
}

void LowFrequencyOscillatorAudioModule::setStateFromXml(const XmlElement& xmlState,
                                                 const juce::String& stateName, bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  /*
  if( wrappedLowFrequencyOscillator == NULL )
    return;
  juce::String name = xmlState.getStringAttribute(T("WaveformName"), T("Sine"));
  if( name == T("Sine") )
  {
    wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::SINE);
  }
  else if( name == T("Triangle") )
  {
    wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::TRIANGLE);
  }
  else if( name == T("Square") )
  {
    wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::SQUARE);
  }
  else if( name == T("Saw") )
  {
    wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::SAW);
  }
  else // ...it must be a custom waveform from a file 
  {
    setWaveformFromFile(name); 
  }
  */
}

void LowFrequencyOscillatorAudioModule::setWaveformFromFile(const juce::String &relativePath)
{
  /*
  juce::String extension = relativePath.fromLastOccurrenceOf(T("."), false, false);
  if( extension == T("flac") || extension == T("wav") )
  {
    AudioSampleBuffer* buffer = 
      AudioFileManager::createAudioSampleBufferFromFile(relativePath, true);
    if( buffer != NULL )
    {
      char* fileNameC = rojue::toZeroTerminatedString(relativePath);
      wrappedLowFrequencyOscillator->setWaveform(buffer->getSampleData(0), 
        buffer->getNumSamples(), fileNameC);
      delete[] fileNameC;
      delete buffer;
    }
    else 
      wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::SILENCE);
  }
  else if( extension == T("xml") )
  {
    juce::String absolutePath = getApplicationDirectory()+File::separatorString+relativePath;
    XmlElement *xmlState = rojue::getXmlFromFile(absolutePath);
    if( xmlState != NULL )
    {
      rosic::BreakpointModulator tmpModulator;
      rosof::breakpointModulatorStateFromXml(&tmpModulator, *xmlState);
      double tmpBuffer[2048];
      tmpModulator.fillBufferWithEnvelope(tmpBuffer, 2048, tmpModulator.getLoopMode() != 0);
      rosic::fitIntoRange(tmpBuffer, 2048, -1.0, 1.0);
      char* fileNameC = rojue::toZeroTerminatedString(relativePath);
      wrappedLowFrequencyOscillator->setWaveform(tmpBuffer, 2048, fileNameC);
      delete[] fileNameC;
      delete xmlState;
    }
  }
  else
    wrappedLowFrequencyOscillator->setWaveform(rosic::WaveTable::SILENCE);
  */
}

//-------------------------------------------------------------------------------------------------
// automation:

void LowFrequencyOscillatorAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedLowFrequencyOscillator == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedLowFrequencyOscillator->setCycleLength(      value);        break;
  case   1: wrappedLowFrequencyOscillator->setTempoSync(        value>=0.5);   break;
  case   2: wrappedLowFrequencyOscillator->setStartPhase(       value);        break;
  case   3: wrappedLowFrequencyOscillator->setStereoPhase(      value);        break;

  //case   4: wrappedLowFrequencyOscillator->setUpwardSlewRate(   value);        break;
  //case   5: wrappedLowFrequencyOscillator->setDownwardSlewRate( value);        break;
  } 
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void LowFrequencyOscillatorAudioModule::initializeAutomatableParameters()
{
  typedef ModulatableParameter Param;
  Param* p;

  std::vector<double> defaultValues;

  p = new Param("CycleLength", 0.125, 1.0, 0.25, Parameter::LINEAR, 0.0125);
  defaultValues.clear();
  defaultValues.push_back(0.125);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.75);
  defaultValues.push_back(1.0);
  defaultValues.push_back(1.5);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(8.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(16.0);
  defaultValues.push_back(24.0);
  defaultValues.push_back(32.0);
  defaultValues.push_back(48.0);
  defaultValues.push_back(64.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 

  p = new Param("TempoSync", 0.0, 1.0, 1.0, Parameter::BOOLEAN, 1.0);
  addObservedParameter(p);

  p = new Param("StartPhase", 0.0, 360.0, 0.0, Parameter::LINEAR, 1.0);
  defaultValues.clear();
  defaultValues.push_back(0.0);
  defaultValues.push_back(45.0);
  defaultValues.push_back(90.0);
  defaultValues.push_back(135.0);
  defaultValues.push_back(180.0);
  defaultValues.push_back(225.0);
  defaultValues.push_back(270.0);
  defaultValues.push_back(315.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 

  p = new Param("StereoPhase", 0.0, 180.0, 0.0, Parameter::LINEAR, 1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 

  // maybe include a rise-time...

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);


  /*
  // the old stuff - todo: copy default-values
  Parameter* p;
  p = new Parameter("CycleLength", 0.125, 64.0, 0.0125, 0.25, Parameter::LINEAR);
  addObservedParameter(p); 

  p = new Parameter("Amount", -200.0, 200.0, 1.0, 100.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.add(-200.0);
  defaultValues.add(-100.0);
  defaultValues.add(0.0);
  defaultValues.add(100.0);
  defaultValues.add(200.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 

  p = new Parameter("StartPhase", 0.0, 360.0, 1.0, 0.0, Parameter::LINEAR);
  defaultValues.clear();
  defaultValues.add(0.0);
  defaultValues.add(45.0);
  defaultValues.add(90.0);
  defaultValues.add(135.0);
  defaultValues.add(180.0);
  defaultValues.add(225.0);
  defaultValues.add(270.0);
  defaultValues.add(315.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 

  p = new Parameter("TempoSync", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  p = new Parameter("UpwardJumpTime", 1.0, 100.0, 0.1, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Parameter("DownwardJumpTime", 1.0, 100.0, 0.1, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Parameter("RiseTime", 1.0, 10000.0, 0.1, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new Parameter("FullWaveWarp", -0.99, 0.99, 0.0001, 0.0, Parameter::LINEAR, -1, false);
  addObservedParameter(p);
  p = new Parameter("HalfWaveWarp", -0.99, 0.99, 0.0001, 0.0, Parameter::LINEAR, -1, false);
  addObservedParameter(p);

  p = new Parameter("SpectralContrast", 0.25, 4.0, 0.0, 1.0, Parameter::EXPONENTIAL, -1, false);
  addObservedParameter(p);
  p = new Parameter("SpectralSlope", -6.0, 6.0, 0.0, 0.0, Parameter::LINEAR, -1, false);
  addObservedParameter(p);
  p = new Parameter("HighestHarmonic", 1.0, 1024.0, 1.0, 1024.0, Parameter::EXPONENTIAL, -1,false);
  addObservedParameter(p);
  p = new Parameter("LowestHarmonic", 1.0, 1024.0, 1.0, 1.0, Parameter::EXPONENTIAL, -1,false);
  addObservedParameter(p);
  p = new Parameter("EvenOddRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR, -1, false);
  addObservedParameter(p);

  //-----------------------------------------------------------------------------------------------
  // phase spectrum related parameters:

  // #018:
  p = new Parameter("EvenOddPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR, -1, false);
  defaultValues.clear();
  defaultValues.add(-90.0);
  defaultValues.add(-45.0);
  defaultValues.add(0.0);
  defaultValues.add(45.0);
  defaultValues.add(90.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #019:
  p = new Parameter("StereoPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR, -1, false);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  // #020:
  p = new Parameter("EvenOddStereoPhaseShift", -180.0, 180.0, 1.0, 0.0, Parameter::LINEAR, -1, 
                    false);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);


  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
    */
}

//=================================================================================================

LowFrequencyOscillatorEditor::LowFrequencyOscillatorEditor(CriticalSection *newPlugInLock, 
  LowFrequencyOscillatorAudioModule* newOscillatorModuleToEdit)                                                       
  : AudioModuleEditor(newOscillatorModuleToEdit)
{
  jassert(newOscillatorModuleToEdit != NULL);  
  oscillatorModuleToEdit = newOscillatorModuleToEdit;
  oscillatorToEdit       = oscillatorModuleToEdit->wrappedLowFrequencyOscillator;
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(BELOW_HEADLINE);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setPresetLabelVisible(false);

  addWidget( cycleLengthSlider = new RSlider("CycleLengthSlider") );
  cycleLengthSlider->assignParameter( moduleToEdit->getParameterByName("CycleLength") );
  cycleLengthSlider->setSliderName(juce::String("Cycle"));
  cycleLengthSlider->setDescription(juce::String("Length of one cycle (in seconds or beats)"));
  cycleLengthSlider->setDescriptionField(infoField);
  cycleLengthSlider->setStringConversionFunction(beatsToStringWithUnit4);

  addWidget( depthSlider = new RSlider("DepthSlider") );
  //depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  //depthSlider->setDescription(juce::String("Depth of the LFO modulation (in percent)"));
  depthSlider->setDescription(juce::String("Not yet implemented"));
  depthSlider->setStringConversionFunction(percentToStringWithUnit1);

  addWidget( tempoSyncButton = new RButton(juce::String("TempoSync")) );
  //tempoSyncButton->assignParameter( moduleToEdit->getParameterByName("TempoSync") );
  tempoSyncButton->setButtonText(juce::String("Sync"));
  //tempoSyncButton->setDescription(juce::String("Toggle tempo synchronization on/off"));
  tempoSyncButton->setDescription(juce::String("Not yet implemented"));
  tempoSyncButton->setDescriptionField(infoField);
  tempoSyncButton->setClickingTogglesState(true);
  tempoSyncButton->addRButtonListener(this);

  waveTableEditor = new WaveTableModuleEditorCompact(lock, 
    oscillatorModuleToEdit->waveTableModule);
  //waveTableEditor->addChangeListener(this);
  addChildEditor(waveTableEditor);

  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void LowFrequencyOscillatorEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{

  /*
  if( oscillatorModuleToEdit == NULL )
  return;
  if( oscillatorModuleToEdit->wrappedLowFrequencyOscillator == NULL )
  return;
  rosic::LowFrequencyOscillator* o = oscillatorModuleToEdit->wrappedLowFrequencyOscillator;

  sendChangeMessage();
  */

}

void LowFrequencyOscillatorEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  waveTableEditor->updateWidgetsAccordingToState();
}

void LowFrequencyOscillatorEditor::resized()
{
  AudioModuleEditor::resized();

  stateWidgetSet->setBounds(4, 8, getWidth()-8, 32);

  int x  = 0;
  int y  = getPresetSectionBottom();
  int h  = getHeight();
  int w  = getWidth();

  cycleLengthSlider->setBounds(x+4, y, w-8, 16);
  y += 14;
  tempoSyncButton->setBounds(cycleLengthSlider->getRight()-40, y, 40, 16);
  y = tempoSyncButton->getBottom();

  y = tempoSyncButton->getBottom()+4;
  //depthSlider->setBounds(x+4, y+4, w-8, 16);
  //y = depthSlider->getBottom()+4;
  waveTableEditor->setBounds(x, y, w, h-y);
}
