
//=================================================================================================
// class StandardWaveformRendererAudioModule:

StandardWaveformRendererAudioModule::StandardWaveformRendererAudioModule(CriticalSection *newPlugInLock,
  rosic::StandardWaveformRenderer *newStandardWaveformRendererToWrap)
  : AudioModule(newPlugInLock)
{
  jassert( newStandardWaveformRendererToWrap != NULL ); // you must pass a valid rosic-object
  wrappedStandardWaveformRenderer = newStandardWaveformRendererToWrap;
  setModuleTypeName("StandardWaveformRenderer");
  initializeAutomatableParameters();
}

void StandardWaveformRendererAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedStandardWaveformRenderer == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: wrappedStandardWaveformRenderer->setWaveform( (int) value); break;
  }

  markStateAsDirty();
}

void StandardWaveformRendererAudioModule::initializeAutomatableParameters()
{
  Parameter* p;

  p = new Parameter(lock, "Shape", 0.0, 3.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue("Sine");
  p->addStringValue("Saw");
  p->addStringValue("Square");
  p->addStringValue("Triangle");
  addObservedParameter(p);

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}


//=================================================================================================
// class WaveformBufferAudioModule:

WaveformBufferAudioModule::WaveformBufferAudioModule(CriticalSection *newPlugInLock, 
  rosic::WaveformBuffer *newWaveformBufferToWrap) : AudioModule(newPlugInLock)
{
  jassert( newWaveformBufferToWrap != NULL ); // you must pass a valid rosic-object
  wrappedWaveformBuffer = newWaveformBufferToWrap;
  setModuleTypeName("WaveformBuffer");
  AudioFileManager::setActiveDirectory(getSupportDirectory() 
    + juce::File::getSeparatorString() + "Samples");
}

XmlElement* WaveformBufferAudioModule::getStateAsXml(const juce::String& stateName,
                                                     bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  if( wrappedWaveformBuffer != NULL )
  {
    juce::String path = juce::String(wrappedWaveformBuffer->getWaveformName());
    xmlState->setAttribute("RelativePath", path);
  }
  return xmlState;
}

void WaveformBufferAudioModule::setStateFromXml(const XmlElement& xmlState,
                                                const juce::String& stateName, bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  juce::String path = xmlState.getStringAttribute("RelativePath");
  setWaveformFromFile(path);
}

bool WaveformBufferAudioModule::loadFile(const juce::File& fileToLoad)
{
  return AudioFileManager::loadFile(fileToLoad);
}

bool WaveformBufferAudioModule::saveToFile(const juce::File& fileToSaveTo)
{
  return AudioFileManager::saveToFile(fileToSaveTo);
}

bool WaveformBufferAudioModule::setAudioData(AudioSampleBuffer* newBuffer,
                                             const juce::File& underlyingFile,
                                             bool markAsClean)
{
  juce::String relativePath = underlyingFile.getRelativePathFrom(rootDirectory);
  char* fileNameC = toZeroTerminatedString(relativePath);
  wrappedWaveformBuffer->setWaveform(newBuffer->getWritePointer(0), newBuffer->getNumSamples(),
    fileNameC);
  delete[] fileNameC;
  markStateAsDirty();
  return true;
}

void WaveformBufferAudioModule::setWaveformFromFile(const juce::String &relativePath)
{
  juce::File fileToLoad = juce::File(rootDirectory.getFullPathName() +
    juce::File::getSeparatorString() + relativePath);
  loadFile(fileToLoad);
  /*
  juce::String extension = relativePath.fromLastOccurrenceOf(T("."), false, false);
  if( extension == T("flac") || extension == T("wav") )
  {
    AudioSampleBuffer* buffer =
      AudioFileManager::createAudioSampleBufferFromFile(relativePath, true);
    if( buffer != NULL )
    {
      char* fileNameC = rojue::toZeroTerminatedString(relativePath);
      wrappedWaveformBuffer->setWaveform(buffer->getSampleData(0),
        buffer->getNumSamples(), fileNameC);
      delete[] fileNameC;
      delete buffer;
    }
    else
      wrappedWaveformBuffer->initWaveform();
  }
  else
    wrappedWaveformBuffer->initWaveform();
  markStateAsDirty();
  */
}

//=================================================================================================
// class WaveformRendererAudioModule:

WaveformRendererAudioModule::WaveformRendererAudioModule(CriticalSection *newPlugInLock,
  rosic::WaveformRenderer *newWaveformRendererToWrap)
  : AudioModule(newPlugInLock)
{
  jassert( newWaveformRendererToWrap != NULL ); // you must pass a valid rosic-object
  wrappedWaveformRenderer = newWaveformRendererToWrap;
  setModuleTypeName("WaveformRenderer");

  // child modules:
  standardRendererModule = new StandardWaveformRendererAudioModule(lock, &wrappedWaveformRenderer->standardRenderer);
  addChildAudioModule(standardRendererModule);

  waveformBufferModule = new WaveformBufferAudioModule(lock, &wrappedWaveformRenderer->waveBuffer);
  addChildAudioModule(waveformBufferModule);

  initializeAutomatableParameters(); // must be called after adding child-modules
}

void WaveformRendererAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedWaveformRenderer == NULL )
    return;

  int    index = getIndexOfParameter(parameterThatHasChanged);
  double value = parameterThatHasChanged->getValue();

  // switch state-save/recall off for inactive child-modules:
  //childModules.getLock().enter();
  for(int c = 0; c < size(childModules); c++)
    childModules[c]->setStateSaveAndRecall(false);
  childModules[(int) value]->setStateSaveAndRecall(true);
  //childModules.getLock().exit();

  switch( index )
  {
  case 0: wrappedWaveformRenderer->setMode((int) value); break;
  }

  markStateAsDirty();
}

void WaveformRendererAudioModule::initializeAutomatableParameters()
{
  Parameter* p;

  p = new Parameter(lock, "Mode", 0.0, 3.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue("StandardWaveform");
  p->addStringValue("AudioFile");
  //p->addStringValue("Algorithm");
  //p->addStringValue("MultiSegment");
  addObservedParameter(p);

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class WaveTableAudioModule:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveTableAudioModule::WaveTableAudioModule(CriticalSection *newPlugInLock,
  rosic::WaveTable *newWaveTableToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newWaveTableToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedWaveTable = newWaveTableToWrap;
  setModuleTypeName("WaveTable");
  initializeAutomatableParameters();
  /*
  audioFileManager.setPermissibleWildcardPatterns(juce::String(T("*.wav;*.flac;")));
  audioFileManager.setActiveDirectory(getApplicationDirectory()
    + File::separatorString + juce::String(T("Samples"))
    + File::separatorString + juce::String(T("SingleCycle"))
    + File::separatorString + juce::String(T("Patterns")) );
  */

  rendererModule = new WaveformRendererAudioModule(lock, &wrappedWaveTable->waveformRenderer);
  //rendererModule->addChangeListener(this);
  addChildAudioModule(rendererModule);
}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* WaveTableAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);
  /*
  if( wrappedWaveTable != NULL )
  {
    juce::String path = juce::String(wrappedWaveTable->getWaveformName());
    xmlState->setAttribute(T("WaveformName"), path);
  }
  */
  return xmlState;
}

void WaveTableAudioModule::setStateFromXml(const XmlElement& xmlState,
                                           const juce::String& stateName, bool markAsClean)
{
  wrappedWaveTable->setAutoUpdateOnParameterChange(false);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  wrappedWaveTable->setAutoUpdateOnParameterChange(true);
  wrappedWaveTable->updateBuffers();


  /*
  if( wrappedWaveTable == NULL )
    return;
  juce::String name = xmlState.getStringAttribute(T("WaveformName"), T("Sine"));
  if( name == T("Sine") )
    wrappedWaveTable->setWaveform(rosic::WaveTable::SINE);
  else if( name == T("Triangle") )
    wrappedWaveTable->setWaveform(rosic::WaveTable::TRIANGLE);
  else if( name == T("Square") )
    wrappedWaveTable->setWaveform(rosic::WaveTable::SQUARE);
  else if( name == T("Saw") )
    wrappedWaveTable->setWaveform(rosic::WaveTable::SAW);
  else // ...it must be a custom waveform from a file
    setWaveformFromFile(name);
  */
}

void WaveTableAudioModule::setWaveformFromFile(const juce::String &relativePath)
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
      wrappedWaveTable->setWaveform(buffer->getSampleData(0),
        buffer->getNumSamples(), fileNameC);
      delete[] fileNameC;
      delete buffer;
    }
    else
      wrappedWaveTable->setWaveform(rosic::WaveTable::SILENCE);
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
      wrappedWaveTable->setWaveform(tmpBuffer, 2048, fileNameC);
      delete[] fileNameC;
      delete xmlState;
    }
  }
  else
    wrappedWaveTable->setWaveform(rosic::WaveTable::SILENCE);
   */
}

void WaveTableAudioModule::markStateAsDirty()
{
  wrappedWaveTable->updateBuffers();
  AudioModule::markStateAsDirty();
}

void WaveTableAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedWaveTable == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  //case   0: wrappedWaveTable->setStartPhase(              value); break;
  //case   1: wrappedWaveTable->setRightChannelPhaseOffset( value); break;

  case 0: wrappedWaveTable->setSmoothAttack(            value); break;
  case 1: wrappedWaveTable->setSmoothRelease(           value); break;
    // more to come, order to change...
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void WaveTableAudioModule::initializeAutomatableParameters()
{
  Parameter* p;
  juce::Array<double> defaultValues;

  p = new Parameter(lock, "SmoothAttack", 0.0, 100.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new Parameter(lock, "SmoothRelease", 0.0, 100.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);


  wrappedWaveTable->setAutoUpdateOnParameterChange(false);
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
  wrappedWaveTable->setAutoUpdateOnParameterChange(true);
  wrappedWaveTable->updateBuffers();
}

//=================================================================================================
// class StandardWaveformEditor:

StandardWaveformEditor::StandardWaveformEditor(CriticalSection *newPlugInLock,
  StandardWaveformRendererAudioModule* newRendererModuleToEdit)
  : AudioModuleEditor(newRendererModuleToEdit)
{
  jassert(newRendererModuleToEdit != NULL);
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  setHeadlineStyle(NO_HEADLINE);

  addWidget( shapeComboBox = new RNamedComboBox("ShapeComboBox", "Shape:") );
  shapeComboBox->setDescription("Select one of the standard waveforms");
  shapeComboBox->assignParameter(moduleToEdit->getParameterByName("Shape"));
  shapeComboBox->registerComboBoxObserver(this);

  updateWidgetsAccordingToState();
}

void StandardWaveformEditor::rComboBoxChanged(RComboBox *rComboBoxChanged)
{
  sendChangeMessage();
}

void StandardWaveformEditor::resized()
{
  AudioModuleEditor::resized();
  shapeComboBox->setBounds(4, 4, getWidth()-8, 16);
}

//=================================================================================================
// class WaveformBufferEditor:

WaveformBufferEditor::WaveformBufferEditor(CriticalSection *newPlugInLock, WaveformBufferAudioModule* newWaveformBufferModuleToEdit)
  : AudioModuleEditor(newWaveformBufferModuleToEdit)
{
  jassert( newWaveformBufferModuleToEdit != NULL );
  waveformBufferModuleToEdit = newWaveformBufferModuleToEdit;
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  setHeadlineStyle(NO_HEADLINE);

  fileSelectionBox = new FileSelectionBox("FileComboBox", waveformBufferModuleToEdit);
  addWidgetSet(fileSelectionBox);
  fileSelectionBox->setDescription("Load a custom single-cycle audiofile");
  fileSelectionBox->setSaveButtonVisible(false);
  fileSelectionBox->setLabelPosition(FileSelectionBox::LABEL_ABOVE);
  fileSelectionBox->setButtonsPosition(FileSelectionBox::BUTTONS_ABOVE);
  fileSelectionBox->addChangeListener(this);

  updateWidgetsAccordingToState();
}

void WaveformBufferEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  if( objectThatHasChanged == fileSelectionBox )
    sendChangeMessage();
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void WaveformBufferEditor::resized()
{
  AudioModuleEditor::resized();
  fileSelectionBox->setBounds(4, 4, getWidth()-8, 32);
}

//=================================================================================================
// class WaveformRendererEditor:

WaveformRendererEditor::WaveformRendererEditor(CriticalSection *newPlugInLock,
  WaveformRendererAudioModule* newWaveformRendererModuleToEdit)
  : AudioModuleEditor(newWaveformRendererModuleToEdit)
{
  jassert(newWaveformRendererModuleToEdit != NULL);
  renderer = newWaveformRendererModuleToEdit->wrappedWaveformRenderer;
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  setHeadlineStyle(NO_HEADLINE);

  addWidget( modeComboBox = new RNamedComboBox("ModeComboBox", "Mode:") );
  modeComboBox->setDescription("Select the mode for raw waveform creation");
  modeComboBox->assignParameter(moduleToEdit->getParameterByName("Mode"));
  modeComboBox->registerComboBoxObserver(this);

  standardEditor = new StandardWaveformEditor(lock,
    newWaveformRendererModuleToEdit->standardRendererModule);
  standardEditor->addChangeListener(this);
  addChildEditor(standardEditor);

  bufferEditor = new WaveformBufferEditor(lock,
    newWaveformRendererModuleToEdit->waveformBufferModule);
  bufferEditor->addChangeListener(this);
  addChildEditor(bufferEditor);

  updateWidgetsAccordingToState();
}

void WaveformRendererEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  sendChangeMessage();
}

void WaveformRendererEditor::rComboBoxChanged(RComboBox *rComboBoxChanged)
{
  updateWidgetVisibility();
  sendChangeMessage();
}

void WaveformRendererEditor::resized()
{
  AudioModuleEditor::resized();
  modeComboBox->setBounds(4, 4, getWidth()-8, 16);
  int y = modeComboBox->getBottom()+4;
  standardEditor->setBounds(0, y, getWidth(), getHeight()-y);
  bufferEditor->setBounds(standardEditor->getBounds());
  //...
}

void WaveformRendererEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  standardEditor->updateWidgetsAccordingToState();
  bufferEditor->updateWidgetsAccordingToState();
  //....

  updateWidgetVisibility();
}

void WaveformRendererEditor::updateWidgetVisibility()
{
  standardEditor->setVisible(false);
  bufferEditor->setVisible(false);
  switch( modeComboBox->getSelectedItemIdentifier() )
  {
  case 0: standardEditor->setVisible(true); break;
  case 1: bufferEditor->setVisible(  true); break;
    //...
  }
}

//=================================================================================================
// class WaveTableModuleEditorPopUp:

WaveTableModuleEditorPopUp::WaveTableModuleEditorPopUp(CriticalSection *newPlugInLock,
  WaveTableAudioModule* newWaveTableModuleToEdit)
  : AudioModuleEditor(newWaveTableModuleToEdit)
{
  jassert(newWaveTableModuleToEdit != NULL);
  waveTableModuleToEdit = newWaveTableModuleToEdit;
  waveTableToEdit       = newWaveTableModuleToEdit->wrappedWaveTable;
  setLinkPosition(INVISIBLE);

  rendererEditor = new WaveformRendererEditor(lock, newWaveTableModuleToEdit->rendererModule);
  rendererEditor->addChangeListener(this);
  addChildEditor(rendererEditor);

  numSamplesInPlot = 0;
  xValues          = NULL;
  yValuesL         = NULL;
  yValuesR         = NULL;
  waveformDisplay = new rsDataPlot(juce::String("Plot"));
  waveformDisplay->setDescription(juce::String("Waveform"));
  waveformDisplay->setAxisLabels(juce::String(""), juce::String(""));
  waveformDisplay->setVerticalCoarseGrid(1.0, false);
  waveformDisplay->setHorizontalCoarseGrid(1.0, false);
  waveformDisplay->setAxisPositionX(rsPlotSettings::INVISIBLE);
  waveformDisplay->setAxisPositionY(rsPlotSettings::INVISIBLE);
  addPlot(waveformDisplay);

  addWidget( closeButton = new RButton(RButton::CLOSE) );
  closeButton->setDescription(juce::String("Closes the LFO popup editor"));
  closeButton->setClickingTogglesState(false);
  // we don't listen to this button ourselves - this is the job of the outlying editor object

  updateWidgetsAccordingToState();
  setSize(180, 536);
}

WaveTableModuleEditorPopUp::~WaveTableModuleEditorPopUp()
{
  delete xValues;
  delete yValuesL;
  delete yValuesR;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void WaveTableModuleEditorPopUp::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  //if( waveTableModuleToEdit != NULL && waveTableToEdit != NULL )
  //  waveTableToEdit->renderWaveform();  // it's not very clean to do this here

  updatePlot();
  sendChangeMessage();
}

void WaveTableModuleEditorPopUp::rComboBoxChanged(RComboBox *rComboBoxChanged)
{
  updatePlot();
  sendChangeMessage();
}

void WaveTableModuleEditorPopUp::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  rendererEditor->updateWidgetsAccordingToState();
  updatePlot();
}

void WaveTableModuleEditorPopUp::updatePlot()
{
  if( waveTableModuleToEdit == NULL || waveTableToEdit == NULL )
    return;

  waveTableToEdit->getWaveform(yValuesL, yValuesR, numSamplesInPlot);
  waveformDisplay->setCurveValues(numSamplesInPlot, xValues, yValuesL);
  // \todo: plot separate curves for left and right channel with their actual phases and slewrates
  // - this may require regular updates (via juce::Timer) due to changing BPM
}

void WaveTableModuleEditorPopUp::resized()
{
  AudioModuleEditor::resized();
  int x  = 0;
  int y  = getPresetSectionBottom();
  int h  = getHeight();
  int w  = getWidth();

  closeButton->setBounds(w-16, 0, 16, 16);

  x = 0;
  y = 24;
  w = getWidth()/2;
  h = 80;

  waveformDisplay->invalidatePointers();
  waveformDisplay->setBounds(x, y, w, h);
  numSamplesInPlot = waveformDisplay->getWidth();
  delete[] xValues;
  delete[] yValuesL;
  delete[] yValuesR;
  xValues  = new double[numSamplesInPlot];
  yValuesL = new double[numSamplesInPlot];
  yValuesR = new double[numSamplesInPlot];

  RAPT::rsArrayTools::fillWithIndex(xValues,  numSamplesInPlot);
  RAPT::rsArrayTools::fillWithZeros(yValuesL, numSamplesInPlot);
  RAPT::rsArrayTools::fillWithZeros(yValuesR, numSamplesInPlot);
  waveformDisplay->setMaximumRange(0.0, numSamplesInPlot, -1.1, 1.1);
  waveformDisplay->setCurrentRange(waveformDisplay->getMaximumRange());
  updatePlot();

  y = waveformDisplay->getBottom();
  h = getHeight()-y;
  rendererEditor->setBounds(x, y, w, h);

  //w = getWidth()/2;
  /*
  waveformLabel->setBounds(0, 4, w, 16);
  editButton->setBounds(w-32-4, y+4, 32, 16);
  y = editButton->getBottom();
  */
}

//=================================================================================================
// class WaveTableModuleEditorCompact:


WaveTableModuleEditorCompact::WaveTableModuleEditorCompact(CriticalSection *newPlugInLock,
  WaveTableAudioModule* newWaveTableModuleToEdit)
  : AudioModuleEditor(newWaveTableModuleToEdit)
{
  jassert(newWaveTableModuleToEdit != NULL);
  waveTableModuleToEdit = newWaveTableModuleToEdit;
  waveTableToEdit       = newWaveTableModuleToEdit->wrappedWaveTable;
  setLinkPosition(INVISIBLE);
  setPresetSectionPosition(INVISIBLE);
  setHeadlineStyle(NO_HEADLINE);

  popUpEditor = new WaveTableModuleEditorPopUp(lock, newWaveTableModuleToEdit);
  popUpEditor->addChangeListener(this);
  popUpEditor->setAlwaysOnTop(true);
  popUpEditor->setOpaque(true);
  popUpEditor->closeButton->addRButtonListener(this);
  addChildEditor(popUpEditor, false, false);
  popUpEditorX = -400;
  popUpEditorW =  400+48;
  popUpEditorY =  16;
  popUpEditorH =  200;

  waveformLabel = new RTextField("Waveform:");
  addWidget(waveformLabel);

  numSamplesInPlot = 0;
  xValues          = NULL;
  yValuesL         = NULL;
  yValuesR         = NULL;
  waveformDisplay = new rsDataPlot(juce::String("Plot"));
  waveformDisplay->setDescription(juce::String("Waveform"));
  waveformDisplay->setAxisLabels(juce::String(""), juce::String(""));
  waveformDisplay->setVerticalCoarseGrid(1.0, false);
  waveformDisplay->setHorizontalCoarseGrid(1.0, false);
  waveformDisplay->setAxisPositionX(rsPlotSettings::INVISIBLE);
  waveformDisplay->setAxisPositionY(rsPlotSettings::INVISIBLE);
  addPlot(waveformDisplay);

  addWidget( editButton = new RButton(juce::String("Edit")) );
  editButton->addRButtonListener(this);
  editButton->setDescription(juce::String("Open/close editor for the wavetable"));
  editButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

WaveTableModuleEditorCompact::~WaveTableModuleEditorCompact()
{
  delete popUpEditor; // this is not a child component
  delete xValues;
  delete yValuesL;
  delete yValuesR;
}

//-------------------------------------------------------------------------------------------------
// setup:

void WaveTableModuleEditorCompact::setPopUpEditorBounds(int x, int y, int w, int h)
{
  popUpEditorX = x;
  popUpEditorW = w;
  popUpEditorY = y;
  popUpEditorH = h;
}

void WaveTableModuleEditorCompact::setHeadlineText(const juce::String& newHeadlineText)
{
  AudioModuleEditor::setHeadlineText(newHeadlineText);
  popUpEditor->setHeadlineText(newHeadlineText);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void WaveTableModuleEditorCompact::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( waveTableToEdit == NULL )
    return;

  if( buttonThatWasClicked == editButton )
  {
    if( editButton->getToggleState() == true )
    {
      int x = editButton->getScreenX();
      int y = editButton->getScreenY();
      popUpEditor->setBounds(x+popUpEditorX, y+popUpEditorY, popUpEditorW, popUpEditorH);
      popUpEditor->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
      popUpEditor->setVisible(true);
      popUpEditor->toFront(true);
    }
    else
      popUpEditor->setVisible(false);
    return;
  }
  else if( buttonThatWasClicked == popUpEditor->closeButton )
  {
    editButton->setToggleState(false, true);
    return;
  }

  moduleToEdit->markStateAsDirty();  // do we need this?
}

void WaveTableModuleEditorCompact::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  if( waveTableToEdit == NULL )
    return;
  if( objectThatHasChanged == popUpEditor )
  {
    moduleToEdit->markStateAsDirty();
    updatePlot();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void WaveTableModuleEditorCompact::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  popUpEditor->updateWidgetsAccordingToState();
  updatePlot();
}

void WaveTableModuleEditorCompact::updatePlot()
{
  if( waveTableModuleToEdit == NULL || waveTableToEdit == NULL )
    return;

  waveTableToEdit->getWaveform(yValuesL, yValuesR, numSamplesInPlot);
  waveformDisplay->setCurveValues(numSamplesInPlot, xValues, yValuesL);
  // \todo: plot separate curves for left and right channel with their actual phases and slewrates
  // - this may require regular updates (via juce::Timer) due to changing BPM
}

void WaveTableModuleEditorCompact::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  waveformLabel->setBounds(0, 4, w, 16);
  editButton->setBounds(w-32-4, y+4, 32, 16);
  y = editButton->getBottom() - RWidget::outlineThickness;

  waveformDisplay->invalidatePointers();
  waveformDisplay->setBounds(x, y, w, h-y);
  numSamplesInPlot = waveformDisplay->getWidth();
  delete[] xValues;
  delete[] yValuesL;
  delete[] yValuesR;
  xValues  = new double[numSamplesInPlot];
  yValuesL = new double[numSamplesInPlot];
  yValuesR = new double[numSamplesInPlot];

  RAPT::rsArrayTools::fillWithIndex(xValues,  numSamplesInPlot);
  RAPT::rsArrayTools::fillWithZeros(yValuesL, numSamplesInPlot);
  RAPT::rsArrayTools::fillWithZeros(yValuesR, numSamplesInPlot);
  waveformDisplay->setMaximumRange(0.0, numSamplesInPlot, -1.1, 1.1);
  waveformDisplay->setCurrentRange(waveformDisplay->getMaximumRange());
  updatePlot();


  //y = waveformDisplay->getBottom();
  /*
  waveformComboBox->setBounds(x, y, w-28, 16);
  x = waveformComboBox->getRight()-2;
  waveMinusButton->setBounds(x, y, 16, 16);
  x += 14;
  wavePlusButton->setBounds(x, y, 16, 16);
  x  = 0;
  y += 16;
  cycleLengthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  tempoSyncButton->setBounds(x+4, y+4, 32, 16);
  triggerButton->setBounds(tempoSyncButton->getRight()+4, y+4, 52, 16);
  y += 24;
  startPhaseSlider->setBounds(x+4, y+4, w/2-8, 16);
  upSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 14;
  stereoPhaseSlider->setBounds(x+4, y+4, w/2-8, 16);
  downSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  */
}
