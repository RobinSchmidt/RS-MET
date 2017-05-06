
//=================================================================================================
// class StandardWaveformRendererAudioModule:

StandardWaveformRendererAudioModule::StandardWaveformRendererAudioModule(CriticalSection *newPlugInLock, 
  rosic::StandardWaveformRenderer *newStandardWaveformRendererToWrap)
  : AudioModule(newPlugInLock)
{
  jassert( newStandardWaveformRendererToWrap != NULL ); // you must pass a valid rosic-object
  wrappedStandardWaveformRenderer = newStandardWaveformRendererToWrap;
  moduleName = juce::String("StandardWaveformRenderer");
  setActiveDirectory(getApplicationDirectory() 
    + juce::File::separatorString + juce::String("StandardWaveformRendererPresets") );
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

WaveformBufferAudioModule::WaveformBufferAudioModule(CriticalSection *newPlugInLock, rosic::WaveformBuffer *newWaveformBufferToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newWaveformBufferToWrap != NULL ); // you must pass a valid rosic-object
  wrappedWaveformBuffer = newWaveformBufferToWrap;
  moduleName = juce::String("WaveformBuffer");
  //setActiveDirectory(getApplicationDirectory() 
  //  + File::separatorString + juce::String(T("WaveformBufferPresets")) );
  AudioFileManager::setActiveDirectory(getApplicationDirectory() + juce::File::separatorString +
    juce::String("Samples") );
  //initializeAutomatableParameters();
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
    juce::File::separatorString + relativePath);
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
  moduleName = juce::String("WaveformRenderer");
  setActiveDirectory(getApplicationDirectory() 
    + juce::File::separatorString + juce::String("WaveformRendererPresets") );

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
  for(int c=0; c<childModules.size(); c++)
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
  p->addStringValue("Algorithm");
  p->addStringValue("MultiSegment");
  addObservedParameter(p); 

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class WaveTableAudioModule:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveTableAudioModule::WaveTableAudioModule(CriticalSection *newPlugInLock, rosic::WaveTable *newWaveTableToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newWaveTableToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedWaveTable = newWaveTableToWrap;
  moduleName = juce::String("WaveTable");
  setActiveDirectory(getApplicationDirectory() 
    + juce::File::separatorString + juce::String("WaveTablePresets") );
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