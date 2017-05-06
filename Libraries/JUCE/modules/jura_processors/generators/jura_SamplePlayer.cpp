

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SamplePlayerAudioModule::SamplePlayerAudioModule(CriticalSection *newPlugInLock, rosic::SamplePlayer *newSamplePlayerToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newSamplePlayerToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedSamplePlayer = newSamplePlayerToWrap;
  moduleName = juce::String("SamplePlayer");
  setActiveDirectory(getApplicationDirectory() + juce::String("/SamplePlayerPresets") );
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void SamplePlayerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedSamplePlayer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedSamplePlayer->setMute(                     value >= 0.5);  break;
  case  1: wrappedSamplePlayer->setSolo(                     value >= 0.5);  break;
  case  2: wrappedSamplePlayer->setLevel(                    value);         break;
  case  3: wrappedSamplePlayer->setLevelByKey(               value);         break;
  case  4: wrappedSamplePlayer->setLevelByVel(               value);         break;
  //case  5: wrappedSamplePlayer->setMidSide(                  value);         break;
  //case  6: wrappedSamplePlayer->setPan(                      value);         break;
  case 10: wrappedSamplePlayer->setRootKey(                  value);         break;
  case  7: wrappedSamplePlayer->setTune(                     value);         break;
  case  8: wrappedSamplePlayer->setTuneByKey(                value);         break;
  case  9: wrappedSamplePlayer->setTuneByVel(                value);         break;
  //case 11: wrappedSamplePlayer->setPlaybackStart(            value);         break;
  case 12: wrappedSamplePlayer->setLoopMode(                 value >= 0.5);  break;
  case 13: wrappedSamplePlayer->setLoopStart(                value);         break;
  case 14: wrappedSamplePlayer->setLoopEnd(                  value);         break;


  //case 15: wrappedSamplePlayer->setLowpassCutoff(            value);         break;
  //case 16: wrappedSamplePlayer->setHighpassCutoff(           value);         break;
  //case 17: wrappedSamplePlayer->setPhaseRandomize(           value >= 0.5);  break;
  //case 18: wrappedSamplePlayer->setPhaseRandomizeSeed( (int) value);         break;
  } // end of switch( parameterIndex )

}

//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* samplePlaybackParametersToXml(SamplePlaybackParameters* parameters, 
  XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the releveant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("PlaybackParameters")); 
  else
    xmlState = xmlElementToStartFrom;

  // general parameters:
  xmlState->setAttribute("AudioFileRelativePath", juce::String(parameters->getSampleName())     );

  // play start/end parameters:
  xmlState->setAttribute("PlayStart",             juce::String(parameters->getPlaybackStart())  );
  if( parameters->getPlaybackStartByKey() != 0.0 )
    xmlState->setAttribute("PlayStartByKey",  juce::String(parameters->getPlaybackStartByKey()) );
  if( parameters->getPlaybackStartByVel() != 0.0 )
    xmlState->setAttribute("PlayStartByVel",  juce::String(parameters->getPlaybackStartByVel()) );
  xmlState->setAttribute("PlayEnd",               juce::String(parameters->getPlaybackEnd())    );

  // loop parameters:
  if( parameters->getLoopMode() == rosic::SamplePlaybackParameters::FORWARD_LOOP )
    xmlState->setAttribute("LoopMode", juce::String("Forward") );
  else
    xmlState->setAttribute("LoopMode", juce::String("NoLoop") );
  xmlState->setAttribute("LoopStart",        juce::String(parameters->getLoopStart())             );
  xmlState->setAttribute("LoopEnd",          juce::String(parameters->getLoopEnd())               );
  xmlState->setAttribute("PitchCyclesInLoop",juce::String(parameters->getNumPitchCyclesInLoop())  );

  // ...crossfade to come

  // key-/vel mapping parameters:
  xmlState->setAttribute("RootKey",     juce::String(parameters->getRootKey())     );
  //xmlState->setAttribute(T("RootDetune"),  juce::String(parameters->getRootDetune())  );
  if( parameters->getLoKey() != 0 )
    xmlState->setAttribute("LoKey",  juce::String(parameters->getLoKey()) );
  if( parameters->getHiKey() != 127 )
    xmlState->setAttribute("HiKey",  juce::String(parameters->getHiKey()) );
  if( parameters->getLoVel() != 0 )
    xmlState->setAttribute("LoVel",  juce::String(parameters->getLoVel()) );
  if( parameters->getHiVel() != 127 )
    xmlState->setAttribute("HiVel",  juce::String(parameters->getHiVel()) );
  //...keyfade to come


  return xmlState;
}

bool samplePlaybackParametersFromXml(SamplePlaybackParameters* parameters, 
  const XmlElement &xmlState)
{
  bool success = true;

  // set up the audio-file name:
  //juce::String samplePath = xmlState.getStringAttribute(T("AudioFileRelativePath"), juce::String::empty);
  //parameters->setSampleName(toZeroTerminatedString(samplePath));

  // play start/end parameters:
  parameters->setPlaybackStart(     xmlState.getDoubleAttribute("PlayStart",        0.0)   );
  parameters->setPlaybackStartByKey(xmlState.getDoubleAttribute("PlayStartByKey",   0.0)   );
  parameters->setPlaybackStartByVel(xmlState.getDoubleAttribute("PlayStartByVel",   0.0)   );
  parameters->setPlaybackEnd(       xmlState.getDoubleAttribute("PlayEnd",          0.0)   );

  // loop parameters:
  if( xmlState.getStringAttribute("LoopMode", "NoLoop") == juce::String("Forward") )
    parameters->setLoopMode( rosic::SamplePlaybackParameters::FORWARD_LOOP );
  else
    parameters->setLoopMode( rosic::SamplePlaybackParameters::NO_LOOP );
  parameters->setLoopStart(xmlState.getDoubleAttribute("LoopStart",    0.0)   );
  parameters->setLoopEnd(  xmlState.getDoubleAttribute("LoopEnd",      0.0)   );
  parameters->setNumPitchCyclesInLoop(  
    xmlState.getDoubleAttribute("PitchCyclesInLoop", 1.0)   );

  // key-/vel mapping parameters:
  parameters->setRootKey(   xmlState.getIntAttribute(   "RootKey",         64)    );
  //parameters->setRootDetune(xmlState.getDoubleAttribute("RootDetune",      0.0)   );
  parameters->setLoKey(     xmlState.getIntAttribute(   "LoKey",           0)     );
  parameters->setHiKey(     xmlState.getIntAttribute(   "HiKey",           127)   );
  parameters->setLoVel(     xmlState.getIntAttribute(   "LoVel",           0)     );
  parameters->setHiVel(     xmlState.getIntAttribute(   "HiVel",           127)   );

  return true;
}

XmlElement* samplePlayerStateToXml(SamplePlayer* player, XmlElement* xmlElementToStartFrom)
{

  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("SamplePlayer")); 
  else
    xmlState = xmlElementToStartFrom;

  XmlElement* playbackParametersState = samplePlaybackParametersToXml(player->parameters, NULL);
  xmlState->addChildElement(playbackParametersState);

  return xmlState;
}

bool samplePlayerStateFromXml(SamplePlayer* player, const XmlElement &xmlState)
{

  XmlElement* playbackParametersState = xmlState.getChildByName("PlaybackParameters");
  if( playbackParametersState == NULL )
    return false;

  // retrieve the audio-file name:
  juce::String samplePath = playbackParametersState->getStringAttribute(
    "AudioFileRelativePath", juce::String::empty);

  // try to load the audio data
  AudioSampleBuffer* buffer = AudioFileManager::createAudioSampleBufferFromFile(samplePath, true);
  if( buffer != NULL )
  {
    // pass the actual audio data:
    float* channelPointers[2];
    channelPointers[0] = buffer->getWritePointer(0, 0);
    if( buffer->getNumChannels() >= 2 )
      channelPointers[1] = buffer->getWritePointer(1, 0);
    else
      channelPointers[1] = buffer->getWritePointer(0, 0);
    player->setSampleData(channelPointers, buffer->getNumSamples(), buffer->getNumChannels() );

    // set up the audio-file name in the SamplePlayer object:
    char* nameC = toZeroTerminatedString(samplePath);
    if( nameC != NULL )
    {
      player->parameters->setSampleName(nameC);
      delete[] nameC;
    }

    // set up the playback parameters:
    samplePlaybackParametersFromXml(player->parameters, *playbackParametersState);

    delete buffer;
    return true;
  }
  else
    return false;
}

XmlElement* SamplePlayerAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // store the parameters of the underlying core object:
  if( wrappedSamplePlayer != NULL )
    xmlState = samplePlayerStateToXml(wrappedSamplePlayer, xmlState);

  return xmlState;
}

void SamplePlayerAudioModule::setStateFromXml(const XmlElement& xmlState,
                                              const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // restore the parameters of the underlying core object:
  if( wrappedSamplePlayer != NULL )
    samplePlayerStateFromXml(wrappedSamplePlayer, xmlState);
}

//-------------------------------------------------------------------------------------------------
// others:

void SamplePlayerAudioModule::setRootKeyFromLoop()
{
  if( wrappedSamplePlayer == NULL )
    return;

  // get the length of the loop in seconds:
  double loopLength = wrappedSamplePlayer->parameters->getLoopLength() / 
    wrappedSamplePlayer->parameters->getRecordingSampleRate() ; 

  // get the corresponding frequency:
  if( wrappedSamplePlayer->parameters->getNumPitchCyclesInLoop() < 1.0 )
    return;
  double freq = wrappedSamplePlayer->parameters->getNumPitchCyclesInLoop() / loopLength;

  // get the corresponding pitch and set it up as our 'RootKey' parameter:
  double rootKey = freqToPitch(freq);
  getParameterByName(juce::String("RootKey"))->setValue(rootKey, true, true);
}

bool SamplePlayerAudioModule::setSampleFromFile(const juce::File &fileToLoadFrom)
{
  if( wrappedSamplePlayer == NULL )
    return false;

  AudioSampleBuffer* buffer 
    = AudioFileManager::createAudioSampleBufferFromFile(fileToLoadFrom, false);

  if( buffer == NULL )
    return false;

  int numSamples  = buffer->getNumSamples();
  int numChannels = buffer->getNumChannels();
  if( buffer != NULL && numChannels > 0 && numSamples > 0 )
  {
    // pass the actual data as well as the filename (the relative path of the file) over to the
    // audio-engine:
    float* channelPointers[2];
    channelPointers[0] = buffer->getWritePointer(0, 0);
    if( numChannels >= 2 )
      channelPointers[1] = buffer->getWritePointer(1, 0);
    else
      channelPointers[1] = buffer->getWritePointer(0, 0);
    wrappedSamplePlayer->setSampleData(channelPointers, numSamples, numChannels);
    delete buffer;

    // update the filename in the audio-engine:
    juce::String relativePath = fileToLoadFrom.getRelativePathFrom(getApplicationDirectory());
    /*
    long  length = relativePath.length();
    char* fileNameC = new char[length+1];
    relativePath.copyToBuffer(fileNameC, length);
    wrappedSamplePlayer->setSampleName(fileNameC);
    */
    char* fileNameC = toZeroTerminatedString(relativePath);
    wrappedSamplePlayer->setSampleName(fileNameC);
    delete[] fileNameC;

    // update the parameter ranges that depend on the length of the sample:
    getParameterByName(juce::String("PlaybackStart"))->setRange(0.0, numSamples);
    getParameterByName(juce::String("LoopStart"))->setRange(0.0, numSamples);
    getParameterByName(juce::String("LoopLength"))->setRange(0.0, numSamples);

    return true;
  }
  else
  {
    delete buffer;
    return false;
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void SamplePlayerAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created parameter-objects:
  AutomatableParameter* p;

  p = new AutomatableParameter(lock, juce::String("Mute"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Solo"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Level"), -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("LevelByKey"), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("LevelByVel"), -12.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("MidSide"),    0.0,  1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Pan"),       -1.0,  1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, juce::String("Tune"), -36.0, 36.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("TuneByKey"), -200.0, 200.0, 0.1, 100.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("TuneByVel"), -200.0, 200.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("RootKey"), 0.0, 127.0, 0.01, 64.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, juce::String("PlaybackStart"), 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Loop"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("LoopStart"), 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("LoopLength"), 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, juce::String("FilterActive"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Lowpass"), 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("Highpass"), 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, juce::String("PhaseRandomize"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, juce::String("PhaseRandomizationSeed"), 0.0, 1000.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

SamplePlayerEditorDisplay::SamplePlayerEditorDisplay(AudioFileBuffer *newBufferToUse) 
  : WaveformDisplay(newBufferToUse) //, InteractiveCoordinateSystem(juce::String(T("SampleDisplay")))
{

  setValueFieldPopup(false);

  //bufferToEdit             = NULL;
  //playbackParametersToEdit = NULL;
  samplePlayerToEdit = NULL;
  locatorBeingDragged      = -1;

  setAxisLabelX(juce::String::empty);
  setAxisLabelY(juce::String::empty);
  setAxisPositionX(CoordinateSystem::INVISIBLE);
  setAxisPositionY(CoordinateSystem::INVISIBLE);

  /*
  lockUsedBufferPointer();
  bufferToUse = &buffer;
  unlockUsedBufferPointer();
  */
  assignAudioFileBuffer(&buffer);

  //buffer.registerUser(this);
  //scaleFactor          = 1.0f;
  //dcOffset             = 0.0f;
}

SamplePlayerEditorDisplay::~SamplePlayerEditorDisplay()
{
  assignAudioFileBuffer(NULL);
  //buffer.deRegisterUser(this);
  //deleteAllChildren();
}


//-------------------------------------------------------------------------------------------------
// 

/*
CoordinateSystemRangeOld SamplePlayerEditorDisplay::getMaximumMeaningfulRange(
double relativeMarginLeft, double relativeMarginRight, 
double relativeMarginTop,  double relativeMarginBottom)
{
return WaveformDisplayOld::getMaximumMeaningfulRange(relativeMarginLeft, relativeMarginRight,
relativeMarginTop, relativeMarginBottom);
}
*/

void SamplePlayerEditorDisplay::setSamplePlayerToEdit(rosic::SamplePlayer *newPlayerToEdit)
{
  samplePlayerToEdit = newPlayerToEdit;

  if( samplePlayerToEdit == NULL )
    return;

  // update the waveform display:
  int numSamples  = samplePlayerToEdit->getNumSamples();
  int numChannels = samplePlayerToEdit->getNumChannels();
  setMaximumRange(0.0, jmax((double)numSamples,1.0), -1.1, 1.1);
  setCurrentRange(0.0, jmax((double)numSamples,1.0), -1.1, 1.1);


  //buffer.setDataToReferTo( samplePlayerToEdit->getSampleData(), numSamples, numChannels );

  /*
  lockUsedBufferPointer();
  if( bufferToUse != NULL )
  {
  juce::String audioFile = getApplicationDirectory() + File::separatorString 
  + juce::String(samplePlayerToEdit->getSampleName());
  bufferToUse->loadAudioDataFromFile(audioFile, true);
  WaveformDisplay::setDirty();
  repaint(); // should not be necesarry?
  }
  unlockUsedBufferPointer();
  */


  juce::String audioFile = getApplicationDirectory() + juce::File::separatorString   
    + juce::String(samplePlayerToEdit->getSampleName());
  setAudioFileToUse(juce::File(audioFile));

  //setWaveform(samplePlayerToEdit->getSampleData(), numSamples, numChannels );
}

bool SamplePlayerEditorDisplay::setAudioFileToUse(const juce::File &newFileToUse)
{
  lockUsedBufferPointer();
  if( bufferToUse != NULL )
  {
    bufferToUse->loadAudioDataFromFile(newFileToUse, false);
    setRangeToBufferLength();

    // pass the data to the audio core-object:
    if( samplePlayerToEdit != NULL )
    {
      bufferToUse->acquireReadLock();

      // todo: retrieve the pointers to the channels and pass it to the rosic-object, this might 
      // also be the place to render samples from modal models etc.

      //bufferToUse->getSampleData(0, 0);
      //samplePlayerToEdit->setSampleData(channelPointers, numSamples numChannles);

      samplePlayerToEdit->parameters->setRecordingSampleRate(bufferToUse->getFileSampleRate());

      bufferToUse->releaseReadLock();
    }

    //setMaximumRangeX(0.0, bufferToUse->getNumSamples()+1);
    //setCurrentRangeX(0.0, bufferToUse->getNumSamples()+1);
    setDirty();
    //updatePlot();
    //repaint(); // should not be necesarry?
  }
  unlockUsedBufferPointer();

  return true; // preliminary
}

/*
void SamplePlayerEditorDisplay::setSamplePlaybackParametersToEdit(
rosic::SamplePlaybackParameters *newParametersToEdit)
{
playbackParametersToEdit = newParametersToEdit;
}
*/

/*
bool SamplePlayerEditorDisplay::setWaveform(double** newWaveformData, int newNumSampleFrames, 
int newNumChannels)
{
bool success = false;

if( modulatorToEdit == NULL )
return false;

// call the setWavefrom-function from the WavefromDisplay base-class - this will result in a 
// copying of the passed data into the (inherited) member 'peakData' 
// (which is pointer to float):
success = WaveformDisplay::setWaveform(newWaveformData, newNumSampleFrames, newNumChannels);

// the float* peakData member now contains the data (plus possibly some padding and decimated 
// data but thsi additional stuff will come later in the array and can therefore be ignored) 
// - we can use this pointer to further pass the data through to the rosic::SampleModulator 
// object which accepts a pointer to float:
if( success == true )
modulatorToEdit->sampleModulator.setSampleData(peakData, newNumSampleFrames);

return success;
}

bool SamplePlayerEditorDisplay::setWaveform(const AudioSampleBuffer& newWaveformBuffer)
{
// similar procedure to the function above:
bool success = false;
if( modulatorToEdit == NULL )
return false;
success = WaveformDisplay::setWaveform(newWaveformBuffer);
if( success == true )
modulatorToEdit->sampleModulator.setSampleData(peakData, numSampleFrames);

return success;
}
*/

void SamplePlayerEditorDisplay::updatePlot(bool resetZoomFactor)
{
  if( samplePlayerToEdit == NULL )
    return;

  /*
  // aquire a lock on the sample data - just in case some other thread corrupts it in between calls
  // getNumSamples/getNumChannels and the actual data-readout via getSampleData():
  samplePlayerToEdit->lockSampleData();

  // update the waveform display:
  int numSamples  = samplePlayerToEdit->getNumSamples();
  int numChannels = samplePlayerToEdit->getNumChannels();

  if( numSamples > 0 && numChannels > 0 )
  {

  // retrieve the filename (the full path) and cut it down to the filename only:
  //juce::String fileName = juce::String(bufferToEdit->getSampleName());
  //if( fileName.contains(T("\\")) )
  //  fileName = fileName.fromLastOccurrenceOf(T("\\"), false, false);
  //if( fileName.contains(T("/")) )

  // setup stuff:
  //sampleDisplay->setCaption(fileName);
  //setWaveform(samplePlayerToEdit->getSampleData(), numSamples, numChannels );
  juce::String audioFile = getApplicationDirectory() + File::separatorString 
  + juce::String(samplePlayerToEdit->getSampleName());
  bufferToUse->loadAudioDataFromFile(audioFile, true);
  setMaximumRange(0.0, numSamples, -1.1, 1.1);
  if( resetZoomFactor == true )
  setCurrentRange(0.0, numSamples, -1.1, 1.1);
  }
  // we may now release the sample data and repaint ourselves:
  samplePlayerToEdit->unlockSampleData();
  repaint();
  */


  /*
  // aquire mutex-locked access to the sample-data (such that no other thread corrupts them while 
  // we are reading them out):
  modulatorToEdit->sampleModulator.suspendAudioProcessing();
  modulatorToEdit->sampleModulator.lockSampleData();

  int      numSamples  = modulatorToEdit->sampleModulator.getNumSamples();
  double*  sampleData  = modulatorToEdit->sampleModulator.getSampleData();
  double** sampleData2 = &sampleData; // we actually need a pointer to a pointer

  WaveformDisplay::setWaveform(sampleData2, numSamples, 1);

  // obtain the name (relative path) of the sample-file
  juce::String fileName = juce::String(modulatorToEdit->sampleModulator.getSampleName());

  // truncate the path-info such that only the file-name will be left:
  fileName = fileName.fromLastOccurrenceOf(T("/"), false, false);
  fileName = fileName.fromLastOccurrenceOf(T("\\"), false, false);

  // ...and show the name of the file in the waveform-display:
  setCaption(fileName);

  // release the mutex-lock:
  modulatorToEdit->sampleModulator.unlockSampleData();
  modulatorToEdit->sampleModulator.resumeAudioProcessing();
  */
}

void SamplePlayerEditorDisplay::paint(juce::Graphics &g)
{
  WaveformDisplay::paint(g); // does not take scale and dc into account

  if( samplePlayerToEdit == NULL )
    return;

  // draw loop-locators:
  double factor    = 1.0 / samplePlayerToEdit->parameters->getRecordingSampleRate();
  double start     = factor * samplePlayerToEdit->parameters->getPlaybackStart();
  double end       = factor * samplePlayerToEdit->parameters->getPlaybackEnd();
  double loopStart = factor * samplePlayerToEdit->parameters->getLoopStart();  
  double loopEnd   = factor * samplePlayerToEdit->parameters->getLoopEnd();
  double xGrayL    = start;
  double xGrayR    = end;
  double y         = 0.0;
  drawCurrentPositionLocator(g,  (float) start);
  drawCurrentPositionLocator(g,  (float) end);
  if( samplePlayerToEdit->parameters->getLoopMode() != rosic::SamplePlaybackParameters::NO_LOOP )
  {
    drawLeftLocator( g,  (float) loopStart);
    drawRightLocator(g,  (float) loopEnd);
    xGrayL = jmin(start, loopStart);
  }


  // gray out the area before the start sample and after the end sample:
  g.setColour( Colours::lightgrey.withAlpha(0.5f) );
  transformToComponentsCoordinates(xGrayL, y);
  g.fillRect(0.f, 0.f, (float) xGrayL, (float) getHeight() );
  y = 0.0;
  transformToComponentsCoordinates(xGrayR, y);
  g.fillRect((float) xGrayR, 0.f, (float) getWidth(), (float) getHeight() );

  //g.fillAll(Colours::grey.withAlpha(0.5f));
}

int SamplePlayerEditorDisplay::whatIsUnderTheMouseCursor(const MouseEvent &e)
{
  if( samplePlayerToEdit == NULL )
    return 0;

  // get the position of the event in components coordinates
  mouseX = e.getMouseDownX();
  mouseY = e.getMouseDownY();

  // a margin for the breakpoint-dots
  double marginInPixels = dotRadius;

  // retrieve the scale and offset variables form the Modulator-object
  double scale  = 1.0;
  double offset = 0.0;

  double factor = 1.0 / samplePlayerToEdit->parameters->getRecordingSampleRate();
  double x1 = factor * samplePlayerToEdit->parameters->getPlaybackStart() ;
  double x2 = factor * samplePlayerToEdit->parameters->getPlaybackEnd() ;
  double x3 = factor * samplePlayerToEdit->parameters->getLoopStart()     ;
  double x4 = factor * samplePlayerToEdit->parameters->getLoopEnd()       ;
  double y1 = 0.0; // dummy
  double y2 = 0.0; // dummy
  double y3 = 0.0; // dummy
  double y4 = 0.0; // dummy
  transformToComponentsCoordinates(x1, y1);
  transformToComponentsCoordinates(x2, y2);
  transformToComponentsCoordinates(x3, y3);
  transformToComponentsCoordinates(x4, y4);


  if( abs(x3+2-mouseX) <= 4.0 && 
    samplePlayerToEdit->parameters->getLoopMode() != rosic::SamplePlaybackParameters::NO_LOOP )
  {
    return LOOP_START_LOCATOR;
  }
  else if( abs(x1+2-mouseX) <= 4.0  )
  {
    return START_LOCATOR;
  }
  else if( abs(x2+2-mouseX) <= 4.0  )
  {
    return END_LOCATOR;
  }
  else if( abs(x4+2-mouseX) <= 4.0  &&
    samplePlayerToEdit->parameters->getLoopMode() != rosic::SamplePlaybackParameters::NO_LOOP )
  {
    return LOOP_END_LOCATOR;
  }

  // no object has been identified to be under the mouse cursor:
  return NO_OBJECT;
}

void SamplePlayerEditorDisplay::mouseDown(const MouseEvent &e)
{
  if( samplePlayerToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  mouseX = e.getMouseDownX();
  mouseY = e.getMouseDownY();

  // check if one of the locators was grabbed:
  if( e.mods.isLeftButtonDown() )
  {
    locatorBeingDragged = whatIsUnderTheMouseCursor(e);

    // don't grab loop-locators when not in loop-mode:
    if( samplePlayerToEdit->parameters->getLoopMode() 
      == rosic::SamplePlaybackParameters::NO_LOOP )
    {
      if( locatorBeingDragged == LOOP_START_LOCATOR ||
        locatorBeingDragged == LOOP_END_LOCATOR )
      {
        locatorBeingDragged = NO_OBJECT;
      }
    }

  }

  // inform our listeners about the change:
  sendChangeMessage();
}

void SamplePlayerEditorDisplay::mouseDrag(const MouseEvent &e)
{
  if( samplePlayerToEdit == NULL )
    return;

  // get the position of the event in components coordinates:
  mouseX = e.getMouseDownX() + e.getDistanceFromDragStartX();
  mouseY = e.getMouseDownY() + e.getDistanceFromDragStartY();

  // get the position of the event in system coordinates (seconds for the x-axis):
  double x = (double) mouseX;
  double y = (double) mouseY;
  transformFromComponentsCoordinates(x, y);

  x *=  samplePlayerToEdit->parameters->getRecordingSampleRate();

  // get the sample-index for the ne loop start or end:
  //x *= playbackParametersToEdit->getRecordingSampleRate();
  //int sampleIndex = (int) floor(x);
  double sampleIndex = x;

  // drag the locator to a new position:
  if ( locatorBeingDragged == START_LOCATOR )
    samplePlayerToEdit->parameters->setPlaybackStart(sampleIndex);
  else if ( locatorBeingDragged == END_LOCATOR )
    samplePlayerToEdit->parameters->setPlaybackEnd(sampleIndex);
  else if ( locatorBeingDragged == LOOP_START_LOCATOR )
    samplePlayerToEdit->setLoopStart(sampleIndex);
  else if ( locatorBeingDragged == LOOP_END_LOCATOR )
    samplePlayerToEdit->setLoopEnd(sampleIndex);

  // trigger a repaint:
  repaint();

  // inform our listeners about the change:
  sendChangeMessage();
}

void SamplePlayerEditorDisplay::mouseMove(const MouseEvent &e)
{
  switch( whatIsUnderTheMouseCursor(e) )
  {
  case NO_OBJECT:       
    currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);          
    break;
  case START_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  case END_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  case LOOP_START_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  case LOOP_END_LOCATOR:    
    currentMouseCursor = MouseCursor(MouseCursor::LeftRightResizeCursor); 
    break;
  }
}

void SamplePlayerEditorDisplay::mouseUp(const MouseEvent &e)
{
  locatorBeingDragged = -1;
  repaint();
  sendChangeMessage();
}

XmlElement* SamplePlayerEditorDisplay::getStateAsXml(
  const juce::String &stateName) const
{
  return InteractiveCoordinateSystem::getStateAsXml(stateName);
}

bool SamplePlayerEditorDisplay::setStateFromXml(const XmlElement &xmlState)
{
  return InteractiveCoordinateSystem::setStateFromXml(xmlState);
}

//=================================================================================================