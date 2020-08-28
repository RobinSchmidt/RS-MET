// Bugs: 
// -SamplePlayerEditorDisplay::paint is called regularly (!) WTF?!
//  ...maybe because i wanted to draw animated position locators - however - do it with 
//  paintOverChildren like in the MultiComp


// ToDo:
// -SamplePlayerEditorDisplay::setAudioFileToUse has beed updated - we now also need to update the
//  sample data in the AudioModule and rosic dsp-core object...
// -check, how the recall of the sample works for the WaveOsc and do it here the same way 
// -maybe factor out a common baseclass to consolidate the code for saving and recalling the 
//  sample-path (SampleBasedAudioModule or something) 
//  ...but there's a proble - the WaveOsc also causes the MipMap to be re-rendered
//  ...maybe we need in rosic a baseclass for the oscillator and SamplePlayer, too that has 
//  functions for setting the sample
// -maybe also have an editor baseclass that contains the sample-load widget set that can be used
//  by SamplePlayer and WaveOsc
// -update parameter creation to new style
// -maybe re-arrange the gui to have some elements above the display, some below and some to the 
//  right
// -let the user switch between waveform and spectrogram view and/or use the spectrogram for the 
//  background and the waveform for the foreground
//  -have an rsHeatMapDisplay or rsColorMapDisplay baseclass and rsSpectrogramPlot subclass
// -maybe change the background of the sample-file text field

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SamplePlayerAudioModule::SamplePlayerAudioModule(
  CriticalSection *newPlugInLock, rosic::SamplePlayer *newSamplePlayerToWrap)
: AudioModule(newPlugInLock)
{
  //jassert( newSamplePlayerToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  if(newSamplePlayerToWrap == nullptr) {
    wrappedSamplePlayer = new rosic::SamplePlayer;
    wrappedSamplePlayerOwned = true;
  }
  else
    wrappedSamplePlayer = newSamplePlayerToWrap;
  setModuleTypeName("SamplePlayer");
  initializeAutomatableParameters();
}

SamplePlayerAudioModule::~SamplePlayerAudioModule()
{
  if(wrappedSamplePlayerOwned)
    delete wrappedSamplePlayer;
}

//-------------------------------------------------------------------------------------------------

AudioModuleEditor* SamplePlayerAudioModule::createEditor(int type)
{
  return new SamplePlayerModuleEditor(lock, this);
}

/*
// get rid of this function - thsi was the odl way of doing it:
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
*/

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
  //bool success = true;

  // set up the audio-file name:
  //juce::String samplePath = xmlState.getStringAttribute(T("AudioFileRelativePath"), juce::String());
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

  auto playbackParametersState = xmlState.getChildByName("PlaybackParameters");
  if( playbackParametersState == NULL )
    return false;

  // retrieve the audio-file name:
  juce::String samplePath = playbackParametersState->getStringAttribute(
    "AudioFileRelativePath", juce::String());

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
  double rootKey = RAPT::rsFreqToPitch(freq);
  getParameterByName("RootKey")->setValue(rootKey, true, true);
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
    getParameterByName("PlaybackStart")->setRange(0.0, numSamples);
    getParameterByName("LoopStart")->setRange(0.0, numSamples);
    getParameterByName("LoopLength")->setRange(0.0, numSamples);

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

  p = new AutomatableParameter(lock, "Mute", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Solo", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Level", -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "LevelByKey", -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "LevelByVel", -12.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "MidSide",    0.0,  1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Pan",       -1.0,  1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "Tune", -36.0, 36.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "TuneByKey", -200.0, 200.0, 0.1, 100.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "TuneByVel", -200.0, 200.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "RootKey", 0.0, 127.0, 0.01, 64.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "PlaybackStart", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Loop", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "LoopStart", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "LoopLength", 0.0, 1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "FilterActive", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Lowpass", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "Highpass", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  p = new AutomatableParameter(lock, "PhaseRandomize", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p = new AutomatableParameter(lock, "PhaseRandomizationSeed", 0.0, 1000.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  //for(int i=0; i < (int) parameters.size(); i++ )
  //  parameterChanged(parameters[i]);
  // i think, it's obsolete, but it also causes an access violation - that should actually not 
  // happen
}

//=================================================================================================

SamplePlayerEditorDisplay::SamplePlayerEditorDisplay(AudioFileBuffer *newBufferToUse)
  : AudioFileBufferUser(newBufferToUse)
  //: WaveformDisplay(newBufferToUse) //, InteractiveCoordinateSystem(juce::String(T("SampleDisplay")))
{
  //setValueFieldPopup(false); // old

  //bufferToEdit             = NULL;
  //playbackParametersToEdit = NULL;
  samplePlayerToEdit = NULL;
  locatorBeingDragged      = -1;

  setAxisLabelX(juce::String());
  setAxisLabelY(juce::String());
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
return rsWaveformPlot::getMaximumMeaningfulRange(relativeMarginLeft, relativeMarginRight,
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
  //int numChannels = samplePlayerToEdit->getNumChannels();
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


  juce::String audioFile = getApplicationDirectory() + juce::File::getSeparatorString()
    + juce::String(samplePlayerToEdit->getSampleName());
  setAudioFileToUse(juce::File(audioFile));

  //setWaveform(samplePlayerToEdit->getSampleData(), numSamples, numChannels );
}

bool SamplePlayerEditorDisplay::setAudioFileToUse(const juce::File &newFileToUse)
{
  lockUsedBufferPointer();
  AudioFileBuffer* buf = bufferToUse; // shorthand
  if( buf != nullptr ) {

    buf->acquireReadLock();  
    // was formerly below load...but it ssems to make more sense to acquire it *before* loading
    // the data
    buf->loadAudioDataFromFile(newFileToUse, false);



    rsWaveformPlot::setWaveform(buf->getSampleData(), buf->getNumSamples(), 
      buf->getNumChannels());
    rsWaveformPlot::setSampleRate(buf->getFileSampleRate());
    // we should probably pass the sampleRate to setWaveform and setWaveform should itself update
    // the current and maximum range

    //setMaximumRangeX(0.0, buf->getNumSamples()+1);
    //setCurrentRangeX(0.0, buf->getNumSamples()+1);
    //setMaximumRangeX(0.0, 5000); // test - probably doesn't work bcs time format is wrong
    //setCurrentRangeX(0.0, 5000); // ..we probably should set the time format to samples




    // pass the data to the underlying rosic-object (todo: just keep a pointer to 
    // SamplePlayerAudioModule and pass the data to *that* which in turn passes it further on to
    // the rosic object):
    if( samplePlayerToEdit != nullptr )  {

      // todo: retrieve the pointers to the channels and pass it to the rosic-object, this might
      // also be the place to render samples from modal models etc.

      //bufferToUse->getSampleData(0, 0);

      // todo:
      //samplePlayerToEdit->setSampleData(
      //  buf->getSampleData(), buf->getNumSamples(), buf->getNumChannels());

      samplePlayerToEdit->parameters->setRecordingSampleRate(buf->getFileSampleRate());
      // RecordingSampleRate means the sample-rate at which the data was recorded
    }
    // hmm...well...actually, it should not be the display's responsibility to update the 
    // sample-data in the AudioModule (and hence, the rosic object) - instead, the editor should be
    // reposnible for this...


    buf->releaseReadLock();
  }
  unlockUsedBufferPointer();


  return buf->isAudioFileValid();
  //return true; // preliminary
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
  rsWaveformPlot::paint(g); // does not take scale and dc into account

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
  //drawCurrentPositionLocator(g,  (float) start); // from InteractiveCoordinateSystem
  //drawCurrentPositionLocator(g,  (float) end);   // from InteractiveCoordinateSystem
  if( samplePlayerToEdit->parameters->getLoopMode() != rosic::SamplePlaybackParameters::NO_LOOP )
  {
    //drawLeftLocator( g,  (float) loopStart);   // from InteractiveCoordinateSystem
    //drawRightLocator(g,  (float) loopEnd);     // from InteractiveCoordinateSystem
    xGrayL = jmin(start, loopStart);
  }


  /*
  // gray out the area before the start sample and after the end sample:
  g.setColour( Colours::lightgrey.withAlpha(0.5f) );
  toPixelCoordinates(xGrayL, y);
  g.fillRect(0.f, 0.f, (float) xGrayL, (float) getHeight() );
  y = 0.0;
  toPixelCoordinates(xGrayR, y);
  g.fillRect((float) xGrayR, 0.f, (float) getWidth(), (float) getHeight() );
  */

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
  //double marginInPixels = dotRadius;

  // retrieve the scale and offset variables form the Modulator-object
  //double scale  = 1.0;
  //double offset = 0.0;

  double factor = 1.0 / samplePlayerToEdit->parameters->getRecordingSampleRate();
  double x1 = factor * samplePlayerToEdit->parameters->getPlaybackStart() ;
  double x2 = factor * samplePlayerToEdit->parameters->getPlaybackEnd() ;
  double x3 = factor * samplePlayerToEdit->parameters->getLoopStart()     ;
  double x4 = factor * samplePlayerToEdit->parameters->getLoopEnd()       ;
  double y1 = 0.0; // dummy
  double y2 = 0.0; // dummy
  double y3 = 0.0; // dummy
  double y4 = 0.0; // dummy
  toPixelCoordinates(x1, y1);
  toPixelCoordinates(x2, y2);
  toPixelCoordinates(x3, y3);
  toPixelCoordinates(x4, y4);


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
  fromPixelCoordinates(x, y);

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
  return nullptr;
  //return InteractiveCoordinateSystem::getStateAsXml(stateName);
}

bool SamplePlayerEditorDisplay::setStateFromXml(const XmlElement &xmlState)
{
  return false;
  //return InteractiveCoordinateSystem::setStateFromXml(xmlState);
}

//=================================================================================================


SamplePlayerEditorContextMenu::SamplePlayerEditorContextMenu(
  SamplePlayerAudioModule* newSamplePlayerModuleToEdit, Component* componentToAttachTo)
  //: ComponentMovementWatcher(componentToAttachTo)
{
  // init the pointer to the samplePlayeror to be edited:
  //jassert( newSamplePlayerModuleToEdit != NULL )
  samplePlayerModuleToEdit = newSamplePlayerModuleToEdit;

  addWidget( ampHeadline = new RTextField( "Amplitude:") );
  ampHeadline->setDescription(juce::String("Manipulations of the amplitude"));

  addWidget( levelSlider = new RSlider("LevelSlider") );
  levelSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Level") );
  levelSlider->setDescription(juce::String("Output level of the samplePlayer"));
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( levelByKeySlider = new RSlider("LevelByKeySlider") );
  levelByKeySlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("LevelByKey") );
  levelByKeySlider->setSliderName(juce::String("K"));
  levelByKeySlider->setDescription(juce::String("Key dependence of samplePlayer's output level"));
  levelByKeySlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( levelByVelSlider = new RSlider("LevelVelSlider") );
  levelByVelSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("LevelByVel") );
  levelByVelSlider->setSliderName(juce::String("V"));
  levelByVelSlider->setDescription(juce::String("Velocity dependence of samplePlayer's output level"));
  levelByVelSlider->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( midSideSlider = new RSlider("MidSideSlider") );
  midSideSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("MidSide") );
  midSideSlider->setSliderName(juce::String("Mid/Side"));
  midSideSlider->setDescription("Mid/side adjustment for stereo samples");
  midSideSlider->setStringConversionFunction(&ratioToString0);

  addWidget( panSlider = new RSlider("PanSlider") );
  panSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Pan") );
  panSlider->setDescription("Panorama position of the samplePlayer");
  panSlider->setStringConversionFunction(&valueToString3);

  addWidget( tuningHeadline = new RTextField( "Tuning:") );
  tuningHeadline->setDescription(juce::String("Manipulations of the tuning/detuning of the samplePlayer"));

  addWidget( tuneSlider = new TuningSlider("TuneSlider") );
  tuneSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Tune") );
  tuneSlider->setDescription(juce::String("Tuning of the samplePlayer in semitones"));
  tuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( tuneByKeySlider = new RSlider("TuneByKeySlider") );
  tuneByKeySlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("TuneByKey") );
  tuneByKeySlider->setSliderName(juce::String("K"));
  tuneByKeySlider->setDescription(juce::String("Key dependence of the tuning"));
  tuneByKeySlider->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( tuneByVelSlider = new RSlider("TuneVelSlider") );
  tuneByVelSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("TuneByVel") );
  tuneByVelSlider->setSliderName("V");
  tuneByVelSlider->setDescription("Velocity dependence of the tuning");
  tuneByVelSlider->setStringConversionFunction(&percentToStringWithUnit2);
  tuneByVelSlider->setVisible(false);

  addWidget( rootKeySlider = new TuningSlider("RootKeySlider") );
  rootKeySlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("RootKey") );
  rootKeySlider->setDescription("Rootkey of the sample");
  rootKeySlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  // todo: add a second slider with the fundamental frequency in Hz and a button for auto detecting it (auto)

  addWidget( timeHeadline = new RTextField( "Time:") );
  timeHeadline->setDescription(juce::String("Time related parameters"));

  addWidget( startSlider = new TuningSlider("StartSlider") ); // TuningSlider?! surely a bug
  //startSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Start") );
  //startSlider->setDescription(juce::String("Start point of playback within the sample"));
  startSlider->setDescription(juce::String("Not yet implemented"));
  startSlider->setStringConversionFunction(&valueToString0);
  //startSlider->addListener(this);

  addWidget( startByVelSlider = new RSlider("StartByVelSlider") );
  //startByVelSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("StartByVel") );
  startByVelSlider->setSliderName("V");
  //startByVelSlider->setDescription(juce::String("Velocity dependence of the start time"));
  startByVelSlider->setDescription("Not yet implemented");
  startByVelSlider->setStringConversionFunction(&valueToString0);

  addWidget( loopButton = new RButton("Loop") );
  loopButton->assignParameter( samplePlayerModuleToEdit->getParameterByName("Loop") );
  loopButton->setDescription("Switch loop on/off");
  //loopButton->setClickingTogglesState(true);
  //loopButton->addRButtonListener(this);

  addWidget( loopStartSlider = new TuningSlider("LoopStartSlider") );
  loopStartSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("LoopStart") );
  loopStartSlider->setDescription("Start point of the loop within the sample");
  loopStartSlider->setStringConversionFunction(&valueToString0);
  //loopStartSlider->addListener(this);

  addWidget( loopLengthSlider = new TuningSlider("LoopLengthSlider") );
  loopLengthSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("LoopLength") );
  loopLengthSlider->setDescription("Length of the loop");
  loopLengthSlider->setStringConversionFunction(&valueToString0);
  //loopLengthSlider->addListener(this);

  addWidget( filterHeadline = new RTextField( "Filter:") );
  filterHeadline->setDescription("Basic filtering of the sample");

  addWidget( lowpassSlider = new TuningSlider("LowpassSlider") );
  lowpassSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Lowpass") );
  lowpassSlider->setDescription("Lowpass cutoff frequency");
  lowpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( highpassSlider = new TuningSlider("HighpassSlider") );
  highpassSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("Highpass") );
  highpassSlider->setDescription("Highpass cutoff frequency");
  highpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( miscHeadline = new RTextField( "Misc:") );
  miscHeadline->setDescription("Miscellaneous sample manipulations");

  addWidget( phaseRandomizeButton = new RButton("PhaseRandomize") );
  phaseRandomizeButton->assignParameter( samplePlayerModuleToEdit->getParameterByName("PhaseRandomize") );
  phaseRandomizeButton->setDescription("Switch phase-randomization on/off");
  //phaseRandomizeButton->addRButtonListener(this);

  addWidget( phaseSeedSlider = new TuningSlider("PhaseRandomizeSeedSlider") );
  //phaseSeedSlider->assignParameter( samplePlayerModuleToEdit->getParameterByName("PhaseSeed") );
  //phaseSeedSlider->setDescription(juce::String("Seed for the random phases"));
  phaseSeedSlider->setDescription("Not yet implemented");
  phaseSeedSlider->setStringConversionFunction(&valueToString0);
  //phaseSeedSlider->addListener(this);


  addWidget( closeButton = new RButton(RButton::CLOSE) );
  closeButton->setDescription("Closes the samplePlayer context menu");
  closeButton->setClickingTogglesState(false);
  // we don't listen to this button ourselves - this is the job of the outlying editor object

  // factor out widget creation int function  createWidgets();


  updateWidgetsAccordingToState();
  setSize(180, 404);
}

SamplePlayerEditorContextMenu::~SamplePlayerEditorContextMenu()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

/*
void SamplePlayerEditorContextMenu::rButtonClicked(RButton *buttonThatWasClicked)
{
if( samplePlayerModuleToEdit == NULL )
return;
if( samplePlayerModuleToEdit->wrappedSamplePlayer == NULL )
return;

rosic::SamplePlayer* sp = samplePlayerModuleToEdit->wrappedSamplePlayer;
if( buttonThatWasClicked == loopButton )
sp->setLoopMode( loopButton->getToggleState() );

sendChangeMessage();
}

void SamplePlayerEditorContextMenu::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
if( samplePlayerModuleToEdit == NULL )
return;
if( samplePlayerModuleToEdit->wrappedSamplePlayer == NULL )
return;
rosic::SamplePlayer* o = samplePlayerModuleToEdit->wrappedSamplePlayer;

sendChangeMessage();
}
*/

/*
void SamplePlayerEditorContextMenu::updateWidgetsAccordingToState()
{
if( samplePlayerModuleToEdit == NULL )
return;
if( samplePlayerModuleToEdit->wrappedSamplePlayer == NULL )
return;

rosic::SamplePlayer* sp = samplePlayerModuleToEdit->wrappedSamplePlayer;
rosic::SamplePlaybackParameters* params = sp->parameters;

// update the widgets:
levelSlider->setValue(                  params->getLevel(),       false);
levelByKeySlider->setValue(             params->getLevelByKey(),  false);
levelByVelSlider->setValue(             params->getLevelByVel(),  false);
//midSideSlider->setValue(                params->getMidSide(),     false);
//panSlider->setValue(                    params->getPan(),         false);

tuneSlider->setValue(      params->getTune(),          false);
tuneSlider->setValue(      params->getTuneByKey(),     false);
tuneSlider->setValue(      params->getTuneByVel(),     false);
rootKeySlider->setValue(   params->getRootKey(),       false);

startSlider->setValue(     params->getPlaybackStart(), false);
loopButton->setToggleState(params->getLoopMode()!=0,   false);
loopStartSlider->setValue( params->getLoopStart(),     false);
loopLengthSlider->setValue(params->getLoopLength(),    false);

// filter and randomize....
}
*/

void SamplePlayerEditorContextMenu::resized()
{
  Component::resized();
  int x  = 0;
  int y  = 0;
  int w  = getWidth();
  int w2 = w/2;
  //int h  = getHeight();

  closeButton->setBounds(w-16, 0, 16, 16);

  int sh   = 16;    // slider height
  int inc  = sh+4;
  int inc2 = inc+8;

  ampHeadline->setBounds(x+4, y+4, w-8-16, sh);
  y += inc;
  levelSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  levelByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  levelByVelSlider->setBounds(w2+4, y+4, w2-8, sh);
  y += inc;
  midSideSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  panSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc2;

  tuningHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  tuneSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  tuneByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  tuneByVelSlider->setBounds(w2+4, y+4, w2-8, sh);
  y += inc;
  rootKeySlider->setBounds(x+4, y+4, w-8, sh);
  y += inc2;

  timeHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  startSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  //startByKeySlider->setBounds(x+4,  y+4, w2-8, sh);
  startByVelSlider->setBounds(w2+4, y+4, w2-8, sh);
  //y += inc;
  y += 8;
  loopButton->setBounds(x+4, y+4, 40, sh);
  y += sh-2;
  loopStartSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  loopLengthSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc2;

  filterHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  lowpassSlider->setBounds(x+4, y+4, w-8, sh);
  y += sh-2;
  highpassSlider->setBounds(x+4, y+4, w-8, sh);
  y += inc2;

  miscHeadline->setBounds(x+4, y+4, w-8, sh);
  y += inc;
  phaseRandomizeButton->setBounds(x+4, y+4, 120, sh);
  y += sh-2;
  phaseSeedSlider->setBounds(x+4, y+4, w-8, sh);
}



//=================================================================================================
// class SamplePlayerModuleEditor:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SamplePlayerModuleEditor::SamplePlayerModuleEditor(CriticalSection *newPlugInLock,
  SamplePlayerAudioModule* newSamplePlayerModuleToEdit)
  : SampleBasedAudioModuleEditor(newSamplePlayerModuleToEdit)
{
  jassert( newSamplePlayerModuleToEdit != nullptr );
  samplePlayerModuleToEdit = newSamplePlayerModuleToEdit;
  //samplePlayerToEdit = newSamplePlayerToEdit->wrappedSamplePlayer; // obsolete?

  setHeadlineText("SamplePlayer");


  addPlot( sampleDisplay = new SamplePlayerEditorDisplay() );
  //addAndMakeVisible( sampleDisplay = new SamplePlayerEditorDisplay() );
  sampleDisplay->setSamplePlayerToEdit(samplePlayerModuleToEdit->wrappedSamplePlayer);
  sampleDisplay->setDescription("Shows the sample waveform data");
  sampleDisplay->addChangeListener(this); // to receive modifications of loop-settings, etc.
  sampleDisplay->setVerticalCoarseGrid(1.0,   false);
  sampleDisplay->setVerticalFineGrid(0.1,     false);
  sampleDisplay->setAxisValuesPositionX(rsPlotSettings::INVISIBLE);
  sampleDisplay->setHorizontalCoarseGrid(1.0, true);  // zero-line and lines at +-1
  sampleDisplay->setHorizontalFineGrid(0.1,   false);
  //sampleDisplay->setSampleRate(1.0);

  // create the zoomer for the sampleDisplay and associate it with the sampleDisplay:
  addChildColourSchemeComponent( sampleDisplayZoomer = new rsPlotZoomer() );
  sampleDisplayZoomer->setRelativeMargins(5.0, 5.0, 10.0, 10.0);
  sampleDisplayZoomer->setCoordinateSystem(sampleDisplay);
  sampleDisplayZoomer->setVerticalMouseWheelMode(
    CoordinateSystemZoomer::horizontalZoomViaVerticalMouseWheel);



  contextMenu = new SamplePlayerEditorContextMenu(samplePlayerModuleToEdit, this);
  contextMenu->addChangeListener(this);
  contextMenu->setAlwaysOnTop(true);
  contextMenu->setOpaque(true);
  contextMenu->closeButton->addRButtonListener(this);
  addChildColourSchemeComponent(contextMenu, false, false);

  addWidget( fileLabel = new RTextField( "File:") );
  fileLabel->setDescription("The currently loaded sample file.");
  fileLabel->setJustification(Justification::centredLeft);

  addWidget( formatLabel = new RTextField( "Format:") );
  formatLabel->setDescription("Data format of currently loaded sample file.");

  addWidget( formatInfoLabel = new RTextField( juce::String()) );
  formatInfoLabel->setDescription(formatLabel->getDescription());

  addWidget( levelSlider = new RSlider("Level") );
  levelSlider->assignParameter( moduleToEdit->getParameterByName("Level") );
  levelSlider->setDescription("Playback level");
  levelSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( tuneSlider = new RSlider("Tune") );
  tuneSlider->assignParameter( moduleToEdit->getParameterByName("Tune") );
  tuneSlider->setDescription("Tuning");
  tuneSlider->setStringConversionFunction(&semitonesToStringWithUnit2);

  addWidget( lowpassSlider = new RSlider("Lowpass") );
  lowpassSlider->assignParameter( moduleToEdit->getParameterByName("Lowpass") );
  lowpassSlider->setSliderName("LPF");
  lowpassSlider->setDescription("Cutoff of the lowpass filter");
  lowpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( highpassSlider = new RSlider("Highpass") );
  highpassSlider->assignParameter( moduleToEdit->getParameterByName("Highpass") );
  highpassSlider->setSliderName("HPF");
  highpassSlider->setDescription("Cutoff of the highpass filter");
  highpassSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( muteButton = new RButton(juce::String("Mute")) );
  muteButton->assignParameter( moduleToEdit->getParameterByName("Mute") );
  muteButton->setDescription("Mute sample player");

  addWidget( soloButton = new RButton("Solo") );
  soloButton->assignParameter( moduleToEdit->getParameterByName("Solo") );
  soloButton->setDescription("Switch sample player to solo mode");

  addWidget( phaseRandomizeButton = new RButton(juce::String("PhaseRandomize")) );
  phaseRandomizeButton->assignParameter( moduleToEdit->getParameterByName("PhaseRandomize") );
  phaseRandomizeButton->setButtonText( "PhRnd" );
  phaseRandomizeButton->setDescription("Turn phase randomization on/off");

  addWidget( loopButton = new RButton("Loop") );
  loopButton->assignParameter( moduleToEdit->getParameterByName("Loop") );
  loopButton->setDescription(juce::String("Switch loop on/off"));
  //loopButton->setClickingTogglesState(true);
  //loopButton->addRButtonListener(this);

  addWidget( moreButton = new RButton(juce::String("More")) );
  moreButton->setDescription("Show context menu with more settings");
  moreButton->setClickingTogglesState(true);
  moreButton->addRButtonListener(this);

  addWidget( fromLoopButton = new RButton("From Loop") );
  fromLoopButton->addRButtonListener(this);
  fromLoopButton->setDescription(
    "Infer root key and detuning from loop-length and number of cycles in loop.");
  fromLoopButton->setClickingTogglesState(true);

  addWidget( loopSnapButton = new RButton("Snap To 0") );
  loopSnapButton->addRButtonListener(this);
  loopSnapButton->setDescription("Snap loop-start and -end to upward zero crossings.");
  loopSnapButton->setClickingTogglesState(true);

  //addWidget( loopEndSnapButton = new RButton(juce::String(T("Snap"))) );
  //loopEndSnapButton->addRButtonListener(this);
  //loopEndSnapButton->setDescription(juce::String(T("Snap loop end to next zero crossing.")));
  //loopEndSnapButton->setClickingTogglesState(false);

  addWidget( loopLengthLockButton = new RButton("Lock Length") );
  loopLengthLockButton->addRButtonListener(this);
  loopLengthLockButton->setDescription("Lock loop length.");
  loopLengthLockButton->setClickingTogglesState(true);

  addWidget( autoNumCyclesButton = new RButton("Auto") );
  autoNumCyclesButton->addRButtonListener(this);
  autoNumCyclesButton->setDescription("Auto-detect the number of pitch cycles in the loop.");
  autoNumCyclesButton->setClickingTogglesState(true);

  addWidget( rootKeySlider = new RSlider("Root") );
  rootKeySlider->addListener(this);
  rootKeySlider->setSliderName("Root:");
  rootKeySlider->setDescription("Root key of the sample.");
  rootKeySlider->setStringConversionFunction(&midiNoteToString);
  rootKeySlider->setRange(0.0, 127.0, 1.0, 64.0);

  addWidget( rootDetuneSlider = new RSlider("Detune") );
  rootDetuneSlider->addListener(this);
  rootDetuneSlider->setSliderName("Detune:");
  rootDetuneSlider->setDescription("Detuning of the fundamental frequency from root key in cents.");
  rootDetuneSlider->setStringConversionFunction(&centsToStringWithUnit2);
  rootDetuneSlider->setRange(-50.0, 50.0, 0.01, 0.0);

  addWidget( startSlider = new RSlider("StartSlider") );
  startSlider->addListener(this);
  startSlider->setSliderName("Play Start:");
  startSlider->setDescription("Playback start point in the sample.");
  startSlider->setStringConversionFunction(&valueToString3);
  startSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( startByVelSlider = new RSlider("StartByVelSlider"));
  startByVelSlider->addListener(this);
  startByVelSlider->setSliderName("V:");
  startByVelSlider->setDescription("Velocity dependence of playback start point.");
  startByVelSlider->setStringConversionFunction(&valueToString3);
  startByVelSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( endSlider = new RSlider("EndSlider") );
  endSlider->addListener(this);
  endSlider->setSliderName("Play End:");
  endSlider->setDescription("Playback end point in the sample.");
  endSlider->setStringConversionFunction(&valueToString3);
  endSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( loopLengthSlider = new RSlider("LoopLengthSlider") );
  loopLengthSlider->addListener(this);
  loopLengthSlider->setSliderName("Length");
  loopLengthSlider->setDescription("Loop length (in samples).");
  loopLengthSlider->setStringConversionFunction(&valueToString3);
  loopLengthSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( loopStartSlider = new RSlider("LoopStartSlider") );
  loopStartSlider->addListener(this);
  loopStartSlider->setSliderName("Loop Start:");
  loopStartSlider->setDescription("Loop start point (in samples).");
  loopStartSlider->setStringConversionFunction(&valueToString3);
  loopStartSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( loopEndSlider = new RSlider("LoopEndSlider") );
  loopEndSlider->addListener(this);
  loopEndSlider->setSliderName("Loop End:");
  loopEndSlider->setDescription("Loop end point (in samples).");
  loopEndSlider->setStringConversionFunction(&valueToString3);
  loopEndSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( loopCrossfadeTimeSlider = new RSlider("LoopCrossfadeTimeSlider") );
  loopCrossfadeTimeSlider->addListener(this);
  loopCrossfadeTimeSlider->setSliderName("X-Fade");
  loopCrossfadeTimeSlider->setDescription("Loop crossfade time (in samples).");
  loopCrossfadeTimeSlider->setStringConversionFunction(&valueToString3);
  loopCrossfadeTimeSlider->setRange(0.0, 1.0, 1.0, 0.0);

  addWidget( loopCrossfadeShapeSlider = new RSlider("LoopCrossfadeShapeSlider") );
  loopCrossfadeShapeSlider->addListener(this);
  loopCrossfadeShapeSlider->setSliderName("Shape");
  loopCrossfadeShapeSlider->setDescription(
    juce::String("Loop crossfade shape (1: constant sum, 2: constant power)."));
  loopCrossfadeShapeSlider->setStringConversionFunction(&valueToString2);
  loopCrossfadeShapeSlider->setRange(1.0, 2.0, 0.01, 0.0);

  addWidget( loopNumCyclesSlider = new RSlider("LoopNumCyclesSlider") );
  loopNumCyclesSlider->addListener(this);
  loopNumCyclesSlider->setSliderName("Cycles:");
  loopNumCyclesSlider->setDescription("Number of pitch cycles in the loop.");
  loopNumCyclesSlider->setStringConversionFunction(&valueToString0);
  loopNumCyclesSlider->setRange(1.0, 100.0, 1.0, 1.0);

  // customize the descriptions for the load/save buttons:
  stateWidgetSet->stateLoadButton->setDescription("Load sample playback settings from file");
  stateWidgetSet->stateSaveButton->setDescription("Save sample playback settings to file");
  stateWidgetSet->statePlusButton->setDescription("Skip to next playback settings file in current directory");
  stateWidgetSet->stateMinusButton->setDescription("Skip to previous playback settings file in current directory");


  // preliminary, for development
  AudioFileManager::setActiveDirectory(getSupportDirectory() + "/Samples/SingleCycle/DampedSinusoids");


  //setSamplePlayerToEdit(newSamplePlayerToEdit->wrappedSamplePlayer);
  updateWidgetsAccordingToState(true);

  setSize(640, 300);
}

SamplePlayerModuleEditor::~SamplePlayerModuleEditor()
{
  delete contextMenu;
}

//-------------------------------------------------------------------------------------------------
// setup:
/*
void SamplePlayerModuleEditor::setSamplePlayerToEdit(rosic::SamplePlayer* newSamplePlayerToEdit)
{
samplePlayerToEdit = newSamplePlayerToEdit;
sampleDisplay->setSamplePlayerToEdit(newSamplePlayerToEdit);
//sampleDisplay->setSampleBufferToEdit(samplePlayerToEdit->theBuffer);
//sampleDisplay->setSamplePlaybackParametersToEdit(samplePlayerToEdit->parameters);
updateWidgetsAccordingToState(true);
}
*/

/*
factored out to AudioFileManager
bool SamplePlayerModuleEditor::setSampleFromFile(const File &fileToLoadFrom)
{
if( samplePlayerModuleToEdit == NULL )
return false;

bool success = samplePlayerModuleToEdit->setSampleFromFile(fileToLoadFrom);
if( success == true )
{
sampleFileNameLabel->setText(fileToLoadFrom.getFileName(), false);
formatInfoLabel->setText(createAudioFileInfoString(fileToLoadFrom), false);
sampleDisplay->setAudioFileToUse(fileToLoadFrom);
return true;
}
else
{
sampleDisplay->setCaption(juce::String(T("Error:")) + fileToLoadFrom.getFileName());
return false;
}
}
*/

//-------------------------------------------------------------------------------------------------
// inquiry:

//-------------------------------------------------------------------------------------------------
// callbacks:

void SamplePlayerModuleEditor::rButtonClicked(RButton *b)
{
  if( samplePlayerModuleToEdit == NULL )
    return;

  else if( b == loopButton )
  {
    /*
    if( loopButton->getToggleState() == true )
    samplePlayerToEdit->parameters->setLoopMode(SamplePlaybackParameters::FORWARD_LOOP);
    else
    samplePlayerToEdit->parameters->setLoopMode(SamplePlaybackParameters::NO_LOOP);
    fromLoopButton->setEnabled(       loopButton->getToggleState() );
    loopSnapButton->setEnabled(       loopButton->getToggleState() );
    loopLengthLockButton->setEnabled( loopButton->getToggleState() );
    loopNumCyclesSlider->setEnabled(  loopButton->getToggleState() );
    loopStartSlider->setEnabled(      loopButton->getToggleState() );
    loopEndSlider->setEnabled(        loopButton->getToggleState() );
    */
  }
  else if( b == fromLoopButton )
  {
    if( fromLoopButton->getToggleState() == true )
      samplePlayerModuleToEdit->setRootKeyFromLoop();
    //setRootKeyAndRootDetuneFromLoop();
  }
  else if( b == loopSnapButton )
  {
    //samplePlayerToEdit->setSnapLoopToZeros( loopSnapButton->getToggleState() );
  }
  else if( b == moreButton )
  {
    if( moreButton->getToggleState() == true )
    {
      //int x = getScreenX() + getWidth();
      int x = moreButton->getScreenX() + moreButton->getWidth();
      int y = getScreenY() - 102; // this is a kludgy thing here with the context menu positioning
      contextMenu->setTopLeftPosition(x, y);
      contextMenu->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
      contextMenu->setVisible(true);
      contextMenu->toFront(true);
    }
    else
      contextMenu->setVisible(false);
    return;
  }
  else if( b == contextMenu->closeButton )
  {
    moreButton->setToggleState(false, true);
    return;
  }
  else
    SampleBasedAudioModuleEditor::rButtonClicked(b);

  //else
  //{
  //  // handle click on the load/save/plus/minus buttons:
  //  Editor::rButtonClicked(buttonThatWasClicked);
  //
  //  // it must have been one of the inherited load/save/plus/minus buttons, so we must update the
  //  // preset field:
  //  updateWidgetsAccordingToState(true);
  //  //updatePresetField();
  //}


  updateWidgetsAccordingToState(true);
  sendChangeMessage();
}

void SamplePlayerModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  if( samplePlayerModuleToEdit == NULL )
    return;

  //startSlider->setValue(samplePlayerToEdit->parameters->getPlaybackStart(),   false);
  //loopStartSlider->setValue(samplePlayerToEdit->parameters->getLoopStart(),   false);
  //loopEndSlider->setValue(samplePlayerToEdit->parameters->getLoopEnd(),       false);
  //loopLengthSlider->setValue(samplePlayerToEdit->parameters->getLoopLength(), false);

  if( objectThatHasChanged == contextMenu )
  {
    moduleToEdit->markStateAsDirty();
    updateWidgetsAccordingToState();
    sendChangeMessage();
  }
  else if( objectThatHasChanged == sampleDisplay )
  {
    if( fromLoopButton->getToggleState() == true )
      samplePlayerModuleToEdit->setRootKeyFromLoop();
    updateWidgetsAccordingToState(false);
    sendChangeMessage();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void SamplePlayerModuleEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  /*
  if( samplePlayerToEdit == NULL )
  return;

  // for sliders, we re-retrieve their value directly after setting up the audio engine because the
  // audio engine may refuse to acept the new value due to consistency constraints:
  if( sliderThatHasChanged == startSlider )
  {
  samplePlayerToEdit->parameters->setPlaybackStart(startSlider->getValue());
  startSlider->setValue(samplePlayerToEdit->parameters->getPlaybackStart(), false);
  sampleDisplay->updatePlot();
  }
  else if( sliderThatHasChanged == endSlider )
  {
  samplePlayerToEdit->parameters->setPlaybackEnd(endSlider->getValue());
  endSlider->setValue(samplePlayerToEdit->parameters->getPlaybackEnd(), false);
  sampleDisplay->updatePlot();
  }
  else if( sliderThatHasChanged == loopStartSlider )
  {
  samplePlayerToEdit->setLoopStart(loopStartSlider->getValue());
  loopStartSlider->setValue(samplePlayerToEdit->parameters->getLoopStart(),   false);
  loopLengthSlider->setValue(samplePlayerToEdit->parameters->getLoopLength(), false);
  sampleDisplay->updatePlot();
  if( fromLoopButton->getToggleState() == true )
  setRootKeyAndRootDetuneFromLoop();
  }
  else if( sliderThatHasChanged == loopEndSlider )
  {
  samplePlayerToEdit->setLoopEnd(loopEndSlider->getValue());
  loopEndSlider->setValue(samplePlayerToEdit->parameters->getLoopEnd(), false);
  loopLengthSlider->setValue(samplePlayerToEdit->parameters->getLoopLength(), false);
  sampleDisplay->updatePlot();
  if( fromLoopButton->getToggleState() == true )
  setRootKeyAndRootDetuneFromLoop();
  }
  else if( sliderThatHasChanged == loopNumCyclesSlider )
  {
  samplePlayerToEdit->parameters->setNumPitchCyclesInLoop(loopNumCyclesSlider->getValue());
  loopNumCyclesSlider->setValue(samplePlayerToEdit->parameters->getNumPitchCyclesInLoop(),
  false);
  if( fromLoopButton->getToggleState() == true )
  setRootKeyAndRootDetuneFromLoop();
  }
  else if( sliderThatHasChanged == loopLengthSlider )
  {
  samplePlayerToEdit->setLoopLength(loopLengthSlider->getValue());
  loopEndSlider->setValue(samplePlayerToEdit->parameters->getLoopEnd(), false);
  loopLengthSlider->setValue(samplePlayerToEdit->parameters->getLoopLength(), false);
  sampleDisplay->updatePlot();
  if( fromLoopButton->getToggleState() == true )
  setRootKeyAndRootDetuneFromLoop();
  }
  else if( sliderThatHasChanged == rootKeySlider )
  {
  samplePlayerToEdit->setRootKey((int) rootKeySlider->getValue());
  fromLoopButton->setToggleState(false, false);
  }

  // although certain sliders would not require an update, most will, so we do the update in any
  // case:
  //sampleDisplay->updatePlot();
  //updateWidgetsAccordingToState(true);
  //setPresetDirty();
  sendChangeMessage();

  */
}

void SamplePlayerModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);

  // draw rectangles for the parameter-groups:

  fillRectWithBilinearGradient(g, fileRectangle, editorColourScheme.topLeft, 
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  g.drawRect(fileRectangle);

  fillRectWithBilinearGradient(g, loopRectangle, editorColourScheme.topLeft, 
    editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  g.drawRect(loopRectangle);
}

void SamplePlayerModuleEditor::resized()
{

  ScopedLock scopedLock(*lock);

  AudioModuleEditor::resized();
  int m = 4; // margin
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  //y = getPresetSectionBottom();
  //y =  getHeadlineBottom();

  y  = getPresetSectionBottom() + m;

  sampleDisplay->setBounds(x, y, w-16, h-y-68-16);
  sampleDisplayZoomer->alignWidgetsToCoordinateSystem();

  y = sampleDisplayZoomer->getBottom();
  y = sampleDisplay->getBottom()+16;  // preliminary until zoomer works


  //sampleDisplay->setVisible(false);
  //sampleDisplayZoomer->setVisible(false); // test

  w = getWidth() / 2;

  fileRectangle.setBounds(0, y, w, getHeight()-y);
  loopRectangle.setBounds(w, y, w, getHeight()-y);

  x = fileRectangle.getX();
  y = fileRectangle.getY();
  w = fileRectangle.getWidth();
  h = fileRectangle.getHeight();

  fileLabel->setBounds(x+4, y+4, 32, 16);
  x = w-84;
  sampleLoadButton->setBounds(x+4, y+4, 40, 16);
  sampleMinusButton->setBounds(sampleLoadButton->getRight()+4,  y+4, 16, 16);
  samplePlusButton->setBounds( sampleMinusButton->getRight()-2, y+4, 16, 16);
  x = fileLabel->getRight();
  w = sampleLoadButton->getX()-8-x;
  sampleFileLabel->setBounds(x+4, y+4, w, 16);

  x = fileRectangle.getX();
  y = sampleFileLabel->getBottom();
  w = fileRectangle.getWidth();

  formatLabel->setBounds(x+4, y+4, 56, 16);
  formatInfoLabel->setBounds(formatLabel->getRight(), y+4, w-formatLabel->getWidth()-8, 16);

  y = formatInfoLabel->getBottom();

  rootKeySlider->setBounds(x+4, y+4, 80, 16);
  x = rootKeySlider->getRight();
  rootDetuneSlider->setBounds(x+4, y+4, 120, 16);
  x = rootDetuneSlider->getRight();
  fromLoopButton->setBounds(x+4, y+4, 72, 16);

  x = loopRectangle.getX();
  y = loopRectangle.getY();
  w = loopRectangle.getWidth();
  h = loopRectangle.getHeight();

  startSlider->setBounds(x+4,   y+4, w/2-8, 16);
  endSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);

  y += 20;
  loopButton->setBounds(x+4, y+4, 40, 16);
  x = loopButton->getRight();
  loopSnapButton->setBounds(x+4, y+4, 80, 16);
  x = loopSnapButton->getRight();
  loopLengthLockButton->setBounds(x+4, y+4, 80, 16);
  x = loopLengthLockButton->getRight();
  loopNumCyclesSlider->setBounds(x+4, y+4, loopRectangle.getRight()-x-8, 16);

  x =  loopRectangle.getX();
  y += 20;
  loopStartSlider->setBounds(x+4,   y+4, w/2-8, 16);
  loopEndSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);

  /*
  loopButton->setBounds(x+4, y+4, 40, 16);
  loopLengthLockButton->setBounds(x+10*w/16-40-4, y+4, 40, 16);
  x = loopButton->getRight();
  w = loopLengthLockButton->getX() - x;
  loopLengthSlider->setBounds(x+4, y+4, w-4, 16);

  x = loopRectangle.getX();
  w = loopRectangle.getWidth();
  y = loopLengthSlider->getBottom();
  loopSnapButton->setBounds(x+10*w/16-40-4, y+4, 40, 16);
  y += 20;
  loopEndSnapButton->setBounds(x+10*w/16-40-4, y+4, 40, 16);
  y -= 20;
  w = loopEndSnapButton->getX() - x;

  loopStartSlider->setBounds(x+4, y+4, w-4, 16);
  y += 20;
  loopEndSlider->setBounds(x+4, y+4, w-4, 16);

  x = loopLengthLockButton->getRight();
  y = loopRectangle.getY();
  w = loopRectangle.getRight()-x;
  loopCrossfadeTimeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  loopCrossfadeShapeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  loopNumCyclesSlider->setBounds(x+4, y+4, w-36-8, 16);
  autoNumCyclesButton->setBounds(loopNumCyclesSlider->getRight(), y+4, 36, 16);
  */

  /*
  x = 0;
  w = getWidth() / 2;

  fundamentalFreqSlider->setBounds(x+4, y+4, w-4-64, 20);
  fromLoopButton->setBounds(fundamentalFreqSlider->getRight()+4, y+4, 64, 20);
  */
}

bool SamplePlayerModuleEditor::setAudioData(AudioSampleBuffer* newBuffer,
  const juce::File& underlyingFile, bool markAsClean)
{
  if( samplePlayerModuleToEdit == NULL )
    return false;
  bool success = samplePlayerModuleToEdit->setSampleFromFile(underlyingFile);
  if( success == true )
  {
    sampleFileLabel->setText(underlyingFile.getFileName());
    formatInfoLabel->setText(createAudioFileInfoString(underlyingFile));

    sampleDisplay->setAudioFileToUse(underlyingFile);
    // hmm - i think, this causes the sampleDisplay to load the same file again - maybe it would be 
    // better if we just pass it the "newBuffer" which already contains the audio data - but then
    // we my need to be more careful about locking - so maybe after all, it's not such a bad idea
    // to let the display load the data independently...hmmm

    return true;
  }
  else
  {
    sampleDisplay->setCaption(juce::String("Error:") + underlyingFile.getFileName());
    return false;
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SamplePlayerModuleEditor::updateWidgetsAccordingToState(bool updateSampleDisplayAlso)
{
  /*
  if( samplePlayerToEdit == NULL )
  return;

  // update the buttons:
  loopButton->setToggleState( samplePlayerToEdit->parameters->getLoopMode() != 0, false);
  loopSnapButton->setToggleState(samplePlayerToEdit->getLoopSnapToZeroMode(),     false);

  // some button should only be visible when loop is on:
  fromLoopButton->setVisible(       loopButton->getToggleState() );
  loopSnapButton->setVisible(       loopButton->getToggleState() );
  loopLengthLockButton->setVisible( loopButton->getToggleState() );
  loopStartSlider->setVisible(      loopButton->getToggleState() );
  loopEndSlider->setVisible(        loopButton->getToggleState() );
  loopNumCyclesSlider->setVisible(  loopButton->getToggleState() );

  // update the sliders:
  rootKeySlider->setValue(      samplePlayerToEdit->parameters->getRootKey(),              false);
  //rootDetuneSlider->setValue(   samplePlayerToEdit->parameters->getRootDetune(),           false);
  startSlider->setValue(        samplePlayerToEdit->parameters->getPlaybackStart(),        false);
  startByVelSlider->setValue(   samplePlayerToEdit->parameters->getPlaybackStartByVel(),   false);
  endSlider->setValue(          samplePlayerToEdit->parameters->getPlaybackEnd(),          false);
  loopStartSlider->setValue(    samplePlayerToEdit->parameters->getLoopStart(),            false);
  loopEndSlider->setValue(      samplePlayerToEdit->parameters->getLoopEnd(),              false);
  loopLengthSlider->setValue(   samplePlayerToEdit->parameters->getLoopLength(),           false);
  loopNumCyclesSlider->setValue(samplePlayerToEdit->parameters->getNumPitchCyclesInLoop(), false);
  */

  AudioModuleEditor::updateWidgetsAccordingToState();
  contextMenu->updateWidgetsAccordingToState();

  // update the waveform display plot if desired:
  if( updateSampleDisplayAlso )
  {
    sampleDisplay->updatePlot(true);
    sampleDisplayZoomer->zoomToAllXY();
  }
}

void SamplePlayerModuleEditor::updateWidgetsAccordingToState()
{
  updateWidgetsAccordingToState(true);
}

//------------------------------------------------------------------------------------------------
// internal methods:

/*
bool SamplePlayerModuleEditor::openSampleLoadingDialog()
{
AudioFileManager::openLoadingDialog();

// AudioFileManager::currentFileFullPath has now be assigned - so we load the file now:
return setSampleFromFile(AudioFileManager::getActiveFile());
}
*/


