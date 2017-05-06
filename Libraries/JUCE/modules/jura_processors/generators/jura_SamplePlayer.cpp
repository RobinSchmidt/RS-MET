

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

