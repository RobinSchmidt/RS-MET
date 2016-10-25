AudioModuleDeletionWatcher::~AudioModuleDeletionWatcher()
{
  for(int i = 0; i < watchedAudioModules.size(); i++ )
    watchedAudioModules[i]->deRegisterDeletionWatcher(this);
}

void AudioModuleDeletionWatcher::addWatchedAudioModule(AudioModule *moduleToBeWatched)
{
  if( moduleToBeWatched != NULL )
  {
    moduleToBeWatched->registerDeletionWatcher(this);
    watchedAudioModules.addIfNotAlreadyThere(moduleToBeWatched);
  }
}

void AudioModuleDeletionWatcher::removeWatchedAudioModule(AudioModule 
  *moduleNotToBeWatchedAnymore)
{
  if( moduleNotToBeWatchedAnymore != NULL )
  {
    moduleNotToBeWatchedAnymore->deRegisterDeletionWatcher(this);
    watchedAudioModules.removeFirstMatchingValue(moduleNotToBeWatchedAnymore);
  }
}

//=================================================================================================

// construction/destruction:

AudioModule::AudioModule()
{
  plugInLock = new CriticalSection;

  ParameterObserver::localAutomationSwitch = true;  // activate automation for this instance
  wantsTempoSyncInfo = true;
  //underlyingRosicInstrument = NULL;  // by default, this does not wrap an instrument
  moduleName = juce::String("AudioModule");
  //versionjuce::String = juce::String(T("1.0"));
  patchFormatIndex = 1;
  triggerInterval = 0.0;
  saveAndRecallState = true;
  initializeAutomatableParameters();
}

AudioModule::~AudioModule()
{
  plugInLock->enter();

  for(int i = 0; i < deletionWatchers.size(); i++)
    deletionWatchers[i]->audioModuleWillBeDeleted(this);

  childModules.getLock().enter();
  AudioModule* childModule = NULL;
  while( childModules.size() > 0 )
  {
    childModule = childModules[childModules.size()-1];
    removeChildAudioModule(childModule, true);
    //childModule->parent = NULL;
    //delete childModules[childModules.size()-1];
    //childModules.removeLast();
  }
  childModules.getLock().exit();

  plugInLock->exit();
  delete plugInLock;
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);

  childModules.getLock().enter();
  for(int i=0; i<childModules.size(); i++)
    childModules[i]->setSampleRate(newSampleRate);
  childModules.getLock().exit();
}

void AudioModule::setBeatsPerMinute(double newBpm)
{
  ScopedLock scopedLock(*plugInLock);

  childModules.getLock().enter();
  for(int i=0; i<childModules.size(); i++)
    childModules[i]->setBeatsPerMinute(newBpm);
  childModules.getLock().exit();
}

void AudioModule::setModuleName(const juce::String &newName)
{
  moduleName = newName;
}

void AudioModule::setModuleNameAppendix(const juce::String& newAppendix)
{
  moduleNameAppendix = newAppendix;
}

void AudioModule::addChildAudioModule(AudioModule* moduleToAdd)
{ 
  ScopedLock scopedLock(*plugInLock);
  childModules.getLock().enter();
  childModules.addIfNotAlreadyThere(moduleToAdd);
  childModules.getLock().exit();
  addChildStateManager(moduleToAdd);
}

void AudioModule::removeChildAudioModule(AudioModule* moduleToRemove, bool deleteObject)
{ 
  ScopedLock scopedLock(*plugInLock);
  childModules.getLock().enter();
  int index = childModules.indexOf(moduleToRemove);
  jassert( index != -1 ); // trying to remove a module which is not a child of this one?
  if( index != -1 )
  { 
    childModules.remove(index);
    removeChildStateManager(moduleToRemove);
    if( deleteObject == true )
      delete moduleToRemove;
  }
  childModules.getLock().exit();
}

bool AudioModule::checkForCrack()
{
  return false;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AudioModule::getIndexAmongNameSakes(AudioModule *child)
{
  ScopedLock scopedLock(*plugInLock);

  juce::String name = child->getModuleName();
  int index = -1;

  childModules.getLock().enter();
  for(int c = 0; c < childModules.size(); c++)
  {
    if( childModules[c]->getModuleName() == name )
    {
      index++;
      if( childModules[c] == child )
        break;
    }
  }
  childModules.getLock().exit();

  return index;
}

juce::String AudioModule::getModuleHeadlineString()
{
  return moduleName + moduleNameAppendix;
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void AudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  markStateAsDirty();
}

XmlElement* AudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  if( !wantsSaveAndRecallState() )
    return NULL;

  // the XmlElement which stores all the relevant state-information:
  XmlElement* xmlState = new XmlElement(moduleName); 

  // store a patch format version (useful when patch-formats change later):
  xmlState->setAttribute("PatchFormat", patchFormatIndex);

  // store controller mappings (if any)
  automatableModuleStateToXml(this, xmlState);

  //// if this module is a polyphonic instrument, we store some global instrument parameters 
  //// (such as number of voices, tuning, etc.):
  //if( underlyingRosicInstrument != NULL )
  //  polyphonicInstrumentStateToXml(underlyingRosicInstrument, xmlState);

  // save the states of all childModules in child-XmlElements:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
  {
    if( childModules[c]->wantsSaveAndRecallState() )
    {
      XmlElement* childState = childModules[c]->getStateAsXml(
        childModules[c]->getStateName(), false);
      xmlState->addChildElement(childState);
    }
  }
  childModules.getLock().exit();
  setStateName(stateName, markAsClean);
  return xmlState;
}

void AudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
  bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  XmlElement convertedState = convertXmlStateIfNecessary(xmlState);
  automatableModuleStateFromXml(this, convertedState);

  //// check, if this module wraps an instrument - if so, we wave some more settings to restore:
  //if( underlyingRosicInstrument != NULL )
  //  polyphonicInstrumentStateFromXml(underlyingRosicInstrument, convertedState);

  juce::String thisName =  this->moduleName;  // for debug

  // if we have child-modules, we try to restore their states by looking for corresponding
  // child XmlElements in the xmlState:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
  {
    childModules[c]->setStateToDefaults();
    int indexAmongNameSakes = getIndexAmongNameSakes(childModules[c]);
    XmlElement* childState = getChildElementByNameAndIndexAmongNameSakes(
      convertedState, childModules[c]->moduleName, indexAmongNameSakes);
    if( childState != NULL )
    {
      //childModules[c]->setStateFromXml(*childState, stateName, markAsClean);
      childModules[c]->setStateFromXml(*childState, juce::String::empty, markAsClean);
    }
  }
  childModules.getLock().exit();

  // this call we need for setting our parent dirty when some embedded sub-module loads a state 
  // - when the call was due to the outer parent loading its state, the subsequent call to 
  // setStateName may set it back to 'clean'  
  if( parent != NULL )
    parent->markStateAsDirty();

  setStateName(stateName, markAsClean);
}

XmlElement AudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*plugInLock);

  int xmlPatchFormatIndex = xmlState.getIntAttribute("PatchFormat", 1);
  jassert( xmlPatchFormatIndex == this->patchFormatIndex ); 
  // You should override the state conversion in your subclass - the idea is that when the stored
  // patch format number in the xml differs from our patchFormatIndex member variable here, a 
  // conversion may be necessary. That way, new versions of the plugin that may use a different
  // patch format can still load patches that were stored by older versions.

  return xmlState;
}
    
void AudioModule::registerDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*plugInLock);
  deletionWatchers.addIfNotAlreadyThere(watcher);
}
           
void AudioModule::deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*plugInLock);
  deletionWatchers.removeFirstMatchingValue(watcher);
}

//-------------------------------------------------------------------------------------------------
// event processing:

void AudioModule::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*plugInLock);
  if( message.isController() )
  {
    int controllerNumber = message.getControllerNumber();
    int controllerValue  = message.getControllerValue();
    setMidiController(controllerNumber, controllerValue);
  }
  else if( message.isNoteOn() )
    noteOn(message.getNoteNumber(), message.getVelocity());
  else if( message.isNoteOff() )
    noteOff(message.getNoteNumber());
  else if( message.isAllNotesOff() )
    allNotesOff();
  else if( message.isPitchWheel() )
    setPitchBend(message.getPitchWheelValue());
}

void AudioModule::noteOn(int noteNumber, int velocity)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->noteOn(noteNumber, velocity);
}

void AudioModule::noteOff(int noteNumber)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->noteOff(noteNumber);
}

void AudioModule::allNotesOff()
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->allNotesOff();
}

void AudioModule::setMidiController(int controllerNumber, int controllerValue)
{
  ScopedLock scopedLock(*plugInLock);
  AutomatableModule::setMidiController(controllerNumber, controllerValue);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->setMidiController(controllerNumber, controllerValue);

  // distribute the controller message to all children (i.e. embedded sub-modules):
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
    childModules[c]->setMidiController(controllerNumber, controllerValue);
  childModules.getLock().exit();
}

void AudioModule::setPitchBend(int pitchBendValue)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //{
  //  double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
  //  underlyingRosicInstrument->setPitchBend(wheelValueMapped);
  //}
}

//-------------------------------------------------------------------------------------------------
// others:

void AudioModule::initializeAutomatableParameters()
{

}

void AudioModule::updateCoreObjectAccordingToParameters()
{
  ScopedLock scopedLock(*plugInLock);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}

//-------------------------------------------------------------------------------------------------
// overrides for juce:AudioProcessor baseclass:

void AudioModule::prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) 
{
  if(getProcessingPrecision() == singlePrecision)
  {
    // The host is going to call the single precision version of the processBlock callback but we
    // need the double precision version to be called internally. So we allocate an internal buffer
    // that is used for back and forth conversion.
    internalAudioBuffer.setSize(numChannels, maximumExpectedSamplesPerBlock,
      false,  // keepExistingContent
      false,  // clearExtraSpace
      true);  // avoidReallocating
  }
  else
  {
    // The host will call the double precision version of the processBlock callback, so we don't 
    // need an internal conversion buffer. We don't want to waste memory, so we request the
    // internal buffer to allocate a zero-sized memory block (i hope, this works - ToDo: check in 
    // the debugger):
    internalAudioBuffer.setSize(0, 0, false, false, false);
  }

  // Maybe we could release the buffer when the host calls releaseResources() - we'll see.

  setSampleRate(sampleRate);  // Subclasses may want to do something in their overriden version
}

void AudioModule::getStateInformation(juce::MemoryBlock& destData)
{
  XmlElement* xmlState = getStateAsXml("StateAsRequestedByHost", false);
  copyXmlToBinary(*xmlState, destData);
  delete xmlState;
}

void AudioModule::setStateInformation(const void* data, int sizeInBytes)
{
  XmlElement* const xmlState = getXmlFromBinary(data, sizeInBytes);
  //ParameterObserver::globalAutomationSwitch = false; // why this - threading problems? -> interferes with total recall in quadrifex
  ParameterObserver::guiAutomationSwitch = false;
  setStateFromXml(*xmlState, "recalled by host", false);
  ParameterObserver::guiAutomationSwitch = true;
  //ParameterObserver::globalAutomationSwitch = true;
  delete xmlState;

  // some hosts (Tracktion) seem to keep and re-use an open GUI when a new project is loaded, so
  // we must make sure, that this re-used GUI is updated according to the new recalled state - we 
  // do this by broadcasting a changeMessage which will be picked up by AudioPlugInEditor
  //sendChangeMessage();
  // (this comment is old and refers to the old way of doing it ...but i want to verify this in 
  // Tracktion someday)
}

// maybe move this helper function to jura_framework/tools:
template<class SourceType, class TargetType>
void convertAudioBuffer(const AudioBuffer<SourceType>& source, AudioBuffer<TargetType>& target)
{
  int numChannels = jmin(source.getNumChannels(), target.getNumChannels());
  int numSamples  = jmin(source.getNumSamples(),  target.getNumSamples());
  const SourceType *sourcePointer;    
        TargetType *targetPointer;
  for(int channel = 0; channel < numChannels; channel++)
  {
    sourcePointer = source.getReadPointer(channel);
    targetPointer = target.getWritePointer(channel);
    for(int sample = 0; sample < numSamples; sample++)
      targetPointer[sample] = (TargetType)sourcePointer[sample];
  }
}

void AudioModule::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
  ScopedLock scopedLock(*plugInLock);

  convertAudioBuffer(buffer, internalAudioBuffer);  // float -> double
  processBlock(internalAudioBuffer, midiMessages);  // process doubles
  convertAudioBuffer(internalAudioBuffer, buffer);  // double -> float
}

void AudioModule::processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages)
{
  ScopedLock scopedLock(*plugInLock);

  int numChannels = buffer.getNumChannels();
  int numSamples  = buffer.getNumSamples();
  double **inOutBuffer = buffer.getArrayOfWritePointers();

  processBlock(inOutBuffer, numChannels, numSamples);
}

bool AudioModule::setPreferredBusArrangement(bool isInput, int bus, 
  const AudioChannelSet& preferredSet)
{
  // I don't really know, if it's necessarry to override this, but i did hit a breakpoint related 
  // to this one day and overriding this made it go away. But it might have been something else to 
  // blame (some things have changed since then). We need some testing/debugging with and without 
  // this override...

  // Reject any bus arrangements that are not compatible with your plugin
  const int numChannels = preferredSet.size();
  if (numChannels != 1 && numChannels != 2)
    return false;
  return AudioProcessor::setPreferredBusArrangement(isInput, bus, preferredSet);
}
