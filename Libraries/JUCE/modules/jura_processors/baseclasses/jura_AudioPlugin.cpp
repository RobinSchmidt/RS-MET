AudioPlugin::AudioPlugin(AudioModule *moduleToWrap)
{
  // this function is from rosic/basics/GlobalFunctions.h, it's currently not available:
  //checkForMemoryLeaksOnExit();
  //juce::String *leakString = new juce::String(T("String that should cause a memory leak")); // to test if leak-detection works

  ScopedLock sl(plugInLock);

  // experimental:
  initialiseJuce_GUI();  // ???

  underlyingAudioModule = moduleToWrap;
  //underlyingAudioModule = nullptr;

  // maybe, here, we could somehow set up the parameter that the plugin exposes to the host and 
  // connect them to the module's internal parameters
}

AudioPlugin::~AudioPlugin()
{    
  ScopedLock sl(plugInLock);
  if( underlyingAudioModule != nullptr )  // not necessary? it's actually safe to delete a nullptr
  {
    delete underlyingAudioModule;
    underlyingAudioModule = nullptr;
  }
}

// mandatory overrides for juce::AudioProcessor:

void AudioPlugin::prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) 
{
  ScopedLock sl(plugInLock);

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

  if( underlyingAudioModule != nullptr )
    underlyingAudioModule->setSampleRate(sampleRate);
}

template<class SourceType, class TargetType>
void convertAudioBuffer(const AudioBuffer<SourceType>& source, AudioBuffer<TargetType>& target)
{
  // maybe move this helper function somewhere else, for example to jura_framework/tools

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
void AudioPlugin::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
  //ScopedLock scopedLock(plugInLock);

  convertAudioBuffer(buffer, internalAudioBuffer);  // float -> double
  processBlock(internalAudioBuffer, midiMessages);  // process doubles
  convertAudioBuffer(internalAudioBuffer, buffer);  // double -> float
}

//bool AudioPlugin::acceptsMidi() const 
//{ 
//  if(dynamic_cast<AudioModuleWithMidiIn*>(underlyingAudioModule) != nullptr)
//    return true;
//  return false;
//  // \todo at some point (when we want to write MIDI processors, i.e. plugins that modify or 
//  // produce MIDI data, we should handle the producesMidi function in a similar way
//} 

AudioProcessorEditor* AudioPlugin::createEditor() 
{ 
  if(underlyingAudioModule == nullptr)
    return nullptr;
  ParameterObserver::guiAutomationSwitch = false;   // don't automate widgets during creation
  AudioModuleEditor *moduleEditor = underlyingAudioModule->createEditor();
  AudioPluginEditor *pluginEditor = new AudioPluginEditor(moduleEditor, this);
  ParameterObserver::guiAutomationSwitch = true;    // now, widgets can be automated again
  return pluginEditor;
}

void AudioPlugin::getStateInformation(juce::MemoryBlock& destData)
{
  if(underlyingAudioModule != nullptr)
  {
    XmlElement* xmlState = underlyingAudioModule->getStateAsXml("StateAsRequestedByHost", false);
    copyXmlToBinary(*xmlState, destData);
    delete xmlState;
  }
}

void AudioPlugin::setStateInformation(const void* data, int sizeInBytes)
{
  if(underlyingAudioModule != nullptr)
  {
    XmlElement* const xmlState = getXmlFromBinary(data, sizeInBytes);
    //ParameterObserver::globalAutomationSwitch = false; // why this - threading problems? -> interferes with total recall in quadrifex
    ParameterObserver::guiAutomationSwitch = false;
    underlyingAudioModule->setStateFromXml(*xmlState, "recalled by host", false);
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
}

// optional overrides for juce::AudioProcessor:

void AudioPlugin::processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages)
{
  ScopedLock scopedLock(plugInLock);

  if(underlyingAudioModule != nullptr)
  {
    int numChannels = buffer.getNumChannels();
    int numSamples  = buffer.getNumSamples();
    double **inOutBuffer = buffer.getArrayOfWritePointers();

    // later, we have to add the MIDI processing here...(and call 
    // underlyingAudioModule->handleMidiMessage accordingly) ...but maybe we should do that in
    // a subclass AudioPluginWithMidiIn

    underlyingAudioModule->processBlock(inOutBuffer, numChannels, numSamples);
  }
  else
    buffer.clear();
}

bool AudioPlugin::setPreferredBusArrangement(bool isInput, int bus, 
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

//=================================================================================================

AudioPluginWithMidiIn::AudioPluginWithMidiIn(AudioModuleWithMidiIn *moduleToWrap)
  : AudioPlugin(moduleToWrap)
{
  wrappedModuleWithMidiIn = moduleToWrap;
}

void AudioPluginWithMidiIn::processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages)
{
  //AudioPlugin::processBlock(buffer, midiMessages); 
  //// preliminary - later we need to handle the midiMessages here....

  ScopedLock sl(plugInLock);

#ifdef DEBUG
  static int callCount = 1;
  if( underlyingAudioModule == NULL )
  {
    buffer.clear();
    jassertfalse;
    return;
  }
  _clearfp();                                        // reset
  unsigned int cw = _controlfp(0, 0);                // get state
  //cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
  cw &= ~(EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
  unsigned int cw0 = _controlfp(cw, _MCW_EM);
#endif

  // request time-info from host and update the bpm-values for the modulators accordingly, if they 
  // are in sync mode (maybe this should be moved up into the baseclass AudioModule):
  int timeToNextTriggerInSamples = -1;
  if( underlyingAudioModule->wantsTempoSyncInfo )
  {
    AudioPlayHead::CurrentPositionInfo info;
    if( getPlayHead() != 0  &&  getPlayHead()->getCurrentPosition(info) )
    {
      if( info.bpm <= 0.0 )
        info.bpm = 120.0;  // fallback value when nothing meaningful is passed
      underlyingAudioModule->setBeatsPerMinute(info.bpm); 
    }

    if( underlyingAudioModule->getTriggerInterval() != 0.0 )
    {
      double timeInBeats = secondsToBeats(info.timeInSeconds, info.bpm); 
      // kludge - will probably not work when tempo changes - use ppqPosition here later....

      double timeToNextTriggerInBeats = underlyingAudioModule->getTriggerInterval()                                 
        - fmod(timeInBeats, underlyingAudioModule->getTriggerInterval());

      timeToNextTriggerInSamples = 
        roundToInt(getSampleRate()*beatsToSeconds(timeToNextTriggerInBeats, info.bpm));

      if( timeToNextTriggerInSamples >= buffer.getNumSamples() )
        timeToNextTriggerInSamples = -1; // indicates that we don't need to trigger in this block
    }
  }

  if(midiMessages.isEmpty() && timeToNextTriggerInSamples == -1)
  {
    //underlyingAudioModule->processBlock(buffer, 0, buffer.getNumSamples()); // old
    AudioPlugin::processBlock(buffer, midiMessages); 
    // no messages to process - baseclass method can handle this
  }
  else
  {
    // some stuff for the input midi-buffer
    MidiBuffer::Iterator midiBufferIterator(midiMessages);
    MidiMessage          currentMidiMessage(0x80,0,0,0.0);
    int                  currentMidiMessageOffset = -1;
    bool                 aMidiMessageWasRetrieved = false;
    int start     = 0;
    int subLength = buffer.getNumSamples();
    while( start < buffer.getNumSamples() )
    {
      // check if there are midi-messages at this instant of time:
      int i = start;
      midiBufferIterator.setNextSamplePosition(i);
      currentMidiMessageOffset = -1; // may be set to another value by function on next line
      aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
        currentMidiMessage, currentMidiMessageOffset);
      while( aMidiMessageWasRetrieved && currentMidiMessageOffset == i ) 
      {
        handleMidiMessage(currentMidiMessage); // respond to the midi-message
        aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(currentMidiMessage, 
          currentMidiMessageOffset);           // get the next message at this sample
      }

      // check if we should spawn a trigger-event at this instant:
      if( start == timeToNextTriggerInSamples )
      {
        underlyingAudioModule->trigger();

        // update the time to the next trigger:
        timeToNextTriggerInSamples = -1;  
        // preliminary - assumes that the next trigger is not inside this buffer (which is assured
        // only when buffersize is much smaller than the retrigger interval - which is probably
        // true in realtime situations but maybe not in rendering)

        int dummy = 0;
      }

      // all messages at this time instant have been processed - now determine the time-offset to 
      // the next event and compute the number of samples that can be processed without 
      // interrupting events (subLength):    
      if( aMidiMessageWasRetrieved ) 
      {
        // ignore clock messages:
        while( aMidiMessageWasRetrieved && currentMidiMessage.isMidiClock() )
        {
          aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
            currentMidiMessage, currentMidiMessageOffset);      
        }
        if( aMidiMessageWasRetrieved )
          subLength = currentMidiMessageOffset - start; 
        else
          subLength = buffer.getNumSamples() - start;
      }
      else if( timeToNextTriggerInSamples != -1 )
        subLength = timeToNextTriggerInSamples - start; 
      else
        subLength = buffer.getNumSamples() - start;

      // process the block of sample-frames before the occurrence of the next event:  
      //underlyingAudioModule->processBlock(buffer, start, subLength); // old
      double* channelPointers[2];   // we assume 2 channels
      channelPointers[0] = buffer.getWritePointer(0) + start;
      channelPointers[1] = buffer.getWritePointer(1) + start;
      underlyingAudioModule->processBlock(channelPointers, 2, subLength);

      start += subLength; // increment start position for next iteration

      //midiMessages.clear(0, start+subLength);
    }
  }

#ifdef DEBUG
  callCount++;
  _controlfp(cw0, _MCW_EM);
#endif
}

void AudioPluginWithMidiIn::handleMidiMessage(MidiMessage message)
{
  ScopedLock sl(plugInLock);
  if( message.isMidiClock() )
    return;
  if( wrappedModuleWithMidiIn != NULL )
    wrappedModuleWithMidiIn->handleMidiMessage(message);
}
