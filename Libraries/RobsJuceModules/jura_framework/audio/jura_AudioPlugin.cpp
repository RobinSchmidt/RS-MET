
void AudioPluginParameter::setValue(float newValue)
{
  if((double)newValue != metaValue)                   // avoid superfluous updates (some DAWs 
    MetaParameter::setMetaValue((double)newValue);    // continuously send constant values)
}

void AudioPluginParameter::parameterChanged(Parameter* p)
{
  beginChangeGesture();
  setValueNotifyingHost((float) p->getNormalizedValue());
  endChangeGesture();
}

//void AudioPluginParameter::setName(const String& newName)
//{
//  name = newName;
//}

//=================================================================================================

AudioPlugin::AudioPlugin(int numParameters)
  : AudioProcessor(BusesProperties().withInput( "Input",  AudioChannelSet::stereo())
                                    .withOutput("Output", AudioChannelSet::stereo()))
{
  ScopedLock sl(plugInLock);

  //setProcessingPrecision(doublePrecision); // nope - this is supposed to be called by th host
  initialiseJuce_GUI();  // why do we need this?
  //configureInsAndOuts();
  smoothingManager.setMutexLock(&plugInLock);
  createHostAutomatableParameters(numParameters);
  metaParaManager.registerObserver(this);
}

AudioPlugin::~AudioPlugin()
{
  ScopedLock sl(plugInLock);
  //metaParaManager.deRegisterObserver(this); // metaParaManager is destroyed anyway, so we may not need to deregister
  if( wrappedAudioModule != nullptr )  // not necessary? it's actually safe to delete a nullptr
  {
    delete wrappedAudioModule;
    wrappedAudioModule = nullptr;
  }
  //delete editorBoundsConstrainer;
}

//void AudioPlugin::configureInsAndOuts()
//{
//  AudioProcessor::BusesLayout layout;
//  setBusesLayout(layout);
//}

void AudioPlugin::setAudioModuleToWrap(AudioModule* moduleToWrap)
{
  jassert(moduleToWrap != nullptr); // you must pass a valid object here
  wrappedAudioModule = moduleToWrap;
  wrappedAudioModule->setSmoothingManager(&smoothingManager);
  wrappedAudioModule->setMetaParameterManager(&metaParaManager);
}

void AudioPlugin::autoAttachMetaParameters()
{
  int metaIndex = 0; // meta-parameter index
  for(int paraIndex = 0; paraIndex < wrappedAudioModule->getNumParameters(); paraIndex++) {
    if(metaIndex >= metaParaManager.getNumMetaParameters())
      break; // end of meta-array reached - we are done
    Parameter* p = wrappedAudioModule->getParameterByIndex(paraIndex);
    MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*>(p);
    if(mcp != nullptr) {
      metaParaManager.setMetaValue(metaIndex, mcp->getNormalizedValue());
      metaParaManager.setMetaName( metaIndex, mcp->getName());
      mcp->attachToMetaParameter(metaIndex);
      metaIndex++;
    }
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

  smoothingManager.setSampleRate(sampleRate);
  if( wrappedAudioModule != nullptr )
    wrappedAudioModule->setSampleRate(sampleRate);
}

template<class SourceType, class TargetType>
void convertAudioBuffer(const AudioBuffer<SourceType>& source, AudioBuffer<TargetType>& target)
{
  // maybe move this helper function somewhere else, for example to jura_framework/tools

  // old:
  //int numChannels = jmin(source.getNumChannels(), target.getNumChannels());
  //int numSamples  = jmin(source.getNumSamples(),  target.getNumSamples());

  // new (especially important for FL Studio with its variable buffer sizes):
  int numChannels = source.getNumChannels();
  int numSamples  = source.getNumSamples();
  target.setSize(numChannels, numSamples, true, false, true);

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
  if(wrappedAudioModule == nullptr)
    return nullptr;
  ParameterObserver::setGuiAutomationSwitch(false);   // don't automate widgets during creation
  AudioModuleEditor *moduleEditor = wrappedAudioModule->createEditor();
  AudioPluginEditor *pluginEditor = new AudioPluginEditor(moduleEditor, this);
  ParameterObserver::setGuiAutomationSwitch(true);    // now, widgets can be automated again
  //pluginEditor->setConstrainer(editorBoundsConstrainer);
  return pluginEditor;
}

//void AudioPlugin::setParameter(int index, float value)
//{
//  //parameters[index] = value;
//  wrappedAudioModule->setMidiController(index, 127.f * value); // preliminary
//}

void AudioPlugin::getStateInformation(juce::MemoryBlock& destData)
{
  if(wrappedAudioModule != nullptr)
  {
    // todo: store values of the MetaParameters

    //XmlElement* xml = wrappedAudioModule->getStateAsXml("StateAsRequestedByHost", false);
    XmlElement* xml = wrappedAudioModule->getStateAsXml(wrappedAudioModule->getStateName(), false);
    xml->setAttribute("EditorWidth",  editorWidth);
    xml->setAttribute("EditorHeight", editorHeight);
    copyXmlToBinary(*xml, destData);
    delete xml;
  }
}

void AudioPlugin::setStateInformation(const void* data, int sizeInBytes)
{
  if(wrappedAudioModule != nullptr)
  {
    // todo: retrieve values of the MetaParameters

    //XmlElement* const xml = getXmlFromBinary(data, sizeInBytes);  // old
    XmlElement* const xml = new XmlElement(*(getXmlFromBinary(data, sizeInBytes).get()));  // new, preliminary

    //ParameterObserver::globalAutomationSwitch = false; // why this - threading problems? -> interferes with total recall in quadrifex
    ParameterObserver::setGuiAutomationSwitch(false);
    wrappedAudioModule->setStateFromXml(*xml, "recalled by host", false);
    editorWidth  = xml->getIntAttribute("EditorWidth",  0);
    editorHeight = xml->getIntAttribute("EditorHeight", 0);
    ParameterObserver::setGuiAutomationSwitch(true);
    //ParameterObserver::globalAutomationSwitch = true;
    delete xml;

    // some hosts (Tracktion) seem to keep and re-use an open GUI when a new project is loaded, so
    // we must make sure, that this re-used GUI is updated according to the new recalled state - we
    // do this by broadcasting a changeMessage which will be picked up by AudioPlugInEditor
    //sendChangeMessage();
    // (this comment is old and refers to the old way of doing it ...but i want to verify this in
    // Tracktion someday)
  }
}

// optional overrides for juce::AudioProcessor:

bool AudioPlugin::isBusesLayoutSupported(const BusesLayout& layout) const
{
  bool r = true;  // result
  int numIns  = layout.getNumChannels(true,  0);
  int numOuts = layout.getNumChannels(false, 0);
  r &= numIns  == 2;
  r &= numOuts == 2;
  return r;
}

// Obsolete - maybe delet or move ot attic:
inline void AudioPlugin::enableFitzdazzing()
{
#if defined(RS_ARCHITECTURE_X64)

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);            // #defined in xmmintrin.h
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);    // #defined in pmmintrin.h
  // see here: https://software.intel.com/en-us/node/523328

#elif defined(RS_ARCHITECTURE_ARM64)


  // see: https://developer.arm.com/documentation/dui0473/m/neon-programming/when-to-use-flush-to-zero-mode-in-neon
  // 
  // todo: implement unit tests for treatment of denormals 
  // -try to feed dernormal inputs into processBlock, this:
  //  https://developer.arm.com/documentation/dui0473/m/pge1423730082105
  //  says that this should trigger an exception. could this become problematic?
  // -see also: https://discourse.julialang.org/t/50x-speed-difference-in-gemv-for-different-values-in-vector/2755/6

#endif

  // See also:
  // https://blog.audio-tk.com/2016/09/20/audio-toolkit-handling-denormals/
}

void AudioPlugin::processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages)
{
  juce::ScopedLock scopedLock(plugInLock);  // Acquire mutex lock
  juce::ScopedNoDenormals scopedDenormals;  // Temoprarily disable denormals.

  // Obsolete:
  // ToDo: int ftzDazState = getFtzDazState();
  //enableFitzdazzing();
  // ToDo: maybe call only once in prepareToPlay() - but no, other plugins in the chain may 
  // potentially modify this between calls
  // Use juce::ScopedDenormals (or something like that) instead

  if(wrappedAudioModule != nullptr)
  {
    int numChannels = buffer.getNumChannels();
    int numSamples  = buffer.getNumSamples();
    double **inOutBuffer = buffer.getArrayOfWritePointers();
    wrappedAudioModule->processBlock(inOutBuffer, numChannels, numSamples);
  }
  else
    buffer.clear();  // Or maybe we should juts leave it as is?


  // setFtzDazState(ftzDazState);
  // ...we should restore the original FTZ/DAZ settings here? other plugins down the chain
  // may want to have denormals...well...that's very unlikely, but still

  // maybe we should have a function that avoids acquiring the lock and setting denormal flags that
  // can be called from the AudioPluginWithMidiIn subclass (assuming that it already holds the lock
  // and has set the denormal flags)
}

/*
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
*/

//// i tried to override it in order to check, if it gets called whe we call updateHostDisplay() in
//// metaNameChanged - but at least the juce host doesn't try to inquire a new name :-(
//String AudioPlugin::getParameterName(int parameterIndex, int maximumStringLength)
//{
//  return AudioProcessor::getParameterName(parameterIndex, maximumStringLength);
//}

void AudioPlugin::createHostAutomatableParameters(int numParameters)
{
  parameters.resize(numParameters);
  for(int i = 0; i < numParameters; i++)
  {
    AudioPluginParameter* p = parameters[i] = new AudioPluginParameter();
    p->setName("Meta " + String(i));
    addParameter(p);
    metaParaManager.addMetaParamater(p);
  }
}

//=================================================================================================

void AudioPluginWithMidiIn::processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages)
{
  juce::ScopedLock scopedLock(plugInLock);  // Acquire mutex lock
  juce::ScopedNoDenormals scopedDenormals;  // Temoprarily disable denormals.

//#ifdef DEBUG
//  static int callCount = 1;
//  if( underlyingAudioModule == NULL )
//  {
//    buffer.clear();
//    jassertfalse;
//    return;
//  }
//  _clearfp();                                        // reset
//  unsigned int cw = _controlfp(0, 0);                // get state
//  //cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
//  cw &= ~(EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
//  unsigned int cw0 = _controlfp(cw, _MCW_EM);
//#endif
// This crashes reaper. What was this supposed to do anyway? Maybe delete.



  /*
  // OLD:
  // Request time-info from host and update the bpm-values for the modulators accordingly, if they
  // are in sync mode (maybe this should be moved up into the baseclass AudioModule):
  int timeToNextTriggerInSamples = -1;
  if( wrappedAudioModule->wantsTempoSyncInfo )
  {
    AudioPlayHead::CurrentPositionInfo info;
    if( getPlayHead() != 0  &&  getPlayHead()->getCurrentPosition(info) )
    {
      if( info.bpm <= 0.0 )
        info.bpm = 120.0;  // Fallback value when nothing meaningful is passed
      wrappedAudioModule->setBeatsPerMinute(info.bpm);
    }

    if( wrappedAudioModule->getTriggerInterval() != 0.0 )
    {
      double timeInBeats = RAPT::rsSecondsToBeats(info.timeInSeconds, info.bpm);
      // kludge - will probably not work when tempo changes - use ppqPosition here later....
      // Also: "info" may contain garbage data in the case when the above conditional statement
      // failed and getPlayHead()->getCurrentPosition(info) didn't do anything

      double timeToNextTriggerInBeats = wrappedAudioModule->getTriggerInterval()
        - fmod(timeInBeats, wrappedAudioModule->getTriggerInterval());

      timeToNextTriggerInSamples =
        roundToInt(getSampleRate()*RAPT::rsBeatsToSeconds(timeToNextTriggerInBeats, info.bpm));

      if( timeToNextTriggerInSamples >= buffer.getNumSamples() )
        timeToNextTriggerInSamples = -1; // indicates that we don't need to trigger in this block
    }
  }
  // ...soo - this whole time-info retrieval code is very quenstionable and should probably be 
  // rewritten entirely. Maybe declare a local variable double bpm = 120; and assign it from the
  // info, if available. Then use that inside the two conditionals. In particular, don't access
  // info.bpm inside the 2nd inner conditional block. Maybe enter 2nd block only, if data was 
  // retrieved
  // ...OK - done - see below. When enough tests were made, this old code may be deleted. We'll
  // just keep it around for a while for reference, if something goes wrong with the new code.
  */

  // NEW - needs more tests - it seems to work fine in standalone, i.e. with the fallback values:
  int timeToNextTriggerInSamples = -1;
  if(wrappedAudioModule->wantsTempoSyncInfo)
  {
    // Initialize the relevant position info variables which we are interested in to their default
    // values which are used as fallback values. For example, they will be used in a standalone 
    // version or when a plugin host for some reason doesn't provide this kind of information:
    double bpm = 140.0;          // Current tempo in beats per minute.
    double timeInSeconds = 0.0;  // Used for periodic retriggering of LFOs and similar stuff.

    // Now try to retrieve the actual position info from host:
    //bool positionInfoRetrieved = false;
    juce::AudioPlayHead* playHead = getPlayHead();
    if(playHead != nullptr)
    {
      AudioPlayHead::CurrentPositionInfo positionInfo;
      if(playHead->getCurrentPosition(positionInfo))
      {
        //positionInfoRetrieved = true;  // This seems to be not used anywhere - maybe get rid
        bpm = positionInfo.bpm;
        timeInSeconds = positionInfo.timeInSeconds;
      }
    }

    // Pass tempo info (either fallback or retrieved value) to wrapped AudioModule:
    wrappedAudioModule->setBeatsPerMinute(bpm);

    // Do the periodic triggering, if desired. This is for syncing internal LFOs and similar stuff 
    // to the host (if I remember correctly - ToDo: check this and document it properly!):
    if( wrappedAudioModule->getTriggerInterval() != 0.0 )
    {
      double timeInBeats = RAPT::rsSecondsToBeats(timeInSeconds, bpm);
      // Kludge! Will probably not work when the tempo changes over time. Maybe use a ppqPosition 
      // variable here later (which also needs to be initialized to a fallback value and possibly 
      // re-assigned from positionInfo, if available)


      double timeToNextTriggerInBeats = wrappedAudioModule->getTriggerInterval()
        - fmod(timeInBeats, wrappedAudioModule->getTriggerInterval());

      timeToNextTriggerInSamples =
        roundToInt(getSampleRate() * RAPT::rsBeatsToSeconds(timeToNextTriggerInBeats, bpm));

      if( timeToNextTriggerInSamples >= buffer.getNumSamples() )
        timeToNextTriggerInSamples = -1; // indicates that we don't need to trigger in this block    
    }
    // ToDo: check, if ToolChainAudioModule handles these retriggering calls correctly and passes 
    // the triggers on to its sub-modules. Well - maybe it doesn't even need to. Not sure. The 
    // modules might obtain their relevant triggers from midi notes. May depend on the module.
  }
  // Maybe wrap into a function updatePositionInfo(). Maybe it should return the 
  // timeToNextTriggerInSamples



  if(midiMessages.isEmpty() && timeToNextTriggerInSamples == -1)
  {
    //underlyingAudioModule->processBlock(buffer, 0, buffer.getNumSamples()); // old
    AudioPlugin::processBlock(buffer, midiMessages);
    // No messages or trigger events to process - the baseclass method can handle this.
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
        wrappedAudioModule->trigger();

        // Update the time to the next trigger:
        timeToNextTriggerInSamples = -1;
        // Preliminary - assumes that the next trigger is not inside this buffer (which is assured
        // only when buffersize is much smaller than the retrigger interval - which is probably
        // true in realtime situations but maybe not in rendering)
      }

      // All messages at this time instant have been processed - now determine the time-offset to
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
      wrappedAudioModule->processBlock(channelPointers, 2, subLength);

      start += subLength; // increment start position for next iteration

      //midiMessages.clear(0, start+subLength); // test
    }
  }

//#ifdef DEBUG
//  callCount++;
//  _controlfp(cw0, _MCW_EM);
//#endif
}

void AudioPluginWithMidiIn::handleMidiMessage(MidiMessage message)
{
  ScopedLock sl(plugInLock); 
  // Is this really needed? I think, it's called only from processBlock which supposedly already 
  // holds the lock. Figure this out and try to get rid of it.

  // We ingore midi clock messages on the highest level (i.e. here) because they are spammy:
  if( message.isMidiClock() )
    return;
  // ToDo: Figure out if we can tell the host that we don't want to receive midi clock messages.
  // The standalone version tends to get bombarded with them (at least, when the JX-305 is 
  // connected via the midi interface. And/or maybe we can pre-filter the midi event buffer in
  // processBlock. But maybe just ignoring the clock events is more efficient than such a 
  // filtering pass.

  // Pass the midi message on to the underlying AudioModule:
  if( wrappedModuleWithMidiIn != nullptr )
    wrappedModuleWithMidiIn->handleMidiMessage(message);
}

void AudioPluginWithMidiIn::setAudioModuleToWrap(AudioModule* moduleToWrap)
{
  AudioPlugin::setAudioModuleToWrap(moduleToWrap);
  wrappedModuleWithMidiIn = dynamic_cast<AudioModuleWithMidiIn*>(moduleToWrap);
  jassert(wrappedModuleWithMidiIn != nullptr); // trying to wrap a non-midi-enabled AudioModule
                                               // into a midi-enabled AudioPlugin?
}

//=================================================================================================

AudioPluginEditor::AudioPluginEditor(AudioModuleEditor* editorToWrap, AudioPlugin* pluginToEdit)
  : AudioProcessorEditor(pluginToEdit)
{
  this->pluginToEdit  = pluginToEdit;
  this->wrappedEditor = editorToWrap;

  // retrieve desired size from values stored in the plugin:
  int w = pluginToEdit->editorWidth;
  int h = pluginToEdit->editorHeight;

  // ...unless there are none stored - then use intitial size of the editor that we wrap:
  if(w == 0)
    w = wrappedEditor->getWidth();
  if(h == 0)
    h = wrappedEditor->getHeight();

  // limit the size according to desired limits:
  //w = RAPT::rsClip(w, pluginToEdit->editorWidthMin,  pluginToEdit->editorWidthMax);
  //h = RAPT::rsClip(h, pluginToEdit->editorHeightMin, pluginToEdit->editorHeightMax);
  // ..but maybe that's not a good idea - we want to allow the gui size to be recalled 
  // programatically, even when it'S not resizable by the user

  // set up resizing limits and initial size:
  bool resizable = pluginToEdit->editorWidthMin  != pluginToEdit->editorWidthMax;
  resizable     &= pluginToEdit->editorHeightMin != pluginToEdit->editorHeightMax;
  setResizable(resizable, resizable);
  setResizeLimits(pluginToEdit->editorWidthMin, pluginToEdit->editorHeightMin, 
    pluginToEdit->editorWidthMax, pluginToEdit->editorHeightMax);// must be called BEFORE setSize
  setSize(w, h);

  addAndMakeVisible(wrappedEditor);
}