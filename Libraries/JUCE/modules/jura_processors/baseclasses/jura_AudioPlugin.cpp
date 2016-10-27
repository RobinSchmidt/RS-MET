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
