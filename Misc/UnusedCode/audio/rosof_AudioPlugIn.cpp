#include "rosof_AudioPlugIn.h"
#include "rosof_AudioPlugInEditor.h"
using namespace rosof;

AudioPlugIn::AudioPlugIn()
{
  checkForMemoryLeaksOnExit();
  //juce::String *leakString = new juce::String(T("String that should cause a memory leak")); // to test if leak-detection works

  ScopedLock sl(plugInLock);

  // experimental:
  initialiseJuce_GUI();

  underlyingAudioModule = NULL;
}

AudioPlugIn::~AudioPlugIn()
{    
  ScopedLock sl(plugInLock);
  if( underlyingAudioModule != NULL )
  {
    delete underlyingAudioModule;
    underlyingAudioModule = NULL;
  }
}

const juce::String AudioPlugIn::getName() const
{
  /*
  return juce::String::empty;
  const void *oldThis = this;
  juce::String result = plugInName;
  const void *newThis = this;
  if( oldThis != newThis )
    jassertfalse; // for some weird reason, the line result = plugInName; messes with the this-pointer
  return result;
  */
  return plugInName;
}

int AudioPlugIn::getNumParameters()
{
  return 0;
}

float AudioPlugIn::getParameter(int index)
{
  return 0.f;
}

void AudioPlugIn::setParameter(int index, float newValue)
{

}

const juce::String AudioPlugIn::getParameterName(int index)
{
  return juce::String::empty;
}

const juce::String AudioPlugIn::getParameterText(int index)
{
  return juce::String::empty;
}

const juce::String AudioPlugIn::getInputChannelName(const int channelIndex) const
{
  if( channelIndex == 0 )
    return juce::String(T("Left Input"));
  else if( channelIndex == 1 )
    return juce::String(T("Right Input"));
  else
    return juce::String(T("Input Channel ")) + juce::String(channelIndex+1);
}

const juce::String AudioPlugIn::getOutputChannelName(const int channelIndex) const
{
  if( channelIndex == 0 )
    return juce::String(T("Left Output"));
  else if( channelIndex == 1 )
    return juce::String(T("Right Output"));
  else
    return juce::String(T("Output Channel ")) + juce::String(channelIndex+1);
}

bool AudioPlugIn::isInputChannelStereoPair(int index) const
{
  return true;
}

bool AudioPlugIn::isOutputChannelStereoPair(int index) const
{
  //return true;
  //return false;

  if( index % 2 == 0 )
    return true;
  else
    return false;
}

bool AudioPlugIn::acceptsMidi() const
{
  return true;
}

bool AudioPlugIn::producesMidi() const
{
  return false;
}

void AudioPlugIn::prepareToPlay(double sampleRate, int samplesPerBlock)
{
  ScopedLock sl(plugInLock);
  if( underlyingAudioModule != NULL )
    underlyingAudioModule->setSampleRate(sampleRate);
}

void AudioPlugIn::releaseResources()
{
  // when playback stops, you can use this as an opportunity to free up any
  // spare memory, etc.
}

void AudioPlugIn::processBlock(AudioSampleBuffer &buffer, MidiBuffer &midiMessages)
{
  ScopedLock sl(plugInLock);

#ifdef DEBUG
  static int callCount = 1;
  if( underlyingAudioModule == NULL )
  {
    buffer.clear();
    jassertfalse;
    return;
  }
  _clearfp();                                                                     // reset
  unsigned int cw = _controlfp(0, 0);                                             // get state
  //cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
  cw &= ~(EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID); // unmask all exceptions
  unsigned int cw0 = _controlfp(cw, _MCW_EM);
#endif



  // request time-info from host and update the bpm-values for the modulators accordingly (if they are in sync mode):
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

  if( midiMessages.isEmpty() && timeToNextTriggerInSamples == -1 )
    underlyingAudioModule->processBlock(buffer, 0, buffer.getNumSamples());
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
      currentMidiMessageOffset = -1; // may be set to another value by the function in the next line
      aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
        currentMidiMessage, currentMidiMessageOffset);
      while( aMidiMessageWasRetrieved && currentMidiMessageOffset == i ) 
      {
        handleMidiMessage(currentMidiMessage);                 
          // respond to the midi-message
        aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(currentMidiMessage, currentMidiMessageOffset);      
          // get the next message at this sample
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
      // the next event and compute the number of samples that can be processed without interrupting
      // events (subLength):    
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
      underlyingAudioModule->processBlock(buffer, start, subLength);
      start += subLength;

      //midiMessages.clear(0, start+subLength);
      //int dummy = 0;
    }
  }


#ifdef DEBUG
  callCount++;
  _controlfp(cw0, _MCW_EM);
#endif
}

void AudioPlugIn::handleMidiMessage(MidiMessage message)
{
  ScopedLock sl(plugInLock);
  if( message.isMidiClock() )
    return;
  if( underlyingAudioModule != NULL )
    underlyingAudioModule->handleMidiMessage(message);
}

//==============================================================================

AudioProcessorEditor* AudioPlugIn::createEditor()
{
  // its not a good idea to update slider during the construction of the GUI, so we temporarily 
  // deactivate automation:
  ParameterObserver::guiAutomationSwitch = false;

  AudioProcessorEditor* editor = new AudioPlugInEditor(this);

  // switch automation for GUI object on again and return the editor:
  ParameterObserver::guiAutomationSwitch = true;
  return editor;
}

void AudioPlugIn::getStateInformation(MemoryBlock& destData)
{
  XmlElement* xmlState = getStateAsXml();
  copyXmlToBinary(*xmlState, destData);
  delete xmlState;
}

void AudioPlugIn::setStateInformation (const void* data, int sizeInBytes)
{
  XmlElement* const xmlState = getXmlFromBinary(data, sizeInBytes);
  //ParameterObserver::globalAutomationSwitch = false; // why this - threading problems? -> interferes with total recall in quadrifex
  ParameterObserver::guiAutomationSwitch = false;
  setStateFromXml(*xmlState, juce::String(T("recalled by host")));
  ParameterObserver::guiAutomationSwitch = true;
  //ParameterObserver::globalAutomationSwitch = true;
  delete xmlState;

  // some hosts (Tracktion) seem to keep and re-use an open GUI when a new project is loaded, so
  // we must make sure, that this re-used GUI is updated according to the new recalled state - we 
  // do this by broadcasting a changeMessage which will be picked up by AudioPlugInEditor
  sendChangeMessage();
}

XmlElement* AudioPlugIn::getStateAsXml(XmlElement* xmlElementToStartFrom) const
{
  ScopedLock sl(plugInLock);
  if( underlyingAudioModule != NULL )
    return underlyingAudioModule->getStateAsXml(juce::String(T("StateAsRequestedByHost")), false);
  else
    return xmlElementToStartFrom;
}

void AudioPlugIn::setStateFromXml(const XmlElement& xmlState, const juce::String& name)
{
  ScopedLock sl(plugInLock);
  if( underlyingAudioModule != NULL )
    underlyingAudioModule->setStateFromXml(xmlState, name, false);
}
