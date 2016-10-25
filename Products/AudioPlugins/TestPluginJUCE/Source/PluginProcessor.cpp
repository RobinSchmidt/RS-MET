#include "PluginProcessor.h"
#include "PluginEditor.h"

#include "../../../../Libraries/RAPT/Code/Library/RAPT.cpp"

TestPluginAudioProcessor::TestPluginAudioProcessor()
{
}

//==============================================================================
const String TestPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool TestPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool TestPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

double TestPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}
int TestPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}
int TestPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}
void TestPluginAudioProcessor::setCurrentProgram (int index)
{
}
const String TestPluginAudioProcessor::getProgramName (int index)
{
    return String();
}
void TestPluginAudioProcessor::changeProgramName (int index, const String& newName)
{
}


void TestPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
  ladder.setSampleRate(sampleRate);
  ladder.reset();  // reset internal buffers
}

void TestPluginAudioProcessor::releaseResources()
{
  // When playback stops, you can use this as an opportunity to free up any  
  // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool TestPluginAudioProcessor::setPreferredBusArrangement (bool isInput, int bus, 
  const AudioChannelSet& preferredSet)
{
    // Reject any bus arrangements that are not compatible with your plugin

    const int numChannels = preferredSet.size();

   #if JucePlugin_IsMidiEffect
    if (numChannels != 0)
        return false;
   #elif JucePlugin_IsSynth
    if (isInput || (numChannels != 1 && numChannels != 2))
        return false;
   #else
    if (numChannels != 1 && numChannels != 2)
        return false;

    if (! AudioProcessor::setPreferredBusArrangement (! isInput, bus, preferredSet))
        return false;
   #endif

    return AudioProcessor::setPreferredBusArrangement (isInput, bus, preferredSet);
}
#endif
// can this function be factored out into a baseclass to avoid boilerplate code for other plugins?

void TestPluginAudioProcessor::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
  if(buffer.getNumChannels() != 2)
  {
    buffer.clear();
    return;
  }
  float *left  = buffer.getWritePointer(0);
  float *right = buffer.getWritePointer(1);
  complex<double> tmp;
  for(int n = 0; n < buffer.getNumSamples(); n++)
  {
    tmp.real((double)left[n]);      // Set up real...
    tmp.imag((double)right[n]);     // ...and imaginary part of temporary variable from input.
    tmp = ladder.getSample(tmp);    // Process one (stereo) sample frame at a time.
    left[n]  = (float)tmp.real();   // Extract real...
    right[n] = (float)tmp.imag();   // ...and imaginary part and write to output.
  }
}

//==============================================================================
bool TestPluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* TestPluginAudioProcessor::createEditor()
{
    return new TestPluginAudioProcessorEditor (*this);
}

void TestPluginAudioProcessor::setParameter(int index, float value)
{
  switch(index)
  {
  case CUTOFF: 
  {
    cutoff = RAPT::rsLinToExp(value, 0, 1, 20, 20000);
    ladder.setCutoff(cutoff);
  } break;
  case RESO: 
  {
    reso = value;
    ladder.setResonance(reso);
  } break;
  }
}
void TestPluginAudioProcessor::setFilterMode(int newMode)
{
  ladder.setMode(newMode);
}

//==============================================================================
void TestPluginAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void TestPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TestPluginAudioProcessor();
}
