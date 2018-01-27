#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED

#include "../JuceLibraryCode/JuceHeader.h"


/** A test plugin implementing a multimode ladder filter  */

class TestPluginAudioProcessor : public AudioProcessor
{
public:

  // enumeration of parameters:
  enum
  {
    CUTOFF = 0,
    RESO
  };

  TestPluginAudioProcessor();

  void prepareToPlay (double sampleRate, int samplesPerBlock) override;
  void releaseResources() override;

#ifndef JucePlugin_PreferredChannelConfigurations
  //bool setPreferredBusArrangement (bool isInput, int bus, const AudioChannelSet& preferredSet) override;
#endif

  void processBlock (AudioSampleBuffer&, MidiBuffer&) override;


  AudioProcessorEditor* createEditor() override;
  bool hasEditor() const override;

  const String getName() const override;

  bool acceptsMidi() const override;
  bool producesMidi() const override;
  double getTailLengthSeconds() const override;

  //==============================================================================
  // factor these program handling functions out into a baseclass - we don't use program
  // handling:
  int getNumPrograms() override;
  int getCurrentProgram() override;
  void setCurrentProgram (int index) override;
  const String getProgramName (int index) override;
  void changeProgramName (int index, const String& newName) override;

  /** Will be called back by the host for automation, the GUI will also call this. */
  void setParameter(int parameterIndex, float newValue) override;

  /** Sets the filter mode. Called from the GUI */
  void setFilterMode(int newMode);

  //==============================================================================
  // state load/save is not yet implemented:
  void getStateInformation (MemoryBlock& destData) override;
  void setStateInformation (const void* data, int sizeInBytes) override;

protected:

  double cutoff, reso;  // mapped parameters - todo: use AudioProcessorParameter class
  rosic::rsLadderFilterStereo ladder;

private:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPluginAudioProcessor)
};

#endif  // PLUGINPROCESSOR_H_INCLUDED