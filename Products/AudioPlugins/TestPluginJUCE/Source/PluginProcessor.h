#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED

#include "../JuceLibraryCode/JuceHeader.h"

#include <complex>    // abused for representing stereo signals
using namespace std;

#include "../../../../Libraries/RAPT/Source/Modules/RAPT.h"
using namespace RAPT;

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
  bool setPreferredBusArrangement (bool isInput, int bus, const AudioChannelSet& preferredSet) override;
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

  LadderFilter<complex<double>, double> ladder;
    // The embedded RAPT LadderFilter object - we use complex numbers for the signal type to 
    // represent stereo signals (real part: left channel, imaginary part: right channel).
    // ToDo: maybe use a different stereo-pair type later, based on a SIMD datatype that represents
    // two doubles - should reduce processing cost by roughly 1/2? Maybe provide such a type in 
    // RAPT - Float64X2 based on this __mm128 (or however it is called) instrinsic, maybe also 
    // Float32X4 (4 single precision floats)

private:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TestPluginAudioProcessor)
};

#endif  // PLUGINPROCESSOR_H_INCLUDED