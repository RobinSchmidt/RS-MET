#ifndef jura_AudioPlugin_h
#define jura_AudioPlugin_h

/** This class wraps a jura::AudioModule into a juce::AudioProcessor, so we can use it as 
plugin.  
\todo - mayb move this pair of h/cpp files into another folder */

class JUCE_API AudioPlugin : public AudioProcessor
{

public:

  AudioPlugin(AudioModule *moduleToWrap);
  virtual ~AudioPlugin();

  //-----------------------------------------------------------------------------------------------
  // mandatory overrides for juce::AudioProcessor baseclass:

  virtual const String getName() const override { return plugInName; }
  virtual void  prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) override;
  virtual void releaseResources() override {}
  virtual double getTailLengthSeconds() const override { return 0.0; }
  virtual bool acceptsMidi() const override 
  { 
    return false; // doesn't seem to get called by JUCE's plugin host ..or maybe only on scan?
  }
  virtual bool producesMidi() const override { return false; }
  virtual bool hasEditor() const override { return true; }
  virtual AudioProcessorEditor* createEditor() override;
  virtual int getNumPrograms() override { return 1; }                // 1, because 0 is not allowed
  virtual int getCurrentProgram() override { return 0; }
  virtual void setCurrentProgram(int index) override {}
  virtual const String getProgramName (int index) override { return String::empty; }
  virtual void changeProgramName(int index, const String& newName) override {}
  virtual void getStateInformation(juce::MemoryBlock& destData) override;
  virtual void setStateInformation(const void* data, int sizeInBytes); 
  virtual void processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages) override;

  //-----------------------------------------------------------------------------------------------
  // optional overrides for juce::AudioProcessor baseclass:

  virtual bool supportsDoublePrecisionProcessing() const override { return true; }
  virtual void processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages) override;
  virtual bool setPreferredBusArrangement(bool isInput, int bus,
    const AudioChannelSet& preferredSet) override;

  //-----------------------------------------------------------------------------------------------
  // data:

  AudioModule *underlyingAudioModule;  
   // the wrapped jura::AudioModule - rename to wrappedAudioModule

  /** Mutex-lock for all accesses to the underlyingAudioModule's member functions - a pointer to 
  the lock is  passed to the embedded AudioModule and should be used there also and the AudioModule
  should also pass this lock on to the GUI Editors. */
  CriticalSection plugInLock; 


protected:

  /** The number of channels that is desired for the in/out buffer that is passed to the 
  processBlock callback. You may set that value in the constructor of your subclass. If the number
  of channels is supposed to change after construction, we may have to make sure that we are in a 
  suspended state before we change that value. I did not yet run into this situation, so i haven't 
  figured it out. */
  int numChannels = 2;

  /** An internal double precision buffer that is used in cases, where the host calls the single
  precision version of the processBlock callback. In such a case, we need to convert back and forth
  between float/double and double/float. That's what this buffer is used for. */
  AudioBuffer<double> internalAudioBuffer; 

  juce::String plugInName;  // assign this in the constructor of your subclass
   // maybe get rid of this and let the wrapper return the name of the wrapped AudioModule

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPlugin)
};

//=================================================================================================

/** A wrapper class that wraps an object of a subclass of AudioModuleWithMidiIn into a 
juce::AudioProcessor. We need a different wrapper class for modules with MIDI than for those 
without MIDI because we must inform the host that we want to receive MIDI events (for this we 
override acceptsMidi here) and furthermore, we actually have to handle the events by passing them 
to the event handler methods of AudioModuleWithMidiIn - which the baseclass AudioModule doesn't 
even have. We will call these event handlers from our overriden processBlock method. */

class JUCE_API AudioPluginWithMidiIn : public AudioPlugin
{

public:

  AudioPluginWithMidiIn(AudioModuleWithMidiIn *moduleToWrap);

  virtual bool acceptsMidi() const override { return true; }

  virtual void processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages) override;

  virtual void handleMidiMessage(MidiMessage message);

  AudioModuleWithMidiIn *wrappedModuleWithMidiIn;  

protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginWithMidiIn)
};



//=================================================================================================

/** A wrapper class that wraps an object of a subclass of AudioModuleEditor into a
juce::AudioProcessorEditor. */

class AudioPluginEditor : public juce::AudioProcessorEditor
{

public:

  AudioPluginEditor(AudioModuleEditor *newContentComponent, AudioProcessor* processorToEdit) 
    : AudioProcessorEditor(processorToEdit)
  {
    contentComponent = newContentComponent;
    setSize(contentComponent->getWidth(), contentComponent->getHeight());
    addAndMakeVisible(contentComponent);
  }

  AudioPluginEditor::~AudioPluginEditor()
  {
    delete contentComponent;  
    // we need to delete it here because the baseclass destructor does not delete it's child 
    // components - is this a change with respect to the old juce?
  }

  virtual void paint(Graphics &g) override {} // we hit a breakpoint if we don't override this

protected:

  AudioModuleEditor *contentComponent;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginEditor)
};

#endif