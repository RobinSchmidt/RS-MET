#ifndef jura_AudioPlugin_h
#define jura_AudioPlugin_h

/** Subclass of juce::AudioProcessorParameter to provide the handling of host automation. It also
derives from jura::MetaParameter in order to provide the "glue" between juce's host automation
handling and jura's MetaParameter handling. Whenever the setValue method, inherited and overriden 
from AudioProcessorParameter, gets called (by the host), we will call MetaParameter's setValue 
method there which in turn will update all the values of the attached MetaControlledParameters. 

\todo override parameterChanged (inherited from MetaParameter) in order to notify host, we need
to call setValueNotifyingHost(float newValue) ---is there a way to only notify the host without
having it internally call setValue()? that would be more convenient bcs otherwise we must somehow
avoid an endless recursion of callbacks


*/

class JUCE_API AudioPluginParameter : public AudioProcessorParameter, public MetaParameter
{

public:

  AudioPluginParameter() {}

  /** Constructor. You n */
  //AudioPluginParameter(const String& parameterName) : name(parameterName) {}

  // mandatory AudioProcessorParameter overrides:
  virtual float getValue() const override { return (float) metaValue; }
  virtual void setValue(float newValue) override;
  virtual float getDefaultValue() const override { return 0.f; }
  virtual String getName(int maximumStringLength) const override { return name; }
  virtual String getLabel() const override { return String::empty; }
  virtual float getValueForText(const String &text) const override { return text.getFloatValue(); }

  // optional AudioProcessorParameter overrides:
  virtual bool isAutomatable() const override { return true; }
  virtual bool isMetaParameter() const override { return true; }

  // overriden from MetaParameter to notify host:
  virtual void parameterChanged(Parameter* p) override;

  /** Sets the name of the parameter that is reported to the host. */
  virtual void setName(const String& newName);
    // move to MetaParameter, rename to setMetaName

protected:

  juce::String name;  

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginParameter)
};


//=================================================================================================

/** This class wraps a jura::AudioModule into a juce::AudioProcessor, so we can use it as 
plugin.  
\todo - mayb move this pair of h/cpp files into another folder */

class JUCE_API AudioPlugin : public AudioProcessor
{

public:

  AudioPlugin(); 

  //AudioPlugin(AudioModule *moduleToWrap); 
  // remove the parameter from the constructor - we should now use setAudioModuleToWrap after the
  // constructor has finished
  

  virtual ~AudioPlugin();

  //-----------------------------------------------------------------------------------------------
  // setup:

  void setPluginName(const String& newName) { plugInName = newName; }

  /** Sets the jura::AudioModule object that is wrapped into an AudioPlugin. */
  virtual void setAudioModuleToWrap(AudioModule* moduleToWrap);

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

  //virtual void setParameter(int parameterIndex, float newValue) override; // preliminary

  virtual int getNumPrograms() override { return 1; }                // 1, because 0 is not allowed
  virtual int getCurrentProgram() override { return 0; }
  virtual void setCurrentProgram(int index) override {}
  virtual const String getProgramName (int index) override { return String::empty; }
  virtual void changeProgramName(int index, const String& newName) override {}
  virtual void getStateInformation(juce::MemoryBlock& destData) override;
  virtual void setStateInformation(const void* data, int sizeInBytes) override; 
  virtual void processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages) override;

  //-----------------------------------------------------------------------------------------------
  // optional overrides for juce::AudioProcessor baseclass:

  virtual bool supportsDoublePrecisionProcessing() const override { return true; }
  virtual void processBlock(AudioBuffer<double>& buffer, MidiBuffer& midiMessages) override;
  virtual bool setPreferredBusArrangement(bool isInput, int bus,
    const AudioChannelSet& preferredSet) override;

  //-----------------------------------------------------------------------------------------------
  // data:

  AudioModule *wrappedAudioModule = nullptr;  // the wrapped jura::AudioModule

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

  // parameter-management:
  //static const int numParameters = 128;
  static const int numParameters = 10;
  AudioPluginParameter* parameters[numParameters];
  MetaParameterManager metaParaManager;

  juce::String plugInName;  // assign this in the constructor of your subclass
   // maybe get rid of this and let the wrapper return the name of the wrapped AudioModule

  int editorWidth  = 0;
  int editorHeight = 0;

  friend class AudioPluginEditor;

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

  AudioPluginWithMidiIn() {}
  //AudioPluginWithMidiIn(AudioModuleWithMidiIn *moduleToWrap);

  virtual bool acceptsMidi() const override { return true; }

  virtual void processBlock(AudioBuffer<double>& buffer, MidiBuffer& midiMessages) override;

  virtual void handleMidiMessage(MidiMessage message);

  virtual void setAudioModuleToWrap(AudioModule* moduleToWrap) override;


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

  AudioPluginEditor(AudioModuleEditor* editorToWrap, AudioPlugin* pluginToEdit) 
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

    setResizeLimits(200, 100, 6000, 3000); // must be called BEFORE setSize
    setSize(w, h);
    addAndMakeVisible(wrappedEditor); 
  }

  AudioPluginEditor::~AudioPluginEditor()
  {
    delete wrappedEditor;  
    // we need to delete it here because the baseclass destructor does not delete its child 
    // components - is this a change with respect to the old juce?
  }

  virtual void paint(Graphics &g) override {} // we hit a breakpoint if we don't override this

  virtual void resized() override
  {
    wrappedEditor->setBounds(0, 0, getWidth(), getHeight());

    // store the size in the plugin, so it can be recalled with the DAW project:
    pluginToEdit->editorWidth  = getWidth();
    pluginToEdit->editorHeight = getHeight();
  }

protected:

  AudioModuleEditor* wrappedEditor;

  AudioPlugin* pluginToEdit;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginEditor)
};

#endif