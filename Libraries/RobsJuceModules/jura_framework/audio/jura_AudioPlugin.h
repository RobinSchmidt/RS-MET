#ifndef jura_AudioPlugin_h
#define jura_AudioPlugin_h

/** Subclass of juce::AudioProcessorParameter to provide the handling of host automation. It also
derives from jura::MetaParameter in order to provide the "glue" between juce's host automation
handling and jura's MetaParameter handling. Whenever the setValue method (inherited and overriden
from AudioProcessorParameter) gets called by the host, we will call MetaParameter's setMetaValue
method there which in turn will update all the values of the attached MetaControlledParameters.
Whenever parameterChanged (overriden from MetaParameter) gets called, which happens when the user
changes a parameter on the plugin's gui, we retrieve the value and call setValueNotifyingHost
(inherited from AudioProcessorParameter) which notifies the host about the change and then calls
our setValue (which in turn updates other Parameters attached to this MetaParameter, if any). */

class JUCE_API AudioPluginParameter : public AudioProcessorParameter, public MetaParameter
{

public:

  AudioPluginParameter() {}

  // mandatory AudioProcessorParameter overrides:
  virtual float getValue() const override { return (float) metaValue; }
  virtual void setValue(float newValue) override;
  virtual float getDefaultValue() const override { return 0.5f; }
  virtual String getName(int maximumStringLength) const override { return name; }
  virtual String getLabel() const override { return String(); }
  virtual float getValueForText(const String &text) const override { return text.getFloatValue(); }

  // optional AudioProcessorParameter overrides:
  virtual bool isAutomatable() const override { return true; }
  virtual bool isMetaParameter() const override { return true; }

  // overriden from MetaParameter to notify host:
  virtual void parameterChanged(Parameter* p) override;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginParameter)
};


//=================================================================================================

/** This class wraps a jura::AudioModule into a juce::AudioProcessor, so we can use it as
plugin. */

class JUCE_API AudioPlugin : public AudioProcessor, public MetaParameterManagerObserver
{

public:

  /** Constructor. You can pass a number of parameters that this plugin will report to the host and
  which will then be available to use as meta-parameters for automation. */
  AudioPlugin(int numParameters = 10);
    // todo: change the default value to 0 - but only after updating Elan's projects to pass in the
    // required number


  virtual ~AudioPlugin();

  //-----------------------------------------------------------------------------------------------
  // setup:

  //void configureInsAndOuts();

  void setPluginName(const String& newName) { plugInName = newName; }

  /** Sets the jura::AudioModule object that is wrapped into an AudioPlugin. */
  virtual void setAudioModuleToWrap(AudioModule* moduleToWrap);

  /** Automatically attaches each of the wrapped AudioModule's parameters to one of our (meta)
  parameters here (up to the minimum of the number of internal and meta parameters.  */
  void autoAttachMetaParameters();

  /** Sets min/max values for the plugin editor width and height. */
  void setEditorSizeLimits(int minWidth, int minHeight, int maxWidth, int maxHeight)
  {
    editorWidthMin  = minWidth;
    editorHeightMin = minHeight;
    editorWidthMax  = maxWidth;
    editorHeightMax = maxHeight;
  }
  // todo: perhaps remove/hide the resizer widget, if min == max for w and h
  // see juce::AudoProcessorEditor::setResizable 
  // https://www.juce.com/doc/classAudioProcessorEditor#a3d36f7385146270fc752ce17418f115a

  //void setEditorBoundsConstrainer(ComponentBoundsConstrainer* newConstrainer)
  //{
  //  editorBoundsConstrainer = newConstrainer;
  //}

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
  virtual const String getProgramName (int index) override { return String(); }
  virtual void changeProgramName(int index, const String& newName) override {}
  virtual void getStateInformation(juce::MemoryBlock& destData) override;
  virtual void setStateInformation(const void* data, int sizeInBytes) override;
  virtual void processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages) override;

  // override for MetaParameterManagerObserver baseclass:
  virtual void metaNameChanged(MetaParameterManager* manager, int index) override
  {
    updateHostDisplay(); // does not yet work
  }

  //-----------------------------------------------------------------------------------------------
  // optional overrides for juce::AudioProcessor baseclass:

  virtual bool supportsDoublePrecisionProcessing() const override { return true; }
  virtual void processBlock(AudioBuffer<double>& buffer, MidiBuffer& midiMessages) override;
  virtual bool isBusesLayoutSupported(const BusesLayout&) const override;
  //virtual bool setPreferredBusArrangement(bool isInput, int bus,
  //  const AudioChannelSet& preferredSet) override;
  //virtual int getBusCount(bool isInput)	const { return 2; } // preliminary
  //virtual String getParameterName(int parameterIndex, int maximumStringLength) override;

  //-----------------------------------------------------------------------------------------------
  // Internals (ToDo: try to move all into the protected section):

  AudioModule *wrappedAudioModule = nullptr;  // the wrapped jura::AudioModule

  /** Mutex-lock for all accesses to the underlyingAudioModule's member functions - a pointer to
  the lock is  passed to the embedded AudioModule and should be used there also and the AudioModule
  should also pass this lock on to the GUI Editors. */
  CriticalSection plugInLock;  // maybe we should use getCallbackLock()


  MetaParameterManager metaParaManager;
  // Needs to be public to be accessible for the AudioModule wrapper functions. 


protected:

  rsSmoothingManager smoothingManager;

  /** Creates the parameters that are reported to the host. Called internally from the
  constructor. */
  void createHostAutomatableParameters(int numParameters);

  /** Sets the flush-to-zero (FTZ) and denormals-are-zero (DAZ) mode. Called from processBlock. */
  inline void enableFitzdazzing();

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
  std::vector<AudioPluginParameter*> parameters;

  juce::String plugInName;  // assign this in the constructor of your subclass
   // maybe get rid of this and let the wrapper return the name of the wrapped AudioModule

  int editorWidth     = 0;
  int editorHeight    = 0;

  //ComponentBoundsConstrainer* editorBoundsConstrainer = nullptr;
  int editorWidthMin  = 1;
  int editorHeightMin = 1;
  int editorWidthMax  = 100000;
  int editorHeightMax = 100000;
  //int editorWidthMax  = INT_MAX;  // breaks with juce 6.0.8, maps to negative value - why?
  //int editorHeightMax = INT_MAX;
  // Maybe replace these with a pointer to a juce::ComponentBoundsConstrainer that is initially a 
  // nullptr an which client code can set via a setEditorBoundsConstrainer function that can be 
  // called in the createPluginFilter function. This object should delete the passed object on
  // destruction ...i tried it - it didn't work



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

  AudioPluginWithMidiIn(int numParameters = 10) : AudioPlugin(numParameters) {}
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

  AudioPluginEditor(AudioModuleEditor* editorToWrap, AudioPlugin* pluginToEdit);

  ~AudioPluginEditor()
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

  /** Returns a pointer to the repaint manager object that is used to trigger periodic calls to 
  repaint for GUI animation. */
  rsRepaintManager* getRepaintManager() { return &repaintManager; }

protected:

  AudioModuleEditor* wrappedEditor;

  AudioPlugin* pluginToEdit;

  rsRepaintManager repaintManager;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioPluginEditor)
};

//=================================================================================================

/** Below are some templates that make it very easy to wrap any subclass of jura::AudioModule into 
a juce::AudioProcessor (actually, a jura::AudioPlugin which is a subclass of juce::AudioProcessor) 
and compile it as plugin. Your plugin project will just need to implement a function 
createPluginFilter(), declare a nullptr of the type of the respective subclass of jura::AudioModule 
and then invoke the appropriate template with that nullptr (there are different templates for 
plugins with and without midi, because we need different wrapper classes). For example, in your 
plugin project, all you need is the following code:

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
jura::Ladder *dummy = nullptr; return createPluginWithoutMidi(dummy);
}

and that will wrap a jura::Ladder module into a plugin (in this case, a plugin without midi 
input). */

template<class AudioModuleType>
AudioPlugin* JUCE_CALLTYPE createPluginWithoutMidi(AudioModuleType *dummy, int numParameters = 10)
{
  // wraps audio module into plugin without midi input
  jura::AudioPlugin *plugIn = new jura::AudioPlugin(numParameters);
  AudioModuleType   *module = new AudioModuleType(&plugIn->plugInLock, &plugIn->metaParaManager);
  module->setSaveAndRecallMetaParameters(true);
  plugIn->setAudioModuleToWrap(module);
  return plugIn;
}

template<class AudioModuleType>
AudioPluginWithMidiIn* JUCE_CALLTYPE createPluginWithMidi(AudioModuleType *dummy, 
  int numParameters = 10)
{
  // wraps audio module into plugin with midi input
  jura::AudioPluginWithMidiIn *plugIn = new jura::AudioPluginWithMidiIn(numParameters);
  AudioModuleType *module = new AudioModuleType(&plugIn->plugInLock, &plugIn->metaParaManager);
  module->setSaveAndRecallMetaParameters(true);
  plugIn->setAudioModuleToWrap(module);
  return plugIn;
  // only the 1st line is different, the others are duplicated -> factor out
}

//// The code below was supposed to dispatch between the createPluginWithoutMidi and 
//// createPluginWithMidi version of the plugin creation code above. Unfortunately, it doesn't work 
//// that way, so we must manually call the appropriate version in out our actual plugin project. 
//// Maybe at some point late, we'll find a solution...
//template<class AudioModuleType>
//AudioProcessor* JUCE_CALLTYPE createPlugin(AudioModuleType *dummy, bool withMidiIn)
//{
//  // dispatcher between the different wrappers above
//  if( withMidiIn == true )
//    return createPluginWithMidi(dummy);
//  else
//    return createPluginWithoutMidi(dummy);
//
//  // todo:
//  // Maybe, we can somehow infer, whether or not we need to create a plugin with or without midi 
//  // from the type of the passed pointer to get rid of the boolean flag "withMidiIn". But 
//  // dynamic_cast and checking agianst a nullptr doesn't work when we actually pass a nullptr 
//  // (which is what we do). Maybe somethign using decltype or type-traits or some other C++11 
//  // feature may help - we'll see.
//}

#endif
