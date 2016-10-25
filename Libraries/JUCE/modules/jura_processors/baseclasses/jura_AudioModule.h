#ifndef jura_AudioModule_h
#define jura_AudioModule_h

class AudioModule;

/** A class that can be informed (via a callback method) when an AudioModule object is going to be 
deleted. Mainly intended as baseclass for GUI elements that keep a pointer to an AudioModule that 
is being edited 

hmm.. - maybe, we don't need that. The documentation of the 
AudioProcessorEditor* AudioProcessor::createEditor() function says:
"It's safe to assume that an editor will be deleted before its filter."

*/

class JUCE_API AudioModuleDeletionWatcher
{

public:

  /** Destructor. */
  virtual ~AudioModuleDeletionWatcher();

  /** Callback function that subclasses should override in order to invalidate their 
  references/pointers to the Audiomodule in question. */
  virtual void audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted) = 0;

  /** Registers ourselves as deletion-watcher with the passed AudioModule such that we will get 
  callbacks to audioModuleWillBeDeleted when the module in question is going to be deleted. */
  virtual void addWatchedAudioModule(AudioModule *moduleToBeWatched);

  /** De-registers ourselves from and AudioModule to which we presumably had previously registered
  via addWatchedAudioModule - if we didn't, this function will do nothing. */
  virtual void removeWatchedAudioModule(AudioModule *moduleNotToBeWatchedAnymore);

protected:

  juce::Array<AudioModule*> watchedAudioModules;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class is the base class for all audio modules. rename to RAudioProcessor or AudioUnit.*/

class JUCE_API AudioModule : public AudioProcessor, public AutomatableModule, 
  public StateFileManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AudioModule();

  /** Destructor. */
  virtual ~AudioModule();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Override this to set the sample-rate for this module. */
  virtual void setSampleRate(double newSampleRate);

  /** Override this to set the tempo in bpm for synced modules. */
  virtual void setBeatsPerMinute(double newBpm);

  /** Sets up the name for this AudioModule. */
  virtual void setModuleName(const juce::String& newName);

  /** Sets up an appendix (like "Demo Version") for the name of this AudioModule. */
  virtual void setModuleNameAppendix(const juce::String& newAppendix);

  /** Adds a child AudioModule to this one. */
  virtual void addChildAudioModule(AudioModule* moduleToAdd);

  /** Removes a child AudioModule from this one and optionally deletes the object (it will only 
  delete it if it is actually in the array of childAudioModules). */
  virtual void removeChildAudioModule(AudioModule* moduleToRemove, bool deleteOject);

  /** Use this function to turn on/off save/recall of the state - the idea is that some module may 
  have child modules which are currently inactive and therefore don't need to save and recall their
  state. This makes preset files more economical. */
  virtual void setStateSaveAndRecall(bool shouldSaveAndRecall)
  {
    saveAndRecallState = shouldSaveAndRecall;
  }

  /** Checks, if this is a cracked version and if so, it sets up the appendix for the headline 
  accordingly. Return value informs also whether or not a cracked version was detected. */
  virtual bool checkForCrack();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the name of this module. */
  virtual juce::String getModuleName() const { return moduleName; }

  /** When we have several child-modules with the same name (member "moduleName"), this function 
  can be used to find the index of the passed child-module among them. It will return 0 when the 
  passed AudioModule is the first (or only one) with that name, 1 for the second and so on. When 
  the passed module is not one of our child modules, it will return -1. */
  virtual int getIndexAmongNameSakes(AudioModule *child);

  /** Returns a string that is to be used as headline for a GUI-editor for this module - this will 
  be typically the name of the module, perhaps appended by some qualifier like "Demo Version", 
  "Cracked By bLaH" or something. */
  virtual juce::String getModuleHeadlineString();

  /** Returns the interval at which the module wants to receive callbacks to trigger(). */
  virtual double getTriggerInterval() const { return triggerInterval; }

  /** Returns, whether this mdoule wants to save/recall its state - the idea is that some module 
  may have child modules which are  currently inactive and therefore don't need to save and recall 
  their state. This makes preset files more economical. */
  bool wantsSaveAndRecallState() const { return saveAndRecallState; }

  //-----------------------------------------------------------------------------------------------
  // automation and state management:

  /** Callback to indicate that a parameter has changed - subclasses should override this and
  update their signal processing accordingly. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  /** Calls a parameterChanged for each of the observed parameters - this should trigger the
  appropriate updating of the signal processing core in the subclasses. */
  virtual void updateCoreObjectAccordingToParameters();

  /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean);

  /** Returns the state (i.e. the settings of all relevant parameters) in form of an XmlElement. */
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  /** Converts a state which might possibly be from an older version to the current patch-format. 
  The baseclass implementation just returns the state as is, but will trigger a debug-break if the 
  patchFormatIndex of the state and the module don't match. Override this function in your subclass
  to do the actual conversion. */
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState);

  //-----------------------------------------------------------------------------------------------
  // Audio processing:

  /** This is the audio callback that your subclass needs to override. */
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) = 0;

  //-----------------------------------------------------------------------------------------------
  // Event processing - move into subclass AudioModuleWithMidi - there, we need to re-implement
  // processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages) callback in order to
  // actually do something with the passed MidiBuffer there - from there, we should call the 
  // individual event-handlers which can int turn be overriden by subclasses of 
  // AudioModuleWithMidi...

  /** Handles a generic MidiMessage. */
  virtual void handleMidiMessage(MidiMessage message);

  /** Triggers a note-on event. */
  virtual void noteOn(int noteNumber, int velocity);

  /** Triggers a note-off event. */
  virtual void noteOff(int noteNumber);

  /** Triggers an all-notes-off event. */
  virtual void allNotesOff();

  /** Overrides setMidiController which is inherited from both base-classes - and we simply we pass
  through the function call to both of them here. */
  virtual void setMidiController(int controllerNumber, int controllerValue);

  /** Triggers a pitch-bend event. */
  virtual void setPitchBend(int pitchBendValue);

  /** Override this and set triggerInterval to some nonzero value if you need to re-trigger 
  something at regular intervals (like LFOs, for example). This function will be called from the 
  process-function at the given triggerInterval (if nonzero) - this value has to be specified in 
  beats. */
  virtual void trigger() {}

  /** Override this to reset this module (audio buffers and such). */
  virtual void reset() {}

  /** Override this to reset the state of this module to defaults (user parameters). */
  virtual void setStateToDefaults() {}

  /** Flag to indicate that this module needs tempo sync information (current BPM). */
  bool wantsTempoSyncInfo;

  // These should all be factored out into a subclass AudioModuleWithMidi -  in this subclass
  // we override the processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages) in order 
  // to actually use the incoming midi data. In this override we take care of sample accurate 
  // handling of midi events...

  //-----------------------------------------------------------------------------------------------
  // mandatory overrides for juce::AudioProcessor baseclass:

  virtual const String getName() const override { return "AudioModule"; }
  virtual void  prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) override;
  virtual void releaseResources() override {}
  virtual double getTailLengthSeconds() const override { return 0.0; }
  virtual bool acceptsMidi()  const override { return false; }
  virtual bool producesMidi() const override { return false; }
  virtual bool hasEditor() const override { return true; }
  virtual AudioProcessorEditor* createEditor() override { return nullptr; } // override in your subclass !!
  virtual int getNumPrograms() override { return 1; }                       // 1, because 0 is not allowed
  virtual int getCurrentProgram() override { return 0; }
  virtual void setCurrentProgram(int index) override {}
  virtual const String getProgramName (int index) override { return String::empty; }
  virtual void changeProgramName(int index, const String& newName) override {}
  virtual void getStateInformation(juce::MemoryBlock& destData) override;
  virtual void setStateInformation(const void* data, int sizeInBytes); 
  virtual void processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages) override;

  // optional AudioProcessor overrides:
  virtual void processBlock(AudioBuffer<double> &buffer, MidiBuffer &midiMessages) override;
  virtual bool supportsDoublePrecisionProcessing() const override { return true; }

  virtual bool setPreferredBusArrangement(bool isInput, int bus,
    const AudioChannelSet& preferredSet) override;

protected:

  /** Must be overriden by subclasses to fill the inherited array of observed parameters. */
  virtual void initializeAutomatableParameters();

  /** Our child modules to which we will distribute MIDI-events and of which we manage the
  states. */
  juce::Array<AudioModule*, CriticalSection> childModules;

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

  CriticalSection *plugInLock;     // mutex to access the wrapped core dsp object - make this a nom-pointer member
  double triggerInterval;          // interval (in beats) for calls to trigger()
  bool saveAndRecallState;         // indicates, that this module wants to save/recall its state
  int patchFormatIndex;            // version of patch format (for backwards compatibility)

  juce::String moduleName;         // name of this AudioModule
  juce::String moduleNameAppendix; // string to be appended to the name on the GUI (such as 
                                   // Demo-Version, etc.) - remove (or factor into some subclass)

  friend class AudioModuleEditor;  // the editor must access our plugInLock member

private:

  /** Registers an AudioModuleDeletionWatcher that will be called back when this object is 
  deleted. */
  virtual void registerDeletionWatcher(AudioModuleDeletionWatcher *watcher);

  /** De-registers a previously registered AudioModuleDeletionWatcher. */
  virtual void deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher);

  juce::Array<AudioModuleDeletionWatcher*> deletionWatchers;
  friend class AudioModuleDeletionWatcher;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

// here goes the subclass AudioModuleWithMidi....



#endif 