#ifndef jura_AudioModule_h
#define jura_AudioModule_h

class AudioModule;
class AudioModuleEditor;

/** A class that can be informed (via a callback method) when an AudioModule object is going to be 
deleted. Mainly intended as baseclass for GUI elements that keep a pointer to an AudioModule that 
is being edited 

hmm.. - maybe, we don't need that. The documentation of the 
AudioProcessorEditor* AudioProcessor::createEditor() function says:
"It's safe to assume that an editor will be deleted before its filter." 

\todo: maybe move these classes down into jura_framework - some plugins may want to use these 
baseclasses without needing any of the dsp-stuff here - for example PrettyScope wouldn't need a 
dependency on jura_processors anymore

*/

class JUCE_API AudioModuleDeletionWatcher
{

public:

  /** Constructor. */
  AudioModuleDeletionWatcher() {}

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

  juce::Array<AudioModule*> watchedAudioModules; // use std::vector

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleDeletionWatcher)
};

//=================================================================================================

/** This class is the base class for all audio modules. 
\todo: maybe move this class into the jura_framework juce-module
*/

class JUCE_API AudioModule : public ParameterManager, public StateFileManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. The caller must pass a pointer to a critical section object that must exist 
  somewhere outside and will be used here to acquire mutually exclusive access from different 
  threads to the underlying dsp object. For example, if the AudioModule is wrapped into a plugin,
  the Criticalsection object could exist on the level of the plugin (an its lifetime be managed 
  there). To enable meta-control, you can also pass
  
  */
  AudioModule(CriticalSection* lockToUse, MetaParameterManager* metaManagerToUse = nullptr);

  /** Destructor. */
  virtual ~AudioModule();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Override this to set the sample-rate for this module. */
  virtual void setSampleRate(double newSampleRate);

  /** Override this to set the tempo in bpm for synced modules. */
  virtual void setBeatsPerMinute(double newBpm);

  /** Sets up the name for this AudioModule. */
  virtual void setModuleName(const juce::String& newName);

  /** Sets up the type name of the module and optionally updates the currently active directory for
  loading and saving presets (which is typically based on that type name). It may also optionally 
  set the module name to the same name. */
  virtual void setModuleTypeName(const juce::String& newName, bool updatePresetDirectory = true, 
    bool setModuleNameAlso = true);

  /** Sets up an appendix (like "Demo Version") for the name of this AudioModule. */
  virtual void setModuleNameAppendix(const juce::String& newAppendix);

  /** Adds a child AudioModule to this one. This module will take ownership over the pased module, 
  i.e. delete it in the destructor (maybe rename to addOwnedChildAudioModule to convey that 
  information at the call site). */
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
  // why virtual? inline...

  /** Call this function, if you want this module to store and recall the *values* of the 
  MetaParameters in our MetaParameterManager member. Initially, this is false for all AudioModules
  but should be activated for top-level modules, i.e. those that are wrapped into an AudioPlugin.
  The function is called in the createPlugin... functions. */
  inline void setSaveAndRecallMetaParameters(bool shouldSaveAndRecall) 
  { 
    saveAndRecallMetas = shouldSaveAndRecall; 
  }

  /** Loads a preset from the given path. The path is assumed to be relative to the preset 
  directory, i.e. the directory that is returned by getPresetDirectory(). */
  virtual void loadPreset(const juce::String& pathFromPresetFolder);

  /** Tries to load the default preset from the location that your subclass can specify by 
  overriding getDefaultPresetLocation (if you don't oveeride this, this function will have no 
  effect). It is called in the Chainer immediately after the MetaParameterManger has been 
  assigned. */
  //virtual void loadDefaultPreset();

  /** Subclasses can override this in order to return a location of a default preset which
  should be loaded when the module is first plugged in. It will be called inside our 
  loadDefaultPreset function, which will try to load the preset from the returned location. */
  //virtual juce::String getDefaultPresetLocation() { return juce::String(); }

  /** Checks, if this is a cracked version and if so, it sets up the appendix for the headline 
  accordingly. Return value informs also whether or not a cracked version was detected. */
  virtual bool checkForCrack();
    // move to another class or get rid

  /** Overrides inherited method to additionaly wire the passed Parameter up to the 
  MetaParameterManager. */
  virtual void addObservedParameter(Parameter *parameterToAdd) override;

  /** Sets up the SmoothingManager and MetaParameterManager for the passed parameter, if 
  applicable (i.e. the parameter is a subclass that supports these things). */
  virtual void setupManagers(Parameter* p);

  //-----------------------------------------------------------------------------------------------
  // \name midi/meta stuff 
  // (midi should go to subclass AudioModuleWithMidiIn):

  /** Assigns a MIDI controller to one of the observed parameters. */
  virtual void assignMidiController(const juce::String& nameOfParameter, int controllerNumber);

  /** Receives MIDI controller messages and dispatches them to the appropriate to one of the
  Parameter object. */
  virtual void setMidiController(int controllerNumber, float controllerValue);

  /** Reverts all observed parameters to their default settings. */
  virtual void revertToDefaultMapping();
    // rename to revertToDefaultMidiMapping or assignDefaultMidiController, move to protected
    // don't acquire lock

  /** Detaches all our parameters from their controlling MetaParameter. */
  virtual void detachMetaParameters();

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry:

  /** Returns the name of this module, for example "Filter1". */
  virtual juce::String getModuleName() const { return moduleName; }

  /** Returns the appendix that should be added to the module name for display on the GUI, such as
  "Demo Version" or whatever (mostly, this is empty). */
  inline juce::String getModuleNameAppendix() const { return moduleNameAppendix; }

  /** Returns the type name of this module, for example "LadderFilter". */
  virtual juce::String getModuleTypeName() const { return moduleTypeName; }

  /** Returns the directory, where the support files such as presets, samples, etc. are supposed to 
  be found. */
  //virtual juce::String getSupportDirectory() const;

  /** Returns the directory, where the presets are supposed to be found. You can switch between 
  factory-presets and user-presets via the "user" parameter (in the baseclass implementation, the 
  same folder will be returned in both cases, but you may want to override this). */
  virtual juce::String getPresetDirectory(bool user = false) const;

  /** When we have several child-modules with the same name (member "moduleName"), this function 
  can be used to find the index of the passed child-module among them. It will return 0 when the 
  passed AudioModule is the first (or only one) with that name, 1 for the second and so on. When 
  the passed module is not one of our child modules, it will return -1. */
  virtual int getIndexAmongNameSakes(AudioModule *child);

  /** Returns a string that is to be used as headline for a GUI-editor for this module - this will 
  be typically the name of the module, perhaps appended by some qualifier like "Demo Version". */
  virtual juce::String getModuleHeadlineString();

  /** Returns the interval at which the module wants to receive callbacks to trigger(). */
  virtual double getTriggerInterval() const { return triggerInterval; }

  /** Returns, whether this mdoule wants to save/recall its state - the idea is that some module 
  may have child modules which are  currently inactive and therefore don't need to save and recall 
  their state. This makes preset files more economical. */
  bool wantsSaveAndRecallState() const { return saveAndRecallState; }

  /** Returns a pointer to the MetaParameterManager that will be used by 
  MetaControlledParameters. */
  MetaParameterManager* getMetaParameterManager() const { return metaParamManager; }

  /** Returns true, if this is a top level AudioModule (i.e. it has no parent module). */
  bool isTopLevelModule() { return parentModule == nullptr; }

  /** Returns the module that is at the top of the parent/children hierarchy, i.e. the outermost 
  AudioModule. */
  AudioModule* getTopLevelModule();

  /** Returns a string that gives the "path" of the AudioModule, for example a string like
  Straightliner.OscSection.Osc2. would be returned from the 2nd oscillator of Straightliner. */
  juce::String getAudioModulePath();
   // todo: get rid of the final dot - i think, it is used to uniquely identify parameters in xml
   // files (for mod-connections, meta-assignements, etc.) so after the dot, there will be the 
   // parameter name

  /** Your subclass may override this to return an object of an appropriate subclass of
  AudioModuleEditor. The baseclass implementation will return a generic editor with sliders, 
  comboboxes and button for all the Parameters of this AudioModule. You can create the editor with
  "new", return it and forget about it - the responsibility for deleting it eventually is on the 
  framework. Via the integer "type" parameter, it is possible to have different kinds of editors
  for the same kind of audio module - if your AudioModule subclass doesn't need to support 
  different editors, just ignore that parameter. */
  virtual AudioModuleEditor* createEditor(int type = 0);

  //-----------------------------------------------------------------------------------------------
  // \name Automation and state management:

  /** Callback to indicate that a parameter has changed - subclasses should override this and
  update their signal processing accordingly. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  /** Calls a parameterChanged for each of the observed parameters - this should trigger the
  appropriate updating of the signal processing core in the subclasses. */
  virtual void updateCoreObjectAccordingToParameters();

  /** Calls the associated callback function for each of our observed Parameters, optionally also
  recursively for all child AudioModules. */
  void callParameterCallbacks(bool recursivelyForChildModules);

  /** Calls the parameterChanged callback for each of our observed parameters, such that each 
  parameter nofies each of its obervers about a parameter update. */
  void notifyParameterObservers(bool recursivelyForChildModules);

  /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

  /** Returns the state (i.e. the settings of all relevant parameters) in form of an XmlElement. */
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  /** Converts a state which might possibly be from an older version to the current patch-format. 
  The baseclass implementation just returns the state as is, but will trigger a debug-break if the 
  patchFormatIndex of the state and the module don't match. Override this function in your subclass
  to do the actual conversion. */
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState);
  // it is a bit inefficient to return the new xml - it has to be created from the old - maybe use
  // XmlElement* convertXmlStateIfNecessary(const XmlElement* oldXml) that just returns the pointer
  // if there's nothing to do. it implies that we have to change setStateFromXml to take a pointer 
  // too - maybe keep both versions - but maybe that's not true since c++11 becuase of return-value 
  // optimization?

  /** Resets all the parameters to their default values. */
  virtual void resetParametersToDefaultValues();

  /** Sets up the MetaParameterManager tha will be used by MetaControlledParameters. This object
  is typically member of some outlying AudioPlugin and the AudioPlugin will set this up in its 
  constructor. */
  void setMetaParameterManager(MetaParameterManager* managerToUse);
    // maybe remove - it should be passed to the constructor

  /** sets up the smoothing manager used for parameter smoothing. @see setMetaParameterManager. */
  void setSmoothingManager(rsSmoothingManager* managerToUse);

  //-----------------------------------------------------------------------------------------------
  // \name Audio processing:

  /** This is the audio callback that your subclass needs to override. */
  //virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) = 0;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) 
  {
    // We provide a baseclass implementation for the important special case of 2 channels here, so
    // subclasses do not necessarily override this (they can get away with just overriding 
    // processStereoFrame). But for optimization purposes, you may want to override this function,
    // anyway.
    if(numChannels == 2)
      for(int n = 0; n < numSamples; n++)
        processStereoFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
    else
      jassertfalse;

  }

  /** A function to override in case you are producing a stereo-signal and need your module to work
  with the modulation system (i.e. uses ModulatableParameters). It is supposed to produce one 
  stereo sample frame at a time and may be called alternatively to processBlock. */
  virtual void processStereoFrame(double *left, double *right) {}

  //-----------------------------------------------------------------------------------------------
  // \name Misc:

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
  bool wantsTempoSyncInfo = true;  // why public?

protected:

  // used internally as sub-tasks in the get/setState functions:
  void parameterToXml(   XmlElement* xmlElementToStartFrom, Parameter* p);
  void parametersToXml(  XmlElement* xmlElementToStartFrom);
  void midiMappingToXml( XmlElement* xmlElementToStartFrom);
  void metaMappingToXml( XmlElement* xmlElementToStartFrom);
  void metaValuesToXml(  XmlElement* xmlElementToStartFrom);
  void childModulesToXml(XmlElement* xmlElementToStartFrom);

  void recallParametersFromXml(  const XmlElement &xmlState);
  void recallChildModulesFromXml(const XmlElement &xmlState, bool markAsClean);
  void recallMidiMappingFromXml( const XmlElement &xmlState);
  void recallMetaMappingFromXml( const XmlElement &xmlState);
  void recallMetaValuesFromXml(  const XmlElement &xmlState);

  /** Our child modules to which we will distribute MIDI-events and of which we manage the
  states. */
  std::vector<AudioModule*> childModules;  // maybe rename to childAudioModules
  AudioModule* parentModule = nullptr; // new - needed for tree traversal from leaves

  rsSmoothingManager* smoothingManager = nullptr;
  MetaParameterManager* metaParamManager = nullptr;

  // ToDo:
  // double dryGain = 0.0;
  // double wetGain = 1.0;
  // to mix between input and output i.e. dry/wet. the idea is that source may modules want to add
  // their output to what's already there (for stacking, layering, parallel connection -> both 
  // gains = 1) whereas filters and effects may want to replace the input with the output (serial 
  // connection: dry = 0, wet = 1). Having these two gains provides a flexible way to facilitate 
  // this and more. These should also be saved and recalled, of course (but saved only when not at 
  // default settings)

  double triggerInterval = 0.0;    // interval (in beats) for calls to trigger()
  bool saveAndRecallState = true;  // indicates, that this module wants to save/recall its state
  bool saveAndRecallMetas = false; // store/recall values of meta-params - should be true for top-level modules only
  int patchFormatIndex = 1;        // version of patch format (for backwards compatibility)


  juce::String moduleTypeName;     // type name like "LadderFilter"
  juce::String moduleName;         // name of this AudioModule, like "Filter1"
  juce::String moduleNameAppendix; // string to be appended to the name on the GUI (such as 
                                   // Demo-Version, etc.) todo: factor into some subclass - only 
                                   // top-level modules actually use this

  friend class AudioModuleEditor;  // the editor must access our plugInLock member

private:

  /** Registers an AudioModuleDeletionWatcher that will be called back when this object is 
  deleted. */
  virtual void registerDeletionWatcher(AudioModuleDeletionWatcher *watcher);

  /** De-registers a previously registered AudioModuleDeletionWatcher. */
  virtual void deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher);

  juce::Array<AudioModuleDeletionWatcher*> deletionWatchers;
  friend class AudioModuleDeletionWatcher;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModule)
};

//=================================================================================================

/** A subclass of AudioModule that can handle modulatable parameters. 
\todo maybe factor this into the regular AudioModule class. */

class JUCE_API ModulatableAudioModule : public AudioModule, public ModulationParticipant
{

public:

  ModulatableAudioModule(CriticalSection* lockToUse, 
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr) 
    : AudioModule(lockToUse, metaManagerToUse), ModulationParticipant(modManagerToUse)
  {
    if(modManager != nullptr)
      modManager->setMetaParameterManager(metaManagerToUse);
  }

  virtual ~ModulatableAudioModule() {}

  /** Overrides inherited method to additionaly wire the passed Parameter up to the 
  ModulationManager, in case it is a ModulatableParameter. */
  virtual void addObservedParameter(Parameter* parameterToAdd) override;

  /** Overriden to additionally set up the ModulationManager. */
  virtual void setupManagers(Parameter* p) override;

  /** Overriden from AudioModule in order to pass the ModulationManger object "modManager" 
  inheridted from ModulationParticipant to the child-module to be added in case, the child is 
  itself a  ModulatableAudioModule. */
  virtual void addChildAudioModule(AudioModule* moduleToAdd) override;

  /** Overriden from ModulationParticipant in order set up the ModulationManager "modMan" in all
  our ModulatableParameters. */
  virtual void setModulationManager(ModulationManager* managerToUse) override;

protected:

  /** Checks, if the passed Parameter is actually a ModulatableParameter and if so, does the 
  required setup work such as registering it with the ModulationManager. */
  void setupForModulationIfModulatable(Parameter* p);

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableAudioModule)
};

//=================================================================================================

/** A subclass of AudioModule that accepts MIDI input. If you derive your effect or instrument from
this baseclass and wrap it into a juce::AudioProcessor (via the wrapper class jura::AudioPlugin),
the plugin will have a MIDI input. Also, you can override the event handler methods in your 
subclass (noteOn, noteOff, setMidiController, etc. - inherited from rsMidiMessageDispatcher) in 
order to respond to incoming MIDI events. You can also override handleMidiMessage and do the 
dispatch yourself. */

class JUCE_API AudioModuleWithMidiIn : public ModulatableAudioModule, 
  public rsMidiMessageDispatcher
{

public:

  //AudioModuleWithMidiIn(CriticalSection *lockToUse) : AudioModule(lockToUse) {}
  AudioModuleWithMidiIn(CriticalSection *lockToUse, 
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr) 
    : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse) {}

  virtual ~AudioModuleWithMidiIn() {}

  //-----------------------------------------------------------------------------------------------
  // \name Event processing:

  /** Handles a generic MidiMessage. */
  void handleMidiMessage(MidiMessage message) override;

  /** Overrides setMidiController */
  void setMidiController(int controllerNumber, float controllerValue) override;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleWithMidiIn)
};

//=================================================================================================

/** Baseclass for polyphonic AudioModules. They have polyphonic audio outputs, so you can chain a 
polyphonic oscillator with a polyphonic filter and drop in some polyphonic modulators in ToolChain.

...under construction....tbc.... */

class JUCE_API AudioModulePoly : public AudioModuleWithMidiIn 
{

public:

  /** Constructor. May set up the managers. */
  AudioModulePoly(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  virtual ~AudioModulePoly() { }

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the voice manager for this module and recursively for all the child modules. This will 
  trigger a call to allocateVoiceResources because the amount of resources needed depends on the
  maxNumVoices value which is a member of rsVoiceManager, so if the voiceManager changes, the
  resources (like states of dsp objects) may need to be (re)allocated. Typically, this just happens
  once, soon after construction. */
  void setVoiceManager(rsVoiceManager* managerToUse);
  // todo: maybe make it optional (true by default) to set it also for the child modules

  /** Overriden in order to set the voice manager for the child modules. */
  void addChildAudioModule(AudioModule* moduleToAdd) override;

  /** Sets the buffer into which the individual voice outputs shall be written. It must be of 
  length 2*maxNumVoice where the factor two comes from the two channels for stereo signals. */
  void setVoiceSignalBuffer(double* buffer) { voicesBuffer = buffer; }

  /** Switches the module into monophonic mode. The idea is that polyphonic modules should be 
  able to optionally behave exactly as if they would be implemented without making any use of the
  polyphony infrastructure. By default i.e. after construrction, an AudioModulePoly will actually 
  be in polyphonic mode, because it would be confusing otherwise (i guess). The feature is mostly 
  meant to facilitate upgrading monophonic module classes to polyphonic ones later without breaking
  any patches. Patches involvong the old monophonic modulators could use the upgraded polyphonic 
  versions the same way as before by just switching them into monophonic mode. */
  virtual void setMonophonic(bool shouldBeMonophonic);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  // int getMaxNumVoices();

  /** Returns a const pointer to our voiceManager which can be used by the caller to inquire 
  things about the current processing state like the number of active voices etc. */
  const rsVoiceManager* getVoiceManager() const { return voiceManager; }

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  /** Overriden to inquires from the voiceManager which voice was used for the most recently 
  triggered note and then calls noteOnForVoice with that. Finalized because subclasses are 
  supposed to override noteOn(int key, int vel, int voice) instead (which is called from this). */
  void noteOn(int key, int vel) override final;

  /** Overriden in order to invoke the per voice callback for the voice that is used for the new 
  note for those polyphonic parameters that are not connected to any modulation sources. For those 
  that are connected, the callbacks will be called in the audio process call by the modulation 
  system but for the disconnected ones, we must do it here. Finalized because subclasses are 
  supposed to override noteOn(int key, int vel, int voice) instead (which is called from this). */
  void noteOnForVoice(int key, int vel, int voice) override final;

  // todo: maybe finalize noteOffForVoice

  /** Can be overriden by subclasses to respond to note-on events for a particular voice. Can be 
  used in audio processors (oscillators, filters, etc.) to reset oscillator phases, filter states,
  etc. and in modulators (envelopes, etc.) to reset the position in the envelope etc. */
  virtual void noteOn(int key, int vel, int voice) {}
  // maybe make this purely virtual

  // todo:
  //virtual void noteOff(int key, int vel, int voice) {}

  /** Supposed to be called on note-on when the given voice should immediately switch to the new 
  pitch due to having just been grabbed from pool of idle voices or was stolen and must now play a 
  new pitch. It may also be called per sample, when the note pitch of a voice is in transition to 
  a new value due to glide. The pitch given is in units of midi-note numbers, but as continuous 
  value. Subclasses can override this in order to implement their required response, like - for 
  example - setting the frequency of an oscillator to the frequency that corresponds to the new 
  pitch. The pitch value already has microtuning baked in, but not pitch-wheel. This is because 
  microtuning is handled globally but the response to pitch-wheel is implemented via the modulation
  system. Pitch-Wheel is available as modulator and can be routed to the parameters at will. */
  virtual void setVoiceKeyPitch(int voice, double pitch) {}
  // obsolete? should this not better be handled by a modulator module? maybe the NotePitch/Freq
  // modulators should have glide already baked in? i think that would make sense

  virtual void setVoiceVelNorm(int voice, double pitch) {}
  // similar to setVoiceKeyPitch but for (normalized) velocity 

  //-----------------------------------------------------------------------------------------------
  // \name Audio processing:

  /** Must be overriden by subclasses to produce a stereo sample frame at a time. */
  virtual void processStereoFrameVoice(double *left, double *right, int voice) = 0;

  /** Overriden in order to first compute all the voice outputs via calling 
  processStereoFramePoly and writing them into our voicesBuffer and them summing them into 
  left/right. */
  void processStereoFrame(double *left, double *right) override;

  /** Fills the given buffer with the stereo outputs of the voices. It just calls 
  processStereoFrameVoice in a loop over the active voices using the buffer for the left/right
  pointer parameters. Stereo signals are stored interleaved - buffer[0], buffer[1] is the pair for
  the 1st voice etc. */
  virtual void processStereoFramePoly(double *buffer, int numActiveVoices);

  // maybe add later for optimization:
  //virtual void processBlockVoice(double **inOutBuffer, int numChannels, int numSamples, int voice) { }
  //virtual void processBlockPoly(double ***inOutBuffer, int numChannels, int numSamples) {}
  //virtual void processStereoFramePoly(double **left, double **right) {}


  //-----------------------------------------------------------------------------------------------
  // \name Misc:

  /** Your subclass should override this to allocate the DSP resources for the voices. When this 
  function gets called, it is assumed that the voiceManager has already been assigned (in fact, it 
  is getting called from setVoiceManager), so in your overriden implementation, you can access the 
  voiceMananger member to inquire things the desired maximum number of voices which will be 
  relevant for knowing how many core DSP objects you need. */
  virtual void allocateVoiceResources() = 0;
  // Maybe it should have a boolean parameter to decide whether it should get called recursively on
  // the child modules. In our setVoiceManager we would pass false because the function calls 
  // itself recursively on the child modules and also calls allocateVoiceResources, so we would end
  // up re-allocating twice. But in other circumstances it may make more sense to call it on the 
  // child modules, so we may need both options.

protected:


  rsVoiceManager* voiceManager = nullptr;

  bool monophonic = false;

  double *voicesBuffer = nullptr; // should be of length 2*maxNumVoices (2 for the 2 channels)
  // hmm - i'm not sure, if that's a good design - the motivation for introducing this buffer is
  // not having to change ToolChain::processBlock such that the realtime code of mono and poly
  // modules can be handled uniformly there...and that requires that we need to write our sample
  // data somewhere else - we cant use the buffers that are passed to processBlock because they
  // are only for monophonic signals. maybe we should use generic processFrame(double* buffer), 
  // processBlock(double* buffer) functions... or introduce processFramePoly already in the
  // AudioModule baseclass and always call that
  // ...or maybe we should actually own the voicesBuffer here - re-using the same buffer for 
  // various modules may be a premature optimization.

};

//=================================================================================================

/** Baseclass for monophonic modulator modules. */

class JUCE_API ModulatorModuleMono : public AudioModuleWithMidiIn, public ModulationSource
{

public:

  using AudioModuleWithMidiIn::AudioModuleWithMidiIn;  // inherit constructors

};

//=================================================================================================

/** Baseclass for polyphonic modulator modules. */

class JUCE_API ModulatorModulePoly : public AudioModulePoly, public ModulationSourcePoly
{
  // ToDo: Maybe experiment with the order of the multiple inheritance. This affects the memory 
  // layout of the class and may therefore affect performance. It should not affect functionality.

public:

  ModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr)
    : AudioModulePoly(lockToUse, metaManagerToUse, modManagerToUse) {}

  /** Must be overriden by subclasses to allocate the DSP sources. It serves the same purpose as 
  AudioModulePoly::allocateVoiceModResources but for modulators. We give it a different name in 
  order to make it purely virtual because we have already overriden the allocateVoiceResources
  here, so we need a different name in order to make it purely virtual again. We do this because
  ModulatorModulePoly needs to always allocate the slots for the voice output signals and we do 
  not want to force subclasses to call the basclass method - that would just add more boilerplate
  and could be forgotten as well. */
  virtual void allocateVoiceModResources() = 0;

  /** Sets the modulator into monophonic mode in which case it should behave in the same way as 
  truly monophonic object that doesn't use the polyphony system would.  */
  void setMonophonic(bool shouldBeMonophonic) override
  {
    AudioModulePoly::setMonophonic(shouldBeMonophonic);
    ModulationSourcePoly::setMonophonic(shouldBeMonophonic);
    // Thats a bit ugly: we inherit a "monophonic" flag and a setter for it from both baseclasses.
    // ToDo: try to find a better design that avoids this
  }

  // Finalize all the audio processing callbacks - they are not suppsoed not be used for anything 
  // anymore and just pass the audio through:
  void processBlock(double** inOutBuffer, int numChannels, int numSamples) override final {}
  void processStereoFrame(double* left, double* right) override final {}
  void processStereoFrameVoice(double* left, double* right, int voice) override final {}
  void processStereoFramePoly(double* buffer, int numActiveVoices) override final {}

  /** Overriden to trigger the allocation of the modulator output signal slots and triggers 
  allocateVoiceModResources which must be overriden by subclasses. */
  void allocateVoiceResources() override final
  {
    ModulationSourcePoly::allocateVoiceOutputPins(voiceManager);
    allocateVoiceModResources();
  }

};

//=================================================================================================

/** Baseclass for GUI editors for AudioModule objects. */

class AudioModuleEditor : public jura::Editor, public ChangeListener, public RDialogBoxListener, 
  public RButtonListener /*, public KeyListener*/
{

public:

  enum positions
  {
    INVISIBLE,
    RIGHT_TO_HEADLINE,
    BELOW_HEADLINE,
    RIGHT_TO_INFOLINE
  };

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  AudioModuleEditor(AudioModule* newModuleToEdit);

  /** Sometimes, the module to edit is irrelevant or unknown to the editor, but we still need to
  access the mutex-lock of the module. In this case, this constructor can be used. */
  AudioModuleEditor(CriticalSection* pluginLockToUse);

  /** Initialization function that is called from the constructors (factored out to consolidate
  stuff that is common to both constructors into one function) */
  void init();

  /** Destructor. */
  virtual ~AudioModuleEditor();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Passes a new AudioModule objcet to be edited. This should be used when the same editor object 
  should be re-used for editing another AudioModule. */
  virtual void setModuleToEdit(AudioModule* newModuleToEdit);

  /** Sets the pointer to the moduleToEdit member to NULL without doing anything else. This should 
  be called whenever the underlying AudioModule was deleted. */
  virtual void invalidateModulePointer();

  /** Makes this a top-level editor meaning that some additional widgets (global preferences 
  button, infoline etc.) should be drawn. */
  virtual void setAsTopLevelEditor(bool isTopLevel) { isTopLevelEditor = isTopLevel; }

  /** Sets the position of the link to the website. @see positions */
  virtual void setLinkPosition(int newPosition) { linkPosition = newPosition; }

  /** Sets the position of preset load/saev section. @see positions */
  virtual void setPresetSectionPosition(int newPosition) { presetSectionPosition = newPosition; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry:

  /** Returns the bottom (in pixels) of the preset section. */
  virtual int getPresetSectionBottom();

  /** Returns a pointer to the AudioModule that is edited by this editor. It's a const pointer 
  because you are not supposed to change the value of the pointer. */
  const AudioModule* getModuleToEdit() { return moduleToEdit; }


  /** Returns a pointer to the parent AudioModuleEditor, if any and a nullptr if none. */
  AudioModuleEditor* getParentAudioModuleEditor()
  {
    return dynamic_cast<AudioModuleEditor*>(getParentComponent());
  }

  /** Returns a pointer to the repaint manager object that is used to trigger periodic calls to 
  repaint for GUI animation. */
  rsRepaintManager* getRepaintManager();

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks:

  //virtual bool keyPressed(const KeyPress &key, Component *originatingComponent) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged) override;
  virtual void rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave) override;
  virtual void rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled) override;
  virtual void rButtonClicked(RButton *buttonThatWasClicked) override;
  virtual void changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged) override;
  virtual void resized() override;

  /** Updates the widgets according to the state of the assignedParameter (if any) and updates the 
  state-widget set. calls updateWidgetEnablement(). */
  virtual void updateWidgetsAccordingToState() override;

  /** Override this if you want to update the enablement of some widgets according to the state
  of the module. Will be called from updateWidgetsAccordingToState(). */
  virtual void updateWidgetEnablement() {}

  //-----------------------------------------------------------------------------------------------
  // \name Public data members:

  StateLoadSaveWidgetSet* stateWidgetSet;  // \todo check, why we have this in the public area?

protected:

  /** Automatically generates a slider for each parameter in the module which is being edited. */
  //virtual void autoGenerateSliders();

  /** Returns a poiner to an RSlider object with the given name or NULL if no such slider exists
  (in our array automatableSliders) */
  //virtual RSlider* getSliderByName(const juce::String& sliderName);

  /** Opens a dialog to adjust the global preferences like the colour-scheme, preset paths etc.
  If your subclass needs some special settings (like, for example, a sample-path), you may override
  this an open a custom dialog in your class. */
  virtual void openPreferencesDialog();

  /** Loads the current colorscheme into a file. @see aveColorSchemeToFile(). */
  virtual void loadPreferencesFromFile();

  /** Saves the current colorscheme into a file. The filename will be given by the name of the 
  underlying AudioModule concatenated with 'Preferences'. Later we may want to store other settings
  there as well (such as preset- and sample-paths etc.) - we may then have to move the function 
  into AudioModule. */
  virtual void savePreferencesToFile();

  // todo: replace loadColorSchemeFromFile()/saveColorSchemeToFile() with loadPreferencesFromFile()/savePreferencesToFile(),
  // introduce methods getPreferencesAsXml/setPreferencesFromXml - these can then be overrided by subclasses

  /** Returns the xml tag-name that should be used for storing the preferences. */
  virtual juce::String getPreferencesTagName();

  /** Returns the xml filename that should be used for storing the preferences. */
  virtual juce::String getPreferencesFileName();

  RTextField*      infoField;             // field for short help texts when mouse is over a widget
  CriticalSection* lock = nullptr;        // pointer to the global plugInLock
  AudioModule*     moduleToEdit = nullptr;
  int presetSectionPosition, linkPosition;

  // clean this up - it's obsolete:
  /** This is an array of the automatable sliders - if add RSlider objects here, they will be
  updated in AudioModuleEditor::updateWidgetsAccordingToState via calls to their
  updateWidgetFromAssignedParameter() methods. In the destructor, this array is cleared first
  without deleting the objects, such that it does not interfere with the deleteAllChildren-function
  (which is supposed to be called in the destructor). */
  //OwnedArray<RSlider,   CriticalSection> sliders;
  //OwnedArray<RButton,   CriticalSection> buttons;
  //OwnedArray<RComboBox, CriticalSection> comboBoxes;

  // factor out into a class TopLevelEditor (or something like that):
  bool   drawGradientsBasedOnOutlines, isTopLevelEditor;
  Colour gradientMidColour;
  int    numHueOffsets; 

  RClickButton             *setupButton;
  RHyperlinkButton         *webLink;
  ColourSchemeSetupDialog  *setupDialog;

  static const bool screenShotsEnabled = true; // maybe disable for release builds

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleEditor)
};

//=================================================================================================

/** Implements a generic editor for AudioModules which just show a simple widget for each of the 
underlying AudioModule's Parameters. Which kind of widget is shown depends on the type of the 
respective Parameter - for continuous numeric parameters, it shows a slider, for multiple choice 
parameters a combobox and for boolean parameters a button. An object of this generic editor class 
is returned in the baseclass implementation AudioModule::createEditor. */

class GenericAudioModuleEditor : public AudioModuleEditor
{

public:

  GenericAudioModuleEditor(AudioModule* newModuleToEdit);
  virtual void resized() override;

  /** Computes the height required to fit all widgets into the editor. */
  virtual int getRequiredHeight() const
  {
    return (widgetHeight+widgetDistance) * (int)parameterWidgets.size() + 28;
  }

protected:

  /** Creates an appropriate widget for each of the parameters in the underlying AudioModule.
  Called from the constructor. */
  virtual void createWidgets();

  std::vector<RWidget*> parameterWidgets; // array of the widgets for the parameters 
  // why do we need this? the baseclass already has an array for them

  int widgetHeight   = 16; 
  //int widgetDistance = 4;
  int widgetDistance = -2;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(GenericAudioModuleEditor)
};

//=================================================================================================

/** A subclass of AudioModuleEditor for sample-based audio modules, like sample-players, 
wavetable-oscillators, etc. It has additional widgets for loading a sample-file and skipping 
forward/backward through sample files in the current directory. */

class SampleBasedAudioModuleEditor : public AudioModuleEditor, public AudioFileManager
{

public:

  SampleBasedAudioModuleEditor(AudioModule* newModuleToEdit);

  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  // do we need to override 
  //virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  // to update the sampleFileLabel?

protected:

  virtual void createWidgets();

  RTextField *sampleFileLabel;
  RButton    *sampleLoadButton, *samplePlusButton, *sampleMinusButton;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SampleBasedAudioModuleEditor)
};

//=================================================================================================


#endif 