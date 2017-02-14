#ifndef jura_AudioModule_h
#define jura_AudioModule_h

class AudioModule;
class AudioModuleEditor;

/** A class that can be informed (via a callback method) when an AudioModule object is going to be 
deleted. Mainly intended as baseclass for GUI elements that keep a pointer to an AudioModule that 
is being edited 

hmm.. - maybe, we don't need that. The documentation of the 
AudioProcessorEditor* AudioProcessor::createEditor() function says:
"It's safe to assume that an editor will be deleted before its filter." */

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

/** This class is the base class for all audio modules. */

class JUCE_API AudioModule : public AutomatableModule, public StateFileManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. The caller must pass a pointer to a critical section object that must exist 
  somewhere outside and will be used here to acquire mutually exclusive access from different 
  threads to the underlying dsp object. For example, if the AudioModule is wrapped into a plugin,
  the Criticalsection object could exist on the level of the plugin (an its lifetime be managed 
  there).  */
  AudioModule(CriticalSection *lockToUse);

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
  be typically the name of the module, perhaps appended by some qualifier like "Demo Version". */
  virtual juce::String getModuleHeadlineString();

  /** Returns the interval at which the module wants to receive callbacks to trigger(). */
  virtual double getTriggerInterval() const { return triggerInterval; }

  /** Returns, whether this mdoule wants to save/recall its state - the idea is that some module 
  may have child modules which are  currently inactive and therefore don't need to save and recall 
  their state. This makes preset files more economical. */
  bool wantsSaveAndRecallState() const { return saveAndRecallState; }

  /** Your subclass may override this to return an object of an appropriate subclass of
  AudioModuleEditor. The baseclass implementation will return a generic editor with sliders, 
  comboboxes and button for all the Parameters of this AudioModule. */
  virtual AudioModuleEditor* createEditor();

  //-----------------------------------------------------------------------------------------------
  // automation and state management:

  /** Callback to indicate that a parameter has changed - subclasses should override this and
  update their signal processing accordingly. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
   // is this obsolete now that we use the callback mechanism? if so, get rid...

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

  /** Resets all the parameters to their default values. */
  virtual void resetParametersToDefaultValues();

  //-----------------------------------------------------------------------------------------------
  // Audio processing:

  /** This is the audio callback that your subclass needs to override. */
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) = 0;

  //-----------------------------------------------------------------------------------------------
  // Misc:

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

protected:

  /** Must be overriden by subclasses to fill the inherited array of observed parameters. */
  virtual void initializeAutomatableParameters();

  /** Our child modules to which we will distribute MIDI-events and of which we manage the
  states. */
  juce::Array<AudioModule*, CriticalSection> childModules;


  CriticalSection *plugInLock;     
  // mutex to access the wrapped core dsp object


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

/** A subclass of AudioModule that accepts MIDI input. If you derive your effect or instrument from
this baseclass and wrap it into a juce::AudioProcessor (via the wrapper class jura::AudioPlugin),
the plugin will have a MIDI input. Also, you can override the event handler methods in your 
subclass (noteOn, noteOff, setMidiController, etc.) in order to respond to incoming MIDI events  */

class JUCE_API AudioModuleWithMidiIn : public AudioModule
{

public:

  AudioModuleWithMidiIn(CriticalSection *lockToUse) : AudioModule(lockToUse) {}

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


protected:


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Baseclass for GUI editors for AudioModule objects. */

class AudioModuleEditor : public jura::Editor, public ChangeListener, public RDialogBoxListener, 
  public RButtonListener
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
  // construction/destruction:

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
  // setup:

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
  // inquiry:

  /** Returns the bottom (in pixels) of the preset section. */
  virtual int getPresetSectionBottom();

  /** Returns a pointer to the AudioModule that is edited by this editor. It's a const pointer 
  because you are not supposed to change the value of the pointer. */
  const AudioModule* getModuleToEdit() { return moduleToEdit; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged);
  virtual void rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave);
  virtual void rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled);
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();

  /** Updates the widgets according to the state of the assignedParameter (if any) and updates the 
  state-widget set. calls updateWidgetEnablement(). */
  virtual void updateWidgetsAccordingToState();

  /** Override this if you want to update the enablement of some widgets according to the state
  of the module. Will be called from updateWidgetsAccordingToState(). */
  virtual void updateWidgetEnablement() {}

  //-----------------------------------------------------------------------------------------------
  // public data members:

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

  RTextField      *infoField;       // field for short help texts, when mouse is over a widget
  CriticalSection *plugInLock;      // pointer to the global plugInLock
  AudioModule     *moduleToEdit;
  int presetSectionPosition, linkPosition;


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

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Implements a generic editor for AudioModules which just show a simple widget for each of the 
underlying AudioModule's Parameters. Which kind of widget is shown depends on the type of the 
respective Parameter - for continuous numeric parameters, it shows a slider, for multiple choice 
parameters a combobox and for boolean parameters a button. An object of this generic editor class 
is returned in the baseclass cimplementation AudioModule::createEditor. */

class GenericAudioModuleEditor : public AudioModuleEditor
{

public:

  GenericAudioModuleEditor(AudioModule* newModuleToEdit);
  virtual void resized() override;


protected:

  /** Creates an appropriate widget for each of the parameters in the underlying AudioModule.
  Called from the constructor. */
  virtual void createWidgets();

  juce::Array<RWidget*> parameterWidgets; // array of the widgets for the paraters

  int widgetHeight   = 16; 
  int widgetDistance = 4;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(GenericAudioModuleEditor)
};

#endif 