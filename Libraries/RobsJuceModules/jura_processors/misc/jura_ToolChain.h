#ifndef jura_ToolChain_h
#define jura_ToolChain_h
  
//=================================================================================================

class JUCE_API ToolChain; // forward declaration

/** Baseclass for objects that must keep track of the state of one (or more) ToolChain 
objects. Observers must override some callback fucntions to take appropriate actions for various 
kinds of state changes. 
\todo: maybe it's sufficient to pass the index in the callbacks - we may not really need to pass 
the pointers along as well
*/


class JUCE_API ToolChainObserver
{

public:
    
  virtual ~ToolChainObserver() {}

  /** Called whenever a module was added to the chain. Your observer subclass may want to keep a 
  pointer to the module to modify it, create an editor, etc. */
  virtual void audioModuleWasAdded(ToolChain *chain, AudioModule *module, int index) = 0;

  /** Called before modules in the chain will be deleted. Your observer subclass will probably want 
  to invalidate any pointers to the module that it keeps, delete editors, etc. */
  virtual void audioModuleWillBeDeleted(ToolChain *chain, AudioModule *module, 
    int index) = 0;

  /** Called whenever a module in the chain was replaced by another module. Note that the old 
  module may also be deleted after being replaced, so you should invalidate all pointers to it that 
  you may have around. */
  virtual void audioModuleWasReplaced(ToolChain *chain, AudioModule *oldModule, 
    AudioModule *newModule, int index) = 0;

};

//=================================================================================================

/** A shell module that can be used to create a chain (i.e. series connection) of some number of
AudioModule objects. 
\todo: 
-have a Dry/Wet and Gain slider
-let Modulators have their own slot-array
-let instruments accumulate their outputs (layering)
-let the plugin have midi-out, so we can make midi effects too
-add more modules
-maybe at some point make an AudioModuleGraph class that allows for free interconnection
 ->especially important for multi I/O modules
-let the plugin switch between chain and graph mode

 */

class JUCE_API ToolChain :  public jura::AudioModuleWithMidiIn 
{

public:

  ToolChain(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr);

  virtual ~ToolChain();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adds an empty slot the end of the chain. */
  void addEmptySlot();

  /** Tries to add a module of the given type at the end of the chain and reports whether this was 
  successful. */
  bool addModule(const juce::String& type);

  /** Adds the passed AudioModule at the end of the chain. This ToolChain object will take over 
  ownership of the module. */
  void addModule(AudioModule* module);

  /** Deletes the module at the given index. */
  void deleteModule(int index);

  /** Removes the last module from the chain. */
  void deleteLastModule();

  /** Replaces the module at the given index with a new module of given type unless the given type 
  matches that of the module which is already there at this position in which case nothing 
  happens. 
  todo: maybe it should return true, if the module was actually replaced, false otherwise */
  void replaceModule(int index, const juce::String& type);

  // todo:
  //void moveModule(int oldIndex, int newIndex);


  //-----------------------------------------------------------------------------------------------

  /** Returns true if the module at the given index matches the type specified by the type 
  string. */
  bool isModuleOfType(int index, const juce::String& type);

  /** Returns the moduel in the chain at the given index. If the index is out of range, it will 
  return a nullptr. */
  AudioModule* getModuleAt(int index);


  /** Ensures that at the end of the module chain, there is exactly one empty slot that can be used
  to insert another module into the chain. If the last slot is not empty, an empty slot will be 
  added, if there are more than one empty slots at the end, the superfluous ones will be 
  deleted. */
  void ensureOneEmptySlotAtEnd();

  // observer stuff:

   /** Adds an observer that will get notified about changes to the state of the chain. */
  void addToolChainObserver(ToolChainObserver *observerToAdd);

  /** Removes an oberver that was previously added by addToolChainObserver. */
  void removeToolChainObserver(ToolChainObserver *observerToRemove);

  /** Called internally, whenever a module was added to the chain. */
  void sendAudioModuleWasAddedNotification(AudioModule *module, int index);

  /** Called internally, whenever a module was removed from the chain. */
  void sendAudioModuleWillBeDeletedNotification(AudioModule *module, int index);

  /** Called internally, whenever a module in the chain was replaced by another module. */
  void sendAudioModuleWasReplacedNotification(AudioModule *oldModule, AudioModule *newModule, 
    int index);


  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor(int type) override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 

  virtual void handleMidiMessage(MidiMessage message) override;
  /*
  virtual void noteOn(int noteNumber, int velocity) override;
  virtual void noteOff(int noteNumber) override;
  virtual void setMidiController(int controllerNumber, float controllerValue) override;
  virtual void setPitchBend(int pitchBendValue) override;
  */

  virtual void reset() override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;



  // temporariliy made public:

  std::vector<AudioModule*> modules;
  // Elan's editor needs to access it, when tabs are moved - todo:
  // provide a function moveModule that can be called from the editor
  // we should better use the inherited childAudioModules array - but there are errors

  AudioModuleFactory moduleFactory;  
  // todo: provide getter - editor subclasses need access to getRegisteredModuleInfos etc.

  int activeSlot = 0;            // slot for which the editor is currently shown 
  // for Elan's subclass



protected:

  void recallSlotsFromXml(      const XmlElement &xmlState, bool markAsClean);
  void recallModulationsFromXml(const XmlElement &xmlState); // move to ModulatbleAudioModule

  /** Sets up the pointers to the SmoothingManager, MetaParameterManager and ModulationManager in
  the passed module. */
  void setupManagers(AudioModule* module);

  // Might become relevant when we want to allow the user to change the maxNumVoices at runtime. 
  // Currently, this is fixed after construction..
  //void allocateVoiceResources(rosic::rsVoiceManager* voiceManager) override;

  /** Checks, if the passed AudioModule can be cast into a ModulationSource and if so, adds it to
  our array of ModulationSources (inherited from ModulationManager). */
  void addToModulatorsIfApplicable(AudioModule* module);
    // maybe factor out into a class ModulatableAudioModule which is subclass of AudioModule and 
    // ModulationManager

  /** Undoes what addToModulatorsIfApplicable does. */
  void removeFromModulatorsIfApplicable(AudioModule* module);

  /** Assigns an appropriate name to the passed ModulationSource which will be used to identify it
  in the modulation setup on the GUI and for state recall. */
  void assignModulationSourceName(ModulationSource* source);

  /** Clears the array of AudioModules which means als to delete all objects. */
  void clearModulesArray();


  void createMidiModSources();

  void deleteMidiModSources();

  /** Just some temporary throwaway code to figure out what is going wrong with the mod-system in
  Elan's SpiralGenerator. */
  void createDebugModSourcesAndTargets();

  /** Populates our AudioModuleFactory with the modules that can be plugged in. */
  void populateModuleFactory();
                     

  //ModulationManager modManager;
  ModulationManagerPoly *modManager;  
  // This is really strange: when switching to using ModulationManagerPoly instead of the baseclass
  // ModulationManager, i also had to change to using a pointer instead of a direct object because
  // otherwise i got heap corruptions. Why is that? Is it because we pass pointers to it other 
  // objects and the address operator does not work properly when it's a subclass? It happened 
  // right after creating the subclass - it literally did not have any additional member variables
  // or functions or overrides.

  // name clash with modManager inherited ModulationParticipant (baseclass of 
  // ModulatableAudioModule) - maybe rename this

  //rsSmoothingManager smoothingManager;
  //std::vector<AudioModule*> modulators;

  double sampleRate = 44100;
  std::vector<ToolChainObserver*> observers;


  rsVoiceManager voiceManager;
  std::vector<double> voiceSignals; // Used to share/re-use a single 
  // buffer for the voice-signals of the modules (using either overwriting or accumulation - 
  // whatever is most appropriate for the particular module). For the time being, each module 
  // allocates it own buffer.

  // The modulation sources that are always available:
  rsConstantOneModulatorModulePoly*  constantModulator;
  rsNotePitchModulatorModulePoly*    notePitchModulator;
  rsNoteFreqModulatorModulePoly*     noteFreqModulator;
  rsNoteVelocityModulatorModulePoly* noteVelocityModulator;



  friend class ToolChainEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ToolChain)
};

//=================================================================================================

/** Implements a GUI editor for the ToolChain.
\todo: 
-Enveloper: set module name of the embedded Modulator module to "Enveloper"
-add bypass switches for each module
-make the infoline work for the selectors
-make it possible to drag the slots up and down to change the order of the modules
 */

class JUCE_API ToolChainEditor : public AudioModuleEditor, public ToolChainObserver,
  public RComboBoxObserver, public ChangeBroadcaster
{

public:

  ToolChainEditor(jura::ToolChain *moduleChainToEdit);
  virtual ~ToolChainEditor();

  /** Returns an editor for the AudioModule in the given slot index. Note that this may return a 
  nullptr in the case when the "modules" and "editors" arrays are empty (this occurs as a 
  transitional situation when recalling the state of the chainer from an xml).*/
  AudioModuleEditor* getEditorForSlot(int index);

  /** Returns the editor for the currently active slot.  */
  inline AudioModuleEditor* getEditorForActiveSlot() 
  { 
    return getEditorForSlot(chain->activeSlot); 
  }

  /** Replaces the module at the given with a new module of given type, if necessary and also 
  replaces the corresponding editor. */
  void replaceModule(int index, const juce::String& type);

  /** Updates our array of selector-widgets (comboboxes) to select the module for each slot. */
  void updateSelectorArray();

  /** Updates our array of AudioModuleEditors to match the number of modules of the edited 
  ToolChain. */
  void updateEditorArray();

  /** Updates this editor to show the module-editor of the currently active slot. This may also 
  cause the GUI to resize itself. */
  void updateActiveEditor();

  /** Updates the highlighting of the selector dropdowns such that the currently active slot is
  highlighted. */
  void updateActiveSelector();

  // overrides:
  virtual void mouseDown(const MouseEvent &e) override;
  virtual void resized() override;
  virtual void paintOverChildren(Graphics& g) override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;
  virtual void changeListenerCallback(ChangeBroadcaster *source) override;
  virtual void audioModuleWasAdded(ToolChain *chain, 
    AudioModule *module, int index) override;
  virtual void audioModuleWillBeDeleted(ToolChain *chain, 
    AudioModule *module, int index) override;
  virtual void audioModuleWasReplaced(ToolChain *chain, 
    AudioModule *oldModule, AudioModule *newModule, int index) override;


protected:

  /** Sends out a change message that we will receive ourselves. On receive, we will call
  updateSelectorArray. This mechanism is used to cause a deferred update of the selectors array 
  from replaceModule. The deferrence is necessarry, because replaceModule is called from 
  rComboBoxChanged - if we would call updateSelectorArray directly in replaceModule, we would 
  possibly delete the combobox that has changed before rComboBoxChanged returns which results
  in a combox trying to update itself with an invalid this-pointer. So, we need a deferred 
  destruction. */
  void scheduleSelectorArrayUpdate();

  /** Deletes the editor at given index in the array. The slot entry will be replaced by 
  nullptr. */
  void deleteEditor(int index);

  /** Deletes all the editors in our array and clears the array itself. */
  void clearEditorArray();

  // Data:
  ToolChain* chain;                    // the edited object
  vector<AudioModuleSelector*> selectors;     // combo-boxes for selecting modules
  vector<AudioModuleEditor*>   editors;       // array of editors for the modules

  AudioModuleEditor* activeEditor = nullptr;  // currently shown editor

  int leftColumnWidth = 160; // for the chainer widgets
  int bottomRowHeight =  16; // for infoline, link, etc.

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ToolChainEditor)
};

#endif 