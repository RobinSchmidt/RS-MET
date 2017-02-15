#ifndef jura_AudioModuleChain_h
#define jura_AudioModuleChain_h
  
/** A do-nothing dummy AudioModule to be used as placeholder. 
todo:
-rename to DummyAudioModule
-maybe move somewhere else in the library
-maybe override create editor to create a special kind of editor that says something like
 "Select Module - Editor will appear here"
*/

class JUCE_API DummyModule : public jura::AudioModule
{
public:
  DummyModule(CriticalSection *lockToUse) : AudioModule(lockToUse) {}
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DummyModule)
};

//=================================================================================================

/** A class for creating objects of various subclasses of AudioModule based on a type string. It 
can also translate back from a given subclass-pointer to the corresponding string and create a list
of all available types. */

class JUCE_API AudioModuleFactory
{

public:

  /** Creates and returns a pointer to an object of some subclass of AudioModule. Which subclass it 
  is, is determined by the passed String parameter. You must also pass the mutex lock object that 
  should be used by the AudioModule. */
  static AudioModule* createModule(const String& type, CriticalSection *lockToUse);

  /** Given a pointer to an object of some subclass of AudioModule, this function returns the
  string that is used to identify the subclass. */
  static String getModuleType(AudioModule *module);

  /** Returns an array of strings with all the available types of AudioModules that can be 
  created. */
  static StringArray getAvailableModuleTypes();

};

//=================================================================================================

/** A widget class for selecting a specific type of AudioModule. */

class JUCE_API AudioModuleSelector : public RComboBox
{
public:
  AudioModuleSelector();
protected:
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleSelector)
};

//=================================================================================================

class JUCE_API AudioModuleChain; // forward declaration

/** Baseclass for objects that must keep track of the state of one (or more) AudioModuleChain 
objects. Observers must override some callback fucntions to take appropriate actions for various 
kinds of state changes. 
\todo: maybe it's sufficient to pass the index in the callbacks - we may not really need to pass 
the pointers along as well
*/

class JUCE_API AudioModuleChainObserver
{

public:

  /** Called whenever a module was added to the chain. Your observer subclass may want to keep a 
  pointer to the module to modify it, create an editor, etc. */
  virtual void audioModuleWasAdded(AudioModuleChain *chain, AudioModule *module, int index) = 0;

  /** Called before modules in the chain will be deleted. Your observer subclass will probably want 
  to invalidate any pointers to the module that it keeps, delete editors, etc. */
  virtual void audioModuleWillBeDeleted(AudioModuleChain *chain, AudioModule *module, 
    int index) = 0;

  /** Called whenever a module in the chain was replaced by another module. Note that the old 
  module may also be deleted after being replaced, so you should invalidate all pointers to it that 
  you may have around. */
  virtual void audioModuleWasBeReplaced(AudioModuleChain *chain, AudioModule *oldModule, 
    AudioModule *newModule, int index) = 0;

};

//=================================================================================================

/** A shell module that can be used to create a chain (i.e. series connection) of some number of
AudioModule objects. 
\todo: 
-keep an array of observers that get notified about state changes
-let the plugin have midi-out, so we can make midi effects too
-organize modules in groups (Generators, Filters, Analyzers, etc.) and use a tree-view for 
 selection
-add more modules
-maybe at some point make an AudioModuleGraph class that allows for free interconnection
 ->especially important for multi I/O modules
-let the plugin switch between chain and graph mode
-rename the plugin to RAPTPlug...or just RAPT
 */

class JUCE_API AudioModuleChain : public jura::AudioModuleWithMidiIn
{

public:

  AudioModuleChain(CriticalSection *lockToUse);
  virtual ~AudioModuleChain();

  /** Adds an empty slot the end of the chain. */
  void addEmptySlot();

  /** Adds a module of the given type at the end of the chain. */
  void addModule(const String& type);

  /** Deletes the module at the given index. */
  void deleteModule(int index);

  /** Removes the last module from the chain. */
  void removeLastModule();
    // rename to deleteLastModule

  /** Replaces the module at the given with a new module of given type unless the given type 
  matches that of the module which is already there at this position in which case nothing 
  happens. Returns true, if the module was replaced, false otherwise. */
  void replaceModule(int index, const String& type);



  /** Returns true if the module at the given index matches the type specified by the type 
  string. */
  bool isModuleOfType(int index, const String& type);

  /** Returns the moduel in the chain at the given index. If the index is out of range, it will 
  return a nullptr. */
  AudioModule* getModuleAt(int index);

  /** Ensures that at the end of the module chain, there is exactly one empty slot that can be used
  to insert another module into the chain. If the last slot is not empty, an empty slot will be 
  added, if there are more than one empty slots at the end, the superfluous ones will be 
  deleted. */
  void ensureOneEmptySlotAtEnd();

   /** Adds an observer that will get notified about changes to the state of the chain. */
  void addAudioModuleChainObserver(AudioModuleChainObserver *observerToAdd);

  /** Removes an oberver that was previously added by addAudioModuleChainObserver. */
  void removeAudioModuleChainObserver(AudioModuleChainObserver *observerToRemove);

  /** Called internally, whenever a module was added to the chain. */
  void sendAudioModuleWasAddedNotification(AudioModule *module, int index);

  /** Called internally, whenever a module was removed from the chain. */
  void sendAudioModuleWillBeDeletedNotification(AudioModule *module, int index);

  /** Called internally, whenever a module in the chain was replaced by another module. */
  void sendAudioModuleWasBeReplacedNotification(AudioModule *oldModule, AudioModule *newModule, 
    int index);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void noteOn(int noteNumber, int velocity) override;
  virtual void noteOff(int noteNumber) override;
  virtual void reset() override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

protected:

  /** Clears the array of AudioModules which means als to delete all objects. */
  void clearModulesArray();
                     
  vector<AudioModule*> modules;  // std::vector better for debugging than juce::Array
  int activeSlot = 0;            // slot for which the editor is currently shown 
  double sampleRate;

  vector<AudioModuleChainObserver*> observers;

  friend class AudioModuleChainEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleChain)
};

//=================================================================================================

/** Implements a GUI editor for the AudioModuleChain.
todo: 
-bug: loading presets from the gui doesn't work (total recall by the host works though, presumably
 because the AudioModuleChain state is recalled before the GUI editor is created) - maybe we need a
 AudioModuleChainOberver baseclass, derive the AudioModuleChainEditor from it and notify the editor 
 about changes ...this may also be used alternatively for informing about module deletion - we 
 could have different callbacks: moduleWillBeDeleted, moduleChainHasChanged, etc.
-when a popup from a combobox is open, the audio throughput is blocked - we need to avoid the lock
-bug: sometimes, we get an access violation when removing a module, the violation occurs in 
 ~AudioModuleEditor() when it calls moduleToEdit->removeStateWatcher(stateWidgetSet); - the 
 moduleToEdit pointer is apparently invalid ....so it seems, we have kept an editor around for a 
 module that has already been deleted - how can this be? we actually delete editors when their
 module gets deleted - also, it does not always happen - strange...
-make it possible to drag the slots up and down to change the order of the modules
-plugin enveloper in 1st slot, plug it out - AudioModuleChainEditor::audioModuleWillBeDeleted is not 
 called - why? bcs we delete the editor already in replaceModule - but it should not be necessary 
 to call it there (see comments in the function)
 */

class JUCE_API AudioModuleChainEditor : public AudioModuleEditor, 
  public AudioModuleDeletionWatcher, public RComboBoxObserver, public ChangeBroadcaster
{

public:

  AudioModuleChainEditor(jura::AudioModuleChain *moduleChainerToEdit);
  virtual ~AudioModuleChainEditor();

  /** Returns an editor for the AudioModule in the given slot index. Note that this may return a 
  nullptr in the case when the "modules" and "editors" arrays are empty (this occurs as a 
  transitional situation when recalling the state of the chainer from an xml).*/
  AudioModuleEditor* getEditorForSlot(int index);

  /** Returns the editor for the currently active slot.  */
  inline AudioModuleEditor* getEditorForActiveSlot() 
  { 
    return getEditorForSlot(chainer->activeSlot); 
  }

  /** Replaces the module at the given with a new module of given type, if necessary and also 
  replaces the corresponding editor. */
  void replaceModule(int index, const String& type);

  /** Updates our array of selector-widgets (comboboxes) to select the module for each slot. */
  void updateSelectorArray();

  /** Updates our array of AudioModuleEditors to match the number of modules of the edited 
  AudioModuleChain. */
  void updateEditorArray();

  /** Updates this editor to show the module-editor of the currently active slot. This may also 
  cause the GUI to resize itself. */
  void updateActiveEditor();

  // overrides:
  virtual void mouseDown(const MouseEvent &e) override;
  virtual void resized() override;
  virtual void paintOverChildren(Graphics& g) override;
  virtual void audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted) override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;
  virtual void changeListenerCallback(ChangeBroadcaster *source) override;

protected:

  /** Sends out a change message that we will receive ourselves. On receive, we will call
  updateSelectorArray. This mechanism is used to cause a deferred update of the selectors array 
  from replaceModule. The deferrence is necessray, because replaceModule is called from 
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
  AudioModuleChain* chainer;                  // the edited object
  vector<AudioModuleSelector*> selectors;     // combo-boxes for selecting modules
  vector<AudioModuleEditor*>   editors;       // array of editors for the modules

  AudioModuleEditor* activeEditor = nullptr;  // currently shown editor

  int leftColumnWidth = 160; // for the chainer widgets
  int bottomRowHeight =  16; // for infoline, link, etc.

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleChainEditor)
};

#endif 