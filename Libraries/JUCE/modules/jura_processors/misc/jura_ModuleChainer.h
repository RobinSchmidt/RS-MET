#ifndef jura_ModuleChainer_h
#define jura_ModuleChainer_h
  
/** A do-nothing dummy AudioModule to be used as placeholder. Maybe move somewhere else... */

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

/** A shell module that can be used to create a chain (i.e. series connection) of some number of
AudioModule objects. */

class JUCE_API ModuleChainer : public jura::AudioModuleWithMidiIn
{

public:

  ModuleChainer(CriticalSection *lockToUse);
  virtual ~ModuleChainer();

  /** Adds a module of the given type at the end of the chain. */
  void addModule(const String& type);
    // Maybe replace this function by addEmptySlot - we just use it to create an empty (bypass)
    // slot. Actual modules are then created by replacing the empty dummy module with some other 
    // type. We should also have a function: removeTrailingEmptySlots which gets invoked from
    // replaceModule. The last slot should always be an empty slot.

  /** Replaces the module at the given with a new module of given type unless the given type 
  matches that of the module which is already there at this position in which case nothing 
  happens. */
  void replaceModule(int index, const String& type);

  /** Returns an editor for the AudioModule in the given slot index. */
  AudioModuleEditor* getEditorForSlot(int index);

  /** Returns the editor for the currently active slot. */
  inline AudioModuleEditor* getEditorForActiveSlot() { return getEditorForSlot(activeSlot); }

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void noteOn(int noteNumber, int velocity) override;
  virtual void noteOff(int noteNumber) override;
  virtual void reset() override;
  // override getStateAsXml, etc...

protected:

  Array<AudioModule*> modules;   // maybe use the inherited childModules array instead?
  Array<AudioModuleEditor*> editors;

  int activeSlot = 0;


  friend class ModuleChainerEditor;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainer)
};

//=================================================================================================

/** Implements a GUI editor for the ModuleChainer. */

class JUCE_API ModuleChainerEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  ModuleChainerEditor(jura::ModuleChainer *moduleChainerToEdit);
  virtual ~ModuleChainerEditor();

  virtual void createWidgets();

  /** Updates this editor. This involves figuring out, which slot is active, retrieving the editor
  for the active slot and adding it as child-editor here. This may also cause the GUI to resize
  itself. */
  virtual void updateEditor();

  // overrides:
  virtual void resized() override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;

protected:

  ModuleChainer* chainer;
  Array<AudioModuleSelector*> selectors;

  AudioModuleEditor* activeEditor = nullptr;

  int leftColumnWidth = 160; // for the chainer widgets
  int bottomRowHeight =  20; // for infoline, link, etc.

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainerEditor)
};

#endif 