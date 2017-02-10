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

/** A shell module that can be used to create a chain (i.e. series connection) of some number of
AudioModule objects. */

class JUCE_API ModuleChainer : public jura::AudioModuleWithMidiIn
{

public:

  ModuleChainer(CriticalSection *lockToUse);

  virtual ~ModuleChainer();


  /** Creates and returns a pointer to an object of some subclass of AudioModule. Which subclass it 
  is, is determined by the passed String parameter. */
  AudioModule* createModule(const String& type);

  /** Adds a module of the given type at the end of the chain. */
  void addModule(const String& type);

  // todo:
  //void addModule(AudioModule *moduleToAdd, int position = -1);
  //void replaceModule(int position, AudioModule replacementModule);
  //void removeModule(int position);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void noteOn(int noteNumber, int velocity) override;
  virtual void noteOff(int noteNumber) override;
  virtual void reset() override;
  // override getStateXml, etc...

protected:

  Array<AudioModule*> modules;
  //OwnedArray<AudioModule*> modules;

  friend class ModuleChainerEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainer)
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

/** Implements a GUI editor for the ModuleChainer. */

class JUCE_API ModuleChainerEditor : public AudioModuleEditor
{

public:

  ModuleChainerEditor(jura::ModuleChainer *moduleChainerToEdit);

  //virtual void createWidgets();
  //virtual void resized() override;

protected:

  ModuleChainer *chainer;
  Array<AudioModuleSelector*> selectors;
  Array<AudioModuleEditor*>   editors;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainerEditor)
};


#endif 