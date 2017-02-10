#ifndef jura_ModuleChainer_h
#define jura_ModuleChainer_h
  
/** A shell module that can be used to create a chain (i.e. series connection) of some number of
AudioModule objects. */

class JUCE_API ModuleChainer : public jura::AudioModuleWithMidiIn
{

public:

  ModuleChainer(CriticalSection *lockToUse);

  //void addModule(AudioModule *moduleToAdd, int position = -1);
  AudioModule* createModule(const String& type);

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

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainer)
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

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModuleChainerEditor)
};


#endif 