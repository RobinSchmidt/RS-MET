#ifndef jura_AudioModuleSlotArray_h
#define jura_AudioModuleSlotArray_h


class JUCE_API AudioModuleSlotArray; // forward declaration


class JUCE_API AudioModuleSlotObserver
{

public:

  virtual ~AudioModuleSlotObserver() {}

  /** Called whenever a module was added to the chain. Your observer subclass may want to keep a 
  pointer to the module to modify it, create an editor, etc. */
  virtual void audioModuleWasAdded(AudioModuleSlotArray *chain, 
    AudioModule *module, int index) = 0;

  /** Called before modules in the chain will be deleted. Your observer subclass will probably want 
  to invalidate any pointers to the module that it keeps, delete editors, etc. */
  virtual void audioModuleWillBeDeleted(AudioModuleSlotArray *chain, AudioModule *module, 
    int index) = 0;

  /** Called whenever a module in the chain was replaced by another module. Note that the old 
  module may also be deleted after being replaced, so you should invalidate all pointers to it that 
  you may have around. */
  virtual void audioModuleWasReplaced(AudioModuleSlotArray *chain, AudioModule *oldModule, 
    AudioModule *newModule, int index) = 0;

};



//=================================================================================================

class JUCE_API AudioModuleSlotArray : public jura::AudioModuleWithMidiIn
{
  // into this class we want to factor out stuff from ToolChain that can be re-used in other
  // slot based modules

public:

  AudioModuleSlotArray(CriticalSection *lockToUse, 
    MetaParameterManager* metaManagerToUse = nullptr);

  virtual ~AudioModuleSlotArray();

  /** Adds an observer that will get notified about changes to the state of the chain. */
  void addModuleSlotObserver(AudioModuleSlotObserver *observerToAdd);

  /** Removes an oberver that was previously added by addToolChainObserver. */
  void removeModuleSlotObserver(AudioModuleSlotObserver *observerToRemove);


protected:

  std::vector<AudioModuleSlotObserver*> observers;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleSlotArray)
};



/*
//=================================================================================================

class JUCE_API AudioModuleSlotArrayEditor
{

public:



protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleSlotArrayEditor)
};

*/


#endif