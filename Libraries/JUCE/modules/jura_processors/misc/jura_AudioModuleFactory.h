#ifndef jura_AudioModuleFactory_h
#define jura_AudioModuleFactory_h



/** A class for creating objects of various subclasses of AudioModule based on a type string. It 
can also translate back from a given subclass-pointer to the corresponding string and create a list
of all available types. 

\todo: maybe make this an abstract factory - that way, the ToolChain could be parameterized
with a factory object and could propagated up into jura_framework (next to AudioModule/AudioPlugin)
without knowing about the actual kinds of AudioModule subclasses that are defined in the 
jura_processors module. The actual Chainer plugin would then somehow need to get an object of a 
subclass of AudioModuleFactory passed - and only that subclass woul know all the different kinds of 
modules defined in jura_processors. the general chaining-logic could be made independent from the 
concrete set of AudioModule types that can be created. */

class JUCE_API AudioModuleFactory
{

public:

  AudioModuleFactory() {} // pass the lock to use

  /** Creates and returns a pointer to an object of some subclass of AudioModule. Which subclass it 
  is, is determined by the passed String parameter. You must also pass the mutex lock object that 
  should be used by the AudioModule. You may also optionally pass a ModulationManager object that 
  will be used for AudioModules with modulatable parameters. 
  ToDo: pass an optional MetaParameterManager the same way as the ModulationManager is passed.
  */
  static AudioModule* createModule(const juce::String& type, CriticalSection* lockToUse, 
    ModulationManager* modManager = nullptr, MetaParameterManager* metaManager = nullptr);


  void registerModuleType(const juce::String& typeName, AudioModule* (*creatorFunction)(), 
    const juce::String& category = String::empty);

  /** Returns an array of strings with all the available types of AudioModules that can be 
  created. */
  //static StringArray getAvailableModuleTypes();

  CriticalSection* lock;

protected:

  /** Structure to represent properties of an AudioModule. */
  struct JUCE_API AudioModuleInfo
  {
    AudioModuleInfo(const juce::String& typeName, AudioModule* (*creatorFunction)(),
      const juce::String& categoryName)
      : type(typeName), createInstance(creatorFunction), category(categoryName)
    {}
    juce::String type;
    juce::String category;
    AudioModule* (*createInstance)();
  };

  /** Our vector of module-infos (todo: maybe use a tree (by category)) */
  std::vector<AudioModuleInfo> moduleInfos;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleFactory)
};

#endif 