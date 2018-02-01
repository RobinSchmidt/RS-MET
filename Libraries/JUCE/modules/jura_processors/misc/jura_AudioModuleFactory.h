#ifndef jura_AudioModuleFactory_h
#define jura_AudioModuleFactory_h

/** Structure to represent properties of an AudioModule. Contains also a factory function to 
create an instance of the module. */

struct JUCE_API AudioModuleInfo
{
  AudioModuleInfo(const juce::String& typeName, 
    std::function<AudioModule*(CriticalSection*)> creatorFunction, 
    const juce::String& categoryName)
    : type(typeName), createInstance(creatorFunction), category(categoryName)
  {}
  juce::String type;
  juce::String category;
  std::function<AudioModule*(CriticalSection*)> createInstance;
};

//=================================================================================================

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

  /** Constructor. You must pass the CriticalSection object that will be passed on to the 
  constructors of the AudioModule objects that we create here. */
  AudioModuleFactory(CriticalSection *lockToUse) : lock(lockToUse) {} 

  /** Registers an AudioModule. You must pass a function that creates and returns the module and 
  you can specify a category (this may be used in a tree-view selector widget). You may also 
  optionally pass the type-name of the module. If you pass nothing, a module will be instatiated 
  and getModuleTypeName will be called and the returned value will be used. If you do pass a 
  string, it must match whatever getModuleTypeName returns (this is for optimization, to avoid
  instantiating all modules just for the purpose of registering). */
  void registerModuleType(std::function<AudioModule*(CriticalSection*)> creatorFunction, 
    const juce::String& category = String::empty, const juce::String& typeName = String::empty); 

  /** Tries to find a module of given type in our registry of AudioModuleInfos and if it finds it,
  it will return an instance of a module of that type. If the type is not found, it returns a 
  "Bypass" module. */
  AudioModule* createModule(const juce::String& type);

  /** Returns a reference to our array of AudioModuleInfos. This is inquired by widgets that are
  used to select a module. */
  const std::vector<AudioModuleInfo>& getRegisteredModuleInfos() const { return moduleInfos; }


protected:

  /** Mutex to be passed to the AudioModule constructors. */
  CriticalSection* lock;

  /** Our vector of module-infos (todo: maybe use a tree (by category)) */
  std::vector<AudioModuleInfo> moduleInfos;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleFactory)
};

#endif 