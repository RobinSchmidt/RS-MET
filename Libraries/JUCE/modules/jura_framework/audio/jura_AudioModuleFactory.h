#ifndef jura_AudioModuleFactory_h
#define jura_AudioModuleFactory_h

/** A special dummy AudioModule type, an instance of which is returned by AudioModuleFactory when a 
module type is requested that is unknown to the factory object. */

class JUCE_API NotFoundAudioModule : public AudioModule
{

public:

  NotFoundAudioModule(CriticalSection *lockToUse, 
    const juce::String& errorText = "ERROR: Unknown module type")
    : AudioModule(lockToUse) 
  {
    setModuleTypeName(errorText);
    setModuleName(errorText);
  }

  // maybe override get/setState function setState just stores the passed xml and getState returns 
  // it - this way, the settings are not lost when a module is not found

  // todo: the error message should appear also in the AudioModuleSelector (currently, it appears
  // only the editor headline)

  // also, we should perhaps save the type-name of the requested module type here (for the state 
  // xml stuff)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NotFoundAudioModule)
};

//=================================================================================================

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

/** A class for creating objects of various subclasses of AudioModule based on a type string. 
Before using it, you must register the AudioModule types that you want to be created. You can do 
this by using registerModuleType. After registering your module types, you can call createModule 
with a string representing the module type. This function will the return an instance of the 
requested type (or an instance of a special NotFoundAudioModule, if the requested type is unknown 
(i.e. was never registered)). */

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
  it will return an instance of a module of that type. If the type is not found, it returns an 
  instance of NotFoundAudioModule. */
  AudioModule* createModule(const juce::String& type);

  /** Returns a reference to our array of AudioModuleInfos. This is inquired by widgets that are
  used to select a module. */
  const std::vector<AudioModuleInfo>& getRegisteredModuleInfos() const { return moduleInfos; }

  /** Returns an array containing the type strings of the registered module types. */
  std::vector<String> getRegisteredModuleTypes() const;

protected:

  /** Mutex to be passed to the AudioModule constructors. */
  CriticalSection* lock;

  /** Our vector of module-infos (todo: maybe use a tree (by category)) */
  std::vector<AudioModuleInfo> moduleInfos;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleFactory)
};

#endif 