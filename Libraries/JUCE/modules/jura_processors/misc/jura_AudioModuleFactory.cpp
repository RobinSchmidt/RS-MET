void AudioModuleFactory::registerModuleType(
  std::function<AudioModule*(CriticalSection*)> creatorFunction, 
  const juce::String& category, const juce::String& typeName)
{
#ifdef JUCE_DEBUG
  AudioModule* m = creatorFunction(lock);
  jassert(typeName == m->getModuleTypeName()); // the name under which the module is registered
  delete m;                                    // must match the module's type name (for recall)
#endif

  if(typeName == String::empty)
  {
    AudioModule* m = creatorFunction(lock);
    moduleInfos.push_back(AudioModuleInfo(m->getModuleTypeName(), creatorFunction, category));
    delete m; 
  }
  else
    moduleInfos.push_back(AudioModuleInfo(typeName, creatorFunction, category));
}

AudioModule* AudioModuleFactory::createModule(const juce::String& type)
{
  for(size_t i = 0; i < moduleInfos.size(); i++)
  {
    if(moduleInfos[i].type == type)
      return moduleInfos[i].createInstance(lock);
  }
  jassertfalse;    // unknown module type requested
  //return nullptr;
  return new BypassAudioModule(lock); // to avoid a crash when loading a messed up xml file
                                      // maybe have a special "Module Not Found" module that shows
                                      // the name of the not found module and error text
  //return new DummyModule(lock); // maybe... 
}