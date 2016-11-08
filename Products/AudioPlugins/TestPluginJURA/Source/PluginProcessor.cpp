#include "../JuceLibraryCode/JuceHeader.h"



template<class AudioModuleType>
AudioProcessor* JUCE_CALLTYPE createPluginWithoutMidi(AudioModuleType *dummy)
{
  // wraps audio module into plugin without midi input
  jura::AudioPlugin *plugIn = new jura::AudioPlugin(nullptr);
  AudioModuleType   *module = new AudioModuleType(&plugIn->plugInLock);
  plugIn->underlyingAudioModule = module;
  return plugIn;
}

template<class AudioModuleType>
AudioProcessor* JUCE_CALLTYPE createPluginWithMidi(AudioModuleType *dummy)
{
  // wraps audio module into plugin with midi input
  jura::AudioPluginWithMidiIn *plugIn = new jura::AudioPluginWithMidiIn(nullptr);
  AudioModuleType *module = new AudioModuleType(&plugIn->plugInLock);
  plugIn->underlyingAudioModule   = module;
  plugIn->wrappedModuleWithMidiIn = module;
  return plugIn;
}

//template<class AudioModuleType>
//AudioProcessor* JUCE_CALLTYPE createPlugin(AudioModuleType *dummy, bool withMidiIn)
//{
//  // dispatcher between the different wrappers above
//  if( withMidiIn == true )
//    return createPluginWithMidi(dummy);
//  else
//    return createPluginWithoutMidi(dummy);
//
//  // todo:
//  // Maybe, we can somehow infer, whether or not we need to create a plugin with or without midi 
//  // from the type of the passed pointer to get rid of the boolean flag "withMidiIn". But 
//  // dynamic_cast and checking agianst a nullptr doesn't work when we actually pass a nullptr 
//  // (which is what we do). Maybe somethign using decltype or type-traits or some other C++11 
//  // feature may help - we'll see.
//}

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  // We just create a dummy pointer here of the subclass of jura::AudioModule that we want to wrap 
  // into a plugin in order to be able to invoke the template above (the dummy is just for the 
  // preprocessor, such that it can infer the type, for which the template should be instantiated):

  //jura::Ladder    *dummy = nullptr; return createPluginWithoutMidi(dummy);
  jura::Enveloper *dummy = nullptr; return createPluginWithMidi(dummy);

  // Now, invoking the template with the dummy pointer will return an object of the appropriate 
  // class (to which the pointer is declared):
  //return createPlugin(dummy, withMidiIn);
  //return createPluginWithoutMidi(dummy);
  //return createPluginWithMidi(dummy);


  // This trick saves us from writing out the following code for each AudioModule subclass, which 
  // we want to wrap into a plugin:
  //  jura::AudioPlugin *plugIn = new jura::AudioPlugin(nullptr);
  //  jura::Ladder      *module = new jura::Ladder(&plugIn->plugInLock);
  //  plugIn->underlyingAudioModule = module;
  //  return plugIn;
  //..for example, if we wanted to create a plugin from the Ladder AudioModule. Instead of writing 
  // out these 4 lines each time and wrapping them into a function (giving another 2 lines), we 
  // just need to change the dummy-pointer declaration above. So, we can use this single project to
  // build many plugins - that's very convenient and maintenance friendly. :-) ...OK - we may have 
  // to manually rename the .dll after the build but that's better than having to maintain a dozen 
  // of projects...
}
