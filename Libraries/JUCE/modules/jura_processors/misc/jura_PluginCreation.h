#ifndef jura_PluginCreation_h
#define jura_PluginCreation_h

/** This file contains some templates that make it very easy to wrap any subclass of 
jura::AudioModule into a juce::AudioProcessor (actually, a jura::AudioPlugin which is a subclass of
juce::AudioProcessor) and compile it as plugin. Your plugin project will just need to implement a 
function createPluginFilter(), declare a nullptr of the type of the respective subclass of 
jura::AudioModule and then invoke the appropriate template with that nullptr (there are different 
templates for plugins with and without midi, because we need different wrapper classes). For 
example, in your plugin project, all you need is the following code:

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  jura::Ladder *dummy = nullptr; return createPluginWithoutMidi(dummy);
}

and that will wrap a jura::Ladder module into a plugin (in this case, a plugin without midi 
input). */

template<class AudioModuleType>
AudioProcessor* JUCE_CALLTYPE createPluginWithoutMidi(AudioModuleType *dummy)
{
  // wraps audio module into plugin without midi input
  jura::AudioPlugin *plugIn = new jura::AudioPlugin(nullptr);
  AudioModuleType   *module = new AudioModuleType(&plugIn->plugInLock);
  plugIn->setAudioModuleToWrap(module);
  return plugIn;
}

template<class AudioModuleType>
AudioProcessor* JUCE_CALLTYPE createPluginWithMidi(AudioModuleType *dummy)
{
  // wraps audio module into plugin with midi input
  jura::AudioPluginWithMidiIn *plugIn = new jura::AudioPluginWithMidiIn();
  AudioModuleType *module = new AudioModuleType(&plugIn->plugInLock);
  plugIn->setAudioModuleToWrap(module);
  return plugIn;

  // only the 1st line is different, the last 3 are duplicated -> factor out
}


//// The code below was supposed to dispatch between the createPluginWithoutMidi and 
//// createPluginWithMidi version of the plugin creation code above. Unfortunately, it doesn't work 
//// that way, so we must manually call the appropriate version in out our actual plugin project. 
//// Maybe at some point late, we'll find a solution...
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


#endif 