#include "../JuceLibraryCode/JuceHeader.h"

template<class AudioModuleType>
AudioProcessor* JUCE_CALLTYPE createPlugin(AudioModuleType *dummy)
{
  jura::AudioPlugin *plugIn = new jura::AudioPlugin(nullptr);
  AudioModuleType   *module = new AudioModuleType(&plugIn->plugInLock);
  plugIn->underlyingAudioModule = module;
  return plugIn;
  // \todo: we need a 2nd version of that to handle plugins with MIDI input - or maybe we can try
  // to cast the dummy pointer to AudioModuleWithMidiIn and if that's successful wrap it into an
  // AudioPluginWithMidiIn, otherwise wrap it into a regular AudioPlugin (without midi).
}

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  // We just create a dummy pointer here of the subclass of jura::AudioModule that we want to wrap 
  // into a plugin in order to be able to invoke the template above (the dummy is just for the 
  // preprocessor, such that it can infer the type, for which the template should be instantiated):
  jura::Ladder *dummy = nullptr;

  // Now, invoking the template with the dummy pointer will return an object of the appropriate 
  // class (to which the pointer is declared):
  return createPlugin(dummy);

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
