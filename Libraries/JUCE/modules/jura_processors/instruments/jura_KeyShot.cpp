
KeyShotAudioModule::KeyShotAudioModule(CriticalSection *newPlugInLock, rosic::KeyShot *keyShotToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, keyShotToWrap)
{
  jassert(keyShotToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedKeyShot = keyShotToWrap;
  underlyingRosicInstrument = keyShotToWrap;
  moduleName = juce::String("KeyShot");
  setActiveDirectory(getApplicationDirectory() + juce::String("/KeyShotPresets") );

  samplePlayerModule = new SamplePlayerAudioModule(lock, 
    &wrappedKeyShot->voiceArray[0].samplePlayer);
  samplePlayerModule->setModuleName(juce::String("SamplePlayer"));
  addChildAudioModule(samplePlayerModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedKeyShot->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}

