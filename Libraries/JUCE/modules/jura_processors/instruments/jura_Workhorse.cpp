
WorkhorseAudioModule::WorkhorseAudioModule(CriticalSection *newPlugInLock, rosic::Workhorse *workhorseToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, workhorseToWrap)
{
  jassert(workhorseToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedWorkhorse        = workhorseToWrap;
  underlyingRosicInstrument = workhorseToWrap;
  moduleName                = juce::String(("Workhorse"));

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String(("/WorkhorsePresets")) );

  vectorMixerModule = new VectorMixerAudioModule(lock, &wrappedWorkhorse->vectorMixer);
  vectorMixerModule->setModuleName(juce::String(("VectorMixer")));
  addChildAudioModule(vectorMixerModule);

  samplePlayerTopLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopLeft);
  samplePlayerTopLeftModule->setModuleName(juce::String(("SamplePlayerTopLeft")));
  addChildAudioModule(samplePlayerTopLeftModule);

  samplePlayerTopRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopRight);
  samplePlayerTopRightModule->setModuleName(juce::String(("SamplePlayerTopRight")));
  addChildAudioModule(samplePlayerTopRightModule);

  samplePlayerBottomLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomLeft);
  samplePlayerBottomLeftModule->setModuleName(juce::String(("SamplePlayerBottomLeft")));
  addChildAudioModule(samplePlayerBottomLeftModule);

  samplePlayerBottomRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomRight);
  samplePlayerBottomRightModule->setModuleName(juce::String(("SamplePlayerBottomRight")));
  addChildAudioModule(samplePlayerBottomRightModule);

  filterModule = new MultiModeFilterAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].filter);
  filterModule->setModuleName(juce::String(("Filter")));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String(("PitchEnvelope")));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String(("FilterEnvelope")));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String(("AmpEnvelope")));
  addChildAudioModule(ampEnvModule);
}

