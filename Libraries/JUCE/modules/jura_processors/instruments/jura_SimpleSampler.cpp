
//-------------------------------------------------------------------------------------------------
// construction/destruction:

SimpleSamplerAudioModule::SimpleSamplerAudioModule(CriticalSection *newPlugInLock, rosic::SimpleSampler *simpleSamplerToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, simpleSamplerToWrap)
{
  jassert(simpleSamplerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedSimpleSampler      = simpleSamplerToWrap;
  underlyingRosicInstrument = simpleSamplerToWrap;
  moduleName = juce::String("SimpleSampler");
  setActiveDirectory(getApplicationDirectory() + juce::String("/SimpleSamplerPresets") );

  samplePlayerModule = new SamplePlayerAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].oscSection.samplePlayer1);
  samplePlayerModule->setModuleName(juce::String("SamplePlayer"));
  addChildAudioModule(samplePlayerModule);

  filterModule = new MultiModeFilterAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].filter);
  filterModule->setModuleName(juce::String("Filter"));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String("PitchEnvelope"));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String("FilterEnvelope"));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}