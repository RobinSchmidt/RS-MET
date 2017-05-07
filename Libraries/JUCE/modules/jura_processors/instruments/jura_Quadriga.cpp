
QuadrigaAudioModule::QuadrigaAudioModule(CriticalSection *newPlugInLock, 
  rosic::Quadriga *quadrigaToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, quadrigaToWrap)
{
  jassert(quadrigaToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedQuadriga = quadrigaToWrap;
  underlyingRosicInstrument = quadrigaToWrap;
  moduleName = juce::String("Quadriga");

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String("/QuadrigaPresets") );

  quadrigenModule = new QuadrigenAudioModule(lock, &wrappedQuadriga->voiceArray[0].quadrigen);
  quadrigenModule->setModuleName(juce::String("Quadrigen"));
  addChildAudioModule(quadrigenModule);

  quadrifexModule = new QuadrifexAudioModule(lock, &wrappedQuadriga->quadrifex);
  quadrifexModule->setModuleName(juce::String("Quadrifex"));
  addChildAudioModule(quadrifexModule);

  equalizerModule = new EqualizerAudioModule(lock, &wrappedQuadriga->equalizer);
  equalizerModule->setModuleName(juce::String("Equalizer"));
  addChildAudioModule(equalizerModule);
}