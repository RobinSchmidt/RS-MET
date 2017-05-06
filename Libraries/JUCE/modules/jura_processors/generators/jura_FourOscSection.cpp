
//-------------------------------------------------------------------------------------------------
// construction/destruction:

FourOscSectionAudioModule::FourOscSectionAudioModule(CriticalSection *newPlugInLock, rosic::FourOscSection *fourOscSectionToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(fourOscSectionToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedFourOscSection = fourOscSectionToWrap;
  moduleName = juce::String("FourOscillatorSection");

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String("/FourOscSectionPresets") );

  osc1Module = new OscillatorStereoAudioModule(lock, &wrappedFourOscSection->osc1);
  osc1Module->setModuleName(juce::String("Osc1"));
  addChildAudioModule(osc1Module);

  osc2Module = new OscillatorStereoAudioModule(lock, &wrappedFourOscSection->osc2);
  osc2Module->setModuleName(juce::String("Osc2"));
  addChildAudioModule(osc2Module);

  osc3Module = new OscillatorStereoAudioModule(lock, &wrappedFourOscSection->osc3);
  osc3Module->setModuleName(juce::String("Osc3"));
  addChildAudioModule(osc3Module);

  osc4Module = new OscillatorStereoAudioModule(lock, &wrappedFourOscSection->osc4);
  osc4Module->setModuleName(juce::String("Osc4"));
  addChildAudioModule(osc4Module);
}
