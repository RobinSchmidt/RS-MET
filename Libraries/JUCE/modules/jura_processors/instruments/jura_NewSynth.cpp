// construction/destruction:

NewSynthAudioModule::NewSynthAudioModule(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse)
{
  setModuleTypeName("NewSynth");
  //createParameters();
}

//=================================================================================================