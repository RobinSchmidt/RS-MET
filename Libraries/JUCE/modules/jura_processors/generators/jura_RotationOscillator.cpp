RotationOscillatorAudioModule::RotationOscillatorAudioModule(CriticalSection *lockToUse)
  : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "RotOsc3D";
  setActiveDirectory(getApplicationDirectory() + "/Presets/RotOsc3D");

  createParameters();
}

void RotationOscillatorAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef RAPT::rsRotationOscillator<double> RO;
  RO* ro = oscCore;

  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("FreqScaleX", 0.0, 10.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<RO>(ro, &RO::setFreqScaleX);



}