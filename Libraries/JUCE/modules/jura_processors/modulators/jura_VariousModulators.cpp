TriSawModulatorModule::TriSawModulatorModule(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  createParameters();
}

void TriSawModulatorModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsTriSawModulator TSM;
  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Attack", 0.1, 1000.0, 0.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setAttackTime);

  p = new Param("Decay", 0.1, 1000.0, 0.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setDecayTime);

}