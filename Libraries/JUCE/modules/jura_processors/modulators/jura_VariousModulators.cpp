TriSawModulatorModule::TriSawModulatorModule(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("TriSawModulator");
  createParameters();
}

void TriSawModulatorModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsTriSawModulator TSM;
  typedef ModulatableParameter Param;
  Param* p;


  p = new Param("TimeScale", 0.01, 100.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setTimeScaler);

  p = new Param("Attack", 0.1, 1000.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setAttackTime);

  p = new Param("Decay", 0.1, 1000.0, 100.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setDecayTime);

  // maybe use another rnage and/or mapping for floor/ceil:
  p = new Param("Floor", -2, +2, -1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setFloor);

  p = new Param("Ceiling", -2, +2, +1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setCeiling);


  p = new Param("AttackBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setAttackBending);

  p = new Param("DecayBending", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setDecayBending);

  p = new Param("AttackSigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setAttackSigmoid);

  p = new Param("DecaySigmoid", -1.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<TSM>(&core, &TSM::setDecaySigmoid);
}