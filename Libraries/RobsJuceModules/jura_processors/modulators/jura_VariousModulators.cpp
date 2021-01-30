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

//=================================================================================================

AttackDecayEnvelopeModule::AttackDecayEnvelopeModule(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("AttackDecayEnvelope");  // maybe use EnvelopeAD or EnvelopeAD_Mono
  createParameters();
}

void AttackDecayEnvelopeModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  //typedef RAPT::rsAttackDecayEnvelope<double> ADE;
  typedef AttackDecayEnvelopeModule ADM;
  typedef ModulatableParameter Param;
  //typedef ModulatableParameterPoly Param;  // for later use
  Param* p;

  p = new Param("Attack", 0.1, 100.0, 10.0, Parameter::EXPONENTIAL); 
  addObservedParameter(p);
  p->setValueChangeCallback<ADM>(this, &ADM::setAttack); 

  p = new Param("Decay", 10., 1000.0, 100.0, Parameter::EXPONENTIAL); 
  addObservedParameter(p);
  p->setValueChangeCallback<ADM>(this, &ADM::setDecay);
}

void AttackDecayEnvelopeModule::setSampleRate(double newSampleRate)
{ 
  sampleRate = newSampleRate;
  core.setAttackSamples(0.001 * attack * sampleRate);
  core.setDecaySamples( 0.001 * decay  * sampleRate);
}

void AttackDecayEnvelopeModule::setAttack(double newAttack)
{
  attack = newAttack;
  core.setAttackSamples(0.001 * attack * sampleRate);
}

void AttackDecayEnvelopeModule::setDecay(double newDecay)
{
  decay = newDecay;
  core.setDecaySamples(0.001 * decay * sampleRate);
}

//=================================================================================================

AttackDecayEnvelopeModulePoly::AttackDecayEnvelopeModulePoly(CriticalSection* lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse) 
  : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("EnvelopeAD");  //
  createParameters();
}

void AttackDecayEnvelopeModulePoly::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef AttackDecayEnvelopeModulePoly ADM;
  typedef ModulatableParameterPoly Param;
  Param* p;

  p = new Param("Attack", 0.1, 100.0, 10.0, Parameter::EXPONENTIAL); 
  addObservedParameter(p);
  p->setValueChangeCallbackPoly([this](double v, int i) { setAttack(v, i); });

  p = new Param("Decay", 10., 1000.0, 100.0, Parameter::EXPONENTIAL); 
  addObservedParameter(p);
  p->setValueChangeCallbackPoly([this](double v, int i) { setDecay(v, i); });
}

void AttackDecayEnvelopeModulePoly::allocateVoiceModResources() 
{
  if(voiceManager)
    cores.resize(voiceManager->getMaxNumVoices());
  else
    cores.resize(1); // monophonic in absence of a voice manager
}

void AttackDecayEnvelopeModulePoly::handleMidiMessage(MidiMessage msg)
{
  if(msg.isNoteOn())
  {
    // We have to trigger the envelope for one of the voices - but which? we need to figure out
    // which voice was assigned to this note


    int dummy = 0;
  }
}

void AttackDecayEnvelopeModulePoly::setAttack(double newAttack, int voice)
{
  jassert(voice < cores.size());
  cores[voice].setAttackSamples(0.001 * newAttack * sampleRate);
}

void AttackDecayEnvelopeModulePoly::setDecay(double newDecay, int voice)
{
  jassert(voice < cores.size());
  cores[voice].setDecaySamples(0.001 * newDecay * sampleRate);
}
