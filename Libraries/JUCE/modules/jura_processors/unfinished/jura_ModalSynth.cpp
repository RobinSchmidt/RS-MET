ModalSynthAudioModule::ModalSynthAudioModule(CriticalSection *lockToUse)
  : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("ModalSynth"); // find a better name - maybe Modacous?
  createParameters();
}

void ModalSynthAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsModalSynth MS;
  typedef Parameter Param;
  typedef ParameterWithKeyVelScaling ParamKV;
  Param* p;
  ParamKV* pkv;

  //int numOptions = 0;  // number of string parameters - can we generally pass 0?

  ratioProfileTopLeft = p = new Param("FreqProfileTopLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile1);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileTopRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile2);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileBottomLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile3);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileBottomRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile4);
  addObservedParameter(p);

  //maxNumModes = p = new Param("MaxNumModes", 10.0, 1024.0, 1024.0, Parameter::EXPONENTIAL, 1.0);
  maxNumModes = p = new Param("MaxNumModes", 10.0, 1024.0, 32.0, Parameter::EXPONENTIAL, 1.0);
  p->setValueChangeCallback<MS>(&core, &MS::setMaxNumPartials);
  addObservedParameter(p);

  freqRatiosX = p = new Param("FreqRatiosX", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioMixX);
  addObservedParameter(p);

  freqRatiosY = p = new Param("FreqRatiosY", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioMixY);
  addObservedParameter(p);

  //attack = pkv = new ParamKV("Attack", 0.1, 100.0, 5.0, Parameter::EXPONENTIAL, 0.0);
  attack = pkv = new ParamKV("Attack", 0.0, 100.0, 5.0, Parameter::LINEAR, 0.0);
  pkv->setKeyScaleRange(-100, 100, 0);
  pkv->setVelScaleRange(-100, 100, 0);
  pkv->setValueChangeCallback<MS>(&core, &MS::setAttack);
  //pkv->setKeyScaleCallback<MS>(&core, &MS::setAttackByKey);
  //pkv->setVelScaleCallback<MS>(&core, &MS::setAttackByVel);
  addObservedParameter(pkv);

  attack = pkv = new ParamKV("Decay", 10.0, 10000.0, 5.0, Parameter::EXPONENTIAL, 0.0);
  pkv->setKeyScaleRange(-100, 100, 0);
  pkv->setVelScaleRange(-100, 100, 0);
  pkv->setValueChangeCallback<MS>(&core, &MS::setDecay);
  //pkv->setKeyScaleCallback<MS>(&core, &MS::setDecayByKey);
  //pkv->setVelScaleCallback<MS>(&core, &MS::setDecayByVel);
  addObservedParameter(pkv);

}

/*
AudioModuleEditor* ModalSynthAudioModule::createEditor(int type)
{
  return new ModalSynthEditor(this);
}
*/

void ModalSynthAudioModule::processBlock(double **buf, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    core.getSampleFrameStereo(&buf[0][n], &buf[1][n]);
}

void ModalSynthAudioModule::processStereoFrame(double *left, double *right)
{
  core.getSampleFrameStereo(left, right);
}

void ModalSynthAudioModule::setSampleRate(double newSampleRate)
{
  core.setSampleRate(newSampleRate);
}

void ModalSynthAudioModule::reset()
{
  core.reset();
}

void ModalSynthAudioModule::populateFreqRatioProfileParam(Parameter* p)
{
  // make sure that the order is the same as in the enumerations in rsModalSynth

  // 0-parametric settings:
  p->addStringValue("All Harmonics");
  p->addStringValue("Odd Harmonics");
  p->addStringValue("Pseudo Harmonic 12-TET");
  p->addStringValue("12-TET");
  p->addStringValue("Rod Free/Free");
  p->addStringValue("Rod Free/Clamped");
}

//=================================================================================================

ModalSynthEditor::ModalSynthEditor(jura::ModalSynthAudioModule *newModalSynthModuleToEdit)
  : AudioModuleEditor(newModalSynthModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  modalModule = newModalSynthModuleToEdit;
  createWidgets();
  setSize(500, 460);
}

void ModalSynthEditor::resized()
{

}

void ModalSynthEditor::createWidgets()
{

}