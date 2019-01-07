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
  Param* p;

  ratioProfileTopLeft = p = new Param("FreqProfileTopLeft", 0, 2, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileTopRight", 0, 2, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileBottomLeft", 0, 2, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("FreqProfileBottomRight", 0, 2, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  addObservedParameter(p);

  freqRatiosX = p = new Param("FreqRatiosX", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(core, &MS::setFreqRatioMixX);
  addObservedParameter(p);

  freqRatiosY = p = new Param("FreqRatiosY", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(core, &MS::setFreqRatioMixY);
  addObservedParameter(p);
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