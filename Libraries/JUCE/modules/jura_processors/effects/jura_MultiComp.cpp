MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "MultiComp";
  setActiveDirectory(getApplicationDirectory() + "/Presets/MultiComp");
  maxNumBands = multiCompCore.getMaxNumberOfBands();
  createParameters();
}

void MultiCompAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsMultiBandCompressor MBC;
  MBC* mbc = &multiCompCore;

  //typedef ModulatableParameter Param;
  typedef Parameter Param;
  Param* p;

  p = new Param("NumBands", 1.0, maxNumBands, 16.0, Parameter::INTEGER); // use 1 as default later
  addObservedParameter(p);
  p->setValueChangeCallback<MBC>(mbc, &MBC::setNumberOfBands);

  p = new Param("SelectedBand", 1.0, maxNumBands, 1.0, Parameter::INTEGER);
  addObservedParameter(p);
  p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::selectBand);

  // SplitMode - steep lowpass, steep highpass, binary tree

  // create per-band parameters:
  for(int i = 0; i < maxNumBands; i++)
  {
    // todo use std::function with lambda functions for callbacks

    juce::String idxStr = juce::String(i+1);

    p = new Param("SplitFrequency" + idxStr, 20.0, 20000.0, 0.0, Parameter::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::setSplitFreq);

    p = new Param("Threshold" + idxStr, -60.0, 0.0, 0.0, Parameter::LINEAR); // in dB
    addObservedParameter(p);
    p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::setThreshold);

    p = new Param("Ratio" + idxStr, 1.0, 100.0, 1.0, Parameter::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::setRatio);

    p = new Param("Attack" + idxStr, 0.1, 1000.0, 10.0, Parameter::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::setAttack);

    p = new Param("Release" + idxStr, 0.1, 1000.0, 100.0, Parameter::EXPONENTIAL);
    addObservedParameter(p);
    p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::setRelease);
  }
}

void MultiCompAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    multiCompCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void MultiCompAudioModule::setSampleRate(double newSampleRate)
{
  multiCompCore.setSampleRate(newSampleRate);
}

void MultiCompAudioModule::reset()
{
  multiCompCore.reset();
}

AudioModuleEditor* MultiCompAudioModule::createEditor()
{
  return new MultiCompModuleEditor(this);
}

//=================================================================================================

MultiCompModuleEditor::MultiCompModuleEditor(MultiCompAudioModule* multiCompModuleToEdit)
  : AudioModuleEditor(multiCompModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  multiCompModule = multiCompModuleToEdit;
  setHeadlineText("MultiComp");
  createWidgets();
  setSize(400, 300);
}

void MultiCompModuleEditor::createWidgets()
{

}

void MultiCompModuleEditor::resized()
{

}

void MultiCompModuleEditor::updateWidgetVisibility()
{

}