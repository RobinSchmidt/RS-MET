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

  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("NumBands", 1.0, maxNumBands, 1.0, Parameter::INTEGER);
  addObservedParameter(p);
  p->setValueChangeCallback<MBC>(mbc, &MBC::setNumberOfBands);



  // SplitMode - steep lowpass, steep highpass, binary tree

  for(int i = 0; i < maxNumBands; i++)
  {
    // let the lowest band formally have a split-frequency of zero, so we can have the same number
    // of splitters as we have compressors ...or maybe not

  }
}

void MultiCompAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  //for(int n = 0; n < numSamples; n++)
  //  multiCompCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
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