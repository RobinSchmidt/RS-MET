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

  p = new Param("SelectedBand", 1.0, maxNumBands, 1.0, Parameter::STRING);
  p->addStringValue("1"); // have a function: addNumericStringValues(int min, int max, int stepSize)
  p->addStringValue("2"); //...
  addObservedParameter(p);
  p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::selectBand);

  p = new Param("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Lowpass");
  p->addStringValue("Steep Highpass");
  p->addStringValue("Binary Tree");
  addObservedParameter(p);
  p->setValueChangeCallback<MBC>(mbc, &MBC::setSplitMode);

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
  RSlider *s;
  RComboBox *c;

  addWidget( s = numBandsSlider = new RSlider );
  s->assignParameter( multiCompModule->getParameterByName("NumBands") );
  s->setSliderName("NumBands");
  s->setDescription("Number of frequency bands");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);

  addWidget( c = bandSelectBox = new RComboBox() );
  c->assignParameter( multiCompModule->getParameterByName("SelectedBand") );
  c->setDescription("Select band to edit");
  c->setDescriptionField(infoField);
  c->registerComboBoxObserver(this);

  addWidget( c = splitModeBox = new RComboBox() );
  c->assignParameter( multiCompModule->getParameterByName("SplitMode") );
  c->setDescription("Mode of the band-splitting");
  c->setDescriptionField(infoField);

  // per band widgets:
  int maxNumBands = multiCompModule->getMaxNumBands();
  splitFreqSliders.resize(maxNumBands);
  thresholdSliders.resize(maxNumBands);
  ratioSliders.resize(maxNumBands);
  attackSliders.resize(maxNumBands);
  releaseSliders.resize(maxNumBands);
  for(int k = 0; k < multiCompModule->getMaxNumBands(); k++)
  {
    juce::String idxStr = juce::String(k+1);

    addWidget( s = splitFreqSliders[k] = new RSlider );
    s->assignParameter( multiCompModule->getParameterByName("SplitFrequency" + idxStr) );

    addWidget( s = thresholdSliders[k] = new RSlider );
    s->assignParameter( multiCompModule->getParameterByName("Threshold" + idxStr) );

    addWidget( s = ratioSliders[k] = new RSlider );
    s->assignParameter( multiCompModule->getParameterByName("Ratio" + idxStr) );

    addWidget( s = attackSliders[k] = new RSlider );
    s->assignParameter( multiCompModule->getParameterByName("Attack" + idxStr) );

    addWidget( s = releaseSliders[k] = new RSlider );
    s->assignParameter( multiCompModule->getParameterByName("Release" + idxStr) );
  }
}

void MultiCompModuleEditor::rComboBoxChanged(RComboBox* box)
{
  updateWidgetVisibility();
}

void MultiCompModuleEditor::resized()
{

}

void MultiCompModuleEditor::updateWidgetVisibility()
{

}