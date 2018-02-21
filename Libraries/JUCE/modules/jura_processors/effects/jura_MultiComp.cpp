MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("MultiComp");
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

  p = new Param("NumBands", 1.0, maxNumBands, 1.0, Parameter::INTEGER); // use 1 as default later
  addObservedParameter(p);
  p->setValueChangeCallback<MBC>(mbc, &MBC::setNumberOfBands);

  p = new Param("SelectedBand", 0.0, maxNumBands-1, 0.0, Parameter::STRING);
  p->addNumericStringValues(1, 16);
  addObservedParameter(p);
  p->setValueChangeCallback<MultiCompAudioModule>(this, &MultiCompAudioModule::selectBand);

  p = new Param("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Lowpass");
  p->addStringValue("Steep Highpass");
  //p->addStringValue("Binary Tree"); // doesn't work yet
  addObservedParameter(p);
  p->setValueChangeCallback<MBC>(mbc, &MBC::setSplitMode);

  // create per-band parameters:
  for(int i = 0; i < maxNumBands; i++)
  {
    // todo use std::function with lambda functions for callbacks

    juce::String idxStr = juce::String(i+1);

    p = new Param("SplitFrequency" + idxStr, 20.0, 20000.0, 0.0, Parameter::EXPONENTIAL);
    p->setValue(mbc->getSplitFrequency(i), false, false);
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

  getParameterByName("SelectedBand")->setValue(0, true, true); // initially select band 1
}

void MultiCompAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    multiCompCore.getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void MultiCompAudioModule::processStereoFrame(double *left, double *right)
{
  multiCompCore.getSampleFrameStereo(left, right);
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

MultiCompPlotEditor::MultiCompPlotEditor(jura::MultiCompAudioModule* multiCompModuleToEdit)
  : multiCompModule(multiCompModuleToEdit)

{
  multiCompCore = multiCompModule->getCore();

  freqRespPlot = new rsFunctionPlot;
  freqRespPlot->setupForDecibelsAgainstLogFrequency(15.625, 32000.0, -48.0, 12.0, 6);
  //for(int i = 0; i < multiCompModule->getMaxNumBands; i++)
  //  freqRespPlot->addFunction([this](double f)->double { return multiCompModule->getDecibelsAt(i, f); });
  freqRespPlot->addMouseListener(this, true);
  addPlot(freqRespPlot);
}

void MultiCompPlotEditor::mouseDown(const MouseEvent& e)
{
  // select band whose rectangle contains the mouse-event
  repaint(); // test
}

void MultiCompPlotEditor::paintOverChildren(Graphics& g)
{
  // draw vertical lines at split frequencies:
  float y1 = 0.f;
  float y2 = (float) getHeight();
  int numBands = multiCompCore->getNumberOfBands(); // for debug 
  g.setColour(Colours::white);
  for(int i = 0; i < multiCompCore->getNumberOfBands(); i++)
  {
    float x = (float)freqRespPlot->toPixelX(multiCompCore->getSplitFrequency(i));
    g.drawLine(x, y1, x, y2, 2.f);
  }

  // highlight rectangle of selected band:

}

void MultiCompPlotEditor::resized()
{
  freqRespPlot->setBounds(0, 0, getWidth(), getHeight());
}

//=================================================================================================

MultiCompModuleEditor::MultiCompModuleEditor(MultiCompAudioModule* multiCompModuleToEdit)
  : AudioModuleEditor(multiCompModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  multiCompModule = multiCompModuleToEdit;
  //setHeadlineText("MultiComp");
  createWidgets();
  updateWidgetVisibility();
  setSize(595, 301);
}

void MultiCompModuleEditor::createWidgets()
{
  RSlider *s;
  RComboBox *c;

  plotEditor = new MultiCompPlotEditor(multiCompModule);
  addChildColourSchemeComponent(plotEditor);

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
  AudioModuleEditor::resized();

  int y = getPresetSectionBottom() + 4;
  plotEditor->setBounds(0, y, getWidth(), getHeight()-110);
  y = plotEditor->getBottom() + 4;

  int x = 4;
  int w = getWidth() / 2 - 8;
  int h = 16;
  int d = h-2;

  numBandsSlider->setBounds(x, y, w, h); y += d;
  bandSelectBox ->setBounds(x, y, w, h); y += d;
  splitModeBox  ->setBounds(x, y, w, h); y += d;

  x = getWidth() / 2 + 4;
  for(int k = 0; k < multiCompModule->getMaxNumBands(); k++)
  {
    //y = getPresetSectionBottom() + 4;
    y = plotEditor->getBottom() + 4;
    splitFreqSliders[k]->setBounds(x, y, w, h); y += d;
    thresholdSliders[k]->setBounds(x, y, w, h); y += d;
    ratioSliders[k]    ->setBounds(x, y, w, h); y += d;
    attackSliders[k]   ->setBounds(x, y, w, h); y += d;
    releaseSliders[k]  ->setBounds(x, y, w, h); y += d;
  }
}

void MultiCompModuleEditor::updateWidgetVisibility()
{
  int k;
  for(k = 0; k < multiCompModule->getMaxNumBands(); k++)
  {
    splitFreqSliders[k]->setVisible(false);
    thresholdSliders[k]->setVisible(false);
    ratioSliders[k]    ->setVisible(false);
    attackSliders[k]   ->setVisible(false);
    releaseSliders[k]  ->setVisible(false);
  }
  k = multiCompModule->getSelectedBand();
  splitFreqSliders[k]->setVisible(true);
  thresholdSliders[k]->setVisible(true);
  ratioSliders[k]    ->setVisible(true);
  attackSliders[k]   ->setVisible(true);
  releaseSliders[k]  ->setVisible(true);
}