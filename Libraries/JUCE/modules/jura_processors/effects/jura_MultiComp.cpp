MultiBandEffect::MultiBandEffect(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  //ScopedLock scopedLock(*lock);
}

void MultiBandEffect::setEffectCore(rosic::rsMultiBandEffect* effectCore)
{
  ScopedLock scopedLock(*lock);
  core = effectCore;
  maxNumBands = core->getMaxNumberOfBands();
  createSplittingParameters();
}

void MultiBandEffect::createSplittingParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsMultiBandEffect MBE;
  MBE* mbe = core;

  typedef Parameter Param;
  Param* p;

  p = new Param("NumBands", 1.0, maxNumBands, 1.0, Parameter::INTEGER); // use 1 as default later
  addObservedParameter(p);
  p->setValueChangeCallback<MultiBandEffect>(this, &MultiBandEffect::setNumBands);

  p = new Param("SelectedBand", 0.0, maxNumBands-1, 0.0, Parameter::STRING);
  p->addNumericStringValues(1, 16);
  addObservedParameter(p);
  p->setValueChangeCallback<MultiBandEffect>(this, &MultiBandEffect::selectBand);

  p = new Param("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Lowpass");
  p->addStringValue("Steep Highpass");
  //p->addStringValue("Binary Tree"); // doesn't work yet
  addObservedParameter(p);
  p->setValueChangeCallback<MBE>(mbe, &MBE::setSplitMode);


  // create per-band parameters:
  for(int i = 0; i < maxNumBands; i++)
  {
    // todo use std::function with lambda functions for callbacks because when doing it like below
    // all parameters use the same callback (and the dispatch is done based on our selectedBand 
    // member) - but this is not suitable for allowing automation

    juce::String idxStr = juce::String(i+1);

    p = new Param("SplitFrequency" + idxStr, 20.0, 20000.0, 0.0, Parameter::EXPONENTIAL);
    p->setValue(mbe->getSplitFrequency(i), false, false);
    addObservedParameter(p);
    p->setValueChangeCallback<MultiBandEffect>(this, &MultiBandEffect::setSplitFreq);

    splitFreqParams.push_back(p);
  }

  getParameterByName("SelectedBand")->setValue(0, true, true); // initially select band 1
}

void MultiBandEffect::parameterChanged(Parameter* p)
{
  ModulatableAudioModule::parameterChanged(p);
  sendChangeMessage();
}

void MultiBandEffect::setNumBands(int newNumBands)
{
  core->setNumberOfBands(newNumBands);
  jassert(areBandsInIncreasingOrder(false)); // for debug
}

void MultiBandEffect::setSplitFreq(int bandIndex, double newFreq)
{
  // we limit the new splitFreq such that bands remain ordered with ascending frequencies
  int numBands = getNumBands();
  if(bandIndex < getNumBands()-1) {
    double freqLimit = core->getSplitFrequency(bandIndex+1); // freq of right neighbour
    if(newFreq > freqLimit)
      getSplitFreqParam(bandIndex)->setValue(freqLimit, true, true);
    else
      core->setSplitFrequency(bandIndex, newFreq);
  }
  else
    core->setSplitFrequency(bandIndex, newFreq); // topmost band has no limit

  // somehow, we must also make sure sure that the topmost band has a freq of 20000 (it's the 
  // highpass band and actually has not lowpass cutoff at all)

  jassert(areBandsInIncreasingOrder(false)); // for debug
}

void MultiBandEffect::setSplitFreq(double newFreq) 
{ 
  setSplitFreq(selectedBand, newFreq);
}

Parameter* MultiBandEffect::getSplitFreqParam(int bandIndex)
{
  return splitFreqParams[bandIndex];
}

int MultiBandEffect::getBandContainingFrequency(double freq)
{
  for(int i = 0; i < core->getNumberOfBands()-1; i++)
    if(core->getSplitFrequency(i) > freq)
      return i;
  return core->getNumberOfBands()-1;
}

bool isLess(double x, double y, bool strictly)
{
  if(strictly)
    return x < y;
  else
    return x <= y;
}

bool MultiBandEffect::areBandsInIncreasingOrder(bool strictly)
{
  for(int i = 1; i < core->getNumberOfBands()-1; i++)
    if(!isLess(core->getSplitFrequency(i-1), core->getSplitFrequency(i), strictly))
      return false;
  return true;
}

//=================================================================================================

MultiBandPlotEditor::MultiBandPlotEditor(jura::MultiBandEffect* moduleToEdit)
  : module(moduleToEdit)

{
  module->addChangeListener(this);
  core = module->getCore();

  freqRespPlot = new rsFunctionPlot;
  freqRespPlot->setupForDecibelsAgainstLogFrequency(15.625, 32000.0, -48.0, 12.0, 6);

  for(int i = 0; i < core->getMaxNumberOfBands(); i++)
    freqRespPlot->addFunction([=](double f)->double { return core->getDecibelsAt(i, f); });

  freqRespPlot->addMouseListener(this, true);
  addPlot(freqRespPlot);
}

MultiBandPlotEditor::~MultiBandPlotEditor()
{
  module->removeChangeListener(this);
}

void MultiBandPlotEditor::changeListenerCallback(ChangeBroadcaster* source)
{
  freqRespPlot->setNumFunctionsToPlot(core->getNumberOfBands());
  repaint();
}

void MultiBandPlotEditor::mouseDown(const MouseEvent& e)
{
  // select band whose rectangle contains the mouse-event:
  double freq = freqRespPlot->fromPixelX(e.x);
  int index = module->getBandContainingFrequency(freq);
  Parameter* p = module->getParameterByName("SelectedBand");
  p->setValue(index, true, true);

  // maybe de-select, if the click was in the currently selected band, so we can have no selection
}

void MultiBandPlotEditor::paintOverChildren(Graphics& g)
{
  int numBands = core->getNumberOfBands();

  // highlight rectangle of selected band:
  int selected = module->getSelectedBand();
  float x1 = 0.f;
  float x2 = (float) getWidth();
  if(selected > 0)
    x1 = (float)freqRespPlot->toPixelX(core->getSplitFrequency(selected-1));
  if(selected < numBands-1)
    x2 = (float)freqRespPlot->toPixelX(core->getSplitFrequency(selected));
  //g.setColour(Colours::red.withAlpha(0.25f));
  g.setColour(Colours::lightblue.withAlpha(0.25f)); // preliminary
  //g.setColour(Colours::magenta.withAlpha(0.25f));
  g.fillRect(x1, 0.f, x2-x1, (float)getHeight());


  // draw vertical lines at split frequencies:
  float y1 = 0.f;
  float y2 = (float) getHeight();
  g.setColour(Colours::white); // preliminary
  for(int i = 0; i < numBands-1; i++)
  {
    //double freq = core->getSplitFrequency(i); // debug
    float x = (float)freqRespPlot->toPixelX(core->getSplitFrequency(i));
    g.drawLine(x, y1, x, y2, 2.f);
  }

  // todo: maybe give all bands a color (going through the rainbow) and just use a higher alpha
  // for the selected band (that may look nice)
}

void MultiBandPlotEditor::resized()
{
  freqRespPlot->setBounds(0, 0, getWidth(), getHeight());
}

//=================================================================================================

MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : MultiBandEffect(lockToUse,  metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("MultiComp");
  MultiBandEffect::setEffectCore(&multiCompCore);
  createCompressionParameters();
}

void MultiCompAudioModule::createCompressionParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsMultiBandCompressor MBC;
  MBC* mbc = &multiCompCore;

  //typedef ModulatableParameter Param;
  typedef Parameter Param;
  Param* p;

  // create per-band parameters:
  for(int i = 0; i < maxNumBands; i++)
  {
    // todo use std::function with lambda functions for callbacks

    juce::String idxStr = juce::String(i+1);

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

  //getParameterByName("SelectedBand")->setValue(0, true, true); // initially select band 1
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
  : MultiBandPlotEditor(multiCompModuleToEdit), multiCompModule(multiCompModuleToEdit)

{

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