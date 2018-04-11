MultiBandEffect::MultiBandEffect(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
  , perBandModuleFactory(lockToUse)
{
  ScopedLock scopedLock(*lock);

  // register the module types that can be used per band (maybe factor out):
  typedef CriticalSection* CS;
  typedef AudioModule* AM;
  CS cs = lock;
  AudioModuleFactory& f = perBandModuleFactory;

  juce::String s = "";
  f.registerModuleType([](CS cs)->AM { return new CompressorAudioModule(cs); }, s, "Compressor");
  //f.registerModuleType([](CS cs)->AM { return new FuncShaperAudioModule(cs); }, s, "FuncShaper");
  // ...

  setEffectType("Compressor");
}

void MultiBandEffect::setEffectCore(rosic::rsMultiBandEffect* effectCore)
{
  ScopedLock scopedLock(*lock);
  core = effectCore;

  // doing this here is a bit dirty:
  Parameter* p = new Parameter("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Lowpass");
  p->addStringValue("Steep Highpass");
  //p->addStringValue("Binary Tree"); // doesn't work yet
  addObservedParameter(p);
  p->setValueChangeCallback<rosic::rsMultiBandEffect>(
    core, &rosic::rsMultiBandEffect::setSplitMode);

  createSplitFreqParams();
}

void MultiBandEffect::setEffectType(const juce::String& typeString)
{
  if(typeString != effectTypeString)
  {
    effectTypeString = typeString;
    // todo: 
    // -replace modules in all pre-band slots with the new selected type
    // -notify gui about that change, so it can replace the editors
  }
}

void MultiBandEffect::parameterChanged(Parameter* p)
{
  ModulatableAudioModule::parameterChanged(p);
  sendChangeMessage();
}

void MultiBandEffect::insertBand(int index, double splitFrequency, bool sendNotification)
{
  core->insertBand(index, splitFrequency);
  addSplitFreqParam(index, splitFrequency); // rename to insertSplitFreqParam
  insertBandEffect(index);
  if(sendNotification)
    sendBandInsertNotification(index);  // notify gui (for creating widgets)
}

void MultiBandEffect::removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification)
{
  if(sendNotification)
    sendBandRemoveNotification(index);  // notify gui (for removing widgets)
  removeSplitFreqParam(index);
  removeBandEffect(index);
  core->removeBand(index);
}

void MultiBandEffect::selectBand(int bandToSelect, bool sendNotification) 
{ 
  selectedBand = bandToSelect; 
  if(sendNotification)
    sendBandSelectNotification(selectedBand);
}

void MultiBandEffect::createSplitFreqParams()
{
  // clearSlpitFreqParams(); // ...maybe later
  size_t numParams = splitFreqParams.size();
  size_t numBands = getNumBands();
  while(numParams < numBands) {
    addSplitFreqParam((int)numParams, getSplitFreq((int)numParams));
    numParams++; }

  jassert(areBandsInIncreasingOrder(false)); // for debug
}

void MultiBandEffect::addSplitFreqParam(int index, double freq)
{
  juce::String idxStr = juce::String(index+1);
  Parameter* p = new Parameter("SplitFrequency" + idxStr, 20.0, 20000.0, 0.0, 
    Parameter::EXPONENTIAL);
  p->setValue(freq, false, false);
  addObservedParameter(p);

  p->setValueChangeCallback<MultiBandEffect>(this, &MultiBandEffect::setSplitFreq);
   // i think, this is wrong - we need a lambda here

  insert(splitFreqParams, p, index);
}

void MultiBandEffect::removeSplitFreqParam(int i)
{
  delete splitFreqParams[i];
  remove(splitFreqParams, i);
}

void MultiBandEffect::setSplitFreq(int bandIndex, double newFreq)
{
  if(bandIndex < 0)  return;  // kludge

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

std::vector<juce::String> MultiBandEffect::getAvailableEffectTypes()
{
  return perBandModuleFactory.getRegisteredModuleTypes();
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

void MultiBandEffect::sendBandInsertNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWasInserted(this, index);
}

void MultiBandEffect::sendBandRemoveNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWillBeRemoved(this, index);
}

void MultiBandEffect::sendBandSelectNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWasSelected(this, index);
}

void MultiBandEffect::insertBandEffect(int i)
{
  AudioModule* m = perBandModuleFactory.createModule(effectTypeString);
  insert(perBandModules, m, i);
}

void MultiBandEffect::removeBandEffect(int i)
{
  delete perBandModules[i];
  remove(perBandModules, i);
}

//=================================================================================================

MultiBandPlotEditor::MultiBandPlotEditor(jura::MultiBandEffect* moduleToEdit)
  : module(moduleToEdit)

{
  module->addChangeListener(this); // obsolete?
  module->registerMultiBandObserver(this);

  core = module->getCore();

  freqRespPlot = new rsFunctionPlot;
  freqRespPlot->setupForDecibelsAgainstLogFrequency(15.625, 32000.0, -48.0, 12.0, 6);

  //for(int i = 0; i < core->getMaxNumberOfBands(); i++)
  //  freqRespPlot->addFunction([=](double f)->double { return core->getDecibelsAt(i, f); });

  freqRespPlot->addMouseListener(this, true);
  addPlot(freqRespPlot);
}

MultiBandPlotEditor::~MultiBandPlotEditor()
{
  module->removeChangeListener(this); // obsolete?
  module->registerMultiBandObserver(this);
  delete bandPopup;
}

void MultiBandPlotEditor::bandWasInserted(MultiBandEffect* mbe, int i)
{
  // insert new graph into plot
  //freqRespPlot->addFunction([=](double f)->double { return core->getDecibelsAt(i, f); });
  //freqRespPlot->insertFunction(i, [=](double f)->double { return core->getDecibelsAt(i, f); });
}

void MultiBandPlotEditor::bandWillBeRemoved(MultiBandEffect* mbe, int i)
{
  //freqRespPlot->removeFunction(i);  // remove graph from plot
}

void MultiBandPlotEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{
  // update selection highlighting
  freqRespPlot->setNumFunctionsToPlot(core->getNumberOfBands());
  repaint();
}

void MultiBandPlotEditor::changeListenerCallback(ChangeBroadcaster* source)
{
  freqRespPlot->setNumFunctionsToPlot(core->getNumberOfBands());
  repaint();
}

void MultiBandPlotEditor::rPopUpMenuChanged(RPopUpMenu* menu)
{
  int index = module->getBandContainingFrequency(freqAtMouse);
  switch(bandPopup->getSelectedIdentifier())
  {
  case ADD_BAND:    module->insertBand(index, freqAtMouse, true);  break;
  case REMOVE_BAND: module->removeBand(index, false, true);        break;
  }
}

void MultiBandPlotEditor::mouseDown(const MouseEvent& e)
{
  freqAtMouse = freqRespPlot->fromPixelX(e.x); 

  if(e.mods.isLeftButtonDown()) {
    // select band whose rectangle contains the mouse-event:
    int index = module->getBandContainingFrequency(freqAtMouse);
    if(index == module->getSelectedBand())
      module->selectBand(-1, true); // de-select, when a band is clicked again
    else
      module->selectBand(index, true);

    // old (we don't have a Parameter for this anymore):
    //Parameter* p = module->getParameterByName("SelectedBand");
    //p->setValue(index, true, true);
  }
  else if(e.mods.isRightButtonDown())
    openRightClickMenu();
}

void MultiBandPlotEditor::paintOverChildren(Graphics& g)
{
  int numBands = core->getNumberOfBands();

  // highlight rectangle of selected band:
  int selected = module->getSelectedBand();
  if(selected >= 0)
  {
    float x1 = 0.f;
    float x2 = (float)getWidth();
    if(selected > 0)
      x1 = (float)freqRespPlot->toPixelX(core->getSplitFrequency(selected-1));
    if(selected < numBands-1)
      x2 = (float)freqRespPlot->toPixelX(core->getSplitFrequency(selected));
    //g.setColour(Colours::red.withAlpha(0.25f));
    g.setColour(Colours::lightblue.withAlpha(0.25f)); // preliminary
    //g.setColour(Colours::magenta.withAlpha(0.25f));
    g.fillRect(x1, 0.f, x2-x1, (float)getHeight());
  }


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

void MultiBandPlotEditor::openRightClickMenu()
{
  // create popup, if necessary:
  if(bandPopup == nullptr)
  {
    bandPopup = new RPopUpMenu(this);
    bandPopup->registerPopUpMenuObserver(this);
    bandPopup->setDismissOnFocusLoss(true);
    bandPopup->addItem(ADD_BAND,    "Add band");
    bandPopup->addItem(REMOVE_BAND, "Remove band");
  }

  // show it:
  int w = bandPopup->getRequiredWidth(true);
  int h = bandPopup->getRequiredHeight(true);
  bandPopup->selectItemByIndex(-1, false);    // select nothing - maybe write a specific fucntion for that
  bandPopup->showAtMousePosition(true, w, h); // showModally = true
}

//=================================================================================================

MultiBandEffectEditor::MultiBandEffectEditor(MultiBandEffect* effect) : AudioModuleEditor(effect)
{
  ScopedLock scopedLock(*lock);
  effectToEdit = effect;
  effectToEdit->registerMultiBandObserver(this);
  //setHeadlineText("MultiBandEffect");
  createWidgets();
  createBandEditors();
  updateEditorVisibility();
  setSize(595, 301);
}

MultiBandEffectEditor::~MultiBandEffectEditor()
{
  effectToEdit->deRegisterMultiBandObserver(this);
  clearBandEditors();
}

void MultiBandEffectEditor::resized()
{

}

void MultiBandEffectEditor::bandWasInserted(MultiBandEffect* mbe, int index)
{

}

void MultiBandEffectEditor::bandWillBeRemoved(MultiBandEffect* mbe, int index)
{

}

void MultiBandEffectEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{

}

void MultiBandEffectEditor::insertBandEditor(int index)
{

}

void MultiBandEffectEditor::removeBandEditor(int index)
{

}

void MultiBandEffectEditor::updateEditorVisibility()
{

}

void MultiBandEffectEditor::createWidgets()
{

}

void MultiBandEffectEditor::createBandEditors()
{

}

void MultiBandEffectEditor::clearBandEditors()
{

}










//=================================================================================================

MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : MultiBandEffect(lockToUse,  metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("MultiComp");
  MultiBandEffect::setEffectCore(&multiCompCore);
  createBandParams();
}

void MultiCompAudioModule::insertBand(int i, double splitFreq, bool sendNotification)
{
  MultiBandEffect::insertBand(i, splitFreq, false);
  addCompressionParams(i);
  if(sendNotification)  
    sendBandInsertNotification(i); 
}

void MultiCompAudioModule::removeBand(int i, bool mergeWithRightNeighbour, bool sendNotification)
{
  MultiBandEffect::removeBand(i, mergeWithRightNeighbour, sendNotification);
  removeCompressionParams(i);
}

void MultiCompAudioModule::createBandParams()
{
  // add a parameter set (freq, threshold, ratio, etc.) for each of the bands 
  createSplitFreqParams(); 
  while(numCompParamSets < getNumBands())
    addCompressionParams((int)numCompParamSets); 
  
  // is the index right? maybe we nee to check of each index, if a parameter set for that index 
  // exists? ...when removing bands the parameter indices are not updated...that is wrong...
}

void MultiCompAudioModule::addCompressionParams(int i)
{
  typedef Parameter Param;
  Param* p;

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

  numCompParamSets++;
}

void MultiCompAudioModule::removeCompressionParams(int i)
{
  juce::String idxStr = juce::String(i+1);
  removeParameter("Threshold" + idxStr, true);
  removeParameter("Ratio"     + idxStr, true);
  removeParameter("Attack"    + idxStr, true);
  removeParameter("Release"   + idxStr, true);

  numCompParamSets++;
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
  multiCompModule->registerMultiBandObserver(this);

  //setHeadlineText("MultiComp");
  createWidgets();
  updateWidgetVisibility();
  setSize(595, 301);
}

MultiCompModuleEditor::~MultiCompModuleEditor()
{
  multiCompModule->deRegisterMultiBandObserver(this);
}

void MultiCompModuleEditor::createWidgets()
{
  RComboBox *c;

  plotEditor = new MultiCompPlotEditor(multiCompModule);
  addChildColourSchemeComponent(plotEditor);

  /*
  addWidget( c = bandSelectBox = new RComboBox() );
  c->assignParameter( multiCompModule->getParameterByName("SelectedBand") );
  c->setDescription("Select band to edit");
  c->setDescriptionField(infoField);
  c->registerComboBoxObserver(this);
  */

  addWidget( c = splitModeBox = new RComboBox() );
  c->assignParameter( multiCompModule->getParameterByName("SplitMode") );
  c->setDescription("Mode of the band-splitting");
  c->setDescriptionField(infoField);

  createBandWidgets();
}

void MultiCompModuleEditor::createBandWidgets()
{
  size_t numWidgetSets = splitFreqSliders.size();
  size_t numBands = multiCompModule->getNumBands();
  while(numWidgetSets < numBands) {
    addBandWidgets((int)numWidgetSets);
    numWidgetSets++; 
  }
}

void MultiCompModuleEditor::rComboBoxChanged(RComboBox* box)
{
  updateWidgetVisibility();
}

void MultiCompModuleEditor::bandWasInserted(MultiBandEffect* mbe, int index)
{
  addBandWidgets(index);
}

void MultiCompModuleEditor::bandWillBeRemoved(MultiBandEffect* mbe, int index)
{
  removeBandWidgets(index);
}

void MultiCompModuleEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{
  updateWidgetVisibility();
}

void MultiCompModuleEditor::addBandWidgets(int i)
{
  juce::String idxStr = juce::String(i+1);
  RSlider *s;

  addWidget(s = new RSlider);
  insert(splitFreqSliders, s, i);
  s->assignParameter( multiCompModule->getParameterByName("SplitFrequency" + idxStr) );

  addWidget( s = new RSlider );
  insert(thresholdSliders, s, i);
  s->assignParameter( multiCompModule->getParameterByName("Threshold" + idxStr) );

  addWidget( s = new RSlider );
  insert(ratioSliders, s, i);
  s->assignParameter( multiCompModule->getParameterByName("Ratio" + idxStr) );

  addWidget( s = new RSlider );
  insert(attackSliders, s, i);
  s->assignParameter( multiCompModule->getParameterByName("Attack" + idxStr) );

  addWidget( s = new RSlider );
  insert(releaseSliders, s, i);
  s->assignParameter( multiCompModule->getParameterByName("Release" + idxStr) );
}

void MultiCompModuleEditor::removeBandWidgets(int i)
{
  remove(splitFreqSliders, i); removeWidget(splitFreqSliders[i], true, true);
  remove(thresholdSliders, i); removeWidget(thresholdSliders[i], true, true);
  remove(ratioSliders,     i); removeWidget(ratioSliders[i],     true, true);
  remove(attackSliders,    i); removeWidget(attackSliders[i],    true, true);
  remove(releaseSliders,   i); removeWidget(releaseSliders[i],   true, true);
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

  //numBandsSlider->setBounds(x, y, w, h); y += d;
  //bandSelectBox ->setBounds(x, y, w, h); y += d;
  splitModeBox  ->setBounds(x, y, w, h); y += d;

  x = getWidth() / 2 + 4;
  for(int k = 0; k < multiCompModule->getNumBands(); k++)
  {
    y = plotEditor->getBottom() + 4;
    splitFreqSliders[k]->setBounds(x, y, w, h); y += d; // maybe position them below the splitModeBox
    thresholdSliders[k]->setBounds(x, y, w, h); y += d;
    ratioSliders[k]    ->setBounds(x, y, w, h); y += d;
    attackSliders[k]   ->setBounds(x, y, w, h); y += d;
    releaseSliders[k]  ->setBounds(x, y, w, h); y += d;
  }
}

void MultiCompModuleEditor::updateWidgetVisibility()
{
  int k;
  for(k = 0; k < multiCompModule->getNumBands(); k++) {
    splitFreqSliders[k]->setVisible(false);
    thresholdSliders[k]->setVisible(false);
    ratioSliders[k]    ->setVisible(false);
    attackSliders[k]   ->setVisible(false);
    releaseSliders[k]  ->setVisible(false);
  }
  k = multiCompModule->getSelectedBand();
  if(k >= 0) {
    splitFreqSliders[k]->setVisible(true);
    thresholdSliders[k]->setVisible(true);
    ratioSliders[k]->setVisible(true);
    attackSliders[k]->setVisible(true);
    releaseSliders[k]->setVisible(true);
  }
}

/*
some thoughts:

MultiBandEffect baseclass:
-each band except the last has a split-freq parameter (for convenience, the last may also have one
 but that is just a dummy - and may need a dummy callback target)
-split-freq parameters are managed in the MultiBandEffect baseclass, this class is also responsible
 for keeping the constraint of ascending split-freqs
-when a band is inserted, all split-freq params above the position of the new band need to be 
 renamed (say from SplitFreq3 to SplitFreq4) and their callbacks must be updated (using a lamda 
 function that calls the setSplitFreq function of the core with an appropriate band index)
-when a band is removed, all split-freq params above the removed need to be renamed (splitFreq4 -> 
 SplitFreq3) and callbacks must be updated
-when a parameter is renamed because of insert/remove actions, its widget should reflect that 
 (respond parameterNameChanged callback?)
-the MultiBandPlot accesses the split-freq parameters (for setting and getting)

MultiComp subclass:
-each band has a set of compression parameters (Threshold, Ratio, Attack, Release)
-when a band is inserted or removed, similar renamings and callback reassignments as for the 
 split-freq parameters have to take place - for this, we override insertBand and removeBand, call 
 the baseclass code there and then do the additional updates
 -maybe these updates can be done in BandCompParameterSet::setBandIndex
 -or maybe a full single-band compressor AudioModule (with its own parameter set and editor) can be 
  created for each band compressor? ...that would be great!
 
-maybe the MultiBandEffect class should hold an array std::vector<AudioModule*> perBandModules
-MultiComp would then override insertBand/removeBand and create a new CompressorAudioModule for 
 each band
-the wrapped rosic::Compressor should not be owned by the jura::Compressor
-a class MultiBandEditor should show the sub-editor of the currently selected per-band module
-the renamings in inserBand/removeBand would then apply to whole modules
-it even could be delegated to MultiBandEffect baseclass: updateNames(String& effectName, int index)
-then, callbacks would not have to re-assigned

OR: templatize the MultiBandEffect class on the type of the effect such that MultiBandCompressor 
 just becomes a template instantiation MultiBandEffect<Compressor>







*/