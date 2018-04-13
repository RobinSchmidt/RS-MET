MultiBandEffect::MultiBandEffect(CriticalSection *lockToUse,
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
  , perBandModuleFactory(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("MultiBandEffect");

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


  // create parameters (factor out):
  Parameter* p = new Parameter("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Lowpass");
  p->addStringValue("Steep Highpass");
  //p->addStringValue("Binary Tree"); // doesn't work yet
  addObservedParameter(p);
  p->setValueChangeCallback<rosic::rsMultiBandEffect>(
    &core, &rosic::rsMultiBandEffect::setSplitMode);
  createSplitFreqParams();

  insertBandEffect(0);
}

MultiBandEffect::~MultiBandEffect()
{
  clearBandEffects();
}

void MultiBandEffect::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    processStereoFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  // todo: optimize to avoid calling processStereoFrame each sample...maybe use internal buffers 
  // and let the splitter/core fill them, then pass the buffers to the per-band modules for 
  // processing
}

void MultiBandEffect::processStereoFrame(double *left, double *right)
{
  jassert(perBandModules.size() == getNumBands());
  core.split(left, right);
  for(int k = 0; k < getNumBands(); k++)  // process individual bands
    perBandModules[k]->processStereoFrame(left, right);
  core.recombine(left, right);
}

void MultiBandEffect::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  core.setSampleRate(newSampleRate);
  for(size_t i = 0; i < perBandModules.size(); i++)
    perBandModules[i]->setSampleRate(newSampleRate);
}

void MultiBandEffect::reset()
{
  ScopedLock scopedLock(*lock);
  core.reset();
  for(size_t i = 0; i < perBandModules.size(); i++)
    perBandModules[i]->reset();;
}

AudioModuleEditor* MultiBandEffect::createEditor()
{
  return new MultiBandEffectEditor(this);
}

void MultiBandEffect::setEffectType(const juce::String& typeString)
{
  ScopedLock scopedLock(*lock);
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
  ScopedLock scopedLock(*lock);
  core.insertBand(index, splitFrequency);
  addSplitFreqParam(index, splitFrequency); // rename to insertSplitFreqParam
  insertBandEffect(index);
  if(sendNotification)
    sendBandInsertNotification(index);  // notify gui (for creating widgets)
}

void MultiBandEffect::removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification)
{
  ScopedLock scopedLock(*lock);

  if(sendNotification)
    sendBandRemovePreNotification(index);  // notify gui (for removing widgets)


  removeBandEffect(index);
  core.removeBand(index);
  removeSplitFreqParam(index); // must be done after core.removeBand


  if(sendNotification)
    sendBandRemovePostNotification(index);
}

void MultiBandEffect::selectBand(int bandToSelect, bool sendNotification) 
{ 
  ScopedLock scopedLock(*lock);
  selectedBand = bandToSelect; 
  if(sendNotification)
    sendBandSelectNotification(selectedBand);
}

void MultiBandEffect::setSplitFreq(int bandIndex, double newFreq)
{
  ScopedLock scopedLock(*lock);

  if(bandIndex < 0)  return;  // kludge

  // preliminary - no limits...but maybe keep it unlimited:
  core.setSplitFrequency(bandIndex, newFreq);
  sendChangeMessage();
  return;


  // we limit the new splitFreq such that bands remain ordered with ascending frequencies
  int numBands = getNumBands();
  if(bandIndex < getNumBands()-1) {

    double freqLimit = core.getSplitFrequency(bandIndex+1); // freq of right neighbour
      // ..that's the upper limit - what about the lower limit?

    //core.setSplitFrequency(bandIndex, newFreq); // test

    if(newFreq > freqLimit)
    {
      //getSplitFreqParam(bandIndex)->setValue(freqLimit, true, true); // stack overflow here?

      getSplitFreqParam(bandIndex)->setValue(freqLimit, true, false); // don't call callbacks
      core.setSplitFrequency(bandIndex, freqLimit);
    }
    else
      core.setSplitFrequency(bandIndex, newFreq);
  }
  else
    core.setSplitFrequency(bandIndex, newFreq); // topmost band has no limit

  // somehow, we must also make sure sure that the topmost band has a freq of 20000 (it's the 
  // highpass band and actually has not lowpass cutoff at all)

  jassert(areBandsInIncreasingOrder(false)); // for debug
}

void MultiBandEffect::setSplitFreq(double newFreq) 
{ 
  ScopedLock scopedLock(*lock);
  setSplitFreq(selectedBand, newFreq);
}

Parameter* MultiBandEffect::getSplitFreqParam(int bandIndex)
{
  ScopedLock scopedLock(*lock);
  return splitFreqParams[bandIndex];
}

int MultiBandEffect::getBandContainingFrequency(double freq)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < core.getNumberOfBands()-1; i++)
    if(core.getSplitFrequency(i) > freq)
      return i;
  return core.getNumberOfBands()-1;
}

std::vector<juce::String> MultiBandEffect::getAvailableEffectTypes() const
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
  ScopedLock scopedLock(*lock);
  for(int i = 1; i < core.getNumberOfBands()-1; i++)
    if(!isLess(core.getSplitFrequency(i-1), core.getSplitFrequency(i), strictly))
      return false;
  return true;
}

// split-freq param handling:

void MultiBandEffect::createSplitFreqParams()
{
  ScopedLock scopedLock(*lock);
  // clearSplitFreqParams(); // ...maybe later
  size_t numParams = splitFreqParams.size();
  size_t numSplits = getNumBands()-1;
  while(numParams < numSplits) {
    addSplitFreqParam((int)numParams, getSplitFreq((int)numParams));
    numParams++; }

  jassert(areBandsInIncreasingOrder(false)); // for debug
}

void MultiBandEffect::clearSplitFreqParams()
{
  // 
}

void MultiBandEffect::addSplitFreqParam(int i, double freq)
{
  ScopedLock scopedLock(*lock);
  Parameter* p = new Parameter("", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL); // name assigned later
  p->setValue(freq, false, false);
  addObservedParameter(p);
  insert(splitFreqParams, p, i);
  updateSplitFreqParamNamesAndCallbacks();
}

void MultiBandEffect::removeSplitFreqParam(int i)
{
  ScopedLock scopedLock(*lock);
  delete splitFreqParams[i];
  remove(splitFreqParams, i);
  updateSplitFreqParamNamesAndCallbacks();
}

void MultiBandEffect::updateSplitFreqParamNamesAndCallbacks()
{
  for(size_t i = 0; i < splitFreqParams.size(); i++){
    Parameter* p = splitFreqParams[i];
    p->setName("SplitFrequency" + String(i+1));
    p->setValueChangeCallback([=](double v){ setSplitFreq((int)i, v); }); 
  }
  // gui needs to be notified to update slider-names...or maybe it can update them inside a 
  // callback that is already being called...
}

void MultiBandEffect::sendBandInsertNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWasInserted(this, index);
}

void MultiBandEffect::sendBandRemovePreNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWillBeRemoved(this, index);
}

void MultiBandEffect::sendBandRemovePostNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWasRemoved(this, index);
}

void MultiBandEffect::sendBandSelectNotification(int index)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->bandWasSelected(this, index);
}

void MultiBandEffect::insertBandEffect(int i)
{
  ScopedLock scopedLock(*lock);
  AudioModule* m = perBandModuleFactory.createModule(effectTypeString);
  insert(perBandModules, m, i);
  updateBandModuleNames();
}

void MultiBandEffect::removeBandEffect(int i)
{
  ScopedLock scopedLock(*lock);
  delete perBandModules[i];
  remove(perBandModules, i);
  updateBandModuleNames();
}

void MultiBandEffect::clearBandEffects()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < perBandModules.size(); i++)
  {
    sendBandRemovePreNotification(i);
    delete perBandModules[i];
    sendBandRemovePostNotification(i);
  }
  perBandModules.clear();
}

void MultiBandEffect::updateBandModuleNames()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < perBandModules.size(); i++)
    perBandModules[i]->setModuleName(effectTypeString + String(i+1));

  // we need to update the splitFreq parameter names, too
}

//=================================================================================================

MultiBandPlotEditor::MultiBandPlotEditor(jura::MultiBandEffect* moduleToEdit)
  : module(moduleToEdit)

{
  module->addChangeListener(this); // obsolete? ..currently used for repainting when freq changes - later connect to the parameters directly
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
  module->deRegisterMultiBandObserver(this);
  delete bandPopup;
}

void MultiBandPlotEditor::bandWasInserted(MultiBandEffect* mbe, int i)
{
  // insert new graph into plot
  //freqRespPlot->addFunction([=](double f)->double { return core->getDecibelsAt(i, f); });
  //freqRespPlot->insertFunction(i, [=](double f)->double { return core->getDecibelsAt(i, f); });
  repaint();
}

void MultiBandPlotEditor::bandWillBeRemoved(MultiBandEffect* mbe, int i)
{
  //freqRespPlot->removeFunction(i);  // remove graph from plot
}

void MultiBandPlotEditor::bandWasRemoved(MultiBandEffect* mbe, int i)
{
  repaint();
}

void MultiBandPlotEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{
  // update selection highlighting
  freqRespPlot->setNumFunctionsToPlot(core->getNumberOfBands());
  repaint();
}

void MultiBandPlotEditor::changeListenerCallback(ChangeBroadcaster* source) // obsolete?
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
  //setSize(595, 301);
  setSize(595, 420);
}

MultiBandEffectEditor::~MultiBandEffectEditor()
{
  effectToEdit->deRegisterMultiBandObserver(this);
  //clearBandEditors();
}

void MultiBandEffectEditor::resized()
{
  AudioModuleEditor::resized();

  int y = getPresetSectionBottom() + 4;
  //plotEditor->setBounds(0, y, getWidth(), getHeight()-110); // 595 x 191 -> ideal
  plotEditor->setBounds(0, y, getWidth(), 191);
  y = plotEditor->getBottom() + 4;

  int x = 4;
  int w = getWidth() / 3 - 8;
  int h = 16;
  int d = h-2;

  effectSelectBox->setBounds(x, y, w, h); y += d;
  splitModeBox   ->setBounds(x, y, w, h); y += d;
  for(size_t i = 0; i < splitFreqSliders.size(); i++) {
    splitFreqSliders[i]->setBounds(x, y, w, h); y += d; }

  positionBandEditors();
  updateSplitSliderPositions();
}

void MultiBandEffectEditor::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
{
  // respond to selection of new per-band effect type
}

void MultiBandEffectEditor::bandWasInserted(MultiBandEffect* mbe, int index)
{
  insertBandEditor(index);
  updateEditorNames();
  updateSplitSliders();
}

void MultiBandEffectEditor::bandWillBeRemoved(MultiBandEffect* mbe, int index)
{
  removeBandEditor(index);
}

void MultiBandEffectEditor::bandWasRemoved(MultiBandEffect* mbe, int i)
{
  updateEditorNames();
  updateSplitSliders();
  updateEditorVisibility();
}

void MultiBandEffectEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{
  updateEditorVisibility();
}

void MultiBandEffectEditor::insertBandEditor(int i)
{
  AudioModuleEditor* e = effectToEdit->getBandEffect(i)->createEditor();
  addChildEditor(e);
  insert(perBandEditors, e, i);
  positionBandEditor(i);
}

void MultiBandEffectEditor::removeBandEditor(int i)
{
  AudioModuleEditor* e = perBandEditors[i];
  remove(perBandEditors, i);
  removeChildEditor(e, true);
}

void MultiBandEffectEditor::updateEditorVisibility()
{
  for(size_t i = 0; i < perBandEditors.size(); i++)
    perBandEditors[i]->setVisible(false);
  int selected = effectToEdit->getSelectedBand();
  jassert(selected < (int)perBandEditors.size());   // no editor for selected band available
  if(selected >= 0)
    perBandEditors[selected]->setVisible(true);
}

void MultiBandEffectEditor::updateEditorNames()
{
  for(int i = 0; i < (int)perBandEditors.size(); i++)
    perBandEditors[i]->setHeadlineText(effectToEdit->getBandEffect(i)->getModuleName());
}

void MultiBandEffectEditor::createWidgets()
{
  RComboBox *c;

  plotEditor = new MultiBandPlotEditor(effectToEdit);
  addChildColourSchemeComponent(plotEditor);

  addWidget( c = effectSelectBox = new RComboBox() );
  c->setDescription("Effect type to be applied to each band");
  c->setDescriptionField(infoField);
  // populate box with available effect types
  //c->addListener(this);

  addWidget( c = splitModeBox = new RComboBox() );
  c->assignParameter( effectToEdit->getParameterByName("SplitMode") );
  c->setDescription("Mode of the band-splitting");
  c->setDescriptionField(infoField);

  updateSplitSliders();
}

void MultiBandEffectEditor::createBandEditors()
{
  for(size_t i = 0; i < effectToEdit->getNumBands(); i++)
    insertBandEditor((int)i);
}

void MultiBandEffectEditor::clearBandEditors()
{
  while(perBandEditors.size() > 0)
    removeBandEditor((int)perBandEditors.size()-1);
}

void MultiBandEffectEditor::positionBandEditors()
{
  for(size_t i = 0; i < perBandEditors.size(); i++)
    positionBandEditor((int)i);
}

void MultiBandEffectEditor::positionBandEditor(int i)
{
  int x = effectSelectBox->getRight() + 4;
  int y = plotEditor->getBottom();
  int w = getWidth()  - x;
  int h = getHeight() - y;
  perBandEditors[i]->setBounds(x, y, w, h);
}

void MultiBandEffectEditor::updateSplitSliders()
{
  RSlider *s;
  int numBands  = effectToEdit->getNumBands();
  int numSplits = numBands-1;

  // make sure that number of split-sliders matches number of split-parameters:
  while(splitFreqSliders.size() > numSplits) {
    s = getAndRemoveLast(splitFreqSliders);
    removeWidget(s, true, true); 
  }
  while(splitFreqSliders.size() < numSplits) {
    s = new RSlider;
    addWidget(s);
    append(splitFreqSliders, s); 
  }

  // (re)assign parameters to sliders and update slider positions:
  for(int i = 0; i < numSplits; i++)
    splitFreqSliders[i]->assignParameter(
      effectToEdit->getParameterByName("SplitFrequency" + String(i+1)));
  updateSplitSliderPositions();
}

void MultiBandEffectEditor::updateSplitSliderPositions()
{
  int x  = splitModeBox->getX();
  int y  = splitModeBox->getBottom() + 4;
  int w  = splitModeBox->getWidth();
  int h  = 16;
  int dy = h-2;
  for(size_t i = 0; i < splitFreqSliders.size(); i++) {
    splitFreqSliders[i]->setBounds(x, y, w, h);
    y += dy;
  }
}









//=================================================================================================

MultiCompAudioModule::MultiCompAudioModule(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse)
  : MultiBandEffect(lockToUse,  metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("MultiComp");
  //MultiBandEffect::setEffectCore(&multiCompCore);
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
  : MultiBandEffectEditor(multiCompModuleToEdit)
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