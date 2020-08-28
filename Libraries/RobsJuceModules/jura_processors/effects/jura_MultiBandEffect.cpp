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
  //f.registerModuleType([](CS cs)->AM { return new GainAudioModule(cs); }, s, "StereoPan"); // good for testing the system
  // Noisifier -> nice effect, Tremolo, Vibrato, WaveShaper

  f.registerModuleType([](CS cs)->AM { return new CompressorAudioModule(cs); }, s, "Compressor");
  //f.registerModuleType([](CS cs)->AM { return new FuncShaperAudioModule(cs); }, s, "FuncShaper");
  // ...
  setEffectType("Compressor");


  // create parameters (factor out):
  Parameter* p = new Parameter("SplitMode", 0.0, 2.0, 0.0, Parameter::STRING);
  p->addStringValue("Steep Highpass");
  p->addStringValue("Steep Lowpass");
  //p->addStringValue("Binary Tree"); // doesn't work yet
  addObservedParameter(p);
  p->setValueChangeCallback<rosic::rsMultiBandEffect>(
    &core, &rosic::rsMultiBandEffect::setSplitMode);


  // maybe remove this and start with 0 bands? ...but maybe not
  core.initBands(1);
  insertBandEffect(0);
  createSplitFreqParams();
  // maybe factor out into init-function (it's called in the same sequence in setState)
}

MultiBandEffect::~MultiBandEffect()
{
  clearBandEffects(true);
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
  double xL = 0, xR = 0, yL = 0, yR = 0;    // for per-band I/O gain measurement
  core.split(left, right);
  for(int k = 0; k < getNumBands(); k++) {  // process individual bands
    // input gain measurement:
    xL = *core.getLeft(k);
    xR = *core.getRight(k);
    inSumOfSquares[k]  += xL*xL + xR*xR;

    // processing:
    perBandModules[k]->processStereoFrame(core.getLeft(k), core.getRight(k));

    // output gain measurement:
    yL = *core.getLeft(k);
    yR = *core.getRight(k);
    outSumOfSquares[k] += yL*yL + yR*yR;
  }
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
    perBandModules[i]->reset();
  // todo: reset gain accumulators
}

AudioModuleEditor* MultiBandEffect::createEditor(int type)
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
    // -or: allow the type of effect only to be changed when the number of bands is 0
    // -maybe have a clear button that completely deletes all bands
  }
}

void MultiBandEffect::setStateFromXml(const XmlElement& xml, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  selectBand(-1, true); 

  sendClearBandsNotification(); // gui will delete all editors and split-freq sliders

  size_t numBands = xml.getIntAttribute("NumBands", 1);
  //effectTypeString = xml.getStringAttribute("EffectType", "Gain");
  effectTypeString = xml.getStringAttribute("EffectType", "Compressor");

  clearSplitFreqParams();
  clearBandEffects(false); // don't notify gui, because it already has deleted the editors
  core.initBands(1);
  insertBandEffect(0);
  createSplitFreqParams();
  for(size_t i = 1; i < numBands; i++) {
    double freq = xml.getDoubleAttribute("SplitFrequency" + String(i), 1000);
    insertBand((int)(i-1), freq, false);
  }

  // restore the states of the per-band modules:
  for(size_t i = 0; i < numBands; i++) {
   auto childXml = xml.getChildByName(effectTypeString + String(i+1));
    if(childXml != nullptr)
      perBandModules[i]->setStateFromXml(*childXml, "", markAsClean);
  }

  sendTotalRefreshNotification(); // gui will create new editors and split-freq sliders
}

XmlElement* MultiBandEffect::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement* xml = ModulatableAudioModule::getStateAsXml(stateName, markAsClean);
  xml->setAttribute("NumBands",   getNumBands());
  xml->setAttribute("EffectType", effectTypeString);
  for(size_t i = 0; i < perBandModules.size(); i++) {
    String name = perBandModules[i]->getModuleName();
    XmlElement* child = perBandModules[i]->getStateAsXml(name, markAsClean);
    xml->addChildElement(child);
  }
  return xml;
}

void MultiBandEffect::parameterChanged(Parameter* p)
{
  ModulatableAudioModule::parameterChanged(p);
  sendChangeMessage();
}

void MultiBandEffect::insertBand(int index, double splitFrequency, bool sendNotification)
{
  ScopedLock scopedLock(*lock);

  selectBand(-1, true); // preliminary - avoid weirdness

  core.insertBand(index, splitFrequency);
  addSplitFreqParam(index, splitFrequency); // rename to insertSplitFreqParam
  insertBandEffect(index);
  if(sendNotification)
    sendBandInsertNotification(index);  // notify gui (for creating widgets)
}

void MultiBandEffect::removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification)
{
  ScopedLock scopedLock(*lock);

  //if(getNumBands() == 1)  return;  // no...only return at 0
  if(getNumBands() == 0)  return;
    

  // preliminary - disallow deletion of last bad (would lead to crash):
  if(index >= getNumBands()-1)
  {
    jassertfalse;
    return;
  }


  selectBand(-1, true); // preliminary...
  // ...perhaps late do something more sophisticated, like:
  //if(selectedBand == index || selectedBand == getNumBands()-1)
  //  selectBand(-1, true);
  //if(index < selectedBand)
  //  selectBand(selectedBand-1, true);


  if(sendNotification)
    sendBandRemovePreNotification(index);  // notify gui (for removing widgets)

  removeBandEffect(index);
  core.removeBand(index);      // crashes whe last band is removed
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
  int numSplits = getNumSplits();
  if(bandIndex < numSplits) {

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

double MultiBandEffect::getBandInOutGain(int i, bool resetAccus)
{
  ScopedLock scopedLock(*lock);
  if(inSumOfSquares[i] < 1.e-10)  // avoid div-by-zero
    return 1.0;
  double g = sqrt(outSumOfSquares[i]/inSumOfSquares[i]);
  if(resetAccus) {
    inSumOfSquares[i]  = 0;
    outSumOfSquares[i] = 0;
  }
  return g; // maybe we need a sample-counter and divide by that?
}

int MultiBandEffect::getNumBands() const 
{ 
  ScopedLock scopedLock(*lock);
  jassert(core.getNumberOfBands() == perBandModules.size()); // if they don't match, something is wrong
  return core.getNumberOfBands(); 
  //return perBandModules.size(); 
}

int MultiBandEffect::getNumSplits() const
{
  ScopedLock scopedLock(*lock);
  return jmax(0, getNumBands()-1);
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

bool MultiBandEffect::isBandRemovable(int index)
{
  return index < getNumBands()-1;
}

// split-freq param handling:

void MultiBandEffect::createSplitFreqParams()
{
  ScopedLock scopedLock(*lock);
  // clearSplitFreqParams(); // ...maybe later
  size_t numParams = splitFreqParams.size();
  size_t numSplits = getNumSplits();

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

void MultiBandEffect::sendClearBandsNotification()
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->allBandsWillBeRemoved(this);
}

void MultiBandEffect::sendTotalRefreshNotification()
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->totalRefreshNeeded(this);
}

void MultiBandEffect::insertBandEffect(int i)
{
  ScopedLock scopedLock(*lock);
  AudioModule* m = perBandModuleFactory.createModule(effectTypeString);

  // todo (maybe): copy effect settings from the current i-th band-effect (or i-1 or i+1 th?)
  // such that the new band starts with the same settings as the old band from which it was split
  // off

  // we don't make it a child module, because in the child-module array, the band-modules may have
  // a messed up order due to the user possibly inserting and removing bands randomly - but we need
  // an ordered array for state storage and recall

  insert(perBandModules, m, i);
  insert(inSumOfSquares,  0.0, i);
  insert(outSumOfSquares, 0.0, i);
  updateBandModuleNames();
}

void MultiBandEffect::removeBandEffect(int i)
{
  ScopedLock scopedLock(*lock);
  delete perBandModules[i];
  remove(perBandModules, i);
  remove(inSumOfSquares,  i);
  remove(outSumOfSquares, i);
  updateBandModuleNames();
}

void MultiBandEffect::clearBandEffects(bool notify)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < perBandModules.size(); i++)
  {
    if(notify) sendBandRemovePreNotification(i);  // maybe generally don't send notifications here
    delete perBandModules[i];
    if(notify) sendBandRemovePostNotification(i);
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
  refreshFunctionsToPlot();

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
  refreshFunctionsToPlot();
  repaint();
}

void MultiBandPlotEditor::bandWillBeRemoved(MultiBandEffect* mbe, int i)
{
  //freqRespPlot->removeFunction(i);  // remove graph from plot
}

void MultiBandPlotEditor::bandWasRemoved(MultiBandEffect* mbe, int i)
{
  refreshFunctionsToPlot();
  repaint();
}

void MultiBandPlotEditor::bandWasSelected(MultiBandEffect* mbe, int index)
{
  // update selection highlighting
  freqRespPlot->setNumFunctionsToPlot(core->getNumberOfBands());
  //refreshFunctionsToPlot(); // for test - may be removed later
  repaint();
}

void MultiBandPlotEditor::allBandsWillBeRemoved(MultiBandEffect* mbe)
{
  repaint(); //?
}

void MultiBandPlotEditor::totalRefreshNeeded(MultiBandEffect* mbe)
{
  refreshFunctionsToPlot();
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
  int index = module->getBandContainingFrequency(freqAtMouse);

  if(e.mods.isLeftButtonDown()) {
    // select band whose rectangle contains the mouse-event:

    if(index == module->getSelectedBand())
      module->selectBand(-1, true); // de-select, when a band is clicked again
    else
      module->selectBand(index, true);

    // old (we don't have a Parameter for this anymore):
    //Parameter* p = module->getParameterByName("SelectedBand");
    //p->setValue(index, true, true);
  }
  else if(e.mods.isRightButtonDown())
    openRightClickMenu(index);
}

void MultiBandPlotEditor::resized()
{
  freqRespPlot->setBounds(0, 0, getWidth(), getHeight());
  //refreshFunctionsToPlot(); // for test - remove later
}

void MultiBandPlotEditor::paintOverChildren(Graphics& g)
{
  paintBandShadings(g);
  paintSplitLines(g);
}

void MultiBandPlotEditor::paintBandShadings(Graphics& g)
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
    if(x2 >= x1)
      g.fillRect(x1, 0.f, x2-x1, (float)getHeight());
  }

  // todo: maybe give all bands a color (going through the rainbow) and just use a higher alpha
  // for the selected band (that may look nice)
}

void MultiBandPlotEditor::paintSplitLines(Graphics& g)
{
  int numBands = core->getNumberOfBands();

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
}

void MultiBandPlotEditor::openRightClickMenu(int bandIndex)
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

  if(module->isBandRemovable(bandIndex))
    bandPopup->setItemEnabled(REMOVE_BAND, true);
  else
    bandPopup->setItemEnabled(REMOVE_BAND, false); // doesn't work - why? bug?

  // show it:
  int w = bandPopup->getRequiredWidth(true);
  int h = bandPopup->getRequiredHeight(true);
  bandPopup->selectItemByIndex(-1, false);    // select nothing - maybe write a specific fucntion for that
  bandPopup->showAtMousePosition(true, w, h); // showModally = true
}

void MultiBandPlotEditor::refreshFunctionsToPlot()
{
  freqRespPlot->init();
  for(int i = 0; i < core->getNumberOfBands(); i++)
    freqRespPlot->addFunction([=](double f)->double { return core->getDecibelsAt(i, f); });
  repaint();
}

//=================================================================================================

MultiBandPlotEditorAnimated::MultiBandPlotEditorAnimated(jura::MultiBandEffect* moduleToEdit)
  : MultiBandPlotEditor(moduleToEdit)
{
  startTimerHz(30); // 30 fps - maybe make adjustable from client code
}

MultiBandPlotEditorAnimated::~MultiBandPlotEditorAnimated()
{

}

void MultiBandPlotEditorAnimated::timerCallback()
{
  repaint();
}

void MultiBandPlotEditorAnimated::paintOverChildren(Graphics& g)
{
  MultiBandPlotEditor::paintOverChildren(g);
  paintInOutGains(g);
}

void MultiBandPlotEditorAnimated::paintInOutGains(Graphics& g)
{
  double freq, gain;           // gain in dB
  float x1, y1, x2 = 0.f, y2;  // rectangle pixel coordinates
  int numBands = core->getNumberOfBands();
  g.setColour(Colours::green.withMultipliedAlpha(0.5f));
  for(int i = 0; i < numBands; i++)
  {
    // get frequency (at rectangle's right) and gain:
    freq = core->getSplitFrequency(i);
    gain = RAPT::rsAmpToDb(module->getBandInOutGain(i));

    // compute rectangle pixel coordinates and draw:
    x1 = x2; // this rect's left is previous rect's right
    x2 = (float)freqRespPlot->toPixelX(freq);
    y1 = (float)freqRespPlot->toPixelY(0.0);
    y2 = (float)freqRespPlot->toPixelY(gain);
    if(y2 < y1)
      swap(y1, y2);
    if(i == numBands-1)
      x2 = (float)getWidth();
    if(x2 >= x1)
      g.fillRect(x1, y1, x2-x1, y2-y1);
  }
}

//=================================================================================================

MultiBandPlotEditorAnimated2::MultiBandPlotEditorAnimated2(jura::MultiBandEffect* moduleToEdit)
  : MultiBandPlotEditorAnimated(moduleToEdit)
{
  freqRespPlot->setVisible(false);
}

void MultiBandPlotEditorAnimated2::paint(Graphics& g)
{
  //MultiBandPlotEditorAnimated::paint(g); return; // preliminary

  //ScopedLock scopedLock(*(module->getCriticalSection())); // do we need this?
  if(backgroundIsDirty)
    updateBackgroundImage();
  g.drawImageAt(background, 0, 0);
  paintInOutGains(g);
}

void MultiBandPlotEditorAnimated2::paintOverChildren(Graphics& g)
{
  // overriden with empty implementation to avoid baseclass method drawing the gain bars a second
  // time (leading to flickering gain bars)
}

void MultiBandPlotEditorAnimated2::resized()
{
  MultiBandPlotEditorAnimated::resized();
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::changeListenerCallback(ChangeBroadcaster* source)
{
  // this function gets called when a split frequency changes
  MultiBandPlotEditorAnimated::changeListenerCallback(source);
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::bandWasInserted(MultiBandEffect* mbe, int index)
{
  MultiBandPlotEditorAnimated::bandWasInserted(mbe, index);
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::bandWasRemoved(MultiBandEffect* mbe, int index)
{
  MultiBandPlotEditorAnimated::bandWasRemoved(mbe, index);
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::bandWasSelected(MultiBandEffect* mbe, int index)
{
  MultiBandPlotEditorAnimated::bandWasSelected(mbe, index);
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::totalRefreshNeeded(MultiBandEffect* mbe)
{
  MultiBandPlotEditorAnimated::totalRefreshNeeded(mbe);
  backgroundIsDirty = true;
}

void MultiBandPlotEditorAnimated2::updateBackgroundImage()
{
  background = Image(Image::PixelFormat::RGB, getWidth(), getHeight(), false); // or maybe ARGB?
  Graphics g(background);
  freqRespPlot->paint(g);
  paintBandShadings(g);
  paintSplitLines(g);
  backgroundIsDirty = false;
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

void MultiBandEffectEditor::allBandsWillBeRemoved(MultiBandEffect* mbe)
{
  clearBandEditors();
  updateEditorVisibility();
}

void MultiBandEffectEditor::totalRefreshNeeded(MultiBandEffect* mbe)
{
  createBandEditors();
  updateSplitSliders();
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

  //plotEditor = new MultiBandPlotEditor(effectToEdit);
  //plotEditor = new MultiBandPlotEditorAnimated(effectToEdit);
  plotEditor = new MultiBandPlotEditorAnimated2(effectToEdit);
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
  //int numBands  = effectToEdit->getNumBands();
  int numSplits = effectToEdit->getNumSplits();

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
