AudioModule* AudioModuleFactory::createModule(const juce::String& type, CriticalSection *lock)
{
  if(type == "None")         return new DummyModule( lock);

  // analysis:
  if(type == "PhaseScope")    return new PhaseScope(              lock);
  //if(type == "PhaseScope2")  return new PhaseScope2( lock);
  if(type == "MultiAnalyzer") return new MultiAnalyzerAudioModule(lock);
  if(type == "TrackMeter")    return new TrackMeterAudioModule(   lock);
  if(type == "MidiMonitor")   return new MidiMonitorAudioModule(  lock);

  // generators:
  if(type == "RayBouncer")    return new RayBouncerAudioModule(lock);

  // filters:
  if(type == "Equalizer")       return new EqualizerAudioModule(      lock);
  if(type == "Ladder")          return new Ladder(                    lock);
  if(type == "PhasorFilter")    return new PhasorFilter(              lock);
  if(type == "EngineersFilter") return new EngineersFilterAudioModule(lock);
  if(type == "CrossOver")       return new CrossOverAudioModule(      lock);

  // effects:
  if(type == "Enveloper")        return new Enveloper(                  lock);
  if(type == "FuncShaper")       return new FuncShaperAudioModule(      lock);
  if(type == "AlgoVerb")         return new AlgoVerbAudioModule(        lock);
  if(type == "EchoLab")          return new EchoLabAudioModule(         lock);
  if(type == "StereoDelay")      return new StereoDelayAudioModule(     lock);
  if(type == "PitchShifter")     return new PitchShifterAudioModule(    lock);
  if(type == "Quadrifex")        return new QuadrifexAudioModule(       lock);
  if(type == "Moduluxury")       return new ModuluxuryAudioModule(      lock);
  if(type == "ChannelMatrix2x2") return new ChannelMatrix2x2AudioModule(lock);
  if(type == "DspWorkbench")     return new DspWorkbenchAudioModule(    lock);

  // instruments:
  if(type == "AciDevil")      return new AciDevilAudioModule(     lock);
  if(type == "Straightliner") return new StraightlinerAudioModule(lock);
  if(type == "MagicCarpet")   return new MagicCarpetAudioModule(  lock);
  if(type == "SimpleSampler") return new SimpleSamplerAudioModule(lock);
  if(type == "KeyShot")       return new KeyShotAudioModule(      lock);
  if(type == "Quadriga")      return new QuadrigaAudioModule(     lock);
  if(type == "Workhorse")     return new WorkhorseAudioModule(    lock);
#ifdef _MSC_VER
  if(type == "Liberty")       return new LibertyAudioModule(      lock);
#endif

  jassertfalse;  // unknown module type requested
  return nullptr;
}

juce::String AudioModuleFactory::getModuleType(AudioModule *m)
{
  if(dynamic_cast<DummyModule*>  (m))              return "None";

  // analysis:
  //if(dynamic_cast<PhaseScope2*>  (m))            return "PhaseScope2"; // always check subclasses before...
  if(dynamic_cast<PhaseScope*>   (m))              return "PhaseScope";  // ...their superclasses
  if(dynamic_cast<MultiAnalyzerAudioModule*> (m))  return "MultiAnalyzer";
  if(dynamic_cast<TrackMeterAudioModule*> (m))     return "TrackMeter";
  if(dynamic_cast<MidiMonitorAudioModule*> (m))    return "MidiMonitor";

  // generators:
  if(dynamic_cast<RayBouncerAudioModule*>(m))      return "RayBouncer";

  // filters:
  if(dynamic_cast<EqualizerAudioModule*>(m))       return "Equalizer";
  if(dynamic_cast<Ladder*>(m))                     return "Ladder";
  if(dynamic_cast<PhasorFilter*>(m))               return "PhasorFilter";
  if(dynamic_cast<EngineersFilterAudioModule*>(m)) return "EngineersFilter";
  if(dynamic_cast<CrossOverAudioModule*>(m))       return "CrossOver";

  // effects:
  if(dynamic_cast<Enveloper*>(m))                   return "Enveloper";
  if(dynamic_cast<FuncShaperAudioModule*>(m))       return "FuncShaper";
  if(dynamic_cast<AlgoVerbAudioModule*>(m))         return "AlgoVerb";
  if(dynamic_cast<EchoLabAudioModule*>(m))          return "EchoLab";
  if(dynamic_cast<StereoDelayAudioModule*>(m))      return "StereoDelay";
  if(dynamic_cast<PitchShifterAudioModule*>(m))     return "PitchShifter";
  if(dynamic_cast<QuadrifexAudioModule*>(m))        return "Quadrifex";
  if(dynamic_cast<ModuluxuryAudioModule*>(m))       return "Moduluxury";
  if(dynamic_cast<ChannelMatrix2x2AudioModule*>(m)) return "ChannelMatrix2x2";
  if(dynamic_cast<DspWorkbenchAudioModule*>(m))     return "DspWorkbench";

  // instruments:
  if(dynamic_cast<AciDevilAudioModule*> (m))       return "AciDevil";
  if(dynamic_cast<StraightlinerAudioModule*> (m))  return "Straightliner";
  if(dynamic_cast<MagicCarpetAudioModule*> (m))    return "MagicCarpet";
  if(dynamic_cast<SimpleSamplerAudioModule*> (m))  return "SimpleSampler";
  if(dynamic_cast<KeyShotAudioModule*> (m))        return "KeyShot";
  if(dynamic_cast<QuadrigaAudioModule*> (m))       return "Quadriga";
  if(dynamic_cast<WorkhorseAudioModule*> (m))      return "Workhorse";

#ifdef _MSC_VER
  if(dynamic_cast<LibertyAudioModule*> (m))        return "Liberty";
#endif

  jassertfalse;  // unknown module type was passed
  return "UnknownType";
}

StringArray AudioModuleFactory::getAvailableModuleTypes()
{
  // here, we can make certain modules temporarily unavailable by simply commenting the
  // corresponding line

  StringArray a;
  a.add("None");        // maybe use "Empty" instead of "None"

  // analysis:
  a.add("PhaseScope");
  //a.add("PhaseScope2");
  a.add("MultiAnalyzer");
  a.add("TrackMeter");
  a.add("MidiMonitor");

  // generators:
  a.add("RayBouncer");

  // filters:
  a.add("Equalizer");
  a.add("Ladder");
  //a.add("PhasorFilter");
  a.add("EngineersFilter");
  //a.add("CrossOver"); // makes actually no sense in 2in/2out plugin - but just for test

  // effects:
  a.add("Enveloper");
  a.add("FuncShaper");
  //a.add("AlgoVerb");  // currently inactive - not yet complete
  //a.add("EchoLab");
  a.add("StereoDelay");
  a.add("PitchShifter");
  //a.add("Quadrifex");
  //a.add("Moduluxury");
  //a.add("ChannelMatrix2x2");
  //a.add("DspWorkbench");

  // instruments:
  a.add("AciDevil");
  a.add("Straightliner");
  //a.add("MagicCarpet");
  //a.add("SimpleSampler");
  //a.add("KeyShot");
  //a.add("Quadriga");
  //a.add("Workhorse");
#ifdef _MSC_VER
  a.add("Liberty"); // not yet available on gcc
#endif

  return a;
}

//=================================================================================================

AudioModuleSelector::AudioModuleSelector() : RComboBox("ModuleSelector")
{
  // old: linear flat array:
  setDescription("Select module type");
  StringArray a = AudioModuleFactory::getAvailableModuleTypes();
  for(int i = 0; i < a.size(); i++)
    addItem(i, a[i]);
  // ...but we want a tree...


//  // ...but that does not yet work:
//  // populate the tree:
//
//  RTreeViewNode *node;
//  int i = 1;           //  the index is actually not used, but we need it as dummy
//
//  node = new RTreeViewNode("Analysis", -1, "Analysis");
//  node->addChildNode(new RTreeViewNode("PhaseScope",      i++));
//  popUpMenu->addTreeNodeItem(node);
//
//  node = new RTreeViewNode("Filters", -1, "Filters");
//  node->addChildNode(new RTreeViewNode("Ladder",          i++));
//  node->addChildNode(new RTreeViewNode("PhasorFilter",    i++));
//  node->addChildNode(new RTreeViewNode("EngineersFilter", i++));
//  popUpMenu->addTreeNodeItem(node);
//
//  node = new RTreeViewNode("Effects", -1, "Effects");
//  node->addChildNode(new RTreeViewNode("Enveloper",  i++));
//  node->addChildNode(new RTreeViewNode("FuncShaper", i++));
//  popUpMenu->addTreeNodeItem(node);
//
//  node = new RTreeViewNode("Instruments", -1, "Instruments");
//  node->addChildNode(new RTreeViewNode("AciDevil", i++));
//#ifdef _MSC_VER
//  node->addChildNode(new RTreeViewNode("Liberty",  i++));
//#endif
//  popUpMenu->addTreeNodeItem(node);
}

//=================================================================================================

AudioModuleChain::AudioModuleChain(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "Chainer";

#ifdef _WIN32
  juce::String presetPath = getApplicationDirectory() + "/ChainerPresets";
#elif __APPLE__
  juce::String presetPath = "/Library/Audio/Presets/RS-MET/Chainer/ChainerPresets";
#elif __linux__
  juce::String presetPath = getApplicationDirectory() + "/ChainerPresets";
#endif
  setActiveDirectory(presetPath);
  //setActiveDirectory(getApplicationDirectory() + "/ChainerPresets");  // old

  addEmptySlot();
}

AudioModuleChain::~AudioModuleChain()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++)
    delete modules[i];
}

void AudioModuleChain::addEmptySlot()
{
  addModule("None");
}

void AudioModuleChain::addModule(const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  AudioModule *m = AudioModuleFactory::createModule(type, lock);
  m->setMetaParameterManager(metaParamManager); // without, we hit jassert(metaParaManager != nullptr) in MetaControlledParameter::attachToMetaParameter
  append(modules, m);
  sendAudioModuleWasAddedNotification(m, size(modules)-1);
}

void AudioModuleChain::deleteModule(int index)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  if(activeSlot == index)
    activeSlot--;
  sendAudioModuleWillBeDeletedNotification(modules[index], index);
  delete modules[index];
  remove(modules, index);
}

void AudioModuleChain::deleteLastModule()
{
  ScopedLock scopedLock(*lock);
  deleteModule(size(modules) - 1);
}

void AudioModuleChain::replaceModule(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  if(!isModuleOfType(index, type)){              // replace only, if new type is different
    AudioModule* oldModule = modules[index];
    AudioModule* newModule = AudioModuleFactory::createModule(type, lock);
    newModule->setMetaParameterManager(metaParamManager);
    newModule->loadDefaultPreset(); // later: either load default preset or recall a stored state
    newModule->setSampleRate(sampleRate);
    modules[index] = newModule;
    sendAudioModuleWasReplacedNotification(oldModule, newModule, index);
    delete oldModule;
    activeSlot = index;
    ensureOneEmptySlotAtEnd();
  }
}

bool AudioModuleChain::isModuleOfType(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(modules)); // index out of range
  return type == AudioModuleFactory::getModuleType(modules[index]);
}

AudioModule* AudioModuleChain::getModuleAt(int index)
{
  if(index < 0 || index >= size(modules))  // no assert, this is supposed to happen
    return nullptr;
  return modules[index];
}

void AudioModuleChain::ensureOneEmptySlotAtEnd()
{
  ScopedLock scopedLock(*lock);

  // if the last module/slot is not the "empty" dummy module, add another slot:
  if(!isModuleOfType(size(modules)-1, "None"))
    addEmptySlot();

  // remove superfluous empty slots at end:
  while(modules.size() > 1 && isModuleOfType(size(modules)-1, "None")
                           && isModuleOfType(size(modules)-2, "None"))
  { // if the last two slots are empty, remove the last
    deleteLastModule();
  }
}

void AudioModuleChain::addAudioModuleChainObserver(AudioModuleChainObserver *observerToAdd)
{
  ScopedLock scopedLock(*lock);
  appendIfNotAlreadyThere(observers, observerToAdd);
}

void AudioModuleChain::removeAudioModuleChainObserver(AudioModuleChainObserver *observerToRemove)
{
  ScopedLock scopedLock(*lock);
  removeFirstOccurrence(observers, observerToRemove);
}

void AudioModuleChain::sendAudioModuleWasAddedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWasAdded(this, module, index);
}

void AudioModuleChain::sendAudioModuleWillBeDeletedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWillBeDeleted(this, module, index);
}

void AudioModuleChain::sendAudioModuleWasReplacedNotification(AudioModule *oldModule,
  AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWasReplaced(this, oldModule, newModule, index);
}

// overrides:

AudioModuleEditor* AudioModuleChain::createEditor()
{
  return new AudioModuleChainEditor(this);
}

void AudioModuleChain::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  ScopedLock scopedLock(*lock);
  jassert(numChannels == 2);
  for(int i = 0; i < size(modules); i++)
    modules[i]->processBlock(inOutBuffer, numChannels, numSamples);
}

void AudioModuleChain::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  sampleRate = newSampleRate;
  for(int i = 0; i < size(modules); i++)
    modules[i]->setSampleRate(sampleRate);
}

void AudioModuleChain::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOn(noteNumber, velocity);
  }
  // todo: maybe let different slots receive MIDI on different channels
  // and/or don't override the noteOn/etc. functions here but rather let the MIDI events also
  // pass through the modules in series. most modules just pass them through, but we can also
  // have MIDI effects such as appregiators and sequencers which modify the sequence and pass
  // the modified sequence to the next module - we could have an appregiator in front of a
  // synth, for example

  // all synthesizer modules should pass through the incoming audio and add their own signal
  // (unless the use it inside for their own signal processing) -> this allows for layering
}

void AudioModuleChain::noteOff(int noteNumber)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOff(noteNumber);
  }
}

void AudioModuleChain::setMidiController(int controllerNumber, float controllerValue)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->setMidiController(controllerNumber, controllerValue);
  }
}

void AudioModuleChain::setPitchBend(int pitchBendValue)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->setPitchBend(pitchBendValue);
  }
}

void AudioModuleChain::reset()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(modules); i++)
    modules[i]->reset();
}

XmlElement* AudioModuleChain::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement *xml = AudioModule::getStateAsXml(stateName, markAsClean);
  xml->setAttribute("ActiveSlot", activeSlot+1);
  for(int i = 0; i < size(modules); i++){
    juce::String typeString = AudioModuleFactory::getModuleType(modules[i]);
    XmlElement *child = new XmlElement("Slot");
    child->setAttribute("Type", typeString);
    //child->setAttribute("Bypass", isSlotBypassed(i)); // add later
    child->addChildElement(modules[i]->getStateAsXml(typeString, markAsClean));
    xml->addChildElement(child);
  }
  return xml;
}

void AudioModuleChain::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean); // actually does nothing?
  activeSlot = -1;  // i think, that's not necessary - should be already -1
  int tmpActiveSlot = xmlState.getIntAttribute("ActiveSlot", 1) - 1;
  clearModulesArray();
  int i = 0;
  forEachXmlChildElementWithTagName(xmlState, slotState, "Slot"){
    juce::String type = slotState->getStringAttribute("Type");
    if(i == tmpActiveSlot)        // hack: we set it before adding the module, so the editor
      activeSlot = tmpActiveSlot; // retrieves the correct value in the moduleAdded callback
    addModule(type);
    XmlElement *moduleState = slotState->getChildElement(0);
    modules[i]->setStateFromXml(*moduleState, "", markAsClean);
    i++;
  }
}

void AudioModuleChain::clearModulesArray()
{
  ScopedLock scopedLock(*lock);
  while(size(modules) > 0)
    deleteLastModule();
}

//=================================================================================================

AudioModuleChainEditor::AudioModuleChainEditor(jura::AudioModuleChain *moduleChainToEdit)
  : AudioModuleEditor(moduleChainToEdit)
{
  ScopedLock scopedLock(*lock);
  chain = moduleChainToEdit;
  setHeadlinePosition(TOP_LEFT);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  updateEditorArray();
  updateSelectorArray();
  updateActiveEditor();
  chain->addAudioModuleChainObserver(this);
  addChangeListener(this); // we listen to ourselves for deferred destruction of selectors
}

AudioModuleChainEditor::~AudioModuleChainEditor()
{
  ScopedLock scopedLock(*lock);
  chain->removeAudioModuleChainObserver(this);
  clearEditorArray();
}

AudioModuleEditor* AudioModuleChainEditor::getEditorForSlot(int index)
{
  ScopedLock scopedLock(*lock);
  if(size(editors) == 0 || index < 0)            // may happen during xml state recall
    return nullptr;
  jassert(index >= 0 && index < size(editors)); // index out of range
  if(editors[index] == nullptr)
    editors[index] = chain->modules[index]->createEditor();
  return editors[index];
}

void AudioModuleChainEditor::replaceModule(int index, const juce::String& type)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(editors));  // index out of range
  if(!chain->isModuleOfType(index, type)){
    chain->replaceModule(index, type);            // will call audioModuleWillBeDeleted
    //AudioModule* m = chain->getModuleAt(index);   // can be 0, if dummy module was placed at end
    updateEditorArray();
    index = chain->activeSlot;
    editors[index] = getEditorForSlot(index);
    updateActiveEditor();
    scheduleSelectorArrayUpdate();                // deferred call to updateSelectorArray
                                                  // may be superfluous now
  }
}

void AudioModuleChainEditor::updateSelectorArray()
{
  ScopedLock scopedLock(*lock);
  int numModules   = size(chain->modules);
  int numSelectors = size(selectors);
  AudioModuleSelector *s;

  // remove superfluous selectors:
  while(numSelectors > numModules){
    s = selectors[numSelectors-1];
    removeWidget(s, true, true);
    remove(selectors, numSelectors-1);
    numSelectors--;
  }

  // add required selectors:
  while(numModules > numSelectors){
    s = new AudioModuleSelector();
    s->setInterceptsMouseClicks(false, false); // we handle them 1st and possibly pass them through
    s->selectItemFromText(
      AudioModuleFactory::getModuleType(chain->modules[numSelectors]), false);
    s->registerComboBoxObserver(this);
    s->setDescriptionField(infoField); // somehow, this doesn't work
    addWidget(s);
    append(selectors, s);
    numSelectors++;
  }

  // without it, the selector for 1st slot is wrong after preset loading:
  if(numSelectors > 0 )
    selectors[0]->selectItemFromText(AudioModuleFactory::getModuleType(chain->modules[0]), false);
}

void AudioModuleChainEditor::updateEditorArray()
{
  ScopedLock scopedLock(*lock);
  int numModules = size(chain->modules);
  int numEditors = size(editors);
  AudioModuleEditor *e;

  // remove superfluous editors:
  while(numEditors > numModules){
    e = editors[numEditors-1];
    if(e == activeEditor)
      activeEditor = nullptr;
    delete e;
    remove(editors, numEditors-1);
    numEditors--;
  }

  // add placeholders for required selectors:
  while(numModules > numEditors){
    append(editors, (AudioModuleEditor*)nullptr);
    numEditors++;
  }
}

void AudioModuleChainEditor::updateActiveEditor()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor* tmpEditor = getEditorForActiveSlot();
  if(tmpEditor != activeEditor){
    removeChildEditor(activeEditor, false);
    addChildEditor(tmpEditor);
    activeEditor = tmpEditor;
    int w = max(240, activeEditor->getWidth());
    int h = max(180, activeEditor->getHeight());
    activeEditor->setBounds(leftColumnWidth, 0, w, h);
    setSize(w + leftColumnWidth, h + bottomRowHeight);
    resized(); // setSize will call resized only if the size actually changes but we need to make
               // sure that it always gets called to arrange the selectors
  }
}

void AudioModuleChainEditor::mouseDown(const MouseEvent &e)
{
  //ScopedLock scopedLock(*lock); // blocks audio when popup is open
  int i = chain->activeSlot;
  juce::Rectangle<int> rect = selectors[i]->getBounds();
  if(rect.contains(e.x, e.y)){
    // click was on active slot selector - pass event through:
    selectors[i]->mouseDown(e.getEventRelativeTo(selectors[i]));
  }
  else{
    for(i = 0; i < size(selectors); i++){
      rect = selectors[i]->getBounds();
      if(rect.contains(e.x, e.y)){
        // click was on inactive slot selector - activate:
        chain->activeSlot = i;
        updateActiveEditor();
        repaint();
      }
    }
  }
}

void AudioModuleChainEditor::resized()
{
  ScopedLock scopedLock(*lock);

  Editor::resized();

  int x, y, w, h, dy, margin;
  margin = 4;
  x = margin;
  y = getHeadlineBottom() + margin;
  w = leftColumnWidth - x - margin;
  h = 16;
  stateWidgetSet->setBounds(x, y, w, 32);

  // arrange selectors:
  y  = getPresetSectionBottom() + margin;
  dy = h-2;
  for(int i = 0; i < size(selectors); i++){
    selectors[i]->setBounds(x, y, w, h);
    y += dy;
  }

  // set up bounds of the editor for the active module:
  if(activeEditor != nullptr){
    y = 0;
    x = leftColumnWidth;
    w = getWidth()  - leftColumnWidth;
    h = getHeight() - bottomRowHeight;
    activeEditor->setBounds(x, y, w, h);
  }

  // set up bounds for inherited widgets:
  x = margin;
  y = getHeight() - bottomRowHeight;
  w = getWidth() - x - margin;
  infoField->setBounds(x, y, w, bottomRowHeight);

  x = getWidth() - 108;
  w = getWidth() - x - margin;
  webLink->setBounds(x, y+3, w, bottomRowHeight);

  //// color setup doesn't work properly yet - the active selector rectangle is drawn on top of it
  //// and when doing other stuff while it's open there can be access violations - so the button
  //// setup code is commented out:
  //int buttonWidth = 40;
  //x = leftColumnWidth - buttonWidth - margin;
  //y = margin;
  //setupButton->setBounds(x, y, buttonWidth, 16);

  // If this is a AudioModuleChain wrapped into an AudioPlugIn, we want to resize the whole parent
  // window as well:
  Component *parent =	getParentComponent();
  if(dynamic_cast<AudioPluginEditor*>(parent))
    parent->setSize(getWidth(), getHeight());
}

void AudioModuleChainEditor::paintOverChildren(Graphics& g)
{
  // highlight active slot by drawing a rectangle around it:
  ScopedLock scopedLock(*lock);
  if(size(selectors) == 0)   // occurs during state recall
    return;
  g.setColour(Colours::black);
  //g.setColour(Colours::darkred);
  //g.setColour(selectors[chain->activeSlot]->getSpecialColour1());
  juce::Rectangle<int> rect = selectors[chain->activeSlot]->getBounds();
  g.drawRect(rect, 2);  // 2nd param: thickness
}

void AudioModuleChainEditor::rComboBoxChanged(RComboBox* box)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < size(selectors); i++){
    if(box == selectors[i]){
      replaceModule(i, box->getSelectedItemText());
    }
  }
}

void AudioModuleChainEditor::changeListenerCallback(ChangeBroadcaster *source)
{
  ScopedLock scopedLock(*lock);
  if(source == this)
  {
    updateSelectorArray();
    resized();  // to arrange selectors
  }
  else
    AudioModuleEditor::changeListenerCallback(source);
}

void AudioModuleChainEditor::audioModuleWasAdded(AudioModuleChain *chain,
  AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  updateEditorArray();
  updateActiveEditor();
  updateSelectorArray();
}

void AudioModuleChainEditor::audioModuleWillBeDeleted(AudioModuleChain *chain,
  AudioModule *module, int index)
{
  ScopedLock scopedLock(*lock);
  deleteEditor(index);
  scheduleSelectorArrayUpdate();
}

void AudioModuleChainEditor::audioModuleWasReplaced(AudioModuleChain *chain,
  AudioModule *oldModule, AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*lock);
  deleteEditor(index);
}

void AudioModuleChainEditor::scheduleSelectorArrayUpdate()
{
  sendChangeMessage();
  // we will receive the message ourselves which causes a call to updateSelectorArray()
}

void AudioModuleChainEditor::deleteEditor(int index)
{
  ScopedLock scopedLock(*lock);
  jassert(index >= 0 && index < size(editors)); // index out of range
  if(activeEditor == editors[index])
    activeEditor = nullptr;
  delete editors[index];
  editors[index] = nullptr;
}

void AudioModuleChainEditor::clearEditorArray()
{
  ScopedLock scopedLock(*lock);
  activeEditor = nullptr;
  for(int i = 0; i < size(editors); i++)
    delete editors[i];
  editors.clear();
}
