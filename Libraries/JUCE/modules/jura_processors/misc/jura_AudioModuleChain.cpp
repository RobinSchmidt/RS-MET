// some little helper/convenience functions to deal with std::vectors (move to RAPT, maybe don't 
// inline them all)

template<class T>
inline int size(const vector<T>& v)
{
  return (int)v.size();
}

template<class T>
inline void append(vector<T>& v, T newElement)
{
  v.push_back(newElement);
}

template<class T>
inline void remove(vector<T>& v, int index)
{
  v.erase(v.begin() + index);
}

template<class T>
inline void removeFirstOccurrence(vector<T>& v, T elementToRemove)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToRemove){
      remove(v, i);
      return;
    }
}

template<class T>
inline bool contains(vector<T>& v, T elementToCheckFor)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToCheckFor)
      return true;
  return false;
}

template<class T>
inline void appendIfNotAlreadyThere(vector<T>& v, T newElement)
{
  if(!contains(v, newElement))
    append(v, newElement);
}


//-------------------------------------------------------------------------------------------------

AudioModule* AudioModuleFactory::createModule(const String& type, CriticalSection *lock)
{
  if(type == "None")         return new DummyModule( lock);
  if(type == "PhaseScope")   return new PhaseScope(  lock);
  //if(type == "PhaseScope2")  return new PhaseScope2( lock);
  if(type == "Enveloper")    return new Enveloper(   lock);
  if(type == "Ladder")       return new Ladder(      lock);
  if(type == "PhasorFilter") return new PhasorFilter(lock);

  jassertfalse;  // unknown module type requested
  return nullptr;
}

String AudioModuleFactory::getModuleType(AudioModule *m)
{
  if(dynamic_cast<DummyModule*>  (m)) return "None";
  //if(dynamic_cast<PhaseScope2*>  (m)) return "PhaseScope2"; // check subclasse before...
  if(dynamic_cast<PhaseScope*>   (m)) return "PhaseScope";  // ...their superclasses
  if(dynamic_cast<Enveloper*>    (m)) return "Enveloper";
  if(dynamic_cast<Ladder*>       (m)) return "Ladder";
  if(dynamic_cast<PhasorFilter*> (m)) return "PhasorFilter";

  jassertfalse;  // unknown module type was passed
  return "UnknownType";
}

StringArray AudioModuleFactory::getAvailableModuleTypes()
{
  StringArray a;
  a.add("None");        // maybe use "Empty" instead of "None"
  a.add("PhaseScope");
  //a.add("PhaseScope2");
  a.add("Enveloper");
  a.add("Ladder");
  a.add("PhasorFilter");
  return a;
}

//=================================================================================================

AudioModuleSelector::AudioModuleSelector() : RComboBox("ModuleSelector") 
{
  setDescription("Select module type");
  StringArray a = AudioModuleFactory::getAvailableModuleTypes();
  for(int i = 0; i < a.size(); i++)
    addItem(i, a[i]); 

  // maybe later we should use a subclass of TreeView and organize the modules in groups: 
  // Analysis, Filter, etc.
}

//=================================================================================================

AudioModuleChain::AudioModuleChain(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "Chainer";
  setActiveDirectory(getApplicationDirectory() + "/ChainerPresets");
  addEmptySlot();
}

AudioModuleChain::~AudioModuleChain()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    delete modules[i];
}

void AudioModuleChain::addEmptySlot()
{
  addModule("None");
}

void AudioModuleChain::addModule(const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  AudioModule *m = AudioModuleFactory::createModule(type, plugInLock);
  append(modules, m);
  sendAudioModuleWasAddedNotification(m, size(modules)-1);
}

void AudioModuleChain::deleteModule(int index)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < modules.size()); // index out of range
  if(activeSlot == index)
    activeSlot--;
  sendAudioModuleWillBeDeletedNotification(modules[index], index);
  delete modules[index];
  remove(modules, index);
}

void AudioModuleChain::deleteLastModule()
{
  ScopedLock scopedLock(*plugInLock);
  deleteModule(size(modules) - 1);
}

void AudioModuleChain::replaceModule(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < modules.size()); // index out of range
  if(!isModuleOfType(index, type)){              // replace only, if new type is different
    AudioModule* oldModule = modules[index];
    AudioModule* newModule = AudioModuleFactory::createModule(type, plugInLock);
    newModule->setSampleRate(sampleRate);
    modules[index] = newModule;
    sendAudioModuleWasReplacedNotification(oldModule, newModule, index);
    delete oldModule;
    activeSlot = index;
    ensureOneEmptySlotAtEnd();
  }
}

bool AudioModuleChain::isModuleOfType(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < modules.size()); // index out of range
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
  ScopedLock scopedLock(*plugInLock);

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
  ScopedLock scopedLock(*plugInLock);
  appendIfNotAlreadyThere(observers, observerToAdd);
}

void AudioModuleChain::removeAudioModuleChainObserver(AudioModuleChainObserver *observerToRemove)
{
  ScopedLock scopedLock(*plugInLock);
  removeFirstOccurrence(observers, observerToRemove);
}

void AudioModuleChain::sendAudioModuleWasAddedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWasAdded(this, module, index);
}

void AudioModuleChain::sendAudioModuleWillBeDeletedNotification(AudioModule *module, int index)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < size(observers); i++)
    observers[i]->audioModuleWillBeDeleted(this, module, index);
}

void AudioModuleChain::sendAudioModuleWasReplacedNotification(AudioModule *oldModule, 
  AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
  jassert(numChannels == 2);
  for(int i = 0; i < modules.size(); i++)
    modules[i]->processBlock(inOutBuffer, numChannels, numSamples);
}

void AudioModuleChain::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);
  sampleRate = newSampleRate;
  for(int i = 0; i < modules.size(); i++)
    modules[i]->setSampleRate(sampleRate);
}

void AudioModuleChain::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++){
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
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOff(noteNumber);
  }
}

void AudioModuleChain::setMidiController(int controllerNumber, float controllerValue)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->setMidiController(controllerNumber, controllerValue);
  }
}

void AudioModuleChain::reset()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    modules[i]->reset();
}

XmlElement* AudioModuleChain::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);
  XmlElement *xml = AudioModule::getStateAsXml(stateName, markAsClean);
  xml->setAttribute("ActiveSlot", activeSlot+1);
  for(int i = 0; i < size(modules); i++){
    String typeString = AudioModuleFactory::getModuleType(modules[i]);
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
  ScopedLock scopedLock(*plugInLock);
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean); // actually does nothing?
  activeSlot = -1;  // i think, that's not necessary - should be already -1
  int tmpActiveSlot = xmlState.getIntAttribute("ActiveSlot", 1) - 1;
  clearModulesArray();
  int i = 0;
  forEachXmlChildElementWithTagName(xmlState, slotState, "Slot"){
    String type = slotState->getStringAttribute("Type");
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
  ScopedLock scopedLock(*plugInLock);
  while(size(modules) > 0)
    deleteLastModule();
}

//=================================================================================================

AudioModuleChainEditor::AudioModuleChainEditor(jura::AudioModuleChain *moduleChainToEdit)
  : AudioModuleEditor(moduleChainToEdit)
{
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
  chain->removeAudioModuleChainObserver(this);
  clearEditorArray();
}

AudioModuleEditor* AudioModuleChainEditor::getEditorForSlot(int index)
{
  ScopedLock scopedLock(*plugInLock);
  if(size(editors) == 0 || index < 0)            // may happen during xml state recall
    return nullptr; 
  jassert(index >= 0 && index < editors.size()); // index out of range
  if(editors[index] == nullptr)
    editors[index] = chain->modules[index]->createEditor();  
  return editors[index];
}

void AudioModuleChainEditor::replaceModule(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < editors.size());  // index out of range
  if(!chain->isModuleOfType(index, type)){
    chain->replaceModule(index, type);            // will call audioModuleWillBeDeleted
    AudioModule* m = chain->getModuleAt(index);   // can be 0, if dummy module was placed at end
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
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
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
  //ScopedLock scopedLock(*plugInLock); // blocks audio when popup is open
  int i = chain->activeSlot;
  Rectangle<int> rect = selectors[i]->getBounds();
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
  ScopedLock scopedLock(*plugInLock);

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
  for(int i = 0; i < selectors.size(); i++){
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
  ScopedLock scopedLock(*plugInLock);
  if(size(selectors) == 0)   // occurs during state recall
    return;
  g.setColour(Colours::black);
  //g.setColour(Colours::darkred);
  //g.setColour(selectors[chain->activeSlot]->getSpecialColour1());
  Rectangle<int> rect = selectors[chain->activeSlot]->getBounds();
  g.drawRect(rect, 2);  // 2nd param: thickness
}

void AudioModuleChainEditor::rComboBoxChanged(RComboBox* box)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < selectors.size(); i++){
    if(box == selectors[i]){
      replaceModule(i, box->getSelectedItemText());
    }
  }
}

void AudioModuleChainEditor::changeListenerCallback(ChangeBroadcaster *source)
{
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
  updateEditorArray();
  updateActiveEditor();
  updateSelectorArray();
}

void AudioModuleChainEditor::audioModuleWillBeDeleted(AudioModuleChain *chain,
  AudioModule *module, int index)
{
  ScopedLock scopedLock(*plugInLock);
  deleteEditor(index);
  scheduleSelectorArrayUpdate();
}

void AudioModuleChainEditor::audioModuleWasReplaced(AudioModuleChain *chain,
  AudioModule *oldModule, AudioModule *newModule, int index)
{
  ScopedLock scopedLock(*plugInLock);
  deleteEditor(index);
}

void AudioModuleChainEditor::scheduleSelectorArrayUpdate()
{
  sendChangeMessage(); 
  // we will receive the message ourselves which causes a call to updateSelectorArray()
}

void AudioModuleChainEditor::deleteEditor(int index)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < editors.size()); // index out of range
  if(activeEditor == editors[index])
    activeEditor = nullptr;
  delete editors[index];
  editors[index] = nullptr;
}

void AudioModuleChainEditor::clearEditorArray()
{
  ScopedLock scopedLock(*plugInLock);
  activeEditor = nullptr;
  for(int i = 0; i < editors.size(); i++)
    delete editors[i];
  editors.clear();
}
