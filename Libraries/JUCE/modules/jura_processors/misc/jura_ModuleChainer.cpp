// some little helper/convenience functions to deal with std::vectors (move to RAPT)

template<class T>
inline int size(const vector<T>& v)
{
  return (int)v.size();
}

template<class T>
inline void append(vector<T>& v, T& newElement)
{
  v.push_back(newElement);
}

template<class T>
inline void remove(vector<T>& v, int index)
{
  v.erase(v.begin() + index);
}

//-------------------------------------------------------------------------------------------------

AudioModule* AudioModuleFactory::createModule(const String& type, CriticalSection *lock)
{
  if(type == "None")         return new DummyModule( lock);
  if(type == "PhaseScope")   return new PhaseScope(  lock);
  if(type == "Enveloper")    return new Enveloper(   lock);
  if(type == "Ladder")       return new Ladder(      lock);
  if(type == "PhasorFilter") return new PhasorFilter(lock);

  jassertfalse;  // unknown module type requested
  return nullptr;
}

String AudioModuleFactory::getModuleType(AudioModule *m)
{
  if(dynamic_cast<DummyModule*>  (m)) return "None";
  if(dynamic_cast<PhaseScope*>   (m)) return "PhaseScope";
  if(dynamic_cast<Enveloper*>    (m)) return "Enveloper";
  if(dynamic_cast<Ladder*>       (m)) return "Ladder";
  if(dynamic_cast<PhasorFilter*> (m)) return "PhasorFilter";

  jassertfalse;  // unknown module type was passed
  return String();
}

StringArray AudioModuleFactory::getAvailableModuleTypes()
{
  StringArray a;
  a.add("None");        // maybe use "Empty" instead of "None"
  a.add("PhaseScope");
  a.add("Enveloper");
  a.add("Ladder");
  a.add("PhasorFilter");
  return a;
}

//=================================================================================================

AudioModuleSelector::AudioModuleSelector() : RComboBox("ModuleSelector") 
{
  StringArray a = AudioModuleFactory::getAvailableModuleTypes();
  for(int i = 0; i < a.size(); i++)
    addItem(i, a[i]); 

  // maybe later we should use a subclass of TreeView and organize the modules in groups: 
  // Analysis, Filter, etc.
}

//=================================================================================================

ModuleChainer::ModuleChainer(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "Chainer";
  setActiveDirectory(getApplicationDirectory() + "/ChainerPresets");
  addEmptySlot();
}

ModuleChainer::~ModuleChainer()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    delete modules[i];
}

void ModuleChainer::addEmptySlot()
{
  addModule("None");
}

void ModuleChainer::addModule(const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  AudioModule *m = AudioModuleFactory::createModule(type, plugInLock);
  append(modules, m);
}

void ModuleChainer::replaceModule(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < modules.size()); // index out of range
  if(!isModuleOfType(index, type)){              // replace only, if new type is different
    delete modules[index];
    modules[index] = AudioModuleFactory::createModule(type, plugInLock);
    modules[index]->setSampleRate(sampleRate);
    activeSlot = index;
    ensureOneEmptySlotAtEnd();
  }
}

void ModuleChainer::removeLastModule()
{
  ScopedLock scopedLock(*plugInLock);
  int index = size(modules) - 1;
  if(activeSlot == index)
    activeSlot--;
  delete modules[index];
  remove(modules, index);
}

bool ModuleChainer::isModuleOfType(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < modules.size()); // index out of range
  return type == AudioModuleFactory::getModuleType(modules[index]);
}

void ModuleChainer::ensureOneEmptySlotAtEnd()
{
  ScopedLock scopedLock(*plugInLock);

  // if the last module/slot is not the "empty" dummy module, add another slot:
  if(!isModuleOfType(size(modules)-1, "None"))
    addEmptySlot();

  // remove superfluous empty slots at end:
  while(modules.size() > 1 && isModuleOfType(size(modules)-1, "None") 
                           && isModuleOfType(size(modules)-2, "None"))
  { // if the last two slots are empty, remove the last
    removeLastModule();
  }
}

// overrides:

AudioModuleEditor* ModuleChainer::createEditor()
{
  return new ModuleChainerEditor(this);
}

void ModuleChainer::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(numChannels == 2);
  for(int i = 0; i < modules.size(); i++)
    modules[i]->processBlock(inOutBuffer, numChannels, numSamples);
}

void ModuleChainer::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);
  sampleRate = newSampleRate;
  for(int i = 0; i < modules.size(); i++)
    modules[i]->setSampleRate(sampleRate);
}

void ModuleChainer::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOn(noteNumber, velocity);
  }
  // todo: maybe let different slots receive MIDI on different channels
  // and/or don't override the noteOn/etc. functions here but rather let the MIDI events also
  // apss through the modules in series. most modules just pass them through, but we can also
  // have MIDI effects such as appregiators and sequencers which modify the sequence and pass
  // the modified sequence to the next module - we could have an appregiator in front of a 
  // synth, for example

  // all synthesizer modules should pass through the incoming audio and add their own signal
  // (unless the use it inside for their own signal processing) -> this allows for layering
}

void ModuleChainer::noteOff(int noteNumber)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++){
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOff(noteNumber);
  }
}

void ModuleChainer::reset()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    modules[i]->reset();
}

XmlElement* ModuleChainer::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  return AudioModule::getStateAsXml(stateName, markAsClean); // preliminary
  // We need to loop over the slots (i.e. the modules array) and for each slot N create a child
  // xml with name SlotN, with attributes like Type="Ladder", Bypass="0", etc and the state
  // of the module as child element.
}

void ModuleChainer::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean); // preliminary
  // We need to clear the modules array, then loop over the child xml elements, infer the Type 
  // attribute and create and add a module of appropriate type to the array and then set up it's
  // state from the child xml-state.
}

//=================================================================================================

ModuleChainerEditor::ModuleChainerEditor(jura::ModuleChainer *moduleChainerToEdit)
  : AudioModuleEditor(moduleChainerToEdit)
{
  ScopedLock scopedLock(*plugInLock);
  chainer = moduleChainerToEdit;
  setHeadlinePosition(TOP_LEFT);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  initEditorArray();
  createSelectorWidgets();
  updateActiveEditor();
}

ModuleChainerEditor::~ModuleChainerEditor()
{
  clearEditorArray();
}

AudioModuleEditor* ModuleChainerEditor::getEditorForSlot(int index)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < editors.size()); // index out of range
  if(editors[index] == nullptr)
    editors.set(index, chainer->modules[index]->createEditor());  
  return editors[index];
}

void ModuleChainerEditor::replaceModule(int index, const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < editors.size()); // index out of range
  if(!chainer->isModuleOfType(index, type)){
    deleteEditor(index);
    chainer->replaceModule(index, type);

    updateEditorArray();
    index = chainer->activeSlot;
    editors.set(index, getEditorForSlot(index));

    updateSelectorArray();
    updateActiveEditor();
  }
}

void ModuleChainerEditor::updateSelectorArray()
{
  ScopedLock scopedLock(*plugInLock);
  int numModules   = size(chainer->modules);
  int numSelectors = selectors.size();
  AudioModuleSelector *s;

  // remove superfluous selectors:
  while(numSelectors > numModules){
    s = selectors[numSelectors-1];
    removeWidget(s, true, true);
    selectors.remove(numSelectors-1);
    numSelectors--;
  }

  // add required selectors:
  while(numModules > numSelectors){
    s = new AudioModuleSelector();
    s->selectItemFromText(
      AudioModuleFactory::getModuleType(chainer->modules[numSelectors]), false);
    s->registerComboBoxObserver(this);
    addWidget(s);
    selectors.add(s);
    numSelectors++;
  }
}

void ModuleChainerEditor::updateEditorArray()
{
  ScopedLock scopedLock(*plugInLock);
  int numModules = size(chainer->modules);
  int numEditors = editors.size();
  AudioModuleEditor *e;

  // remove superfluous editors:
  while(numEditors > numModules){
    e = editors[numEditors-1];
    if(e == activeEditor)
      activeEditor = nullptr;
    delete e;
    editors.remove(numEditors-1);
    numEditors--;
  }

  // add placeholders for required selectors:
  while(numModules > numEditors){
    editors.add(nullptr);
    numEditors++;
  }
}

void ModuleChainerEditor::updateActiveEditor()
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

void ModuleChainerEditor::resized()
{
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

  // maybe, we could have bypass switches for each plugin
  // arrange setup button for color scheme, infoline, link, etc.

  // If this is a ModuleChainer wrapped into an AudioPlugIn, we want to resize the whole parent 
  // window as well:
  Component *parent =	getParentComponent();
  if(dynamic_cast<AudioPluginEditor*>(parent))
    parent->setSize(getWidth(), getHeight());
}

void ModuleChainerEditor::rComboBoxChanged(RComboBox* box)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < selectors.size(); i++){
    if(box == selectors[i]){
      replaceModule(i, box->getSelectedItemText());
    }
  }
}

void ModuleChainerEditor::deleteEditor(int index)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(index >= 0 && index < editors.size()); // index out of range
  if(activeEditor == editors[index])
    activeEditor = nullptr;
  delete editors[index];
  editors.set(index, nullptr);
}

void ModuleChainerEditor::clearEditorArray()
{
  ScopedLock scopedLock(*plugInLock);
  activeEditor = nullptr;
  for(int i = 0; i < editors.size(); i++)
    delete editors[i];
  editors.clear();
}

void ModuleChainerEditor::initEditorArray()
{
  ScopedLock scopedLock(*plugInLock);
  clearEditorArray();
  for(int i = 0; i < chainer->modules.size(); i++)
    editors.add(nullptr);
}

void ModuleChainerEditor::createSelectorWidgets()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < chainer->modules.size(); i++){
    AudioModuleSelector *s = new AudioModuleSelector();
    s->selectItemFromText(AudioModuleFactory::getModuleType(chainer->modules[i]), false);
    s->registerComboBoxObserver(this);
    addWidget(s);
    selectors.add(s);
  }

  // this method may be superfluous now - we can use updateSelectorArray()
}