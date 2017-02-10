
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
  a.add("None");
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

  addModule("None"); // always have at least one dummy module in the chain
  //addModule("Ladder"); // for test
}

ModuleChainer::~ModuleChainer()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    delete modules[i];
}

void ModuleChainer::addModule(const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  AudioModule *m = AudioModuleFactory::createModule(type, plugInLock);
  modules.add(m);
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
  for(int i = 0; i < modules.size(); i++)
    modules[i]->setSampleRate(newSampleRate);
}

void ModuleChainer::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
  {
    AudioModuleWithMidiIn *m = dynamic_cast<AudioModuleWithMidiIn*> (modules[i]);
    if(m != nullptr)
      m->noteOn(noteNumber, velocity);
  }
  // todo: maybe let different slots receive MIDI on different channels
}

void ModuleChainer::noteOff(int noteNumber)
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
  {
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

//=================================================================================================

ModuleChainerEditor::ModuleChainerEditor(jura::ModuleChainer *moduleChainerToEdit)
  : AudioModuleEditor(moduleChainerToEdit)
{
  ScopedLock scopedLock(*plugInLock);
  chainer = moduleChainerToEdit;
  setHeadlinePosition(TOP_LEFT);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  createWidgets();
  setSize(200, 100);
}

ModuleChainerEditor::~ModuleChainerEditor()
{
  // delete editors 
}

void ModuleChainerEditor::createWidgets()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < chainer->modules.size(); i++)
  {
    AudioModuleSelector *s = new AudioModuleSelector();
    s->selectItemFromText(AudioModuleFactory::getModuleType(chainer->modules[i]), false);
    addWidget(s);
    selectors.add(s);
  }
}

void ModuleChainerEditor::resized()
{
  Editor::resized();

  int x, y, w, h, dy, margin;
  margin = 4;
  x = margin;
  y = getHeadlineBottom() + margin;
  w = 160;
  h = 16;
  stateWidgetSet->setBounds(x, y, w, 32);

  // arrange selectors:
  y  = getPresetSectionBottom() + margin;
  dy = h;
  for(int i = 0; i < selectors.size(); i++)
  {
    selectors[i]->setBounds(x, y, w, h);
    y += dy;
  }

  // maybe, we could have bypass switches for each plugin
  // arrange setup button for color scheme, infoline, link, etc.
}
