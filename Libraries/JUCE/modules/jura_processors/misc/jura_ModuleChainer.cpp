
ModuleChainer::ModuleChainer(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "Chainer";
  setActiveDirectory(getApplicationDirectory() + "/ChainerPresets");

  addModule("None"); // always have at least one dummy module in the chain
}

ModuleChainer::~ModuleChainer()
{
  ScopedLock scopedLock(*plugInLock);
  for(int i = 0; i < modules.size(); i++)
    delete modules[i];
}

AudioModule* ModuleChainer::createModule(const String& type)
{
  ScopedLock scopedLock(*plugInLock);

  if(type == "None")         return new DummyModule( plugInLock);
  if(type == "PhaseScope")   return new PhaseScope(  plugInLock);
  if(type == "Enveloper")    return new Enveloper(   plugInLock);
  if(type == "Ladder")       return new Ladder(      plugInLock);
  if(type == "PhasorFilter") return new PhasorFilter(plugInLock);
  // add more module types here...

  jassertfalse;  // unknown module type requested
  return nullptr;
}

void ModuleChainer::addModule(const String& type)
{
  ScopedLock scopedLock(*plugInLock);
  AudioModule *m = createModule(type);
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

AudioModuleSelector::AudioModuleSelector() : RComboBox("ModuleSelector") 
{
  // add menu items for the different module types:
  int i = 0;
  addItem(i, "None");         i++;
  addItem(i, "Ladder");       i++;
  addItem(i, "PhasorFilter"); i++;
  addItem(i, "PhaseScope");   i++;
  addItem(i, "Enveloper");    i++;

  // maybe later we should use a subclass of TreeView and have groups: Analysis, Filter, etc.
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
  for(int i = 0; i < chainer->modules.size(); i++)
  {
    AudioModuleSelector *s = new AudioModuleSelector();

    // set up the selector to reflect the correct module type
    // ...

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
