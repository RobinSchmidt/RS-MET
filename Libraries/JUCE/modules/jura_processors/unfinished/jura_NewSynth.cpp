// construction/destruction:

NewSynthAudioModule::NewSynthAudioModule(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse), sourceFactory(lockToUse), filterFactory(lockToUse)
  , modulatorFactory(lockToUse)
{
  setModuleTypeName("NewSynth");
  addChildAudioModule(sourceModule = new QuadSourceAudioModule(lock/*, &synthCore.source*/));
  addChildAudioModule(filterModule = new DualFilterAudioModule(lock/*, &synthCore.filter*/));
  addChildAudioModule(modulatorsModule = new PolyModulatorsAudioModule(lock));

  sourceModule->setModuleFactory(&sourceFactory);
  filterModule->setModuleFactory(&filterFactory);
  modulatorsModule->setModuleFactory(&modulatorFactory);

  populateModuleFactories();
  //createParameters();
}

AudioModuleEditor* NewSynthAudioModule::createEditor()
{
  return new jura::NewSynthEditor(lock, this); // get rid of passing the lock
}

void NewSynthAudioModule::populateModuleFactories()
{
  // something to do...
}

//=================================================================================================

NewSynthEditor::NewSynthEditor(CriticalSection* lockToUse, NewSynthAudioModule* synthToEdit)
  : AudioModuleEditor(lockToUse), synthModule(synthToEdit)
{
  ScopedLock scopedLock(*lock);
  addChildEditor(sourceEditor = new QuadSourceEditor(lock, synthModule->sourceModule));
  addChildEditor(filterEditor = new DualFilterEditor(lock, synthModule->filterModule));
  addChildEditor(modulatorsEditor = new PolyModulatorsEditor(lock, synthModule->modulatorsModule));
  setSize(720, 720+26); // 26 = bottom of preset section + 4
}

void NewSynthEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  h = h/3;
  int dy = h-2;
  sourceEditor->setBounds(x, y, w, h); y += dy;
  filterEditor->setBounds(x, y, w, h); y += dy;
  h = getHeight() - y;
  modulatorsEditor->setBounds(x, y, w, h);
}