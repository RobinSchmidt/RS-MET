// construction/destruction:

NewSynthAudioModule::NewSynthAudioModule(CriticalSection *lockToUse) 
  : AudioModuleWithMidiIn(lockToUse)
{
  setModuleTypeName("NewSynth");
  //createParameters();

  addChildAudioModule(sourceModule = new QuadSourceAudioModule(lock, &synthCore.source));
  addChildAudioModule(filterModule = new DualFilterAudioModule(lock, &synthCore.filter));
  addChildAudioModule(modulatorsModule = new PolyModulatorsAudioModule(lock));
}

//=================================================================================================

NewSynthEditor::NewSynthEditor(CriticalSection* lockToUse, NewSynthAudioModule* synthToEdit)
  : AudioModuleEditor(lockToUse), synthModule(synthToEdit)
{
  /*
  ScopedLock scopedLock(*lock);
  addChildEditor( sourceEditor = new QuadSourceEditor(lock) );
  addChildEditor( filterEditor = new DualFilterEditor(lock) );
  addChildEditor( modulatorsEditor = new PolyModulatorsEditor(lock) );
  setSize(720, 720);
  */
}

void NewSynthEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  //int y = getPresetSectionBottom();
  int y = 0;
  int w = getWidth();
  int h = getHeight() - y;

  h = h/3;

  sourceEditor->setBounds(x, y, w, h);
  filterEditor->setBounds(x, y, w, h);
  h = getHeight() - filterEditor->getBottom();
  modulatorsEditor->setBounds(x, y, w, h);
}