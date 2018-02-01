DualFilterAudioModule::DualFilterAudioModule(
  CriticalSection *lockToUse/*, rosic::rsDualFilter* coreToUse*/)
  : AudioModule(lockToUse)/*, filterCore(coreToUse)*/
{
  addChildAudioModule(leftModule  = new PolySlotAudioModule(lock));
  addChildAudioModule(rightModule = new PolySlotAudioModule(lock));
  addChildAudioModule(vectorMixerModule = new PolyAudioModule(lock));
}

AudioModuleEditor* DualFilterAudioModule::createEditor()
{
  return new jura::DualFilterEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

DualFilterEditor::DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit)
  : AudioModuleEditor(lockToUse), filterModule(filterToEdit)
{
  ScopedLock scopedLock(*lock);
  addChildEditor(leftEditor  = new PolySlotEditor(lock, filterModule->leftModule));
  addChildEditor(rightEditor = new PolySlotEditor(lock, filterModule->rightModule));
  addChildEditor(vectorMixerEditor = new AudioModuleEditor(filterModule->vectorMixerModule));
}

void DualFilterEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  w = h;
  x = (getWidth() - w) / 2;
  vectorMixerEditor->setBounds(x, y, w, h);
}