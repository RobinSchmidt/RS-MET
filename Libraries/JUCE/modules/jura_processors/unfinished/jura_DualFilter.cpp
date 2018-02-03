DualFilterAudioModule::DualFilterAudioModule(
  CriticalSection *lockToUse/*, rosic::rsDualFilter* coreToUse*/)
  : AudioModule(lockToUse)/*, filterCore(coreToUse)*/
{
  setModuleTypeName("DualFilter");
  setModuleName("Filters");
  addChildAudioModule(leftModule  = new PolySlotAudioModule(lock));
  addChildAudioModule(rightModule = new PolySlotAudioModule(lock));
  addChildAudioModule(vectorMixerModule = new PolyAudioModule(lock));
}

AudioModuleEditor* DualFilterAudioModule::createEditor()
{
  return new jura::DualFilterEditor(lock, this); // get rid of passing the lock
}

void DualFilterAudioModule::setModuleFactory(AudioModuleFactory* newFactory)
{
  leftModule->setModuleFactory(newFactory);
  rightModule->setModuleFactory(newFactory);
}

//=================================================================================================

DualFilterEditor::DualFilterEditor(CriticalSection* lockToUse, DualFilterAudioModule* filterToEdit)
  : AudioModuleEditor(lockToUse), filterModule(filterToEdit)
{
  ScopedLock scopedLock(*lock);
  addChildEditor(leftEditor  = new PolySlotEditor(lock, filterModule->leftModule));
  addChildEditor(rightEditor = new PolySlotEditor(lock, filterModule->rightModule));
  addChildEditor(vectorMixerEditor = new AudioModuleEditor(filterModule->vectorMixerModule));
  //setPresetSectionPosition(INVISIBLE);
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

  w = vectorMixerEditor->getX() + 2;
  leftEditor->setBounds(0, y, w, h);

  x = vectorMixerEditor->getRight() - 2;
  w = getWidth() - x;
  rightEditor->setBounds(x, y, w, h);
}