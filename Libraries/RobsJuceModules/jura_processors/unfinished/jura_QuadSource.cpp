QuadSourceAudioModule::QuadSourceAudioModule(
  CriticalSection *lockToUse/*, rosic::rsQuadSource* coreToUse*/)
  : AudioModule(lockToUse)/*, sourceCore(coreToUse)*/
{
  setModuleTypeName("QuadSource");
  setModuleName("Sources");
  addChildAudioModule(topLeftModule     = new PolySlotAudioModule(lock));
  addChildAudioModule(topRightModule    = new PolySlotAudioModule(lock));
  addChildAudioModule(bottomLeftModule  = new PolySlotAudioModule(lock));
  addChildAudioModule(bottomRightModule = new PolySlotAudioModule(lock));
  addChildAudioModule(vectorMixerModule = new PolyAudioModule(lock));
}

AudioModuleEditor* QuadSourceAudioModule::createEditor(int type)
{
  return new jura::QuadSourceEditor(lock, this); // get rid of passing the lock
}

void QuadSourceAudioModule::setModuleFactory(AudioModuleFactory* newFactory)
{
  topLeftModule->setModuleFactory(newFactory);
  topRightModule->setModuleFactory(newFactory);
  bottomLeftModule->setModuleFactory(newFactory);
  bottomRightModule->setModuleFactory(newFactory);
}

//=================================================================================================

QuadSourceEditor::QuadSourceEditor(CriticalSection* lockToUse, QuadSourceAudioModule* sourceToEdit)
  : AudioModuleEditor(sourceToEdit), sourceModule(sourceToEdit)
{
  ScopedLock scopedLock(*lock);
  addChildEditor(topLeftEditor     = new PolySlotEditor(lock, sourceModule->topLeftModule));
  addChildEditor(topRightEditor    = new PolySlotEditor(lock, sourceModule->topRightModule));
  addChildEditor(bottomLeftEditor  = new PolySlotEditor(lock, sourceModule->bottomLeftModule));
  addChildEditor(bottomRightEditor = new PolySlotEditor(lock, sourceModule->bottomRightModule));
  addChildEditor(vectorMixerEditor = new AudioModuleEditor(sourceModule->vectorMixerModule));
}

void QuadSourceEditor::resized()
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

  h = h/2 + 1;
  w = vectorMixerEditor->getX() + 2;
  topLeftEditor->setBounds(0, y, w, h);

  x = vectorMixerEditor->getRight() - 2;
  w = getWidth() - x;
  topRightEditor->setBounds(x, y, w, h);

  y = topLeftEditor->getBottom() - 2;
  w = topLeftEditor->getWidth();
  h = getHeight() - y;
  bottomLeftEditor->setBounds(0, y, w, h);

  x = topRightEditor->getX();
  w = topRightEditor->getWidth();
  bottomRightEditor->setBounds(x, y, w, h);
}