BlepOscArrayModule::BlepOscArrayModule(CriticalSection *lockToUse)
  : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("BlepOscArray");
  setModuleName("OscArray");
}

AudioModuleEditor* BlepOscArrayModule::createEditor(int type)
{
  return new jura::BlepOscArrayEditor(this); // get rid of passing the lock
}


//=================================================================================================

BlepOscArrayEditor::BlepOscArrayEditor(BlepOscArrayModule* oscArrayToEdit)
  : AudioModuleEditor(oscArrayToEdit->lock), oscArrayModule(oscArrayToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...
}

void BlepOscArrayEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}