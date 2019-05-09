FooBarModule::FooBarModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("FooBar");
  setModuleName("FooBar1");
}

AudioModuleEditor* FooBarModule::createEditor(int type)
{
  return new jura::FooBarEditor(this);
}

//=================================================================================================

FooBarEditor::FooBarEditor(FooBarModule* fooBarToEdit)
  : AudioModuleEditor(fooBarToEdit->lock), fooBarModule(fooBarToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...
}

void FooBarEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}