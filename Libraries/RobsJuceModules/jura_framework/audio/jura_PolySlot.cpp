AudioModuleEditor* PolySlotAudioModule::createEditor(int type)
{
  return new jura::PolySlotEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

PolySlotEditor::PolySlotEditor(CriticalSection* lockToUse, PolySlotAudioModule* slotToEdit)
  : AudioModuleEditor(lockToUse), slotModule(slotToEdit)
{
  setPresetSectionPosition(INVISIBLE);
  addWidget(moduleSelector = new AudioModuleSelector(slotModule->moduleFactory));
}

//void PolySlotEditor::paint(Graphics& g)
//{
//  g.fillAll(Colours::black); // just for test
//}

void PolySlotEditor::resized() 
{
  moduleSelector->setBounds(4, 4, getWidth()-8, 16);
  if(slotInsertEditor) {
    int y = moduleSelector->getBottom() + 4;
    slotInsertEditor->setBounds(0, y, getWidth(), getHeight()-y);
  }
  // maybe use the full size for the slotInsertEditor and intercept mouse clicks and open the
  // selector popup on (right) click
}