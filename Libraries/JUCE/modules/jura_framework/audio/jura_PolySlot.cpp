AudioModuleEditor* PolySlotAudioModule::createEditor()
{
  return new jura::PolySlotEditor(lock, this); // get rid of passing the lock
}

//=================================================================================================

PolySlotEditor::PolySlotEditor(CriticalSection* lockToUse, PolySlotAudioModule* slotToEdit)
  : AudioModuleEditor(lockToUse), slotModule(slotToEdit)
{

}