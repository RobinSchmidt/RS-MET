
AudioModuleSlotArray::AudioModuleSlotArray(CriticalSection* lockToUse,
  MetaParameterManager* metaManagerToUse)
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse)
{


  int dummy = 0;
}

AudioModuleSlotArray::~AudioModuleSlotArray()
{

}

void AudioModuleSlotArray::addModuleSlotObserver(AudioModuleSlotObserver *observerToAdd)
{
  ScopedLock scopedLock(*lock);
  appendIfNotAlreadyThere(observers, observerToAdd);
}

void AudioModuleSlotArray::removeModuleSlotObserver(AudioModuleSlotObserver *observerToRemove)
{
  ScopedLock scopedLock(*lock);
  removeFirstOccurrence(observers, observerToRemove);
}