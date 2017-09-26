
void rsDeletionRequester::requestDeletion()
{
  deletor->deleteObject(this);
}

//-------------------------------------------------------------------------------------------------

rsGarbageCollector::rsGarbageCollector(int _cleanUpInterval)
{
  cleanUpInterval = _cleanUpInterval;
}

rsGarbageCollector::~rsGarbageCollector()
{
  deleteGarbage();
}

void rsGarbageCollector::setCleanUpInterval(int newInterval)
{
  cleanUpInterval = newInterval;
  if(isTimerRunning())
    startTimer(cleanUpInterval);
}

void rsGarbageCollector::deleteObject(rsDeletionRequester* objectToDisposeOf)
{
  garbage.push_back(objectToDisposeOf);
  if(!isTimerRunning()) // bcs we don't want to reset it, if it is already running
    startTimer(cleanUpInterval);
}

void rsGarbageCollector::timerCallback()
{
  deleteGarbage();
}

void rsGarbageCollector::deleteGarbage()
{
  for(size_t i = 0; i < garbage.size(); i++)
    delete garbage[i];
  garbage.clear();
  stopTimer(); // no more timer callbacks needed until new garbage is added
}