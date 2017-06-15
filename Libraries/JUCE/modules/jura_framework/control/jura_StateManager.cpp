
StateWatcher::StateWatcher()
{
  watchedStateManager = NULL;
}

StateWatcher::~StateWatcher()
{
  if( watchedStateManager != NULL )
  {
    watchedStateManager->removeStateWatcher(this);
    watchedStateManager = NULL;
  }
}

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StateManager::StateManager()
{
  parent              = NULL;
  //stateName           = NULL;
  stateName           = String::empty;
  stateIsDirty        = false;
  ignoreDirtification = false;
  //setStateName("Init", true);
}

StateManager::~StateManager()
{
  removeAllChildStateManagers();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void StateManager::setStateName(const String& newStateName, bool markStateAsClean_)
{
  stateName = newStateName;
  if( markStateAsClean_ == true )
  {
    markStateAsClean();
  }
}

void StateManager::markStateAsClean()
{
  stateIsDirty = false;
  int i;
  watchers.getLock().enter();
  for(i=0; i<watchers.size(); i++)
    watchers[i]->stateDirtyFlagChanged(this);
  watchers.getLock().exit();

  children.getLock().enter();
  for(i=0; i<children.size(); i++)
    children[i]->markStateAsClean();
  children.getLock().exit();
}

void StateManager::markStateAsDirty()
{
  if( ignoreDirtification )
    return;

  stateIsDirty = true;
  if( parent != NULL )
    parent->markStateAsDirty();

  watchers.getLock().enter();
  for(int i=0; i<watchers.size(); i++)
    watchers[i]->stateDirtyFlagChanged(this);
  watchers.getLock().exit();
}

void StateManager::addChildStateManager(StateManager *newChild)
{
  children.getLock().enter();
  children.addIfNotAlreadyThere(newChild);
  children.getLock().exit();
  newChild->parent = this;
}

void StateManager::removeChildStateManager(StateManager *childToRemove)
{
  children.getLock().enter();
  children.removeFirstMatchingValue(childToRemove);
  children.getLock().exit();
  childToRemove->parent = NULL;
}

void StateManager::removeAllChildStateManagers()
{
  children.getLock().enter();
  for(int i=0; i<children.size(); i++)
  {
    //StateManager* child = children[i];
    children[i]->parent = NULL;
  }
  children.getLock().exit();
  children.clear();
}

void StateManager::addStateWatcher(StateWatcher *watcherToAdd)
{
  watchers.getLock().enter();
  watchers.addIfNotAlreadyThere(watcherToAdd);
  watcherToAdd->setStateManagerToWatch(this);
  watchers.getLock().exit();
}

bool StateManager::removeStateWatcher(StateWatcher *watcherToRemove)
{
  watchers.getLock().enter();
  bool result = watchers.indexOf(watcherToRemove) != -1;
  jassert( result == true ); // trying to remove a watcher that is not watching this object?
  watchers.removeFirstMatchingValue(watcherToRemove);
  watcherToRemove->setStateManagerToWatch(NULL);
  watchers.getLock().exit();
  return result;
}

void StateManager::removeAllStateWatchers()
{
  watchers.getLock().enter();
  while( watchers.size() > 0 )
    removeStateWatcher( watchers[watchers.size()-1]  );
  watchers.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

const String& StateManager::getStateName() const
{
  return stateName;
}

const String StateManager::getStateNameWithStarIfDirty() const
{
  if( stateIsDirty == true && stateName != String::empty )
    return stateName + String("*");
  else
    return stateName;
}

bool StateManager::isStateDirty() const
{
  return stateIsDirty;
}

