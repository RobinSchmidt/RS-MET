#ifndef jura_StateManager_h
#define jura_StateManager_h

class StateManager;

class JUCE_API StateWatcher
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StateWatcher();

  /** Destructor. */
  virtual ~StateWatcher();

  /** The callback that you must implement to respond to events that change the state dirty
  flag. */
  virtual void stateDirtyFlagChanged(StateManager* stateManager) = 0;

  /** Set up a pointer to the watched StateManager. */
  virtual void setStateManagerToWatch(StateManager* stateManagerToWatch)
  {
    watchedStateManager = stateManagerToWatch;
  }



  StateManager* watchedStateManager;

protected:

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class can remember the name of a current (module) preset and inform whether or not the
current parameter settings still reflect the preset - for this to work properly, subclasses should 
call the member markPresetAsDirty whenever some (preset relevant) parameter was changed. This will 
then trigger a callback to all the attached StateWatchers and their parents. */

class JUCE_API StateManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StateManager();

  /** Destructor. */
  virtual ~StateManager();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the name (i.e. the relative path) of the current preset) - this function is only
  for conviniently handling preset management in a plugIn-context. */
  virtual void setStateName(const juce::String& newStateName, bool markStateAsClean = true);

  /** Marks the state as clean (the current parameter settings reflect the preset. */
  virtual void markStateAsClean();

  /** Marks the state as dirty (the current parameter settings do not reflect the preset
  anymore, probably due to some tweaking). */
  virtual void markStateAsDirty();

  /** Causes this object to ignore calls to markStateAsDirty. At times it might be desirable to 
  turn the response to these calls temporarily off. */
  virtual void setIgnoreDirtification(bool shouldIgnore) { ignoreDirtification = shouldIgnore; }

  /** Should be overriden by subclasses and is intended to recall the state of the StateManager 
  (i.e. the settings of all relevant
  parameters) from an XmlElement.  */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) = 0;

  /** Adds a child StateManager to this one and set its parent to this. */
  virtual void addChildStateManager(StateManager *newChild);

  /** Removes a child StateManager from this one and set its parent to NULL. */
  virtual void removeChildStateManager(StateManager *childToRemove);

  /** Removes all the child-StateManagers (if any) frm this one and set theri parent to NULL. */
  virtual void removeAllChildStateManagers();

  /** Adds a watcher to this StateManager which will be notified when the dirty-flag of the state
  changes. */
  virtual void addStateWatcher(StateWatcher *watcherToAdd);

  /** Removes a watcher from this StateManager and returns true if there was actually something
  to remove (i.e. the to-removed StateWatcher was watching this object). */
  virtual bool removeStateWatcher(StateWatcher *watcherToRemove);

  /** Removes all watchers from this StateManager. */
  virtual void removeAllStateWatchers();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the name of the state. */
  virtual const juce::String& getStateName() const;

  /** Returns the name of the state with an added star if the state is dirty. */
  virtual const juce::String getStateNameWithStarIfDirty() const;

  /** Should be overriden by subclasses and is intended to return the state of the StateManager
  (i.e. the settings of all relevant parameters) in form of an XmlElement. The caller must take
  care to delete the XmlElement when it is not needed anymore.  */
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) = 0;

  /** Returns the state without changing its name or setting the dirty-flag. */
  virtual XmlElement* getStateAsXml() { return getStateAsXml(getStateName(), false); }

  /** Returns whether or not some parameter was changed since the last call to setPresetName()
  with ture as last argument. */
  virtual bool isStateDirty() const;

protected:

  bool          ignoreDirtification;
  bool          stateIsDirty;
  juce::String  stateName;
  StateManager* parent;
  juce::Array<StateManager*, CriticalSection> children;
  juce::Array<StateWatcher*, CriticalSection> watchers;

  juce_UseDebuggingNewOperator;
};

#endif 
