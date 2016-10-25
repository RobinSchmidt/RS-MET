#include "rosof_AudioModule.h"
using namespace rosof;


//=========================================================================================================================================
// class AudioModuleDeletionWatcher

AudioModuleDeletionWatcher::~AudioModuleDeletionWatcher()
{
  for(int i = 0; i < watchedAudioModules.size(); i++ )
    watchedAudioModules[i]->deRegisterDeletionWatcher(this);
}

void AudioModuleDeletionWatcher::addWatchedAudioModule(AudioModule *moduleToBeWatched)
{
  if( moduleToBeWatched != NULL )
  {
    moduleToBeWatched->registerDeletionWatcher(this);
    watchedAudioModules.addIfNotAlreadyThere(moduleToBeWatched);
  }
}

void AudioModuleDeletionWatcher::removeWatchedAudioModule(AudioModule *moduleNotToBeWatchedAnymore)
{
  if( moduleNotToBeWatchedAnymore != NULL )
  {
    moduleNotToBeWatchedAnymore->deRegisterDeletionWatcher(this);
    watchedAudioModules.removeValue(moduleNotToBeWatchedAnymore);
  }
}


//=========================================================================================================================================
// class AudioModule

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

AudioModule::AudioModule(CriticalSection *newPlugInLock)
{
  plugInLock = newPlugInLock;

  ParameterObserver::localAutomationSwitch = true;  // activate automation for this instance
  wantsTempoSyncInfo                       = true;
  underlyingRosicInstrument                = NULL;  // by default, this does not wrap an instrument
  moduleName                               = juce::String(T("AudioModule"));
  //versionjuce::String                            = juce::String(T("1.0"));
  patchFormatIndex                         = 1;
  triggerInterval                          = 0.0;
  saveAndRecallState                       = true;
  initializeAutomatableParameters();
}

AudioModule::~AudioModule()
{
  ScopedLock scopedLock(*plugInLock);

  for(int i = 0; i < deletionWatchers.size(); i++)
    deletionWatchers[i]->audioModuleWillBeDeleted(this);

  childModules.getLock().enter();
  AudioModule* childModule = NULL;
  while( childModules.size() > 0 )
  {
    childModule = childModules[childModules.size()-1];
    removeChildAudioModule(childModule, true);
    //childModule->parent = NULL;
    //delete childModules[childModules.size()-1];
    //childModules.removeLast();
  }
  childModules.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);

  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->setSampleRate(newSampleRate);

  childModules.getLock().enter();
  for(int i=0; i<childModules.size(); i++)
    childModules[i]->setSampleRate(newSampleRate);
  childModules.getLock().exit();
}

void AudioModule::setBeatsPerMinute(double newBpm)
{
  ScopedLock scopedLock(*plugInLock);

  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->setBeatsPerMinute(newBpm);

  childModules.getLock().enter();
  for(int i=0; i<childModules.size(); i++)
    childModules[i]->setBeatsPerMinute(newBpm);
  childModules.getLock().exit();
}

void AudioModule::setModuleName(const juce::String &newName)
{
  moduleName = newName;
}

void AudioModule::setModuleNameAppendix(const juce::String& newAppendix)
{
  moduleNameAppendix = newAppendix;
}

void AudioModule::addChildAudioModule(AudioModule* moduleToAdd)
{ 
  ScopedLock scopedLock(*plugInLock);
  childModules.getLock().enter();
  childModules.addIfNotAlreadyThere(moduleToAdd);
  childModules.getLock().exit();
  addChildStateManager(moduleToAdd);
}

void AudioModule::removeChildAudioModule(AudioModule* moduleToRemove, bool deleteObject)
{ 
  ScopedLock scopedLock(*plugInLock);
  childModules.getLock().enter();
  int index = childModules.indexOf(moduleToRemove);
  jassert( index != -1 ); // trying to remove a module which is not a child of this one?
  if( index != -1 )
  { 
    childModules.remove(index);
    removeChildStateManager(moduleToRemove);
    if( deleteObject == true )
      delete moduleToRemove;
  }
  childModules.getLock().exit();
}

bool AudioModule::checkForCrack()
{
  return false;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AudioModule::getIndexAmongNameSakes(AudioModule *child)
{
  ScopedLock scopedLock(*plugInLock);

  juce::String name = child->getModuleName();
  int index = -1;

  childModules.getLock().enter();
  for(int c = 0; c < childModules.size(); c++)
  {
    if( childModules[c]->getModuleName() == name )
    {
      index++;
      if( childModules[c] == child )
        break;
    }
  }
  childModules.getLock().exit();

  return index;
}

juce::String AudioModule::getModuleHeadlineString()
{
  return moduleName + moduleNameAppendix;
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void AudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  markStateAsDirty();
}

XmlElement* AudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  if( !wantsSaveAndRecallState() )
    return NULL;

  // the XmlElement which stores all the relevant state-information:
  XmlElement* xmlState = new XmlElement(moduleName); 

  // store a version string - useful when patch-formats change later:
  //xmlState->setAttribute(T("Version"), versionjuce::String);
  xmlState->setAttribute(T("PatchFormat"), patchFormatIndex);

  // store controller mappings (if any)
  automatableModuleStateToXml(this, xmlState);

  // if this module is a polyphonic instruments, we store some global instrument parameters (such 
  // as number of voices, tuning, etc.):
  if( underlyingRosicInstrument != NULL )
    polyphonicInstrumentStateToXml(underlyingRosicInstrument, xmlState);

  // save the states of all childModules in child-XmlElements:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
  {
    if( childModules[c]->wantsSaveAndRecallState() )
    {
      XmlElement* childState = childModules[c]->getStateAsXml(childModules[c]->getStateName(),false);
      xmlState->addChildElement(childState);
    }
  }
  childModules.getLock().exit();

  setStateName(stateName, markAsClean);

  return xmlState;
}

void AudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  XmlElement convertedState = convertXmlStateIfNecessary(xmlState);
  automatableModuleStateFromXml(this, convertedState);

  // check, if this module wraps an instrument - if so, we wave some more settings to restore:
  if( underlyingRosicInstrument != NULL )
    polyphonicInstrumentStateFromXml(underlyingRosicInstrument, convertedState);

  juce::String thisName =  this->moduleName;  // for debug

  // if we have child-modules, we try to restore their states by looking for corresponding
  // child XmlElements in the xmlState:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
  {
    childModules[c]->setStateToDefaults();
    int indexAmongNameSakes = getIndexAmongNameSakes(childModules[c]);
    XmlElement* childState  = rojue::getChildElementByNameAndIndexAmongNameSakes(convertedState, childModules[c]->moduleName, 
                                                                                 indexAmongNameSakes);
    if( childState != NULL )
    {
      //childModules[c]->setStateFromXml(*childState, stateName, markAsClean);
      childModules[c]->setStateFromXml(*childState, juce::String::empty, markAsClean);
    }
  }
  childModules.getLock().exit();


  // this call we need for setting our parent dirty when some embedded sub-module loads a state 
  // - when the call was due to the outer parent loading its state, the subsequent call to 
  // setStateName may set it back to 'clean'  
  if( parent != NULL )
    parent->markStateAsDirty();

  setStateName(stateName, markAsClean);
}

XmlElement AudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*plugInLock);

  int xmlPatchFormatIndex = xmlState.getIntAttribute(T("PatchFormat"), 1);
  jassert( xmlPatchFormatIndex == this->patchFormatIndex ); // you should override the state conversion in your subclass

  return xmlState;
}
    
void AudioModule::registerDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*plugInLock);
  deletionWatchers.addIfNotAlreadyThere(watcher);
}
           
void AudioModule::deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*plugInLock);
  deletionWatchers.removeValue(watcher);
}

//-------------------------------------------------------------------------------------------------
// event processing:

void AudioModule::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*plugInLock);
  if( message.isController() )
  {
    int controllerNumber = message.getControllerNumber();
    int controllerValue  = message.getControllerValue();
    setMidiController(controllerNumber, controllerValue);
  }
  else if( message.isNoteOn() )
    noteOn(message.getNoteNumber(), message.getVelocity());
  else if( message.isNoteOff() )
    noteOff(message.getNoteNumber());
  else if( message.isAllNotesOff() )
    allNotesOff();
  else if( message.isPitchWheel() )
    setPitchBend(message.getPitchWheelValue());
}

void AudioModule::noteOn(int noteNumber, int velocity)
{
  ScopedLock scopedLock(*plugInLock);
  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->noteOn(noteNumber, velocity);
}

void AudioModule::noteOff(int noteNumber)
{
  ScopedLock scopedLock(*plugInLock);
  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->noteOff(noteNumber);
}

void AudioModule::allNotesOff()
{
  ScopedLock scopedLock(*plugInLock);
  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->allNotesOff();
}

void AudioModule::setMidiController(int controllerNumber, int controllerValue)
{
  ScopedLock scopedLock(*plugInLock);
  AutomatableModule::setMidiController(controllerNumber, controllerValue);
  if( underlyingRosicInstrument != NULL )
    underlyingRosicInstrument->setMidiController(controllerNumber, controllerValue);

  int dummy = 0;

  // distribute the controller message to all children (i.e. embedded sub-modules):
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
    childModules[c]->setMidiController(controllerNumber, controllerValue);
  childModules.getLock().exit();
}

void AudioModule::setPitchBend(int pitchBendValue)
{
  ScopedLock scopedLock(*plugInLock);
  if( underlyingRosicInstrument != NULL )
  {
    double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
    underlyingRosicInstrument->setPitchBend(wheelValueMapped);
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void AudioModule::initializeAutomatableParameters()
{

}

void AudioModule::updateCoreObjectAccordingToParameters()
{
  ScopedLock scopedLock(*plugInLock);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) observedParameters.size(); i++ )
    parameterChanged(observedParameters[i]);
}