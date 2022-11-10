AudioModuleDeletionWatcher::~AudioModuleDeletionWatcher()
{
  for(int i = 0; i < watchedAudioModules.size(); i++ )
    watchedAudioModules[i]->deRegisterDeletionWatcher(this);
}

void AudioModuleDeletionWatcher::addWatchedAudioModule(AudioModule *moduleToBeWatched)
{
  if( moduleToBeWatched != nullptr ) {
    moduleToBeWatched->registerDeletionWatcher(this);
    watchedAudioModules.addIfNotAlreadyThere(moduleToBeWatched); }
}

void AudioModuleDeletionWatcher::removeWatchedAudioModule(AudioModule
  *moduleNotToBeWatchedAnymore)
{
  if( moduleNotToBeWatchedAnymore != nullptr ) {
    moduleNotToBeWatchedAnymore->deRegisterDeletionWatcher(this);
    watchedAudioModules.removeFirstMatchingValue(moduleNotToBeWatchedAnymore); }
}

//=================================================================================================
// class AudioModule

// construction/destruction:

AudioModule::AudioModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse) 
  : ParameterManager(lockToUse), metaParamManager(metaManagerToUse)
{
  ParameterObserver::setLocalAutomationSwitch(true);  // activate automation for this instance
  moduleName = juce::String("AudioModule");

  // initialized in header already:
  //wantsTempoSyncInfo = true;
  //patchFormatIndex = 1;
  //triggerInterval = 0.0;
  //saveAndRecallState = true;
}

AudioModule::~AudioModule()
{
  ScopedPointerLock scopedLock(lock);

  for(int i = 0; i < deletionWatchers.size(); i++)
    deletionWatchers[i]->audioModuleWillBeDeleted(this);

  AudioModule* childModule = nullptr;
  while( childModules.size() > 0 )
  {
    childModule = childModules[childModules.size()-1];
    removeChildAudioModule(childModule, true);
    //childModule->parent = NULL;
    //delete childModules[childModules.size()-1];
    //childModules.removeLast();
  }
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < (int)childModules.size(); i++)
    childModules[i]->setSampleRate(newSampleRate);
}

void AudioModule::setBeatsPerMinute(double newBpm)
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < (int)childModules.size(); i++)
    childModules[i]->setBeatsPerMinute(newBpm);
}

void AudioModule::setModuleName(const juce::String &newName)
{
  moduleName = newName;
}

void AudioModule::setModuleTypeName(const juce::String& newName, bool updatePresetDirectory,
  bool setModuleNameAlso)
{
  moduleTypeName = newName;
  if(setModuleNameAlso)
    moduleName = newName;
  if(updatePresetDirectory)
    setActiveDirectory(getPresetDirectory());
}

void AudioModule::setModuleNameAppendix(const juce::String& newAppendix)
{
  moduleNameAppendix = newAppendix;
}

void AudioModule::addChildAudioModule(AudioModule* moduleToAdd)
{
  ScopedLock scopedLock(*lock);
  appendIfNotAlreadyThere(childModules, moduleToAdd);
  moduleToAdd->parentModule = this;
  moduleToAdd->setSmoothingManager(smoothingManager);
  moduleToAdd->setMetaParameterManager(metaParamManager);
  addChildStateManager(moduleToAdd);
}

void AudioModule::removeChildAudioModule(AudioModule* moduleToRemove, bool deleteObject)
{
  ScopedLock scopedLock(*lock);
  int index = find(childModules, moduleToRemove);
  jassert( index != -1 ); // trying to remove a module which is not a child of this one?
  if( index != -1 ) {
    remove(childModules, index);
    removeChildStateManager(moduleToRemove);
    moduleToRemove->parentModule = nullptr;
    moduleToRemove->setMetaParameterManager(nullptr);
    if( deleteObject == true )
      delete moduleToRemove; }
}

void AudioModule::loadPreset(const juce::String& pathFromPresetFolder)
{
  juce::String path = getPresetDirectory() + "/" + pathFromPresetFolder;
  jassert(loadFile(File(path))); // such a preset file doesn't exist
}

/*
void AudioModule::loadDefaultPreset()
{
  juce::String fileName = getDefaultPresetLocation();
  juce::File presetFile = juce::File(fileName);
  if( presetFile.existsAsFile() )
  {
    XmlDocument myDocument(presetFile);
    XmlElement* xmlState = myDocument.getDocumentElement();
    if( xmlState != nullptr )
    {
      setStateFromXml(*xmlState, presetFile.getFileNameWithoutExtension(), true);
      setActiveFile(presetFile);
      delete xmlState;
    }
  }
}
*/

bool AudioModule::checkForCrack()
{
  return false;
}

void AudioModule::addObservedParameter(Parameter *p)
{
  ParameterManager::addObservedParameter(p);
  setupManagers(p);
}

void AudioModule::setupManagers(Parameter* p)
{
  MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*> (p);
  if(mcp != nullptr)
    mcp->setMetaParameterManager(getMetaParameterManager());

  rsSmoothableParameter* sp = dynamic_cast<rsSmoothableParameter*> (p);
  if(sp != nullptr)
    sp->setSmoothingManager(smoothingManager);
}

//-------------------------------------------------------------------------------------------------
// MIDI controller stuff:

void AudioModule::assignMidiController(const String& nameOfParameter, int controllerNumber)
{
  ScopedLock scopedLock(*lock);
  Parameter* p = getParameterByName(nameOfParameter);
  if( p != nullptr ) {
    AutomatableParameter* ap = dynamic_cast<AutomatableParameter*>(p);
    if( ap != nullptr )
      ap->assignMidiController(controllerNumber); }
}

void AudioModule::setMidiController(int controllerNumber, float controllerValue)
{
  // loop through all the observed parameters and pass the controller value to them - the
  // parameters themselves will take care to respond only to controller-numbers which are assigned
  // to them:
  ScopedLock scopedLock(*lock);
  AutomatableParameter *ap;
  for(int i = 0; i < size(parameters); i++) {
    ap = dynamic_cast<AutomatableParameter*>(parameters[i]);
    if( ap != nullptr )
      ap->setMidiController(controllerNumber, controllerValue); }
}

void AudioModule::revertToDefaultMapping()
{
  ScopedLock scopedLock(*lock);
  AutomatableParameter *ap;
  for(int i = 0; i < size(parameters); i++) {
    ap = dynamic_cast<AutomatableParameter*>(parameters[i]);
    if( ap != nullptr )
      ap->revertToDefaults(false, false, false); }
}

void AudioModule::detachMetaParameters()
{
  ScopedLock scopedLock(*lock);
  MetaControlledParameter *mcp;
  for(int i = 0; i < size(parameters); i++) {
    mcp = dynamic_cast<MetaControlledParameter*>(parameters[i]);
    if( mcp != nullptr )
      mcp->detachFromMetaParameter(); }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

/*
juce::String AudioModule::getSupportDirectory() const
{
  juce::String supportDir = getUserAppDataDirectory() + File::getSeparatorString() + "RS-MET";
  return supportDir;

#ifdef _WIN32
  return getApplicationDirectory(); // preliminary - use user-domcuments folder + "/RS-MET"
#elif __APPLE__
  return "/Library/Audio/RS-MET";
#elif __linux__
  return getApplicationDirectory();
#else
  return getApplicationDirectory(); 
#endif
}
*/

juce::String AudioModule::getPresetDirectory(bool user) const
{
  juce::String presetDir = getSupportDirectory() + "/Presets/" + getModuleTypeName();
  juce::File presetDirAsFile(presetDir);
  if(!presetDirAsFile.exists())
  {
    return getApplicationDirectory();
    //presetDirAsFile.createDirectory();
  }
  return presetDir;
}

int AudioModule::getIndexAmongNameSakes(AudioModule *child)
{
  ScopedLock scopedLock(*lock);
  juce::String name = child->getModuleName();
  int index = -1;
  for(int c = 0; c < (int)childModules.size(); c++) {
    if( childModules[c]->getModuleName() == name ) {
      index++;
      if( childModules[c] == child )
        break; }}
  return index;
}

juce::String AudioModule::getModuleHeadlineString()
{
  return moduleName + moduleNameAppendix;
}

AudioModule* AudioModule::getTopLevelModule()
{
  if(isTopLevelModule())
    return this;
  else
    return parentModule->getTopLevelModule();
}

juce::String AudioModule::getAudioModulePath()
{
  // We use a period "." instead of a slash "/" as path separator because periods are allowed in 
  // xml tag names whereas slashes are not
  if(isTopLevelModule())
    return moduleName + ".";
  else
    return parentModule->getAudioModulePath() + moduleName + ".";
}

AudioModuleEditor* AudioModule::createEditor(int type)
{
  return new GenericAudioModuleEditor(this);
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void AudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  markStateAsDirty();
}

void AudioModule::parameterToXml(XmlElement* xml, Parameter* p)
{
  if(p == nullptr) return; // why do we need this?
  p->saveToXml(xml);
}

void AudioModule::parametersToXml(XmlElement* xml)
{
  for(int i = 0; i < getNumParameters(); i++) 
    parameterToXml(xml, getParameterByIndex(i));
}

void AudioModule::midiMappingToXml(XmlElement* xml)
{
  // child element will be added when there are any relevant controller-mappings to be stored:
  XmlElement* xmlMapping = new XmlElement("MidiMapping");
  for(int i = 0; i < getNumParameters(); i++) {
    AutomatableParameter* ap = dynamic_cast<AutomatableParameter*>(getParameterByIndex(i));
    if( ap != nullptr ) {
      if( ap->isInDefaultState() == false ) { // store only non-default mappings
        XmlElement* xmlParameterSetup = new XmlElement(ap->getName());
        xmlParameterSetup->setAttribute("MidiCC", ap->getAssignedMidiController());
        xmlParameterSetup->setAttribute("Min",    ap->getLowerAutomationLimit()  );
        xmlParameterSetup->setAttribute("Max",    ap->getUpperAutomationLimit()  );
        xmlMapping->addChildElement(xmlParameterSetup); }}}
  if( xmlMapping->getNumChildElements() != 0 )
    xml->addChildElement(xmlMapping);
  else
    delete xmlMapping;
}

void AudioModule::metaMappingToXml(XmlElement* xml)
{
  XmlElement* xmlMapping = new XmlElement("MetaMapping");
  for(int i = 0; i < getNumParameters(); i++) {
    MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*>(getParameterByIndex(i));
    if( mcp != nullptr ) {
      if( mcp->getMetaParameterIndex() != -1 ) { // store only meta-assigned parameters
        XmlElement* xmlParameterSetup = new XmlElement(mcp->getName());
        xmlParameterSetup->setAttribute("MetaIndex", mcp->getMetaParameterIndex());
        // todo: store the mapping function
        xmlMapping->addChildElement(xmlParameterSetup); }}}
  if( xmlMapping->getNumChildElements() != 0 )
    xml->addChildElement(xmlMapping);
  else
    delete xmlMapping;
}

void AudioModule::metaValuesToXml(XmlElement* xml)
{
  if(saveAndRecallMetas == true && metaParamManager != nullptr) {
    bool allMetasAtDefault = true;
    XmlElement* xmlValues = new XmlElement("MetaParameterValues");
    for(int i = 0; i < metaParamManager->getNumMetaParameters(); i++) {
      MetaParameter* mp = metaParamManager->getMetaParameter(i);
      double value = mp->getMetaValue();
      if(value != 0.5) {     // .5 is the default, we store only values different from that
        allMetasAtDefault = false;
        xmlValues->setAttribute("M" + String(i), value); }}
    if( !allMetasAtDefault ) // no child xml needs to be stored
      xml->addChildElement(xmlValues);
    else
      delete xmlValues;
  }
}

void AudioModule::childModulesToXml(XmlElement* xml)
{
  for(int c = 0; c < (int)childModules.size(); c++) {
    if( childModules[c]->wantsSaveAndRecallState() ) {
      XmlElement* childState = childModules[c]->getStateAsXml(
        childModules[c]->getStateName(), false);
      xml->addChildElement(childState); }}
}

XmlElement* AudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  if( !wantsSaveAndRecallState() )
    return nullptr;
  XmlElement* xml = new XmlElement(moduleName);
  xml->setAttribute("PatchFormat", patchFormatIndex); // useful for later patch-format changes
  midiMappingToXml( xml);
  metaMappingToXml( xml);
  metaValuesToXml(  xml);
  parametersToXml(  xml);
  childModulesToXml(xml);
  setStateName(stateName, markAsClean);
  return xml;
}

void AudioModule::recallParametersFromXml(const XmlElement &xml)
{
  for(int i = 0; i < size(parameters); i++) 
    parameters[i]->recallFromXml(xml);
}

void AudioModule::recallChildModulesFromXml(const XmlElement &xml, bool markAsClean)
{
  for(int c = 0; c < (int)childModules.size(); c++)
  {
    childModules[c]->setStateToDefaults();
    int indexAmongNameSakes = getIndexAmongNameSakes(childModules[c]);
    XmlElement* childState = getChildElementByNameAndIndexAmongNameSakes(
      xml, childModules[c]->moduleName, indexAmongNameSakes);
    if( childState != NULL )
    {
      //childModules[c]->setStateFromXml(*childState, stateName, markAsClean);
      childModules[c]->setStateFromXml(*childState, juce::String(), markAsClean);
    }
  }
}

void AudioModule::recallMidiMappingFromXml(const XmlElement &xml)
{
  revertToDefaultMapping(); // rename to revertToDefaultMidiMapping
  XmlElement* xmlMapping = xml.getChildByName("MidiMapping");
  if( xmlMapping == nullptr )
    return; // no mapping stored, nothing to do
  forEachXmlChildElement(*xmlMapping, xmlParamSetup) {
    Parameter* p = getParameterByName(xmlParamSetup->getTagName());
    AutomatableParameter *ap = dynamic_cast<AutomatableParameter*>(p);
    if( ap != nullptr ) {
      ap->assignMidiController(   xmlParamSetup->getIntAttribute("MidiCC", -1));
      ap->setLowerAutomationLimit(xmlParamSetup->getDoubleAttribute("Min", ap->getMinValue()));
      ap->setUpperAutomationLimit(xmlParamSetup->getDoubleAttribute("Max", ap->getMaxValue())); }}
}

void AudioModule::recallMetaMappingFromXml(const XmlElement &xml)
{
  detachMetaParameters();
  XmlElement* xmlMapping = xml.getChildByName("MetaMapping");
  if( xmlMapping == nullptr )
    return; // no mapping stored, nothing to do
  forEachXmlChildElement(*xmlMapping, xmlParamSetup) {
    Parameter* p = getParameterByName(xmlParamSetup->getTagName());
    MetaControlledParameter *mcp = dynamic_cast<MetaControlledParameter*>(p);
    if(mcp != nullptr) {
      mcp->attachToMetaParameter(xmlParamSetup->getIntAttribute("MetaIndex", -1)); }}
      // todo: retrieve mapping function/curve
}

void AudioModule::recallMetaValuesFromXml(const XmlElement &xml)
{
  if(saveAndRecallMetas == true && metaParamManager != nullptr) {
    metaParamManager->resetAllToDefaults();
    XmlElement* xmlValues = xml.getChildByName("MetaParameterValues");
    if(xmlValues == nullptr)
      return;
    for(int i = 0; i < xmlValues->getNumAttributes(); i++) {
      String tmp = xmlValues->getAttributeName(i);
      tmp = tmp.fromLastOccurrenceOf("M", false, false);
      int metaIndex = tmp.getIntValue();
      tmp = xmlValues->getAttributeValue(i);
      double metaValue = tmp.getDoubleValue();
      metaParamManager->setMetaValue(metaIndex, metaValue); }}
}

void AudioModule::setStateFromXml(const XmlElement& xml, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);

#ifdef JUCE_DEBUG
  juce::String xmlString = xml.toString(); // to look at the xml content in the debugger
#endif

  bool smoothingIsBypassed = true;
  if (smoothingManager != nullptr)
  {
    smoothingIsBypassed = smoothingManager->isSmoothingBypassed();
    smoothingManager->setBypassSmoothing(true);
  }

  XmlElement convertedState = convertXmlStateIfNecessary(xml);
  recallParametersFromXml(convertedState);
  recallMidiMappingFromXml(convertedState);
  if(metaParamManager != nullptr)
  {
    recallMetaMappingFromXml(convertedState);
    recallMetaValuesFromXml(convertedState);
    // maybe this should be done before retrieving the parameter values - because when you change
    // a parameter range in the code and then load an older preset, a parameter that is connected
    // to a meta parameter will get the wrong value.
  }

  recallChildModulesFromXml(convertedState, markAsClean);

  // this call we need for setting our parent dirty when some embedded sub-module loads a state
  // - when the call was due to the outer parent loading its state, the subsequent call to
  // setStateName may set it back to 'clean'
  if( parent != NULL )
    parent->markStateAsDirty();

  setStateName(stateName, markAsClean);

  if (smoothingManager != nullptr)
    smoothingManager->setBypassSmoothing(smoothingIsBypassed); // restore old bypass state
}

XmlElement AudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);

  int xmlPatchFormatIndex = xmlState.getIntAttribute("PatchFormat", 1);
  jassert( xmlPatchFormatIndex == this->patchFormatIndex );
  // You should override the state conversion in your subclass - the idea is that when the stored
  // patch format number in the xml differs from our patchFormatIndex member variable here, a
  // conversion may be necessary. That way, new versions of the plugin that may use a different
  // patch format can still load patches that were stored by older versions.

  return xmlState;
}

void AudioModule::resetParametersToDefaultValues()
{
  ScopedLock scopedLock(*lock);
  for(int i = 0; i < (int)parameters.size(); i++)
    parameters[i]->resetToDefaultValue(true, true);
}

void AudioModule::setMetaParameterManager(MetaParameterManager* managerToUse)
{
  ScopedLock scopedLock(*lock);
  metaParamManager = managerToUse;
  unsigned int i;

  for(i = 0; i < childModules.size(); i++)
    childModules[i]->setMetaParameterManager(metaParamManager);

  for(i = 0; i < parameters.size(); i++)
  {
    MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*>(parameters[i]);
    if(mcp != nullptr)
      mcp->setMetaParameterManager(metaParamManager);
  }
}

void AudioModule::setSmoothingManager(rsSmoothingManager* managerToUse)
{
  ScopedLock scopedLock(*lock);
  smoothingManager = managerToUse;
  unsigned int i;

  for(i = 0; i < childModules.size(); i++)
    childModules[i]->setSmoothingManager(smoothingManager);

  for(i = 0; i < parameters.size(); i++)
  {
    rsSmoothableParameter* sp = dynamic_cast<rsSmoothableParameter*>(parameters[i]);
    if(sp != nullptr)
      sp->setSmoothingManager(smoothingManager);
  }
}

void AudioModule::registerDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*lock);
  deletionWatchers.addIfNotAlreadyThere(watcher);
}

void AudioModule::deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher)
{
  ScopedLock scopedLock(*lock);
  deletionWatchers.removeFirstMatchingValue(watcher);
}

//-------------------------------------------------------------------------------------------------
// others:

void AudioModule::updateCoreObjectAccordingToParameters()
{
  ScopedLock scopedLock(*lock);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

void AudioModule::callParameterCallbacks(bool recursivelyForChildModules)
{
  //ScopedLock scopedLock(*lock); // wasn't there but should be?
  int i;
  for(i = 0; i < (int)parameters.size(); i++)
    parameters[i]->callValueChangeCallbacks(parameters[i]->getValue());
  if(recursivelyForChildModules)
    for(i = 0; i < (int)childModules.size(); i++)
      childModules[i]->callParameterCallbacks(true);
}

void AudioModule::notifyParameterObservers(bool recursivelyForChildModules)
{
  //ScopedLock scopedLock(*lock); // wasn't there but should be?
  int i;
  for(i = 0; i < (int)parameters.size(); i++)
    parameters[i]->notifyObservers();
  if(recursivelyForChildModules)
    for(i = 0; i < (int)childModules.size(); i++)
      childModules[i]->notifyParameterObservers(true);
}

//=================================================================================================
// class ModulatableAudioModule

void ModulatableAudioModule::addObservedParameter(Parameter* p)
{
  AudioModule::addObservedParameter(p);
  setupForModulationIfModulatable(p);
}

void ModulatableAudioModule::setupManagers(Parameter* p)
{
  AudioModule::setupManagers(p);
  setupForModulationIfModulatable(p);
}

void ModulatableAudioModule::addChildAudioModule(AudioModule* moduleToAdd)
{
  AudioModule::addChildAudioModule(moduleToAdd);
  ModulatableAudioModule* ma = dynamic_cast<ModulatableAudioModule*> (moduleToAdd);
  if(ma != nullptr)
    ma->setModulationManager(modManager);
}

void ModulatableAudioModule::setModulationManager(ModulationManager* managerToUse)
{
  ModulationParticipant::setModulationManager(managerToUse);
  for(int i = 0; i < size(childModules); i++) 
  {
    ModulatableAudioModule* mm = dynamic_cast<ModulatableAudioModule*>(childModules[i]);
    if(mm)
      mm->setModulationManager(managerToUse); 
  }
  for(size_t i = 0; i < parameters.size(); i++)
    setupForModulationIfModulatable(parameters[i]);
}

void ModulatableAudioModule::setupForModulationIfModulatable(Parameter* p)
{
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (p);
  if(mp != nullptr)
  {
    registerModulationTarget(mp);
    mp->setOwnerAudioModule(this);
  }
}

//=================================================================================================
// class AudioModuleWithMidiIn

//-------------------------------------------------------------------------------------------------
// event processing (must be moved to subclass):

void AudioModuleWithMidiIn::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*lock);
  rsMidiMessageDispatcher::handleMidiMessage(message);
}

void AudioModuleWithMidiIn::setMidiController(int controllerNumber, float controllerValue)
{
  //ScopedLock scopedLock(*lock); // lock is already held
  //AudioModule::setMidiController(controllerNumber, controllerValue);
  ModulatableAudioModule::setMidiController(controllerNumber, controllerValue);
  for(int c = 0; c < (int)childModules.size(); c++)
    childModules[c]->setMidiController(controllerNumber, controllerValue);
}

/*
void AudioModuleWithMidiIn::setPitchBend(int pitchBendValue)
{
  //int dummy = 0;
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //{
  //  double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
  //  underlyingRosicInstrument->setPitchBend(wheelValueMapped);
  //}
}
*/

//=================================================================================================
// class AudioModulePoly

AudioModulePoly::AudioModulePoly(CriticalSection *lockToUse, 
  MetaParameterManager* metaManagerToUse, ModulationManager* modManagerToUse) 
  : AudioModuleWithMidiIn(lockToUse, metaManagerToUse, modManagerToUse) 
{
  //setVoiceManager(voiceManagerToUse);
  // This call to setVoiceManager will not call the subclass implementation of
  // allocateVoiceResources, it will resolve to the empty baseclass version, see here:
  // https://stackoverflow.com/questions/14552412/is-it-possible-to-use-the-template-method-pattern-in-the-constructor
  // "A very good thing to keep in mind in these situations is the old chant "during construction, 
  // virtual methods aren't". Because in the constructor of the base class, the object is still
  // considered to be of that base type, and so calls to virtual functions will resolve to the
  // version implemented for that base class. This is seldom what you want, and If it is pure 
  // virtual in the base class, it is definitely not what you want."

  // Then, maybe the voiceManager should not be passed to the constructor. Instead, client code
  // should be forced to call setVoiceManager manually, such that the right version of 
  // allocateVoiceResources gets called.
}

void AudioModulePoly::setVoiceManager(rsVoiceManager* managerToUse)
{
  voiceManager = managerToUse;
  allocateVoiceResources();
  for(int i = 0; i < size(childModules); i++) {
    AudioModulePoly* pm = dynamic_cast<AudioModulePoly*>(childModules[i]);
    if(pm)
      pm->setVoiceManager(voiceManager); }
}

void AudioModulePoly::addChildAudioModule(AudioModule* moduleToAdd)
{
  AudioModuleWithMidiIn::addChildAudioModule(moduleToAdd);
  AudioModulePoly* pm = dynamic_cast<AudioModulePoly*> (moduleToAdd);
  if(pm != nullptr)
    pm->setVoiceManager(voiceManager);
}

void AudioModulePoly::setMonophonic(bool shouldBeMonophonic)
{
  monophonic = shouldBeMonophonic;
  for(int i = 0; i < getNumParameters(); i++) {
    Parameter* p = getParameterByIndex(i);            // maybe factor out this sequence of 3 calls
    ModulatableParameterPoly* mp;                     // its used in the same form in
    mp = dynamic_cast<ModulatableParameterPoly*>(p);  // noteOnForVoice
    if(mp)
      mp->setMonophonic(shouldBeMonophonic); }
  // ToDo: maybe set the child-modules monophonic, too - but maybe not? Maybe a module should be
  // able to have any mix of monophonic and polyphonic modules as child modules? maybe have a 
  // boolean parameter that lets the module optionally set up the child-modules
}

void AudioModulePoly::noteOn(int key, int vel)
{
  // When we receive a noteOn, we assume that the voiceManager has just immediately before that 
  // also has received the same noteOn. The main module that drives this module (for example, 
  // ToolChain) is supposed to first set up the voice manager and then pass the event on to the 
  // child modules. 
  if(!voiceManager) return; // maybe assert that it's not null
  int voice = voiceManager->getNewestActiveVoice();
  jassert(voice >= 0); 
  // it returns -1 only when there's no active voice but that should never be the case immediately 
  // after a note-on
  noteOnForVoice(key, vel, voice);
}

void AudioModulePoly::noteOnForVoice(int key, int vel, int voice)
{
  // For all modulatable parameters that have no modulation connections set up, we call the 
  // callback for the voice that was allocated for the note, which we assume to be the newest note
  // in the voiceManager.
  for(int i = 0; i < getNumParameters(); i++) {
    Parameter* p = getParameterByIndex(i);
    ModulatableParameterPoly* mp;
    mp = dynamic_cast<ModulatableParameterPoly*>(p);
    if(mp != nullptr && !mp->hasConnectedSources())
      mp->callCallbackForVoice(voice); }
  noteOn(key, vel, voice);
}

void AudioModulePoly::processStereoFrame(double* left, double* right)
{
  if(!voiceManager) return;  // maybe assert that it's not null



  // maybe for optimization, have a boolean flag "needsMonophonicMixdown" (member) to suppress the
  // mixing of the voices because if several polyphonic modules are chained, the later modules 
  // actually do not use the mixed signals...or do they? hmmm...well, i guess they may but most 
  // will not

  double sumL = 0.0, sumR = 0.0;  // accumulators, todo: use rsFloat64x2
  if(monophonic) {
    // In mono mode, use only output of most recently triggered voice that is still active:
    int k = voiceManager->getNewestActiveVoice();
    if(k >= 0)  
      processStereoFrameVoice(&sumL, &sumR, k); }
  else {
    // In poly mode, accumulate outputs of all active voices:
    int numActiveVoices = voiceManager->getNumActiveVoices();
    processStereoFramePoly(voicesBuffer, numActiveVoices);
    //sumL = sumR = 0.0;
    for(int i = 0; i < numActiveVoices; i += 1) {
      sumL += voicesBuffer[2*i];
      sumR += voicesBuffer[2*i+1]; }}

  *left  = sumL; // todo: use *left = thruGain * *left + outGain * sumL
  *right = sumR; // ditto
}

void AudioModulePoly::processStereoFramePoly(double *buffer, int numActiveVoices)
{
  jassert(buffer); 
  // Must be a valid pointer, length should be at least 2*numActiveVoices, more typically, it will
  // be 2*voiceManager->getNumActiveVoices()...this is a bit dangerous: we rely on the assumption 
  // that the buffer is large enough, but don't check that. Maybe provide a means to make that 
  // safer.

  for(int i = 0; i < numActiveVoices; i+=1) {
    int k = voiceManager->getActiveVoiceIndex(i);
    processStereoFrameVoice(&buffer[2*i], &buffer[2*i+1], k); }
    // we use an interleaved format for easier interfacing with rsFloat64x2
}


//=================================================================================================
// class AudioModuleEditor

AudioModuleEditor::AudioModuleEditor(AudioModule* newModuleToEdit)
{
  moduleToEdit = newModuleToEdit;
  if(moduleToEdit != nullptr)
    lock = moduleToEdit->lock;
  init();
}

AudioModuleEditor::AudioModuleEditor(CriticalSection* pluginLockToUse)
{
  lock = pluginLockToUse;
  moduleToEdit = nullptr;
  init();

  //jassertfalse;
  // we need to factor out all the initialization code from the other constructor above (that
  // takes a pointer to the AudioModule and call this init-function it from there and from here...
  // -> done -> test it
}

void AudioModuleEditor::init()
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);

  if( moduleToEdit != NULL )
  {
    moduleToEdit->checkForCrack(); // todo: move this to somewhere else - maybe triggered by a Timer-thread
    setHeadlineText(moduleToEdit->getModuleHeadlineString());
  }

  addWidget( infoField = new RTextField() );
  infoField->setNoBackgroundAndOutline(true);
  infoField->setDescription(juce::String("Description of GUI elements will appear here"));
  setDescriptionField(infoField, true);

  addWidget( webLink =
    new RHyperlinkButton("www.rs-met.com", URL("http://www.rs-met.com")) );
  webLink->setNoBackgroundAndOutline(true);
  webLink->setDescription(juce::String("Visit www.rs-met.com in the web"));
  webLink->setDescriptionField(infoField);

  stateWidgetSet = new StateLoadSaveWidgetSet();
  addWidgetSet(stateWidgetSet, true, true);
  //addChildColourSchemeComponent( stateWidgetSet ); // old
  if( moduleToEdit != NULL )
    moduleToEdit->addStateWatcher(stateWidgetSet);
  stateWidgetSet->setDescriptionField(infoField);
  stateWidgetSet->stateLabel->setText(juce::String("Preset"));
  stateWidgetSet->addChangeListener(this);

  addWidget( setupButton = new RClickButton(juce::String("Setup")) );
  setupButton->setDescription(juce::String("Opens a dialog for the general settings"));
  setupButton->setDescriptionField(infoField);
  setupButton->setClickingTogglesState(false);
  setupButton->addRButtonListener(this);

  setupDialog = NULL; // we create it only when needed the first time - i.e. 'lazy initialization'

  isTopLevelEditor      = false;
  presetSectionPosition = RIGHT_TO_HEADLINE;
  linkPosition          = RIGHT_TO_INFOLINE;
  numHueOffsets         = 0;

  loadPreferencesFromFile();
  updateWidgetsAccordingToState();

  //setWantsKeyboardFocus(true);
}

AudioModuleEditor::~AudioModuleEditor()
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);
  stateWidgetSet->removeChangeListener(this);
  if( moduleToEdit != NULL )
    moduleToEdit->removeStateWatcher(stateWidgetSet);
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModuleEditor::setModuleToEdit(AudioModule *newModuleToEdit)
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);
  if( moduleToEdit != NULL )
    moduleToEdit->removeStateWatcher(stateWidgetSet);
  moduleToEdit = newModuleToEdit;
  if( moduleToEdit != NULL )
  {
    setHeadlineText(moduleToEdit->getModuleHeadlineString());
    moduleToEdit->addStateWatcher(stateWidgetSet);
  }
}

void AudioModuleEditor::invalidateModulePointer()
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);
  moduleToEdit = NULL;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AudioModuleEditor::getPresetSectionBottom()
{
  return stateWidgetSet->getBottom();
}

rsRepaintManager* AudioModuleEditor::getRepaintManager()
{
  AudioModuleEditor* parentEditor = getParentAudioModuleEditor();
  if(parentEditor != nullptr)
    return parentEditor->getRepaintManager();

  AudioPluginEditor* pluginEditor = dynamic_cast<AudioPluginEditor*>(getParentComponent());
  if(pluginEditor != nullptr)
    return pluginEditor->getRepaintManager();

  jassertfalse; // an AudioModuleEditor should be child of another AudioModuleEditor or of an AudioPluginEditor
  return nullptr;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

// was meant for copy/paste of the preset but it doesn't work because we don't receive keypresses
//bool AudioModuleEditor::keyPressed(const KeyPress &key, Component *originatingComponent)
//{
//
//  return true; // true indicates, that the event was consumed
//}

void AudioModuleEditor::mouseDown(const MouseEvent& e)
{
#if JUCE_DEBUG
  if(e.mods.isAltDown() && e.mods.isCommandDown() && screenShotsEnabled)
  {
    Image im = createComponentSnapshot(Rectangle<int>(0, 0, getWidth(), getHeight()));
    String name = moduleToEdit->getModuleTypeName();

    // factor out into saveToPngFile(const Image& im, const String& name):
    File file = File(getApplicationDirectory() + "/" + name + ".png");
    PNGImageFormat pngFormat;
    FileOutputStream fileStream(file);
    bool success = pngFormat.writeImageToStream(im, fileStream);
    // we also need to either overwrite the file if it exists or use a variation of the filename,
    // otherwise the data just gets appended
  }
#endif
}

void AudioModuleEditor::rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged)
{
  copyColourSettingsFrom(setupDialog);
}

void AudioModuleEditor::rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave)
{
  copyColourSettingsFrom(setupDialog);
  setupDialog->setVisible(false);
  savePreferencesToFile();
}

void AudioModuleEditor::rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled)
{
  copyColourSettingsFrom(setupDialog);
  setupDialog->setVisible(false);
}

void AudioModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == setupButton )
    openPreferencesDialog();
}

void AudioModuleEditor::changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged)
{
  if( objectThatHasChanged == stateWidgetSet )
    updateWidgetsAccordingToState();
}

void AudioModuleEditor::updateWidgetsAccordingToState()
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);
  Editor::updateWidgetsAccordingToState();
  if( moduleToEdit != NULL )
    stateWidgetSet->stateFileNameLabel->setText(moduleToEdit->getStateNameWithStarIfDirty());
  updateWidgetEnablement();
}

void AudioModuleEditor::resized()
{
  Editor::resized();

  if( isTopLevelEditor )
  {
    setupButton->setVisible(true);
    setupButton->setBounds(getWidth()-44, 4, 40, 16);


    infoField->setVisible(true);
    infoField->setBounds(0, getHeight()-16, getWidth(), 16);

    //int webLinkWidth = 60;
    webLink->setVisible(true);
    //int webLinkWidth = boldFont10px.getTextPixelWidth(webLink->getButtonText(), 1);
    int webLinkWidth = BitmapFontRoundedBoldA10D0::instance.getTextPixelWidth(webLink->getButtonText(), 1);

    if( linkPosition == RIGHT_TO_HEADLINE )
      webLink->setBounds(getWidth()-webLinkWidth-6, 0, webLinkWidth, 20);
    else if( linkPosition == RIGHT_TO_INFOLINE )
      webLink->setBounds(getWidth()-webLinkWidth-6, getHeight()-16+3, webLinkWidth, 16);
    else
      webLink->setVisible(false);
  }
  else
  {
    setupButton->setBounds(getWidth(), 4, 40, 16); // shifts it out to the right
    setupButton->setVisible(false);
    infoField->setVisible(false);

    infoField->setBounds(0, getHeight(), getWidth(), 16);
    // despite being invisible, it needs well-defined bounds anyway because some subclasses
    // use them to arrange their widgets - we postion the infoField just below the actual
    // component in this case

    webLink->setVisible(false);
  }

  if( presetSectionPosition != INVISIBLE )
  {
    int w = setupButton->getX();
    stateWidgetSet->setVisible(true);
    if( Editor::headlineStyle == NO_HEADLINE )
      stateWidgetSet->setBounds(0, 4, w, 16);
    else
    {
      int offset = 0;
      if( closeButton != NULL )
        offset = 20;
      if( presetSectionPosition == RIGHT_TO_HEADLINE )
        stateWidgetSet->setBounds(getHeadlineRight(), 6, w-getHeadlineRight()-offset, 16);
      else if( presetSectionPosition == BELOW_HEADLINE )
        stateWidgetSet->setBounds(0, getHeadlineBottom()+4, w, 16);
    }
  }
  else
  {
    stateWidgetSet->setVisible(false);
    stateWidgetSet->setBounds(0, 0, 0, 0);
  }
}

void AudioModuleEditor::openPreferencesDialog()
{
  if( setupDialog == nullptr )
  {
    addChildComponent( setupDialog = new ColourSchemeSetupDialog(this, numHueOffsets) );
    setupDialog->setDescriptionField(infoField, true);
    setupDialog->addListener(this);
  }
  setupDialog->setCentreRelative(0.5f, 0.5f);
  setupDialog->toFront(true);
  setupDialog->setVisible(true);
}

void AudioModuleEditor::loadPreferencesFromFile()
{
  XmlElement *xmlPreferences = getXmlFromFile( getPreferencesFileName() );
  if( xmlPreferences == nullptr )
    return;
  XmlElement *xmlColors = xmlPreferences->getChildByName(juce::String("ColorScheme"));
  if(xmlColors == nullptr)
  {
    delete xmlPreferences; // new - needs test
    return;
  }
  setColourSchemeFromXml(*xmlColors);
  delete xmlPreferences;
}

void AudioModuleEditor::savePreferencesToFile()
{
  XmlElement *xmlPreferences = new XmlElement( getPreferencesTagName() );
  XmlElement *xmlColors      = getColourSchemeAsXml();
  xmlPreferences->addChildElement(xmlColors);
  saveXmlToFile(*xmlPreferences, File(getPreferencesFileName()), false);
  delete xmlPreferences; // will also delete xmlColors
}

juce::String AudioModuleEditor::getPreferencesTagName()
{
  //ScopedLock scopedLock(*lock);
  ScopedPointerLock spl(lock);
  if( moduleToEdit != nullptr )
    return moduleToEdit->getModuleName() + juce::String("Preferences");
  else
    return juce::String("Preferences");
}

juce::String AudioModuleEditor::getPreferencesFileName()
{
  return getSupportDirectory() + File::getSeparatorString() + getPreferencesTagName()
    + juce::String(".xml");
}

//void AudioModuleEditor::autoGenerateSliders()
//{
//  //ScopedLock scopedLock(*lock);
//  ScopedPointerLock spl(lock);
//  if( moduleToEdit == NULL )
//    return;
//
//  sliders.getLock().enter();
//
//  Parameter* p;
//  RSlider*   s;
//  for(int i=0; i < moduleToEdit->getNumParameters(); i++)
//  {
//    // retrieve data of the parameter:
//    p           = moduleToEdit->getParameterByIndex(i);
//    juce::String name = juce::String(p->getName());
//
//    s = new RSlider(name + juce::String(T("Slider")));
//    addWidget(s);
//    s->setRange(p->getLowerLimit(), p->getUpperLimit(), p->getInterval(), p->getDefaultValue());
//    s->assignParameter(p);
//    s->setSliderName(name);
//    s->setDescriptionField(infoField);
//    s->setStringConversionFunction(&valueToString);
//    sliders.addIfNotAlreadyThere(s);
//  }
//
//  sliders.getLock().exit();
//}
//
//RSlider* AudioModuleEditor::getSliderByName(const juce::String &sliderName)
//{
//  automatableSliders.getLock().enter();
//  for(int i=0; i<automatableSliders.size(); i++)
//  {
//    if( automatableSliders[i]->getSliderName() == sliderName )
//    {
//      automatableSliders.getLock().exit();
//      return automatableSliders[i];
//    }
//  }
//  automatableSliders.getLock().exit();
//  return NULL;
//}

//=================================================================================================

GenericAudioModuleEditor::GenericAudioModuleEditor(AudioModule* newModuleToEdit)
  : AudioModuleEditor(newModuleToEdit)
{
  ScopedPointerLock spl(lock);
  setPresetSectionPosition(RIGHT_TO_HEADLINE);
  createWidgets();
  setSize(360, getRequiredHeight());
}

void GenericAudioModuleEditor::resized()
{
  ScopedPointerLock spl(lock);
  AudioModuleEditor::resized();

  // preliminary - this arrangement is still ugly:
  int x  = 0;
  int y  = getPresetSectionBottom() + 4;
  int w  = getWidth();
  //int h  = getHeight();
  for(int i = 0; i < parameterWidgets.size(); i++)
  {
    parameterWidgets[i]->setBounds(x+4, y, w-8, widgetHeight);
    y += widgetHeight + widgetDistance;
  }
}

void GenericAudioModuleEditor::createWidgets()
{
  ScopedPointerLock spl(lock);
  jassert(moduleToEdit != nullptr);

  // for each of the module's parameter, create an appropriate widget and add it to this editor
  // using the inherited addWidget method (this will add the widget to our inherited widgets
  // array)

  typedef rsAutomatableButton   Btn;
  typedef rsAutomatableComboBox Cmb;
  typedef rsModulatableSlider   Sld;
  //typedef rsAutomatableSlider   Sld; // doesn't show modulation connection

  Parameter *p;
  Sld *s;
  Cmb *c;
  Btn *b;
  for(int i = 0; i < moduleToEdit->getNumParameters(); i++)
  {
    p = moduleToEdit->getParameterByIndex(i);
    juce::String name = p->getName();

    if(p->getScaling() == Parameter::BOOLEAN)
    {
      // on/off parameter - create button:
      b = new Btn(name);
      b->assignParameter(p);
      b->setButtonText(name);
      b->setDescriptionField(infoField);
      addWidget(b);
      parameterWidgets.push_back(b);
    }
    else if(p->getScaling() == Parameter::STRING)
    {
      // multiple-choice parameter - create combobox:
      c = new Cmb();
      c->assignParameter(p);
      c->setDescriptionField(infoField);
      addWidget(c);
      parameterWidgets.push_back(c);
    }
    else
    {
      // numeric parameter - create slider:
      s = new Sld();
      s->setRange(p->getMinValue(), p->getMaxValue(), p->getInterval(), p->getDefaultValue());
      s->assignParameter(p);
      s->setSliderName(name);
      s->setDescriptionField(infoField);
      s->setStringConversionFunction(&valueToString);
      addWidget(s);
      parameterWidgets.push_back(s);
    }
  }
}

//=================================================================================================

SampleBasedAudioModuleEditor::SampleBasedAudioModuleEditor(AudioModule* newModuleToEdit) 
  : AudioModuleEditor(newModuleToEdit)
{
  createWidgets();

  // initialize the current directory for sample loading:
  AudioFileManager::setActiveDirectory(getSupportDirectory() + "/Samples");
}

void SampleBasedAudioModuleEditor::rButtonClicked(RButton *b)
{
  if(     b == samplePlusButton)  AudioFileManager::loadNextFile();
  else if(b == sampleMinusButton) AudioFileManager::loadPreviousFile();
  else if(b == sampleLoadButton)  AudioFileManager::openLoadingDialog();
}

void SampleBasedAudioModuleEditor::createWidgets()
{
  addWidget( sampleFileLabel = new RTextField() );
  sampleFileLabel->setNoBackgroundAndOutline(true);
  sampleFileLabel->setJustification(Justification::centredBottom);
  sampleFileLabel->setDescription("Name of the currently loaded audio file");

  addWidget( sampleLoadButton = new RButton("Load") );
  sampleLoadButton->addRButtonListener(this);
  sampleLoadButton->setDescription("Load audio file");
  sampleLoadButton->setClickingTogglesState(false);

  addWidget( samplePlusButton = new RButton(RButton::ARROW_RIGHT) );
  samplePlusButton->addRButtonListener(this);
  samplePlusButton->setDescription("Load next audio file in current directory");
  samplePlusButton->setClickingTogglesState(false);

  addWidget( sampleMinusButton = new RButton(RButton::ARROW_LEFT) );
  sampleMinusButton->addRButtonListener(this);
  sampleMinusButton->setDescription("Load previous audio file in current directory");
  sampleMinusButton->setClickingTogglesState(false);
}