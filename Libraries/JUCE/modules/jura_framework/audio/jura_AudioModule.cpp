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

AudioModule::AudioModule(CriticalSection *lockToUse) : ParameterManager(lockToUse)
{
  ParameterObserver::localAutomationSwitch = true;  // activate automation for this instance
  moduleName = juce::String("AudioModule");

  // initialized in header already:
  //wantsTempoSyncInfo = true;
  //patchFormatIndex = 1;
  //triggerInterval = 0.0;
  //saveAndRecallState = true;
}

AudioModule::~AudioModule()
{
  ScopedLock scopedLock(*lock);

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

void AudioModule::setModuleNameAppendix(const juce::String& newAppendix)
{
  moduleNameAppendix = newAppendix;
}

void AudioModule::addChildAudioModule(AudioModule* moduleToAdd)
{
  ScopedLock scopedLock(*lock);
  appendIfNotAlreadyThere(childModules, moduleToAdd);
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
    moduleToRemove->setMetaParameterManager(nullptr);
    if( deleteObject == true )
      delete moduleToRemove; }
}

bool AudioModule::checkForCrack()
{
  return false;
}

void AudioModule::addObservedParameter(Parameter *parameterToAdd)
{
  ParameterManager::addObservedParameter(parameterToAdd);
  MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*> (parameterToAdd);
  if(mcp != nullptr)
    mcp->setMetaParameterManager(getMetaParameterManager());
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

AudioModuleEditor* AudioModule::createEditor()
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

void AudioModule::midiMappingToXml(XmlElement* xmlState)
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
    xmlState->addChildElement(xmlMapping);
  else
    delete xmlMapping;
}

void AudioModule::metaMappingToXml(XmlElement* xmlState)
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
    xmlState->addChildElement(xmlMapping);
  else
    delete xmlMapping;
}

void AudioModule::metaValuesToXml(XmlElement* xmlState)
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
      xmlState->addChildElement(xmlValues);
    else
      delete xmlValues;
  }
}

XmlElement* AudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  if( !wantsSaveAndRecallState() )
    return nullptr;

  // the XmlElement which stores all the relevant state-information:
  XmlElement* xmlState = new XmlElement(moduleName);

  // store a patch format version (useful when patch-formats change later):
  xmlState->setAttribute("PatchFormat", patchFormatIndex);

  // store midi and/or meta mappings (if any):
  midiMappingToXml(xmlState);
  metaMappingToXml(xmlState);
  metaValuesToXml(xmlState);

  // maybe store smoothing values here

  // store current parameter values:
  for(int i = 0; i < getNumParameters(); i++) {
    Parameter* p = getParameterByIndex(i);
    if( p != nullptr ) {  // do we need this?
      if( p->shouldBeSavedAndRecalled() && !p->isCurrentValueDefaultValue() ) {
        if( p->isStringParameter() )
          xmlState->setAttribute(p->getName(), p->getStringValue());
        else
          xmlState->setAttribute(p->getName(), juce::String(p->getValue()) ); }}}

  // save the states of all childModules in child-XmlElements:
  for(int c = 0; c < (int)childModules.size(); c++) {
    if( childModules[c]->wantsSaveAndRecallState() ) {
      XmlElement* childState = childModules[c]->getStateAsXml(
        childModules[c]->getStateName(), false);
      xmlState->addChildElement(childState); }}

  setStateName(stateName, markAsClean);
  return xmlState;
}

void AudioModule::midiMappingFromXml(const XmlElement &xmlState)
{
  revertToDefaultMapping(); // rename to revertToDefaultMidiMapping
  XmlElement* xmlMapping = xmlState.getChildByName("MidiMapping");
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

void AudioModule::metaMappingFromXml(const XmlElement &xmlState)
{
  detachMetaParameters();
  XmlElement* xmlMapping = xmlState.getChildByName("MetaMapping");
  if( xmlMapping == nullptr )
    return; // no mapping stored, nothing to do
  forEachXmlChildElement(*xmlMapping, xmlParamSetup) {
    Parameter* p = getParameterByName(xmlParamSetup->getTagName());
    MetaControlledParameter *mcp = dynamic_cast<MetaControlledParameter*>(p);
    if(mcp != nullptr) {
      mcp->attachToMetaParameter(xmlParamSetup->getIntAttribute("MetaIndex", -1)); }}
      // todo: retrieve mapping function/curve
}

void AudioModule::metaValuesFromXml(const XmlElement &xmlState)
{
  if(saveAndRecallMetas == true && metaParamManager != nullptr) {
    metaParamManager->resetAllToDefaults();
    XmlElement* xmlValues = xmlState.getChildByName("MetaParameterValues");
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

void AudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ScopedLock scopedLock(*lock);

  XmlElement convertedState = convertXmlStateIfNecessary(xmlState);

  // retrieve parameter values:
  juce::String name;
  Parameter* p;
  for(int i = 0; i < size(parameters); i++) {
    p = parameters[i];
    name = p->getName();
    if(p->shouldBeSavedAndRecalled()) {
      if(p->isStringParameter())
        p->setStringValue(xmlState.getStringAttribute(
          name, p->getDefaultStringValue()), true, true);
      else
        p->setValue(xmlState.getDoubleAttribute(name, p->getDefaultValue()), true, true); }}

  // retrieve midi- and/or meta mapping:
  midiMappingFromXml(convertedState);
  metaMappingFromXml(convertedState);
  metaValuesFromXml(convertedState);


  // if we have child-modules, we try to restore their states by looking for corresponding
  // child XmlElements in the xmlState:
  for(int c = 0; c < (int)childModules.size(); c++)
  {
    childModules[c]->setStateToDefaults();
    int indexAmongNameSakes = getIndexAmongNameSakes(childModules[c]);
    XmlElement* childState = getChildElementByNameAndIndexAmongNameSakes(
      convertedState, childModules[c]->moduleName, indexAmongNameSakes);
    if( childState != NULL )
    {
      //childModules[c]->setStateFromXml(*childState, stateName, markAsClean);
      childModules[c]->setStateFromXml(*childState, juce::String::empty, markAsClean);
    }
  }

  // this call we need for setting our parent dirty when some embedded sub-module loads a state
  // - when the call was due to the outer parent loading its state, the subsequent call to
  // setStateName may set it back to 'clean'
  if( parent != NULL )
    parent->markStateAsDirty();

  setStateName(stateName, markAsClean);
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
  int i;

  for(i = 0; i < (int)childModules.size(); i++)
    childModules[i]->setMetaParameterManager(metaParamManager);

  // maybe factor out - we may need to do this in other places as well:
  for(i = 0; i < (int)parameters.size(); i++)
  {
    MetaControlledParameter* mcp = dynamic_cast<MetaControlledParameter*>(parameters[i]);
    if(mcp != nullptr)
      mcp->setMetaParameterManager(metaParamManager);
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
    parameters[i]->callValueChangeCallbacks();
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
// class AudioModuleWithMidiIn

//-------------------------------------------------------------------------------------------------
// event processing (must be moved to subclass):

void AudioModuleWithMidiIn::handleMidiMessage(MidiMessage message)
{
  ScopedLock scopedLock(*lock);
  if( message.isController() )
  {
    int controllerNumber = message.getControllerNumber();
    int controllerValue  = message.getControllerValue();
    setMidiController(controllerNumber, (float) controllerValue);
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

void AudioModuleWithMidiIn::noteOn(int noteNumber, int velocity)
{

}

void AudioModuleWithMidiIn::noteOff(int noteNumber)
{

}

void AudioModuleWithMidiIn::allNotesOff()
{

}

void AudioModuleWithMidiIn::setMidiController(int controllerNumber, float controllerValue)
{
  ScopedLock scopedLock(*lock);
  AudioModule::setMidiController(controllerNumber, controllerValue);
  for(int c = 0; c < (int)childModules.size(); c++)
    childModules[c]->setMidiController(controllerNumber, controllerValue);
}

void AudioModuleWithMidiIn::setPitchBend(int pitchBendValue)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //{
  //  double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
  //  underlyingRosicInstrument->setPitchBend(wheelValueMapped);
  //}
}

//=================================================================================================
// class AudioModuleEditor

AudioModuleEditor::AudioModuleEditor(AudioModule* newModuleToEdit)
{
  lock = newModuleToEdit->lock;
  moduleToEdit = newModuleToEdit;
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
  ScopedLock scopedLock(*lock);

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
  addChildColourSchemeComponent( stateWidgetSet );
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
}

AudioModuleEditor::~AudioModuleEditor()
{
  ScopedLock scopedLock(*lock);
  stateWidgetSet->removeChangeListener(this);
  if( moduleToEdit != NULL )
    moduleToEdit->removeStateWatcher(stateWidgetSet);
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModuleEditor::setModuleToEdit(AudioModule *newModuleToEdit)
{
  ScopedLock scopedLock(*lock);
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
  ScopedLock scopedLock(*lock);
  moduleToEdit = NULL;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AudioModuleEditor::getPresetSectionBottom()
{
  return stateWidgetSet->getBottom();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

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
  ScopedLock scopedLock(*lock);
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
    int w = setupButton->getX();;
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
  if( setupDialog == NULL )
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
  if( xmlPreferences == NULL )
    return;
  XmlElement *xmlColors = xmlPreferences->getChildByName(juce::String("ColorScheme"));
  if( xmlColors == NULL )
    return;
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
  ScopedLock scopedLock(*lock);
  if( moduleToEdit != NULL )
    return moduleToEdit->getModuleName() + juce::String("Preferences");
  else
    return juce::String("Preferences");
}

juce::String AudioModuleEditor::getPreferencesFileName()
{
  return getApplicationDirectory() + File::separatorString + getPreferencesTagName()
    + juce::String(".xml");
}

//void AudioModuleEditor::autoGenerateSliders()
//{
//  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*lock);

  setPresetSectionPosition(RIGHT_TO_HEADLINE);
  createWidgets();

  int height = (widgetHeight+widgetDistance) * parameterWidgets.size() + 28;
  setSize(360, height);
}

void GenericAudioModuleEditor::resized()
{
  ScopedLock scopedLock(*lock);
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
  ScopedLock scopedLock(*lock);
  jassert(moduleToEdit != nullptr);

  // for each of the module's parameter, create an appropriate widget and add it to this editor
  // using the inherited addWidget method (this will add the widget to out inherited widgets
  // array)

  Parameter *p;
  RSlider   *s;
  RComboBox *c;
  //RButton   *b;
  for(int i = 0; i < moduleToEdit->getNumParameters(); i++)
  {
    p = moduleToEdit->getParameterByIndex(i);
    juce::String name = juce::String(p->getName());

    if(p->getScaling() == Parameter::BOOLEAN)
    {
      // on/off parameter - create button:
      jassertfalse; // not yet implemented
    }
    else if(p->getScaling() == Parameter::STRING)
    {
      // multiple-choice parameter - create combobox:
      c = new RComboBox(name + "Menu");
      c->assignParameter(p);
      c->setDescriptionField(infoField);
      addWidget(c);
      parameterWidgets.add(c);
    }
    else
    {
      // numeric parameter - create slider:
      s = new RSlider(name + "Slider");
      s->setRange(p->getMinValue(), p->getMaxValue(), p->getInterval(), p->getDefaultValue());
      s->assignParameter(p);
      s->setSliderName(name);
      s->setDescriptionField(infoField);
      s->setStringConversionFunction(&valueToString);
      addWidget(s);
      parameterWidgets.add(s);
    }
  }
}
