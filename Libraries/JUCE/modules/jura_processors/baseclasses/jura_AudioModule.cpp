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

void AudioModuleDeletionWatcher::removeWatchedAudioModule(AudioModule 
  *moduleNotToBeWatchedAnymore)
{
  if( moduleNotToBeWatchedAnymore != NULL )
  {
    moduleNotToBeWatchedAnymore->deRegisterDeletionWatcher(this);
    watchedAudioModules.removeFirstMatchingValue(moduleNotToBeWatchedAnymore);
  }
}

//=================================================================================================
// class AudioModule

// construction/destruction:

AudioModule::AudioModule(CriticalSection *lockToUse)
{
  //plugInLock = new CriticalSection;

  plugInLock = lockToUse;

  ParameterObserver::localAutomationSwitch = true;  // activate automation for this instance
  wantsTempoSyncInfo = true;
  //underlyingRosicInstrument = NULL;  // by default, this does not wrap an instrument
  moduleName = juce::String("AudioModule");
  //versionjuce::String = juce::String(T("1.0"));
  patchFormatIndex = 1;
  triggerInterval = 0.0;
  saveAndRecallState = true;
  initializeAutomatableParameters();
}

AudioModule::~AudioModule()
{
  ScopedLock scopedLock(*plugInLock);
  //plugInLock->enter();

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

  //plugInLock->exit();
  //delete plugInLock;
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);

  childModules.getLock().enter();
  for(int i=0; i<childModules.size(); i++)
    childModules[i]->setSampleRate(newSampleRate);
  childModules.getLock().exit();
}

void AudioModule::setBeatsPerMinute(double newBpm)
{
  ScopedLock scopedLock(*plugInLock);

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

  // store a patch format version (useful when patch-formats change later):
  xmlState->setAttribute("PatchFormat", patchFormatIndex);

  // store controller mappings (if any)
  automatableModuleStateToXml(this, xmlState);

  //// if this module is a polyphonic instrument, we store some global instrument parameters 
  //// (such as number of voices, tuning, etc.):
  //if( underlyingRosicInstrument != NULL )
  //  polyphonicInstrumentStateToXml(underlyingRosicInstrument, xmlState);

  // save the states of all childModules in child-XmlElements:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
  {
    if( childModules[c]->wantsSaveAndRecallState() )
    {
      XmlElement* childState = childModules[c]->getStateAsXml(
        childModules[c]->getStateName(), false);
      xmlState->addChildElement(childState);
    }
  }
  childModules.getLock().exit();
  setStateName(stateName, markAsClean);
  return xmlState;
}

void AudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
  bool markAsClean)
{
  ScopedLock scopedLock(*plugInLock);

  XmlElement convertedState = convertXmlStateIfNecessary(xmlState);
  automatableModuleStateFromXml(this, convertedState);

  //// check, if this module wraps an instrument - if so, we wave some more settings to restore:
  //if( underlyingRosicInstrument != NULL )
  //  polyphonicInstrumentStateFromXml(underlyingRosicInstrument, convertedState);

  juce::String thisName =  this->moduleName;  // for debug

  // if we have child-modules, we try to restore their states by looking for corresponding
  // child XmlElements in the xmlState:
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
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

  int xmlPatchFormatIndex = xmlState.getIntAttribute("PatchFormat", 1);
  jassert( xmlPatchFormatIndex == this->patchFormatIndex ); 
  // You should override the state conversion in your subclass - the idea is that when the stored
  // patch format number in the xml differs from our patchFormatIndex member variable here, a 
  // conversion may be necessary. That way, new versions of the plugin that may use a different
  // patch format can still load patches that were stored by older versions.

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
  deletionWatchers.removeFirstMatchingValue(watcher);
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

//=================================================================================================
// class AudioModuleWithMidiIn

//-------------------------------------------------------------------------------------------------
// event processing (must be moved to subclass):

void AudioModuleWithMidiIn::handleMidiMessage(MidiMessage message)
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

void AudioModuleWithMidiIn::noteOn(int noteNumber, int velocity)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->noteOn(noteNumber, velocity);
}

void AudioModuleWithMidiIn::noteOff(int noteNumber)
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->noteOff(noteNumber);
}

void AudioModuleWithMidiIn::allNotesOff()
{
  //ScopedLock scopedLock(*plugInLock);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->allNotesOff();
}

void AudioModuleWithMidiIn::setMidiController(int controllerNumber, int controllerValue)
{
  ScopedLock scopedLock(*plugInLock);
  AutomatableModule::setMidiController(controllerNumber, controllerValue);
  //if( underlyingRosicInstrument != NULL )
  //  underlyingRosicInstrument->setMidiController(controllerNumber, controllerValue);

  // distribute the controller message to all children (i.e. embedded sub-modules):
  childModules.getLock().enter();
  for(int c=0; c<childModules.size(); c++)
    childModules[c]->setMidiController(controllerNumber, controllerValue);
  childModules.getLock().exit();
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
  plugInLock   = newModuleToEdit->plugInLock;
  moduleToEdit = newModuleToEdit;
  init();
}

AudioModuleEditor::AudioModuleEditor(CriticalSection* pluginLockToUse)
{
  plugInLock   = pluginLockToUse;
  moduleToEdit = nullptr;
  init();

  //jassertfalse;
  // we need to factor out all the initialization code from the other constructor above (that 
  // takes a pointer to the AudioModule and call this init-function it from there and from here...
  // -> done -> test it
}

void AudioModuleEditor::init()
{
  ScopedLock scopedLock(*plugInLock);

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

  setSize(400, 300); // do we need this?
}

AudioModuleEditor::~AudioModuleEditor()
{
  ScopedLock scopedLock(*plugInLock);
  stateWidgetSet->removeChangeListener(this);
  if( moduleToEdit != NULL )
    moduleToEdit->removeStateWatcher(stateWidgetSet);
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioModuleEditor::setModuleToEdit(AudioModule *newModuleToEdit)
{
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
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
  ScopedLock scopedLock(*plugInLock);
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