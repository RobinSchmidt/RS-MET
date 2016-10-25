
AudioModuleEditor::AudioModuleEditor(AudioModule* newModuleToEdit)
{
  plugInLock = newModuleToEdit->plugInLock;

  ScopedLock scopedLock(*plugInLock);

  moduleToEdit = newModuleToEdit;

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

/*
void AudioModuleEditor::autoGenerateSliders()
{
  ScopedLock scopedLock(*plugInLock);
  if( moduleToEdit == NULL )
    return;

  sliders.getLock().enter();

  Parameter* p;
  RSlider*   s;
  for(int i=0; i < moduleToEdit->getNumParameters(); i++)
  {
    // retrieve data of the parameter:
    p           = moduleToEdit->getParameterByIndex(i);
    juce::String name = juce::String(p->getName());

    s = new RSlider(name + juce::String(T("Slider")));
    addWidget(s);
    s->setRange(p->getLowerLimit(), p->getUpperLimit(), p->getInterval(), p->getDefaultValue());
    s->assignParameter(p);
    s->setSliderName(name);
    s->setDescriptionField(infoField);
    s->setStringConversionFunction(&valueToString);
    sliders.addIfNotAlreadyThere(s);
  }

  sliders.getLock().exit();
}
*/

/*
RSlider* AudioModuleEditor::getSliderByName(const juce::String &sliderName)
{
  automatableSliders.getLock().enter();
  for(int i=0; i<automatableSliders.size(); i++)
  {
    if( automatableSliders[i]->getSliderName() == sliderName )
    {
      automatableSliders.getLock().exit();
      return automatableSliders[i];
    }
  }
  automatableSliders.getLock().exit();
  return NULL;
}
*/

