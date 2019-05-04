
KeyShotAudioModule::KeyShotAudioModule(CriticalSection *newPlugInLock)
  : PolyphonicInstrumentAudioModule(newPlugInLock)
{
  wrappedKeyShot = new rosic::KeyShot;
  setInstrumentToWrap(wrappedKeyShot);
  setModuleTypeName("KeyShot");

  samplePlayerModule = new SamplePlayerAudioModule(lock, 
    &wrappedKeyShot->voiceArray[0].samplePlayer);
  samplePlayerModule->setModuleName(juce::String("SamplePlayer"));
  addChildAudioModule(samplePlayerModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedKeyShot->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}

KeyShotAudioModule::~KeyShotAudioModule()
{
  delete wrappedKeyShot;
}

AudioModuleEditor* KeyShotAudioModule::createEditor(int type)
{
  return new KeyShotModuleEditor(lock, this);
}

//=================================================================================================

// construction/destruction:

KeyShotModuleEditor::KeyShotModuleEditor(CriticalSection *newPlugInLock, 
  KeyShotAudioModule* newKeyShotAudioModule) 
  : PolyphonicInstrumentEditor(newPlugInLock, newKeyShotAudioModule)
{
  // set the plugIn-headline:
  //setHeadlineText( juce::String(T("Simple Sampler")) );
  setHeadlineStyle(AudioModuleEditor::MAIN_HEADLINE);
  setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);

  // assign the pointer to the rosic::KeyShot object to be used as aduio engine:
  jassert(newKeyShotAudioModule != NULL ); // you must pass a valid module here
  keyShotAudioModule = newKeyShotAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the browser:

  // create WildcardFileFilter component that will be used by the browser:
  fileFilter = new WildcardFileFilter(juce::String("*.xml"), juce::String(), 
    juce::String("XML Files"));

  /*
  // create FilePreviewComponent component that will be used by the browser:
  browser = new FileBrowserComponent(FileBrowserComponent::chooseDirectoryMode, 
  FileManager::getApplicationDirectory(), fileFilter, NULL, false, true);
  //browser->addListener(this);
  addAndMakeVisible(browser);
  */

  /*
  //---------------------------------------------------------------------------
  // create and setup the button to switch between the GUI pages:

  // create the tab-selector and assign the embedded module editors:
  moduleEditorTabButtonBar = new TabbedButtonBar(TabbedButtonBar::TabsAtTop);
  moduleEditorTabButtonBar->addTab(juce::String(T("Performance")),        Colours::white,  0);
  moduleEditorTabButtonBar->addTab(juce::String(T("Sample")),             Colours::white,  1);
  moduleEditorTabButtonBar->addTab(juce::String(T("Pitch Envelope")),     Colours::white,  2);
  moduleEditorTabButtonBar->addTab(juce::String(T("Filter")),             Colours::white,  3);
  moduleEditorTabButtonBar->addTab(juce::String(T("Filter Envelope")),    Colours::white,  4);
  moduleEditorTabButtonBar->addTab(juce::String(T("Amplitude Envelope")), Colours::white,  5);
  moduleEditorTabButtonBar->setCurrentTabIndex(0);
  moduleEditorTabButtonBar->addChangeListener(this);
  addAndMakeVisible(moduleEditorTabButtonBar);
  */

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  /*
  samplePlayerEditor = new SamplePlayerModuleEditor(
  keyShotAudioModule->samplePlayerModule);
  addAndMakeVisible( samplePlayerEditor );
  samplePlayerEditor->addChangeListener(this);
  samplePlayerEditor->setDescriptionField(infoField, true);
  samplePlayerEditor->setActiveDirectory(getApplicationDirectory() + juce::String(T("/Samples")) );
  */

  addAndMakeVisible( waveformDisplay = new rsWaveformPlot() );  
  // todo: addPlot

  //ampEnvEditor = new BreakpointModulatorEditor(keyShotAudioModule->ampEnvModule);
  //ampEnvEditor->setHeadlineText(juce::String(T("Amp Env")));
  //ampEnvEditor->setDescription(juce::String(T("This is the modulation generator for the amplitude")));
  //addAndMakeVisible( ampEnvEditor );
  //ampEnvEditor->addChangeListener(this);
  //ampEnvEditor->setDescriptionField(infoField, true);
  

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(600, 400);
}

/*
KeyShotModuleEditor::~KeyShotModuleEditor()
{
deleteAllChildren();
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

void KeyShotModuleEditor::changeListenerCallback(ChangeBroadcaster* objectThatHasChanged)
{
  PolyphonicInstrumentEditor::changeListenerCallback(objectThatHasChanged);
}

void KeyShotModuleEditor::updateWidgetsAccordingToState()
{
  if( keyShotAudioModule == NULL )
    return;
  if( keyShotAudioModule->wrappedKeyShot == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = keyShotAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  //samplePlayerEditor->updateWidgetsAccordingToState();
  //ampEnvEditor->updateWidgetsAccordingToState();

  // update waveform display....

  // pretend a click on the selected tab in order to update the visibility of the editors:
  //changeListenerCallback(moduleEditorTabButtonBar);

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )
    keyShotAudioModule->markStateAsDirty();
  else
    keyShotAudioModule->markStateAsClean();
}

void KeyShotModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  x = getWidth()/2;
  y = getPresetSectionBottom();
  w = getWidth()-x;
  //h = (infoField->getY()-y)/2;

  //samplePlayerEditor->setBounds(x, y, w, h);
  //ampEnvEditor->setBounds(      x, y, w, h);

  //y = getHeadlineBottom();
  //w = getWidth()-x;
  h = (getHeight()-y-20)/3;

  waveformDisplay->setBounds(x, y, w, h);
  y = waveformDisplay->getBottom();
  h = infoField->getY()-y;
  //ampEnvEditor->setBounds(   x, y, w, h);


  //stateWidgetSet->stateLabel->setBounds(x+4, y+4, 52, 20);
  //stateWidgetSet->statePlusButton->setBounds(x+w/2-4-20, y+4, 20, 20);
  //stateWidgetSet->stateMinusButton->setBounds(stateWidgetSet->statePlusButton->getX()-20, y+4, 20, 20);
  //stateWidgetSet->stateLoadButton->setBounds(stateWidgetSet->stateMinusButton->getX()-40-4, y+4, 40, 20);
  //stateWidgetSet->stateSaveButton->setBounds(stateWidgetSet->stateLoadButton->getX()-40-4, y+4, 40, 20);

  //y = stateWidgetSet->stateLabel->getBottom();
  //stateWidgetSet->stateFileNameLabel->setBounds(stateWidgetSet->stateLabel->getX(), y+4, w/2-8, 20);

  //x += w/2;
  //levelSlider->setBounds(x+4, y+4, w/4-8, 20);
  //y += 24;
  //levelByVelSlider->setBounds(x+4, y+4, w/4-8, 20);
  //y -= 24;
  //x += w/4;
  numVoicesSlider->setBounds(x+4, y+4, w/4-8, 20);
}
