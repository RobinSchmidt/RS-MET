
//-------------------------------------------------------------------------------------------------
// construction/destruction:

SimpleSamplerAudioModule::SimpleSamplerAudioModule(CriticalSection *newPlugInLock)
: PolyphonicInstrumentAudioModule(newPlugInLock)
{
  wrappedSimpleSampler = new rosic::SimpleSampler;
  setInstrumentToWrap(wrappedSimpleSampler);

  setModuleTypeName("SimpleSampler");

  samplePlayerModule = new SamplePlayerAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].oscSection.samplePlayer1);
  samplePlayerModule->setModuleName(juce::String("SamplePlayer"));
  addChildAudioModule(samplePlayerModule);

  filterModule = new MultiModeFilterAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].filter);
  filterModule->setModuleName(juce::String("Filter"));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String("PitchEnvelope"));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String("FilterEnvelope"));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedSimpleSampler->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String("AmpEnvelope"));
  addChildAudioModule(ampEnvModule);
}

SimpleSamplerAudioModule::~SimpleSamplerAudioModule()
{
  delete wrappedSimpleSampler;
}

AudioModuleEditor* SimpleSamplerAudioModule::createEditor(int type)
{
  return new SimpleSamplerModuleEditor(lock, this);
}

//=================================================================================================


SimpleSamplerModuleEditor::SimpleSamplerModuleEditor(CriticalSection *newPlugInLock, 
  SimpleSamplerAudioModule* newSimpleSamplerAudioModule) 
  : PolyphonicInstrumentEditor(newPlugInLock, newSimpleSamplerAudioModule)
{
  // assign the pointer to the rosic::SimpleSampler object to be used as aduio engine:
  jassert(newSimpleSamplerAudioModule != NULL ); // you must pass a valid module here
  simpleSamplerAudioModule = newSimpleSamplerAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the browser:

  // create WildcardFileFilter component that will be used by the browser:
  fileFilter = new WildcardFileFilter(juce::String("*.xml"), 
    juce::String(), juce::String("XML Files"));

  /*
  // create FilePreviewComponent component that will be used by the browser:
  browser = new FileBrowserComponent(FileBrowserComponent::chooseDirectoryMode, 
  FileManager::getApplicationDirectory(), fileFilter, NULL, false, true);
  //browser->addListener(this);
  addAndMakeVisible(browser);
  */

  //---------------------------------------------------------------------------
  // create and setup the button to switch between the GUI pages:

  // create the tab-selector and assign the embedded module editors:
  moduleEditorTabButtonBar = new TabbedButtonBar(TabbedButtonBar::TabsAtTop);
  moduleEditorTabButtonBar->addTab(juce::String("Performance"),        Colours::white,  0);
  moduleEditorTabButtonBar->addTab(juce::String("Sample"),             Colours::white,  1);
  moduleEditorTabButtonBar->addTab(juce::String("Pitch Envelope"),     Colours::white,  2);
  moduleEditorTabButtonBar->addTab(juce::String("Filter"),             Colours::white,  3);
  moduleEditorTabButtonBar->addTab(juce::String("Filter Envelope"),    Colours::white,  4);
  moduleEditorTabButtonBar->addTab(juce::String("Amplitude Envelope"), Colours::white,  5);
  moduleEditorTabButtonBar->setCurrentTabIndex(0);
  moduleEditorTabButtonBar->addChangeListener(this);
  addAndMakeVisible(moduleEditorTabButtonBar);

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  samplePlayerEditor = new SamplePlayerModuleEditor(lock, simpleSamplerAudioModule->samplePlayerModule);
  addAndMakeVisible( samplePlayerEditor );
  samplePlayerEditor->addChangeListener(this);
  /*
  samplePlayerEditor->setPresetRemembererToEdit( 
  &(simpleSamplerAudioModule->wrappedSimpleSampler->voiceArray[0].oscSection.samplePlayer1) );
  */
  samplePlayerEditor->setDescriptionField(infoField, true);

  filterEditor = new MultiModeFilterModuleEditor(lock, simpleSamplerAudioModule->filterModule);
  addAndMakeVisible( filterEditor );
  //filterEditor->addChangeListener(this);  
  /*
  filterEditor->setPresetRemembererToEdit( 
  &(simpleSamplerAudioModule->wrappedSimpleSampler->voiceArray[0].filter) );
  */
  filterEditor->setDescriptionField(infoField, true);

  pitchEnvEditor = new BreakpointModulatorEditor(lock, simpleSamplerAudioModule->pitchEnvModule);
  pitchEnvEditor->setHeadlineText(juce::String("Pitch Env"));
  pitchEnvEditor->setDescription(juce::String("This is the modulation generator for the pitch"));
  addAndMakeVisible( pitchEnvEditor );
  pitchEnvEditor->addChangeListener(this);
  /*
  pitchEnvEditor->setPresetRemembererToEdit(
  &(simpleSamplerAudioModule->wrappedSimpleSampler->voiceArray[0].pitchEnv) );
  */
  pitchEnvEditor->setDescriptionField(infoField, true);

  ampEnvEditor = new BreakpointModulatorEditor(lock, simpleSamplerAudioModule->ampEnvModule);
  ampEnvEditor->setHeadlineText(juce::String("Amp Env"));
  ampEnvEditor->setDescription(juce::String("This is the modulation generator for the amplitude"));
  addAndMakeVisible( ampEnvEditor );
  ampEnvEditor->addChangeListener(this);
  /*
  ampEnvEditor->setPresetRemembererToEdit(
  &(simpleSamplerAudioModule->wrappedSimpleSampler->voiceArray[0].ampEnv) );
  */
  ampEnvEditor->setDescriptionField(infoField, true);

  filterEnvEditor = new BreakpointModulatorEditor(lock, simpleSamplerAudioModule->filterEnvModule);
  filterEnvEditor->setHeadlineText(juce::String("Filter Env"));
  filterEnvEditor->setDescription(juce::String("This is the modulation generator for the filter frequency"));
  addAndMakeVisible( filterEnvEditor );
  filterEnvEditor->addChangeListener(this);
  /*
  filterEnvEditor->setPresetRemembererToEdit(
  &(simpleSamplerAudioModule->wrappedSimpleSampler->voiceArray[0].filterEnv) );
  */
  filterEnvEditor->setDescriptionField(infoField, true);

  samplePlayerEditor->setActiveDirectory(getApplicationDirectory() + juce::String("/Samples") );

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(600, 400);
}

/*
SimpleSamplerModuleEditor::~SimpleSamplerModuleEditor()
{
deleteAllChildren();
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

void SimpleSamplerModuleEditor::changeListenerCallback(ChangeBroadcaster* objectThatHasChanged)
{
  if( objectThatHasChanged == moduleEditorTabButtonBar )
  {
    stateWidgetSet->stateLabel->setVisible(false);
    stateWidgetSet->statePlusButton->setVisible(false);
    stateWidgetSet->stateMinusButton->setVisible(false);
    stateWidgetSet->stateLoadButton->setVisible(false);
    stateWidgetSet->stateSaveButton->setVisible(false);
    stateWidgetSet->stateFileNameLabel->setVisible(false);
    levelSlider->setVisible(false);
    levelByVelSlider->setVisible(false);
    numVoicesSlider->setVisible(false);
    compSlider->setVisible(false);
    tuningLabel->setVisible(false);
    tuningPlusButton->setVisible(false);
    tuningMinusButton->setVisible(false);
    tuningLoadButton->setVisible(false);
    tuningFileNameLabel->setVisible(false);
    glideButton->setVisible(false);
    glideTimeSlider->setVisible(false);
    masterTuneSlider->setVisible(false);
    wheelRangeSlider->setVisible(false);
    samplePlayerEditor->setVisible(false);
    pitchEnvEditor->setVisible(false);
    filterEditor->setVisible(false);
    filterEnvEditor->setVisible(false);
    ampEnvEditor->setVisible(false);
    //poleZeroModelEditor->setVisible(false);

    switch( moduleEditorTabButtonBar->getCurrentTabIndex() )
    {
    case 0: // performance page
    {
      stateWidgetSet->stateLabel->setVisible(true);
      stateWidgetSet->statePlusButton->setVisible(true);
      stateWidgetSet->stateMinusButton->setVisible(true);
      stateWidgetSet->stateLoadButton->setVisible(true);
      stateWidgetSet->stateSaveButton->setVisible(true);
      stateWidgetSet->stateFileNameLabel->setVisible(true);
      levelSlider->setVisible(true);
      levelByVelSlider->setVisible(true);
      numVoicesSlider->setVisible(true);
      compSlider->setVisible(true);
      tuningLabel->setVisible(true);
      tuningPlusButton->setVisible(true);
      tuningMinusButton->setVisible(true);
      tuningLoadButton->setVisible(true);
      tuningFileNameLabel->setVisible(true);
      glideButton->setVisible(true);
      glideTimeSlider->setVisible(true);
      masterTuneSlider->setVisible(true);
      wheelRangeSlider->setVisible(true);
    }
    break;
    case 1: samplePlayerEditor->setVisible(true); break;
    case 2: pitchEnvEditor->setVisible(true);     break;
    case 3: filterEditor->setVisible(true);       break;
    case 4: filterEnvEditor->setVisible(true);    break;
    case 5: ampEnvEditor->setVisible(true);       break;
    }
  }
}

void SimpleSamplerModuleEditor::updateWidgetsAccordingToState()
{
  if( simpleSamplerAudioModule == NULL )
    return;
  if( simpleSamplerAudioModule->wrappedSimpleSampler == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = simpleSamplerAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  samplePlayerEditor->updateWidgetsAccordingToState();
  filterEditor->updateWidgetsAccordingToState();
  pitchEnvEditor->updateWidgetsAccordingToState();
  filterEnvEditor->updateWidgetsAccordingToState();
  ampEnvEditor->updateWidgetsAccordingToState();

  // pretend a click on the selected tab in order to update the visibility of the editors:
  changeListenerCallback(moduleEditorTabButtonBar);

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )
    simpleSamplerAudioModule->markStateAsDirty();
  else
    simpleSamplerAudioModule->markStateAsClean();
}

void SimpleSamplerModuleEditor::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  webLink->setBounds(w-112, 0, 112-4, 20);

  x = 0;
  y = getHeadlineBottom();
  w = getWidth()/4;
  h = getHeight()-y-20;

  //browser->setBounds(x, y-4, w, h+8);

  x = 0;
  //x = browser->getRight();
  w = getWidth()-x;
  //y = browser->getY()+4;
  y = getHeadlineBottom();

  moduleEditorTabButtonBar->setBounds(x, y, w, 28);


  w = getWidth()-x;
  y = moduleEditorTabButtonBar->getBottom();
  h = getHeight()-y-20;


  samplePlayerEditor->setBounds(x, y, w, h);
  pitchEnvEditor->setBounds(    x, y, w, h);
  filterEditor->setBounds(      x, y, w, h);
  ampEnvEditor->setBounds(      x, y, w, h);
  filterEnvEditor->setBounds(   x, y, w, h);


  stateWidgetSet->stateLabel->setBounds(x+4, y+4, 52, 20);

  stateWidgetSet->statePlusButton->setBounds(x+w/2-4-20, y+4, 20, 20);
  stateWidgetSet->stateMinusButton->setBounds(stateWidgetSet->statePlusButton->getX()-20, y+4, 20, 20);
  stateWidgetSet->stateLoadButton->setBounds(stateWidgetSet->stateMinusButton->getX()-40-4, y+4, 40, 20);
  stateWidgetSet->stateSaveButton->setBounds(stateWidgetSet->stateLoadButton->getX()-40-4, y+4, 40, 20);

  y = stateWidgetSet->stateLabel->getBottom();
  stateWidgetSet->stateFileNameLabel->setBounds(stateWidgetSet->stateLabel->getX(), y+4, w/2-8, 20);

  x += w/2;
  //y = globalButton->getBottom();
  moduleEditorTabButtonBar->getBottom();
  levelSlider->setBounds(x+4, y+4, w/4-8, 20);
  y += 24;
  levelByVelSlider->setBounds(x+4, y+4, w/4-8, 20);
  y -= 24;
  x += w/4;
  numVoicesSlider->setBounds(x+4, y+4, w/4-8, 20);
  y += 24;
  compSlider->setBounds(x+4, y+4, w/4-8, 20);

  //x = browser->getRight();
  x = 0;
  y = compSlider->getBottom()+4;

  tuningLabel->setBounds(x+4, y+4, 52, 20);

  tuningPlusButton->setBounds(x+w/2-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-20, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40-4, y+4, 40, 20);

  y = tuningLabel->getBottom();
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y+4, w/2-8, 20);


  x += w/2;
  y  = compSlider->getBottom()+4;


  glideButton->setBounds(x+4, y+4, 56, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w/2-glideButton->getWidth()-12, 20);

  //w /= 2;
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w/4-8, 20);
  x += w/4;
  wheelRangeSlider->setBounds(x+4, y+4, w/4-8, 20);

  //x = 0;
  //y = presetFileNameLabel->getBottom();
  //w = getWidth()/2;
  //h = 304;

  //infoLabel->setBounds(0, getHeight()-20, 40, 20);
  //infoField->setBounds(infoLabel->getRight(), getHeight()-20, getWidth()-infoLabel->getRight(),20);
}
