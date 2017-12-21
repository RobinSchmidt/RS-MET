
QuadrigaAudioModule::QuadrigaAudioModule(CriticalSection *newPlugInLock)
  : PolyphonicInstrumentAudioModule(newPlugInLock)
{
  //jassert(quadrigaToWrap != NULL); // you must pass a valid rosic-object to the constructor

  //wrappedQuadriga = quadrigaToWrap;
  //underlyingRosicInstrument = quadrigaToWrap;

  wrappedQuadriga = new rosic::Quadriga;
  setInstrumentToWrap(wrappedQuadriga);
  setModuleTypeName("Quadriga");

  quadrigenModule = new QuadrigenAudioModule(lock, &wrappedQuadriga->voiceArray[0].quadrigen);
  quadrigenModule->setModuleName(juce::String("Quadrigen"));
  addChildAudioModule(quadrigenModule);

  quadrifexModule = new QuadrifexAudioModule(lock, nullptr, nullptr, &wrappedQuadriga->quadrifex); // nullptrs preliminary
  quadrifexModule->setModuleName(juce::String("Quadrifex"));
  addChildAudioModule(quadrifexModule);

  equalizerModule = new EqualizerAudioModule(lock, &wrappedQuadriga->equalizer);
  equalizerModule->setModuleName(juce::String("Equalizer"));
  addChildAudioModule(equalizerModule);
}

QuadrigaAudioModule::~QuadrigaAudioModule()
{
  delete wrappedQuadriga;
}

AudioModuleEditor* QuadrigaAudioModule::createEditor()
{
  return new QuadrigaModuleEditor(lock, this);
}

//=================================================================================================


QuadrigaModuleEditor::QuadrigaModuleEditor(CriticalSection *newPlugInLock, 
  QuadrigaAudioModule* newQuadrigaAudioModule) 
  : PolyphonicInstrumentEditor(newPlugInLock, newQuadrigaAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);
  setHeadlineText( juce::String(("Quadriga")) );

  setHeadlinePosition(AudioModuleEditor::TOP_LEFT);
  setPresetSectionPosition(AudioModuleEditor::BELOW_HEADLINE);
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

  // assign the pointer to the rosic::Quadriga object to be used as aduio engine:
  jassert(newQuadrigaAudioModule != NULL ); // you must pass a valid module here
  quadrigaAudioModule = newQuadrigaAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  // \todo: make radio-button work again

  quadrigenEditor = new QuadrigenModuleEditor(lock, quadrigaAudioModule->quadrigenModule);
  quadrigenEditor->setDescription(juce::String(("Settings for the sources/generators section")));
  addChildEditor(quadrigenEditor, true, false);

  quadrifexEditor = new QuadrifexModuleEditor(lock, quadrigaAudioModule->quadrifexModule);
  quadrifexEditor->setDescription(juce::String(("Settings for the effects section")));
  addChildEditor(quadrifexEditor, true, false);

  addWidget( performanceButton = new RButton(juce::String(("Performance"))) );
  performanceButton->setDescription(juce::String(("Switch to the performance page")));
  performanceButton->setDescriptionField(infoField);
  //performanceButton->setRadioGroupId(1);
  performanceButton->addRButtonListener(this);

  addWidget( generatorsButton = new RButton(juce::String(("Generators"))) );
  generatorsButton->setDescription(juce::String(("Switch to the sound-generators page")));
  //generatorsButton->setRadioGroupId(1);
  generatorsButton->addRButtonListener(this);

  addWidget( filtersButton = new RButton(juce::String(("Filters"))) );
  filtersButton->setDescription(juce::String(("Switch to the filters page")));
  //filtersButton->setRadioGroupId(1);
  filtersButton->addRButtonListener(this);

  addWidget( modulatorsButton = new RButton(juce::String(("Modulators"))) );
  modulatorsButton->setDescription(juce::String(("Switch to the modulators page")));
  //modulatorsButton->setRadioGroupId(1);
  modulatorsButton->addRButtonListener(this);

  addWidget( effectsButton = new RButton(juce::String(("Effects"))) );
  effectsButton->setDescription(juce::String(("Switch to the effects page")));
  //effectsButton->setRadioGroupId(1);
  effectsButton->addRButtonListener(this);

  // set up the widgets:
  updateWidgetsAccordingToState();

  generatorsButton->setToggleState(true, true);  // \todo: save the currently visible page int the state and recall it on construction

  setSize(800, 600);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void QuadrigaModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  PolyphonicInstrumentEditor::rButtonClicked(buttonThatWasClicked);
  updateSubEditorVisibility();
}

void QuadrigaModuleEditor::changeListenerCallback(ChangeBroadcaster* objectThatHasChanged)
{

}

void QuadrigaModuleEditor::updateWidgetsAccordingToState()
{
  if( quadrigaAudioModule == NULL )
    return;
  if( quadrigaAudioModule->wrappedQuadriga == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = quadrigaAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  quadrigenEditor->updateWidgetsAccordingToState();
  quadrifexEditor->updateWidgetsAccordingToState();

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )  
    quadrigaAudioModule->markStateAsDirty();
  else
    quadrigaAudioModule->markStateAsClean();
}

void QuadrigaModuleEditor::updateSubEditorVisibility()
{
  quadrigenEditor->setVisible(false);
  quadrifexEditor->setVisible(false);
  if( generatorsButton->getToggleState() == true )
    quadrigenEditor->setVisible(true);
  else if( effectsButton->getToggleState() == true )
    quadrifexEditor->setVisible(true);
}

void QuadrigaModuleEditor::resized()
{
  //PolyphonicInstrumentEditor::resized();
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  w = getWidth()/4;
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+8, w-8, 40);

  x = stateWidgetSet->getRight();
  y = stateWidgetSet->getBottom() - 18;
  w = 100;
  performanceButton->setBounds(x+4, y, w-4, 20);
  x += w+8;
  generatorsButton->setBounds(x+4, y, w-4, 20);
  x += w;
  filtersButton->setBounds(x+4, y, w-4, 20);
  x += w;
  modulatorsButton->setBounds(x+4, y, w-4, 20);
  x += w+8;
  effectsButton->setBounds(x+4, y, w-4, 20);

  x = 0;
  y = stateWidgetSet->getBottom();
  w = getWidth();
  h = getHeight()-y-20;
  quadrigenEditor->setBounds(x, y, w, h);
  quadrifexEditor->setBounds(quadrigenEditor->getBounds());
}


