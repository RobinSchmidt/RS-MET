
WorkhorseAudioModule::WorkhorseAudioModule(CriticalSection *newPlugInLock, rosic::Workhorse *workhorseToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, workhorseToWrap)
{
  jassert(workhorseToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedWorkhorse        = workhorseToWrap;
  underlyingRosicInstrument = workhorseToWrap;
  moduleName                = juce::String(("Workhorse"));

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String(("/WorkhorsePresets")) );

  vectorMixerModule = new VectorMixerAudioModule(lock, &wrappedWorkhorse->vectorMixer);
  vectorMixerModule->setModuleName(juce::String(("VectorMixer")));
  addChildAudioModule(vectorMixerModule);

  samplePlayerTopLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopLeft);
  samplePlayerTopLeftModule->setModuleName(juce::String(("SamplePlayerTopLeft")));
  addChildAudioModule(samplePlayerTopLeftModule);

  samplePlayerTopRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopRight);
  samplePlayerTopRightModule->setModuleName(juce::String(("SamplePlayerTopRight")));
  addChildAudioModule(samplePlayerTopRightModule);

  samplePlayerBottomLeftModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomLeft);
  samplePlayerBottomLeftModule->setModuleName(juce::String(("SamplePlayerBottomLeft")));
  addChildAudioModule(samplePlayerBottomLeftModule);

  samplePlayerBottomRightModule = new SamplePlayerAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomRight);
  samplePlayerBottomRightModule->setModuleName(juce::String(("SamplePlayerBottomRight")));
  addChildAudioModule(samplePlayerBottomRightModule);

  filterModule = new MultiModeFilterAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].filter);
  filterModule->setModuleName(juce::String(("Filter")));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String(("PitchEnvelope")));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String(("FilterEnvelope")));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(lock, 
    &wrappedWorkhorse->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String(("AmpEnvelope")));
  addChildAudioModule(ampEnvModule);
}

//=================================================================================================

// construction/destruction:

WorkhorseModuleEditor::WorkhorseModuleEditor(CriticalSection *newPlugInLock, 
  WorkhorseAudioModule* newWorkhorseAudioModule) 
  : PolyphonicInstrumentEditor(newPlugInLock, newWorkhorseAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);
  setHeadlineText( juce::String(("Workhorse")) );

  // assign the pointer to the rosic::Workhorse object to be used as aduio engine:
  jassert(newWorkhorseAudioModule != NULL ); // you must pass a valid module here
  workhorseAudioModule = newWorkhorseAudioModule;

  //---------------------------------------------------------------------------
  // create and setup the sub-module editors:

  samplePlayerTopLeftEditor = new SamplePlayerModuleEditor(lock, 
    workhorseAudioModule->samplePlayerTopLeftModule);
  addAndMakeVisible( samplePlayerTopLeftEditor );
  samplePlayerTopLeftEditor->addChangeListener(this);
  samplePlayerTopLeftEditor->setDescriptionField(infoField, true);

  samplePlayerTopRightEditor = new SamplePlayerModuleEditor(lock, 
    workhorseAudioModule->samplePlayerTopRightModule);
  addChildComponent( samplePlayerTopRightEditor );
  samplePlayerTopRightEditor->addChangeListener(this);
  samplePlayerTopRightEditor->setDescriptionField(infoField, true);

  samplePlayerBottomLeftEditor = new SamplePlayerModuleEditor(lock, 
    workhorseAudioModule->samplePlayerBottomLeftModule);
  addChildComponent( samplePlayerBottomLeftEditor );
  samplePlayerBottomLeftEditor->addChangeListener(this);
  samplePlayerBottomLeftEditor->setDescriptionField(infoField, true);

  samplePlayerBottomRightEditor = new SamplePlayerModuleEditor(lock, 
    workhorseAudioModule->samplePlayerBottomRightModule);
  addChildComponent( samplePlayerBottomRightEditor );
  samplePlayerBottomRightEditor->addChangeListener(this);
  samplePlayerBottomRightEditor->setDescriptionField(infoField, true);

  filterEditor = new MultiModeFilterModuleEditor(lock, 
    workhorseAudioModule->filterModule);
  addChildComponent( filterEditor );
  filterEditor->setDescriptionField(infoField, true);

  pitchEnvEditor = new BreakpointModulatorEditor(lock, 
    workhorseAudioModule->pitchEnvModule);
  pitchEnvEditor->setHeadlineText(juce::String(("Pitch Env")));
  pitchEnvEditor->setDescription(juce::String(("This is the modulation generator for the pitch")));
  addChildComponent( pitchEnvEditor );
  pitchEnvEditor->addChangeListener(this);
  pitchEnvEditor->setDescriptionField(infoField, true);

  ampEnvEditor = new BreakpointModulatorEditor(lock, 
    workhorseAudioModule->ampEnvModule);
  ampEnvEditor->setHeadlineText(juce::String(("Amp Env")));
  ampEnvEditor->setDescription(juce::String(("This is the modulation generator for the amplitude")));
  addChildComponent( ampEnvEditor );
  ampEnvEditor->addChangeListener(this);
  ampEnvEditor->setDescriptionField(infoField, true);

  filterEnvEditor = new BreakpointModulatorEditor(lock, 
    workhorseAudioModule->filterEnvModule);
  filterEnvEditor->setHeadlineText(juce::String(("Filter Env")));
  filterEnvEditor->setDescription(juce::String(("This is the modulation generator for the filter frequency")));
  addChildComponent( filterEnvEditor );
  filterEnvEditor->addChangeListener(this);
  filterEnvEditor->setDescriptionField(infoField, true);

  addAndMakeVisible( vectorMixerPad = new VectorMixerModuleEditor(lock, 
    workhorseAudioModule->vectorMixerModule) );

  //---------------------------------------------------------------------------
  // create and setup the buttons to switch between the GUI pages:

  // \todo: make the radio-group functionality work again

  addAndMakeVisible( samplePlayerTopLeftButton = new PlotPreviewButton(
    juce::String("TopLeftSamplePlayer"), NULL) );
  samplePlayerTopLeftButton->addRButtonListener(this);
  //samplePlayerTopLeftButton->setRadioGroupId(1);
  samplePlayerTopLeftButton->setDescription(
    juce::String(("Show editor for the top left sample-player")));
  samplePlayerTopLeftButton->setDescriptionField(infoField);
  samplePlayerTopLeftButton->setClickingTogglesState(true);
  samplePlayerTopLeftButton->setToggleState(true, false);

  addAndMakeVisible( samplePlayerTopRightButton = new PlotPreviewButton(
    juce::String("TopRightSamplePlayer")) );
  samplePlayerTopRightButton->addRButtonListener(this);
  //samplePlayerTopRightButton->setRadioGroupId(1);
  samplePlayerTopRightButton->setDescription(
    juce::String(("Show editor for the top right sample-player")));
  samplePlayerTopRightButton->setDescriptionField(infoField);
  samplePlayerTopRightButton->setClickingTogglesState(true);
  samplePlayerTopRightButton->setToggleState(false, false);

  addAndMakeVisible( samplePlayerBottomLeftButton = new PlotPreviewButton(
    juce::String("BottomLeftSamplePlayer")) );
  samplePlayerBottomLeftButton->addRButtonListener(this);
  //samplePlayerBottomLeftButton->setRadioGroupId(1);
  samplePlayerBottomLeftButton->setDescription(
    juce::String(("Show editor for the bottom left sample-player")));
  samplePlayerBottomLeftButton->setDescriptionField(infoField);
  samplePlayerBottomLeftButton->setClickingTogglesState(true);
  samplePlayerBottomLeftButton->setToggleState(false, false);

  addAndMakeVisible( samplePlayerBottomRightButton = new PlotPreviewButton(
    juce::String("BottomRightSamplePlayer")) );
  samplePlayerBottomRightButton->addRButtonListener(this);
  //samplePlayerBottomRightButton->setRadioGroupId(1);
  samplePlayerBottomRightButton->setDescription(
    juce::String(("Show editor for the bottom right sample-player")));
  samplePlayerBottomRightButton->setDescriptionField(infoField);
  samplePlayerBottomRightButton->setClickingTogglesState(true);
  samplePlayerBottomRightButton->setToggleState(false, false);

  addAndMakeVisible( lowFreqOscXButton = new PlotPreviewButton(juce::String("X-LFO")) );
  lowFreqOscXButton->addRButtonListener(this);
  //lowFreqOscXButton->setRadioGroupId(1);
  lowFreqOscXButton->setDescription(
    juce::String(("Show editor for the LFO for the x-coordinate")));
  lowFreqOscXButton->setDescriptionField(infoField);
  lowFreqOscXButton->setClickingTogglesState(true);
  lowFreqOscXButton->setToggleState(false, false);

  addAndMakeVisible( lowFreqOscYButton = new PlotPreviewButton(juce::String("Y-LFO")) );
  lowFreqOscYButton->addRButtonListener(this);
  //lowFreqOscYButton->setRadioGroupId(1);
  lowFreqOscYButton->setDescription(
    juce::String(("Show editor for the LFO for the y-coordinate")));
  lowFreqOscYButton->setDescriptionField(infoField);
  lowFreqOscYButton->setClickingTogglesState(true);
  lowFreqOscYButton->setToggleState(false, false);

  addAndMakeVisible( filterButton = new PlotPreviewButton(juce::String("Filter")) );
  filterButton->addRButtonListener(this);
  //filterButton->setRadioGroupId(1);
  filterButton->setDescription(juce::String(("Show editor for the filter")));
  filterButton->setDescriptionField(infoField);
  filterButton->setClickingTogglesState(true);
  filterButton->setToggleState(false, false);

  addAndMakeVisible( pitchEnvButton = new PlotPreviewButton(juce::String("PitchEnv")) );
  pitchEnvButton->addRButtonListener(this);
  //pitchEnvButton->setRadioGroupId(1);
  pitchEnvButton->setDescription(juce::String(("Show editor for the pitch envelope")));
  pitchEnvButton->setDescriptionField(infoField);
  pitchEnvButton->setClickingTogglesState(true);
  pitchEnvButton->setToggleState(false, false);

  addAndMakeVisible( filterEnvButton = new PlotPreviewButton(juce::String("FilterEnv")) );
  filterEnvButton->addRButtonListener(this);
  //filterEnvButton->setRadioGroupId(1);
  filterEnvButton->setDescription(juce::String(("Show editor for the filter envelope")));
  filterEnvButton->setDescriptionField(infoField);
  filterEnvButton->setClickingTogglesState(true);
  filterEnvButton->setToggleState(false, false);

  addAndMakeVisible( ampEnvButton = new PlotPreviewButton(juce::String("AmpEnv")) );
  ampEnvButton->addRButtonListener(this);
  //ampEnvButton->setRadioGroupId(1);
  ampEnvButton->setDescription(juce::String(("Show editor for the amplitude envelope")));
  ampEnvButton->setDescriptionField(infoField);
  ampEnvButton->setClickingTogglesState(true);
  ampEnvButton->setToggleState(false, false);

  addAndMakeVisible( masterFilterButton = new PlotPreviewButton(juce::String("MasterFilter")) );
  masterFilterButton->addRButtonListener(this);
  //masterFilterButton->setRadioGroupId(1);
  masterFilterButton->setDescription(juce::String(("Show editor for the master filter")));
  masterFilterButton->setDescriptionField(infoField);
  masterFilterButton->setClickingTogglesState(true);
  masterFilterButton->setToggleState(false, false);

  addAndMakeVisible( masterFilterEnvButton = new PlotPreviewButton(juce::String("MasterFilterEnv")) );
  masterFilterEnvButton->addRButtonListener(this);
  //masterFilterEnvButton->setRadioGroupId(1);
  masterFilterEnvButton->setDescription(
    juce::String(("Show editor for the master-filters envelope")));
  masterFilterEnvButton->setDescriptionField(infoField);
  masterFilterEnvButton->setClickingTogglesState(true);
  masterFilterEnvButton->setToggleState(false, false);

  addAndMakeVisible( masterAmpEnvButton = new PlotPreviewButton(juce::String("MasterAmpEnv")) );
  masterAmpEnvButton->addRButtonListener(this);
  //masterAmpEnvButton->setRadioGroupId(1);
  masterAmpEnvButton->setDescription(
    juce::String(("Show editor for the master amplitude envelope/gate")));
  masterAmpEnvButton->setDescriptionField(infoField);
  masterAmpEnvButton->setClickingTogglesState(true);
  masterAmpEnvButton->setToggleState(false, false);

  addAndMakeVisible( equalizerButton = new PlotPreviewButton(juce::String("Equalizer")) );
  equalizerButton->addRButtonListener(this);
  //equalizerButton->setRadioGroupId(1);
  equalizerButton->setDescription(juce::String(("Show editor for the equalizer")));
  equalizerButton->setDescriptionField(infoField);
  equalizerButton->setClickingTogglesState(true);
  equalizerButton->setToggleState(false, false);

  addAndMakeVisible( delayButton = new PlotPreviewButton(juce::String("Delay")) );
  delayButton->addRButtonListener(this);
  //delayButton->setRadioGroupId(1);
  delayButton->setDescription(juce::String(("Show editor for the delay")));
  delayButton->setDescriptionField(infoField);
  delayButton->setClickingTogglesState(true);
  delayButton->setToggleState(false, false);

  addAndMakeVisible( reverbButton = new PlotPreviewButton(juce::String("Reverb")) );
  reverbButton->addRButtonListener(this);
  //reverbButton->setRadioGroupId(1);
  reverbButton->setDescription(juce::String(("Show editor for the reverb")));
  reverbButton->setDescriptionField(infoField);
  reverbButton->setClickingTogglesState(true);
  reverbButton->setToggleState(false, false);

  samplePlayerTopLeftEditor->setActiveDirectory(
    getApplicationDirectory() + juce::String(("/Samples")) );

  // set up the widgets:
  updateWidgetsAccordingToState();
}

/*
WorkhorseModuleEditor::~WorkhorseModuleEditor()
{
deleteAllChildren();
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

void WorkhorseModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == samplePlayerTopLeftButton )
  {
    makeSubEditorsInvisible();
    samplePlayerTopLeftEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == samplePlayerTopRightButton )
  {
    makeSubEditorsInvisible();
    samplePlayerTopRightEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == samplePlayerBottomLeftButton )
  {
    makeSubEditorsInvisible();
    samplePlayerBottomLeftEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == samplePlayerBottomRightButton )
  {
    makeSubEditorsInvisible();
    samplePlayerBottomRightEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == lowFreqOscXButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == lowFreqOscYButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == filterButton )
  {
    makeSubEditorsInvisible();
    filterEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == pitchEnvButton )
  {
    makeSubEditorsInvisible();
    pitchEnvEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == filterEnvButton )
  {
    makeSubEditorsInvisible();
    filterEnvEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == ampEnvButton )
  {
    makeSubEditorsInvisible();
    ampEnvEditor->setVisible(true);
  }
  else if( buttonThatWasClicked == masterFilterButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == masterFilterEnvButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == masterAmpEnvButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == equalizerButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == delayButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else if( buttonThatWasClicked == reverbButton )
  {
    makeSubEditorsInvisible();
    //...
  }
  else
    PolyphonicInstrumentEditor::rButtonClicked(buttonThatWasClicked);


  repaint();
}

void WorkhorseModuleEditor::changeListenerCallback(ChangeBroadcaster* objectThatHasChanged)
{

}

void WorkhorseModuleEditor::updateWidgetsAccordingToState()
{
  if( workhorseAudioModule == NULL )
    return;
  if( workhorseAudioModule->wrappedWorkhorse == NULL )
    return;

  // remember if the preset was clean or dirty before making a few calls that may lead to a 
  // dirtification of the preset-state:
  bool presetIsDirty = workhorseAudioModule->isStateDirty();
  const MessageManagerLock mmLock;     
  // the event loop will now be locked so it's safe to make a few calls..

  // update the global widgets and automatable sliders:
  PolyphonicInstrumentEditor::updateWidgetsAccordingToState();

  // update the sub-editors:
  samplePlayerTopLeftEditor->updateWidgetsAccordingToState();
  samplePlayerTopRightEditor->updateWidgetsAccordingToState();
  samplePlayerBottomLeftEditor->updateWidgetsAccordingToState();
  samplePlayerBottomRightEditor->updateWidgetsAccordingToState();
  filterEditor->updateWidgetsAccordingToState();
  pitchEnvEditor->updateWidgetsAccordingToState();
  filterEnvEditor->updateWidgetsAccordingToState();
  ampEnvEditor->updateWidgetsAccordingToState();
  // todo: update the XY Pad

  // preserve the clean/dirty state of the preset regardless of any parameter changes that may take
  // place that may take place - note that not sending a change notification from the widgets is 
  // not enough to make that sure because some of the have AutomatableParameters associated with 
  // them which themselves may dirtify the preset:
  if( presetIsDirty )  
    workhorseAudioModule->markStateAsDirty();
  else
    workhorseAudioModule->markStateAsClean();
}

void WorkhorseModuleEditor::makeSubEditorsInvisible()
{
  samplePlayerTopLeftEditor->setVisible(false);
  samplePlayerTopRightEditor->setVisible(false);
  samplePlayerBottomLeftEditor->setVisible(false);
  samplePlayerBottomRightEditor->setVisible(false);
  filterEditor->setVisible(false);
  pitchEnvEditor->setVisible(false);
  filterEnvEditor->setVisible(false);
  ampEnvEditor->setVisible(false);
}

void WorkhorseModuleEditor::paint(Graphics &g)
{
  PolyphonicInstrumentEditor::paint(g);

  g.setColour(Colours::blue);
  g.drawRect(filterEditor->getBounds(), 3);

  int x1, x2, y1, y2;

  if( samplePlayerTopLeftButton->getToggleState() == true )
  {
    x1 = x2 = samplePlayerTopLeftButton->getX();
    y1 =      samplePlayerTopLeftButton->getY();
  }
  else if( samplePlayerTopRightButton->getToggleState() == true )
  {
    x1 = x2 = samplePlayerTopRightButton->getRight();
    y1 =      samplePlayerTopRightButton->getY();
  }
  else if( samplePlayerBottomLeftButton->getToggleState() == true )
  {
    x1 = x2 = samplePlayerBottomLeftButton->getX();
    y1 =      samplePlayerBottomLeftButton->getY();
  }
  else if( samplePlayerBottomRightButton->getToggleState() == true )
  {
    x1 = x2 = samplePlayerBottomRightButton->getRight();
    y1 =      samplePlayerBottomRightButton->getY();
  }

  else if( lowFreqOscXButton->getToggleState() == true )
  {
    x1 = x2 = lowFreqOscXButton->getX();
    y1 =      lowFreqOscXButton->getY();
  }
  else if( lowFreqOscYButton->getToggleState() == true )
  {
    x1 = x2 = lowFreqOscYButton->getX();
    y1 =      lowFreqOscYButton->getY();
  }
  else if( filterButton->getToggleState() == true )
  {
    x1 = x2 = filterButton->getX();
    y1 =      filterButton->getY();
  }
  else if( pitchEnvButton->getToggleState() == true )
  {
    x1 = x2 = pitchEnvButton->getRight();
    y1 =      pitchEnvButton->getY();
  }
  else if( filterEnvButton->getToggleState() == true )
  {
    x1 = x2 = filterEnvButton->getRight();
    y1 =      filterEnvButton->getY();
  }
  else if( ampEnvButton->getToggleState() == true )
  {
    x1 = x2 = ampEnvButton->getRight();
    y1 =      ampEnvButton->getY();
  }
  else if( masterFilterButton->getToggleState() == true )
  {
    x1 = x2 = masterFilterButton->getX();
    y1 =      masterFilterButton->getY();
  }
  else if( masterFilterEnvButton->getToggleState() == true )
  {
    x1 = x2 = masterFilterEnvButton->getX();
    y1 =      masterFilterEnvButton->getY();
  }
  else if( masterAmpEnvButton->getToggleState() == true )
  {
    x1 = x2 = masterAmpEnvButton->getX();
    y1 =      masterAmpEnvButton->getY();
  }
  else if( equalizerButton->getToggleState() == true )
  {
    x1 = x2 = equalizerButton->getRight();
    y1 =      equalizerButton->getY();
  }
  else if( delayButton->getToggleState() == true )
  {
    x1 = x2 = delayButton->getRight();
    y1 =      delayButton->getY();
  }
  else if( reverbButton->getToggleState() == true )
  {
    x1 = x2 = reverbButton->getRight();
    y1 =      reverbButton->getY();
  }



  else
  {
    x1 = x2 = 0;
    y1 = y2 = 0;
  }

  y2 = filterEditor->getY();
  g.drawLine((float) (x1-0), (float) y1, (float) (x2-0), (float) y2, 8.f);
  g.drawRect(filterEditor->getX()-4, filterEditor->getY()-4, 
    filterEditor->getWidth()+8, filterEditor->getHeight()+8, 4);

}


void WorkhorseModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();


  //webLink->setBounds(w-112, 0, 112-4, 20);

  x = 0;
  y = getHeadlineBottom();
  w = getWidth()/4;

  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+4, w-8, 40);
  /*
  stateWidgetSet->stateLabel->setBounds(x+4, 4, 52, 20);
  stateWidgetSet->stateFileNameLabel->setBounds(stateWidgetSet->stateLabel->getX(), 4, w-8, 20);
  stateWidgetSet->statePlusButton->setBounds(w-4-20, 24, 20, 20);
  stateWidgetSet->stateMinusButton->setBounds(stateWidgetSet->statePlusButton->getX()-20, 24, 20, 20);
  stateWidgetSet->stateLoadButton->setBounds(stateWidgetSet->stateMinusButton->getX()-40-4, 24, 40, 20);
  stateWidgetSet->stateSaveButton->setBounds(stateWidgetSet->stateLoadButton->getX()-40-4, 24, 40, 20);
  //y = stateWidgetSet->stateLabel->getBottom();
  */

  x  = w;
  //w /= 2;
  y = getHeadlineBottom()+4;
  levelSlider->setBounds(       x+4,     y, w-8,   16);
  y += 16;
  levelByKeySlider->setBounds(  x+4,     y, w/2-8, 12);
  levelByVelSlider->setBounds(  x+w/2+4, y, w/2-8, 12);

  y += 16;
  w /= 3;
  midSideRatioSlider->setBounds(      x+4, y, w-8, 12);
  x = midSideRatioSlider->getRight();
  numVoicesSlider->setBounds(         x+4, y, w-8, 12);
  x = numVoicesSlider->getRight();
  compSlider->setBounds(              x+4, y, w,   12);

  x = getWidth()/2;
  y = getHeadlineBottom();
  w = getWidth()/4;

  tuningLabel->setBounds(x+4, y+4, 80, 20);
  tuningPlusButton->setBounds(x+w-4-20, y+4, 20, 20);
  tuningMinusButton->setBounds(tuningPlusButton->getX()-20, y+4, 20, 20);
  tuningLoadButton->setBounds(tuningMinusButton->getX()-40-4, y+4, 40, 20);
  y = tuningLabel->getBottom()-2;
  tuningFileNameLabel->setBounds(tuningLabel->getX(), y, w-8, 20);

  y = getHeadlineBottom();
  x = x+w;
  glideButton->setBounds(x+4, y+4, 56, 20);
  glideTimeSlider->setBounds(glideButton->getRight()+4, y+4, w-glideButton->getWidth()-12, 20);

  w /= 2;
  y += 24;
  masterTuneSlider->setBounds(x+4, y+4, w-8, 20);
  x += w;
  wheelRangeSlider->setBounds(x+4, y+4, w-8, 20);

  // the preview components for the editors (preliminary as mostly as buttons):

  x = 4;
  //y = stateWidgetSet->stateFileNameLabel->getBottom()+8;
  //y = getPresetSectionBottom();
  y = numVoicesSlider->getBottom()+8;

  int xyPadEdgeLength = 180; // should be divisible by 2 and by 3
  int halfEdgeLength  = xyPadEdgeLength/2;
  int thirdEdgeLength = xyPadEdgeLength/3;

  int samplerPreviewWidth  = 120;
  samplePlayerTopLeftButton->setBounds(   x, y,                samplerPreviewWidth, halfEdgeLength);
  samplePlayerBottomLeftButton->setBounds(x, y+halfEdgeLength, samplerPreviewWidth, halfEdgeLength);
  x = samplePlayerTopLeftButton->getRight();
  y = samplePlayerTopLeftButton->getY();
  vectorMixerPad->setBounds(x, y, xyPadEdgeLength, xyPadEdgeLength);
  x = vectorMixerPad->getRight();
  //x += xyPadEdgeLength;
  samplePlayerTopRightButton->setBounds(   x, y,                samplerPreviewWidth, halfEdgeLength);
  samplePlayerBottomRightButton->setBounds(x, y+halfEdgeLength, samplerPreviewWidth, halfEdgeLength);

  int filterPreviewWidth  = 120;
  x = samplePlayerTopRightButton->getRight()+8;
  lowFreqOscXButton->setBounds(x, y,                   filterPreviewWidth, thirdEdgeLength);
  lowFreqOscYButton->setBounds(x, y+thirdEdgeLength,   filterPreviewWidth, thirdEdgeLength);
  filterButton->setBounds(     x, y+2*thirdEdgeLength, filterPreviewWidth, thirdEdgeLength);

  int envPreviewWidth  = 120;
  x = lowFreqOscXButton->getRight();
  pitchEnvButton->setBounds( x, y,                   envPreviewWidth, thirdEdgeLength);

  filterEnvButton->setBounds(x, y+thirdEdgeLength,   envPreviewWidth, thirdEdgeLength);
  ampEnvButton->setBounds(   x, y+2*thirdEdgeLength, envPreviewWidth, thirdEdgeLength);

  int effectsPreviewWidth = 120;
  x = pitchEnvButton->getRight()+8;
  masterFilterButton->setBounds(   x, y,                   effectsPreviewWidth, thirdEdgeLength);
  masterFilterEnvButton->setBounds(x, y+thirdEdgeLength,   effectsPreviewWidth, thirdEdgeLength);
  masterAmpEnvButton->setBounds(   x, y+2*thirdEdgeLength, effectsPreviewWidth, thirdEdgeLength);
  x += effectsPreviewWidth;
  equalizerButton->setBounds(      x, y,                   effectsPreviewWidth, thirdEdgeLength);
  delayButton->setBounds(          x, y+thirdEdgeLength,   effectsPreviewWidth, thirdEdgeLength);
  reverbButton->setBounds(         x, y+2*thirdEdgeLength, effectsPreviewWidth, thirdEdgeLength);

  //infoLabel->setBounds(0, getHeight()-20, 40, 20);
  //infoField->setBounds(infoLabel->getRight(), getHeight()-20, getWidth()-infoLabel->getRight(),20);

  x = 4;
  y = reverbButton->getBottom()+8;
  w = getWidth()-8;
  h = infoField->getY()-y-4;
  samplePlayerTopLeftEditor->setBounds(x, y, w, h);
  samplePlayerTopRightEditor->setBounds(x, y, w, h);
  samplePlayerBottomLeftEditor->setBounds(x, y, w, h);
  samplePlayerBottomRightEditor->setBounds(x, y, w, h);
  //...
  filterEditor->setBounds(x, y, w, h);
  pitchEnvEditor->setBounds(x, y, w, h);
  filterEnvEditor->setBounds(x, y, w, h);
  ampEnvEditor->setBounds(x, y, w, h);
}