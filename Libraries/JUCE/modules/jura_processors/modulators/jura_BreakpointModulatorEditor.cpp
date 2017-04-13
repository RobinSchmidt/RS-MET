//-------------------------------------------------------------------------------------------------
// construction/destruction:

BreakpointModulatorGlobalEditor::BreakpointModulatorGlobalEditor(CriticalSection *newPlugInLock,                                                                
  BreakpointModulatorAudioModule* newModulatorToEdit) 
  : AudioModuleEditor(newModulatorToEdit)
  //: AudioModuleEditor(newPlugInLock, newModulatorToEdit)
{
  layout = 0;

  // init the pointer to the modulator to be edited:
  modulatorToEdit = newModulatorToEdit;
  jassert(modulatorToEdit != NULL);

  addWidget( timeScaleSlider = new AutomatableSlider );
  timeScaleSlider->setSliderName("Time Scale");
  timeScaleSlider->setDescription("Scales overall time duration by a factor");
  timeScaleSlider->setStringConversionFunction(&valueToString4);
  timeScaleSlider->setRange(0.0625, 16.0, 0.01, 1.0);
  timeScaleSlider->setScaling(Parameter::EXPONENTIAL);

  addWidget( timeScaleByKeySlider = new AutomatableSlider );
  timeScaleByKeySlider->setSliderName("Key");
  timeScaleByKeySlider->setDescription("Key dependence of the the overall time duration");
  timeScaleByKeySlider->setStringConversionFunction(&percentToStringWithUnit0);
  timeScaleByKeySlider->setRange(-150.0, 150.0, 1.0, 0.0);
  timeScaleByKeySlider->setScaling(Parameter::LINEAR_BIPOLAR);

  addWidget( timeScaleByVelSlider = new AutomatableSlider );
  timeScaleByVelSlider->setSliderName("Vel");
  timeScaleByVelSlider->setDescription("Velocity dependence of the the overall time duration");
  timeScaleByVelSlider->setStringConversionFunction(&percentToStringWithUnit0);
  timeScaleByVelSlider->setRange(-150.0, 150.0, 1.0, 0.0);
  timeScaleByVelSlider->setScaling(Parameter::LINEAR_BIPOLAR);

  addWidget( depthSlider = new AutomatableSlider );
  depthSlider->setSliderName("Depth");
  depthSlider->setDescription("Depth of the modulation");
  depthSlider->setStringConversionFunction(&valueToString2);
  depthSlider->setRange(0.0, 4.0, 0.01, 1.0);

  addWidget( depthByKeySlider = new AutomatableSlider );
  depthByKeySlider->setSliderName("Key");
  depthByKeySlider->setDescription("Key dependence of the modulation depth");
  depthByKeySlider->setStringConversionFunction(&percentToStringWithUnit0);
  depthByKeySlider->setRange(-150.0, 150.0, 1.0, 0.0);
  depthByKeySlider->setScaling(Parameter::LINEAR_BIPOLAR);

  addWidget( depthByVelSlider = new AutomatableSlider );
  depthByVelSlider->setSliderName("Vel");
  depthByVelSlider->setDescription("Velocity dependence of the modulation depth");
  depthByVelSlider->setStringConversionFunction(&percentToStringWithUnit0);
  depthByVelSlider->setRange(-150.0, 150.0, 1.0, 0.0);
  depthByVelSlider->setScaling(Parameter::LINEAR_BIPOLAR);

  addWidget( loopButton = new AutomatableButton("Loop") );
  loopButton->setDescription("Toggle sustain loop on/off");
  loopButton->setClickingTogglesState(true);
  loopButton->addRButtonListener(this);

  addWidget( syncButton = new AutomatableButton("Sync") );
  syncButton->setDescription("Toggle sync on/off. Time unit is beats in sync-mode, seconds otherwise");
  syncButton->setClickingTogglesState(true);
  syncButton->addRButtonListener(this);

  addWidget( editButton = new RButton("Edit") );
  editButton->setDescription("Selects the envelope for editing");

  setHeadlineStyle(AudioModuleEditor::NO_HEADLINE);
  webLink->setVisible(false);
  infoField->setVisible(false);

  setModulatorToEdit(modulatorToEdit); // will also call updateWidgetsAccordingToState
}

//-------------------------------------------------------------------------------------------------
// setup:

void BreakpointModulatorGlobalEditor::setModulatorToEdit(
  BreakpointModulatorAudioModule* newModulatorToEdit)
{
  modulatorToEdit = newModulatorToEdit;

  if( modulatorToEdit == NULL )
    return;

  timeScaleSlider->assignParameter(modulatorToEdit->getParameterByName("TimeScale"));
  timeScaleByKeySlider->assignParameter(modulatorToEdit->getParameterByName("TimeScaleByKey"));
  timeScaleByKeySlider->setSliderName("Key");
  timeScaleByVelSlider->assignParameter(modulatorToEdit->getParameterByName("TimeScaleByVel"));
  timeScaleByVelSlider->setSliderName("Vel");

  depthSlider->assignParameter(modulatorToEdit->getParameterByName("Depth"));
  depthByKeySlider->assignParameter(modulatorToEdit->getParameterByName("DepthByKey"));
  depthByKeySlider->setSliderName("Key");
  depthByVelSlider->assignParameter(modulatorToEdit->getParameterByName("DepthByVel"));
  depthByVelSlider->setSliderName("Vel");

  setLinkPosition(INVISIBLE);

  updateWidgetsAccordingToState();
}

void BreakpointModulatorGlobalEditor::setLayout(int newLayout)
{
  layout = newLayout;
  resized();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void BreakpointModulatorGlobalEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  if( buttonThatWasClicked == loopButton )
    modulatorToEdit->wrappedBreakpointModulator->setLoopMode( loopButton->getToggleState() );  
  else if( buttonThatWasClicked == syncButton )
    modulatorToEdit->wrappedBreakpointModulator->setSyncMode( syncButton->getToggleState() );  

  moduleToEdit->markStateAsDirty();
}

void BreakpointModulatorGlobalEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  moduleToEdit->markStateAsDirty();
}

void BreakpointModulatorGlobalEditor::resized()
{
  int x, y, w, h;
  if( layout == 0 ) // bottom section of a single-modulator editor
  {
    AudioModuleEditor::resized(); // do we need this?

    //// old:
    //stateWidgetSet->stateLabel->setVisible(true);
    //stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);

    // new:
    stateWidgetSet->setVisible(false);
    stateWidgetSet->stateLabel->setVisible(false);

    x = 0;
    y = 2;
    w = (getWidth()-60) / 2;
    timeScaleSlider->setBounds(x+4, y+4, w-8, 16);
    y += 14;
    w /= 2;
    timeScaleByKeySlider->setBounds(x+4,   y+4, w-8, 16);
    x = x+w+4;
    w = timeScaleSlider->getRight()-x;
    timeScaleByVelSlider->setBounds(x, y+4, w, 16);

    x = 0;
    y = 2;
    w = (getWidth()-60) / 2;
    x = getWidth()-w;
    depthSlider->setBounds(x+4, y+4, w-8, 16);
    y += 14;
    w /= 2;
    depthByKeySlider->setBounds(x+4,   y+4, w-8, 16);
    x = x+w+4;
    w = depthSlider->getRight()-x;
    depthByVelSlider->setBounds(x, y+4, w, 16);

    x  = 0;
    y  = +4;
    w  = getWidth() / 3;
    x += getWidth() / 2;
    loopButton->setBounds(x-24, y, 48, 16);
    y += 20;
    syncButton->setBounds(x-24, y, 48, 16);
  }
  else if( layout == 1 )
  {
    stateWidgetSet->setVisible(true);
    stateWidgetSet->stateLabel->setVisible(true);
    stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
    //stateWidgetSet->setLayout(StateLoadSaveWidgetSet::layouts::ONE_LINE); // for debug
    AudioModuleEditor::resized();

    int right = getWidth();
    w = 40;
    x = right-40-4;
    y = 4;
    h = 16;
    editButton->setBounds(x, y, w, 20);
    y += 30;
    loopButton->setBounds(x, y, w, h);
    y += h+4;
    syncButton->setBounds(x, y, w, h);

    x = 4;
    y = 4;
    w = loopButton->getX()-x-4;
    h = 32;

    stateWidgetSet->setBounds(x, y, w, h);

    x = stateWidgetSet->getX();
    y = stateWidgetSet->getBottom()+4;
    w = stateWidgetSet->getRight()-x;
    h = 16;
    int w2 = w/2;
    timeScaleSlider->setBounds(x, y, w2-4, 16);
    int w3 = timeScaleSlider->getWidth();
    y += 14;
    timeScaleByKeySlider->setBounds(x,        y, w3/2-2, 16);
    timeScaleByVelSlider->setBounds(x+w3/2+2, y, w3/2-2, 16);

    y  = timeScaleSlider->getY();
    x += w2+4;
    depthSlider->setBounds(x, y, w2-4, 16);
    w3 = depthSlider->getWidth();
    y += 14;
    depthByKeySlider->setBounds(x,        y, w3/2-2, 16);
    depthByVelSlider->setBounds(x+w3/2+2, y, w3/2-2, 16);
  }
}

void BreakpointModulatorGlobalEditor::updateWidgetsAccordingToState()
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  // restore the slider settings:
  timeScaleSlider->setValue(     m->getTimeScale(),        false);
  timeScaleByKeySlider->setValue(m->getTimeScaleByKey(),   false);
  timeScaleByVelSlider->setValue(m->getTimeScaleByVel(),   false);
  depthSlider->setValue(         m->getDepth(),            false);
  depthByKeySlider->setValue(    m->getDepthByKey(),       false);
  depthByVelSlider->setValue(    m->getDepthByVel(),       false);

  // restore the button settings:
  loopButton->setToggleState((   m->getLoopMode() != 0),   false);
  syncButton->setToggleState(    m->isInSyncMode(),        false);

  stateWidgetSet->updateStateNameField();
}

//=================================================================================================
// class BreakpointParameterEditor:

BreakpointParameterEditor::BreakpointParameterEditor(CriticalSection *newPlugInLock)
//: AudioModuleEditor(newPlugInLock, NULL) // mmm...will this 'NULL' cause problems?
: AudioModuleEditor(newPlugInLock)
{
  modulatorToEdit         = NULL;
  selectedBreakpointIndex = -1;

  addWidget( indexLabel = new RTextField("Breakpoint") );
  indexLabel->setDescription("Index of selected breakpoint");
  indexLabel->setNoBackgroundAndOutline(true);
  indexLabel->setJustification(Justification::centred);

  addWidget( indexValueLabel = new RTextField() );
  indexValueLabel->setDescription("Index of selected breakpoint");
  indexValueLabel->setNoBackgroundAndOutline(true);
  indexValueLabel->setJustification(Justification::centred);

  addWidget( timeSlider = new RSlider("TimeSlider") );
  timeSlider->addListener(this);
  timeSlider->setSliderName("Time");
  timeSlider->setDescription("Time stamp of slected breakpoint (in seconds or beats)");
  timeSlider->setStringConversionFunction(&secondsToStringWithUnitTotal4);
  timeSlider->setRange(0.0, 5.0, 0.0001, 0.0);
  timeSlider->setLayout(RSlider::NAME_ABOVE);

  addWidget( levelSlider = new RSlider("LevelSlider") );
  levelSlider->addListener(this);
  levelSlider->setSliderName("Level");
  levelSlider->setDescription("Level of selected breakpoint");
  levelSlider->setStringConversionFunction(&valueToString3);
  //levelSlider->setRange(-2.0, 2.0, 0.001, 1.0);
  levelSlider->setRange(0.0, 4.0, 0.001, 1.0);
  levelSlider->setLayout(RSlider::NAME_ABOVE);

  std::vector<double> defaultValues;
  defaultValues.push_back(0.0);
  defaultValues.push_back(0.25);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.75);
  defaultValues.push_back(1.0);
  defaultValues.push_back(1.5);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  levelSlider->setDefaultValues(defaultValues);

  addWidget( shapeComboBox = new RNamedComboBox("ShapeComboBox", "Shape") );
  shapeComboBox->setDescription("Shape of the curve approaching the selected breakpoint");
  shapeComboBox->setNameLabelPosition(RNamedComboBox::ABOVE_BOX);
  shapeComboBox->registerComboBoxObserver(this);
  shapeComboBox->addItem(0, "Stairstep"  );
  shapeComboBox->addItem(1, "Linear"     );
  shapeComboBox->addItem(2, "Smooth"     );
  shapeComboBox->addItem(3, "Analog"     );
  shapeComboBox->addItem(4, "AntiAnalog" );
  shapeComboBox->addItem(5, "Sigmoid"    );
  shapeComboBox->addItem(6, "Spikey"     );
  shapeComboBox->addItem(7, "Sine 1"     );
  shapeComboBox->addItem(8, "Sine 2"     );

  addWidget( shapeSlider = new RSlider("ShapeSlider") );
  shapeSlider->addListener(this);
  shapeSlider->setSliderName(juce::String::empty);
  shapeSlider->setDescription("Amount of the shape (shapiness)");
  shapeSlider->setStringConversionFunction(&valueToString2);
  shapeSlider->setRange(0.1, 10.0, 0.01, 1.0);
  shapeSlider->setScaling(Parameter::EXPONENTIAL);

  addWidget( shapeToAllButton = new RButton("ToAll") );
  shapeToAllButton->setDescription("Apply the current shape setting to all breakpoints");
  shapeToAllButton->setClickingTogglesState(true);
  shapeToAllButton->addRButtonListener(this);
}

//-------------------------------------------------------------------------------------------------
// setup:

void BreakpointParameterEditor::setModulatorToEdit(
  BreakpointModulatorAudioModule* newModulatorToEdit)
{
  modulatorToEdit = newModulatorToEdit;
  deSelectBreakpoint();
  updateWidgetsAccordingToState();
}

void BreakpointParameterEditor::selectBreakpoint(int index)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  if( index == -1 )
    deSelectBreakpoint();
  else
  {
    jassert( index >= 0 && index < m->getNumBreakpoints() );
    if( index < 0 || index >= m->getNumBreakpoints() )
      return;

    selectedBreakpointIndex = index;
    indexLabel->setVisible(true);
    indexValueLabel->setVisible(true);
    timeSlider->setVisible(true);
    levelSlider->setVisible(true);
    shapeComboBox->setVisible(true);
    shapeSlider->setVisible(true);
    shapeToAllButton->setVisible(true);
    updateWidgetsAccordingToState();
  }
}

void BreakpointParameterEditor::deSelectBreakpoint()
{
  selectedBreakpointIndex = -1;
  indexLabel->setVisible(false);  
  indexValueLabel->setVisible(false);
  timeSlider->setVisible(false);
  levelSlider->setVisible(false);
  shapeComboBox->setVisible(false);
  shapeSlider->setVisible(false);
  shapeToAllButton->setVisible(false);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void BreakpointParameterEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  if( buttonThatWasClicked == shapeToAllButton )
  {
    if( shapeToAllButton->getToggleState() == true )
    {
      for(int p=0; p<=m->lastBreakpointIndex(); p++)
      {
        m->setBreakpointShapeAmount(p, shapeSlider->getValue());
        m->setBreakpointShape(      p, shapeComboBox->getSelectedItemIdentifier());
      }
    }
  }

  modulatorToEdit->markStateAsDirty();
  sendChangeMessage();
}

void BreakpointParameterEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  if( rComboBoxThatHasChanged == shapeComboBox )
  {
    int newShape = shapeComboBox->getSelectedItemIdentifier();

    // enable the shape-slider only where it applies:
    shapeSlider->setEnabled(false);
    if( newShape >= 3 && newShape <= 6 ) // 3: analog, 4: anti-analog, 5: sigmoid, 6: spikey
      shapeSlider->setEnabled(true);

    if( shapeToAllButton->getToggleState() == false )
      m->setBreakpointShape(selectedBreakpointIndex, newShape);
    else
    {
      for(int p=0; p<=m->lastBreakpointIndex(); p++)
        m->setBreakpointShape(p, newShape);
    }
  }

  modulatorToEdit->markStateAsDirty(); 
  sendChangeMessage();
}

void BreakpointParameterEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  if( sliderThatHasChanged == timeSlider )
    m->setBreakpointTime(selectedBreakpointIndex, timeSlider->getValue());
  else if( sliderThatHasChanged == levelSlider )
    m->setBreakpointLevel(selectedBreakpointIndex, levelSlider->getValue());
  else if( sliderThatHasChanged == shapeSlider )
  {
    if( shapeToAllButton->getToggleState() == false )
      m->setBreakpointShapeAmount( selectedBreakpointIndex, shapeSlider->getValue());
    else
    {
      for(int p=0; p<=m->lastBreakpointIndex(); p++)
        m->setBreakpointShapeAmount(p, shapeSlider->getValue());
    }
  }

  modulatorToEdit->markStateAsDirty();  
  sendChangeMessage();
}

void BreakpointParameterEditor::resized()
{
  int x = 4;
  int y = 4;
  int w = getWidth()-8;
  indexLabel->setBounds(      x-4, y, w, 16);
  y += 16;
  indexValueLabel->setBounds( x,   y, w, 16);
  y += 32;
  timeSlider->setBounds(      x,   y, w, 32);
  y += 36;
  levelSlider->setBounds(     x,   y, w, 32);
  y += 48;
  shapeComboBox->setBounds(   x,     y, w,    32);
  shapeToAllButton->setBounds(shapeComboBox->getRight()-44, y+2, 44, 16);
  y = shapeComboBox->getBottom()-2;
  shapeSlider->setBounds(     x,   y, w, 16);
}

void BreakpointParameterEditor::updateWidgetsAccordingToState()
{
  if( modulatorToEdit == NULL )
    return;
  if( modulatorToEdit->wrappedBreakpointModulator == NULL )
    return;

  RAPT::rsBreakpointModulator<double>* m = modulatorToEdit->wrappedBreakpointModulator;

  indexValueLabel->setText(valueToString0(selectedBreakpointIndex+1) + juce::String("/") 
    + valueToString0(m->getNumBreakpoints()));

  double minTime = m->getBreakpointMinTime(selectedBreakpointIndex);
  double maxTime = m->getBreakpointMaxTime(selectedBreakpointIndex);
  if( minTime >= maxTime ) // minTime == maxTime can occur for simultaneous breakpoints
    maxTime = minTime + 0.00001;
  timeSlider->setRange(minTime, maxTime, 0.001, 0.0, false);
  timeSlider->setValue(m->getBreakpointTime(selectedBreakpointIndex), false);

  double minLevel = jmin(0.0, m->getMinLevel());
  double maxLevel = m->getMaxLevel();
  //levelSlider->setRange(minLevel, maxLevel, 0.001, 1.0, false);
  levelSlider->setValue(m->getBreakpointLevel(selectedBreakpointIndex), false);

  shapeSlider->setValue(m->getBreakpointShapeAmount(selectedBreakpointIndex),  false);
  //int shapeIndex = rmax(m->getBreakpointShape(selectedBreakpointIndex), 
  //  (int) rsModBreakpoint::STAIRSTEP);
  int shapeIndex = jmax(m->getBreakpointShape(selectedBreakpointIndex), 
    (int) rsModBreakpoint<double>::STAIRSTEP);
  shapeComboBox->selectItemByIndex(shapeIndex, false);
  shapeSlider->setEnabled(false);
  if( shapeIndex >= 3 && shapeIndex <= 6 ) 
    shapeSlider->setEnabled(true);
}

//=================================================================================================
// class BreakpointModulatorEditor:

BreakpointModulatorEditor::BreakpointModulatorEditor(CriticalSection *newPlugInLock, 
  BreakpointModulatorAudioModule* newBreakpointModulatorAudioModule)                                                    
: AudioModuleEditor(newBreakpointModulatorAudioModule)
//  : AudioModuleEditor(newPlugInLock, newBreakpointModulatorAudioModule)
{
  modulatorToEdit = NULL; // ? old and obsolete ?
  jassert(newBreakpointModulatorAudioModule != NULL);
  setLinkPosition(INVISIBLE);

  // change the headline from the default "Sub-Editor" to "Modulator-Editor":
  setHeadlineText("Modulator");
  setDescription("This is a multi-breakpoint modulation generator");

  // create the breakpoint-editor:
  breakpointEditor = new ModulatorCurveEditor("PlotEditor");
  breakpointEditor->addChangeListener(this);
  addPlot(breakpointEditor);

  // create the zoomer for the breakpointEditor and associate it with the breakpointEditor:
  breakpointZoomer = new CoordinateSystemZoomerOld();
  breakpointZoomer->setRelativeMargins(5.0, 5.0, 10.0, 10.0);
  breakpointZoomer->setVerticalMouseWheelMode(CoordinateSystemZoomerOld
    ::horizontalZoomViaVerticalMouseWheel);
  addChildColourSchemeComponent(breakpointZoomer);
  breakpointZoomer->setCoordinateSystem(breakpointEditor);

  globalEditor = new BreakpointModulatorGlobalEditor(lock, newBreakpointModulatorAudioModule);
  globalEditor->loopButton->addRButtonListener(this); 
  addChildEditor( globalEditor );

  breakpointParameterEditor = new BreakpointParameterEditor(lock);
  breakpointParameterEditor->addChangeListener(this);
  breakpointParameterEditor->setModulatorToEdit(newBreakpointModulatorAudioModule);
  addChildEditor( breakpointParameterEditor );

  addWidget( snapXButton = new RButton("#X:") );
  snapXButton->setDescription("Toggle time-quantization on/off.");
  snapXButton->setClickingTogglesState(true);
  snapXButton->addRButtonListener(this);

  addWidget( snapXComboBox = new RComboBox("SnapXComboBox") );
  snapXComboBox->registerComboBoxObserver(this);
  snapXComboBox->setDescription("Select spacing of the vertical grid lines");
  snapXComboBox->addItem(0, "1/2"   );
  snapXComboBox->addItem(1, "1/4"   );
  snapXComboBox->addItem(2, "1/8"   );
  snapXComboBox->addItem(3, "0.1"   );
  snapXComboBox->addItem(4, "1/16"  );
  snapXComboBox->addItem(5, "1/32"  );
  snapXComboBox->addItem(6, "1/64"  );
  snapXComboBox->addItem(7, "0.01"  );
  snapXComboBox->addItem(8, "1/128" );
  snapXComboBox->selectItemByIndex(2, false);

  addWidget( snapYButton = new RButton("#Y:") );
  snapYButton->setDescription("Toggle level-quantization on/off.");
  snapYButton->setClickingTogglesState(true);
  snapYButton->addRButtonListener(this);

  addWidget( snapYComboBox = new RComboBox("SnapYComboBox") );
  snapYComboBox->registerComboBoxObserver(this);
  snapYComboBox->setDescription("Select spacing of the horizontal grid lines");
  snapYComboBox->addItem(0, "1/2"   );
  snapYComboBox->addItem(1, "1/4"   );
  snapYComboBox->addItem(2, "1/8"   );
  snapYComboBox->addItem(3, "0.1"   );
  snapYComboBox->addItem(4, "1/16"  );
  snapYComboBox->addItem(5, "1/32"  );
  snapYComboBox->addItem(6, "1/64"  );
  snapYComboBox->addItem(7, "0.01"  );
  snapYComboBox->addItem(8, "1/128" );
  snapYComboBox->selectItemByIndex(2, false);
   // get rid of this code duplication...

  addWidget(closeButton = new RButton(RButton::CLOSE), true, false); // invisible by default
  closeButton->setDescription("Closes the modulator editor");
  closeButton->setClickingTogglesState(false);
  // we don't listen to this button ourselves - this is the job of the outlying editor object

  // customize the descriptions for the load/save buttons:
  stateWidgetSet->stateLoadButton->setDescription("Load modulator settings from file");
  stateWidgetSet->stateSaveButton->setDescription("Save modulator settings to file");
  stateWidgetSet->statePlusButton->setDescription("Skip to next modulation curve in current directory");
  stateWidgetSet->stateMinusButton->setDescription("Skip to previous modulation curve in current directory");
  stateWidgetSet->stateFileNameLabel->setDescription("Name of current preset for the breakpoint modulator (if any)");

  setModulatorToEdit(newBreakpointModulatorAudioModule->wrappedBreakpointModulator);
     // this will also set up the widgets according to the state of the modulator
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void BreakpointModulatorEditor::setModulatorToEdit(RAPT::rsBreakpointModulator<double>* newModulatorToEdit)
{
  modulatorToEdit = newModulatorToEdit;
  //setPresetRemembererToEdit( modulatorToEdit );
  //updatePresetField();
  breakpointEditor->setModulatorToEdit(newModulatorToEdit);
  //breakpointEditor->updateMaximumRange(true);
  //breakpointEditor->setMaximumRange(breakpointEditor->getMaximumMeaningfulRange(10.0, 5.0, 5.0, 10.0));
  breakpointZoomer->zoomToAllXY();

  updateWidgetsAccordingToState(true);
}

/*
void BreakpointModulatorEditor::setDescriptionField(RLabel *newDescriptionField)
{
  Editor::setDescriptionField(newDescriptionField);

  stateWidgetSet->stateLoadButton->setDescriptionField(newDescriptionField);
  stateWidgetSet->stateSaveButton->setDescriptionField(newDescriptionField);
  stateWidgetSet->statePlusButton->setDescriptionField(newDescriptionField);
  stateWidgetSet->stateMinusButton->setDescriptionField(newDescriptionField);
  stateWidgetSet->stateFileNameLabel->setDescriptionField(newDescriptionField);

  breakpointEditor->setDescriptionField(newDescriptionField);
  breakpointZoomer->setWidgetDescriptionField(newDescriptionField);

  breakpointParameterEditor->setDescriptionField(newDescriptionField);
  globalEditor->setDescriptionField(newDescriptionField);

  snapXButton->setDescriptionField(newDescriptionField);
  snapXComboBox->setDescriptionField(newDescriptionField);
  snapYButton->setDescriptionField(newDescriptionField);
  snapYComboBox->setDescriptionField(newDescriptionField);
}
*/

//-------------------------------------------------------------------------------------------------
// callbacks:

void BreakpointModulatorEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( modulatorToEdit == NULL )
    return;

  /*
  if( buttonThatWasClicked == breakpointParameterEditor->shapeToAllButton )
  {
    if( breakpointParameterEditor->shapeToAllButton->getToggleState() == true )
    {
      breakpointEditor->updatePlotImage();
    }
  }
  */
  if( buttonThatWasClicked == globalEditor->loopButton )
    breakpointEditor->updatePlotImage();
  else if( buttonThatWasClicked == snapXButton )
  {
    breakpointEditor->setSnapToFineGridX(snapXButton->getToggleState());
    breakpointEditor->setVerticalFineGridVisible(snapXButton->getToggleState());
  }
  else if( buttonThatWasClicked == snapYButton )
  {
    breakpointEditor->setSnapToFineGridY(snapYButton->getToggleState());
    breakpointEditor->setHorizontalFineGridVisible(snapYButton->getToggleState());
  }
  moduleToEdit->markStateAsDirty();
}

void BreakpointModulatorEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;

  if( objectThatHasChanged == stateWidgetSet )
  {
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);  
    sendChangeMessage();
    return;
  }
  else if( objectThatHasChanged == breakpointEditor )
  {
    breakpointParameterEditor->selectBreakpoint(breakpointEditor->getSelectedBreakpointIndex());
    breakpointZoomer->updateScrollbars();
  }
  else if( objectThatHasChanged == breakpointParameterEditor )
    breakpointEditor->updatePlotCurveData();

  moduleToEdit->markStateAsDirty();
  sendChangeMessage();
}

void BreakpointModulatorEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;

    /*
  // todo: move this into class BreapointParameterEditor:
  if( rComboBoxThatHasChanged == breakpointParameterEditor->shapeComboBox )
  {
    breakpointEditor->updatePlotImage();

    int newShape = breakpointParameterEditor->shapeComboBox->getSelectedId();

    // enable the shape-slider only where it applies:
    breakpointParameterEditor->shapeSlider->setEnabled(false);
    if( newShape >= 4 && newShape <= 7 ) // 4: analog, 5: anti-analog, 6: sigmoid, 7: spikey
      breakpointParameterEditor->shapeSlider->setEnabled(true);

    if( breakpointParameterEditor->shapeToAllButton->getToggleState() == false )
      breakpointEditor->setSelectedBreakpointShape(newShape);
    else
      breakpointEditor->setAllBreakpointShapes(newShape);
  }
    */

  if( rComboBoxThatHasChanged == snapXComboBox )
  {
    int newGridIntervalIndex = snapXComboBox->getSelectedItemIdentifier();
    breakpointEditor->setVerticalFineGrid(gridIntervalFromIndex(newGridIntervalIndex), true);
    breakpointEditor->repaint();
  }
  else if( rComboBoxThatHasChanged == snapYComboBox )
  {
    int newGridIntervalIndex = snapYComboBox->getSelectedItemIdentifier();
    breakpointEditor->setHorizontalFineGrid(gridIntervalFromIndex(newGridIntervalIndex), true);
    breakpointEditor->repaint();
  }

  moduleToEdit->markStateAsDirty();
}

/*
void BreakpointModulatorEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( modulatorToEdit == NULL )
    return;

  // todo: move this into class BreakpointParameterEditor:
  if( sliderThatHasChanged == breakpointParameterEditor->timeSlider )
  {
    breakpointEditor->setSelectedBreakpointTime(breakpointParameterEditor->timeSlider->getValue());
    //breakpointEditor->updateModulationCurve();
  }
  else if( sliderThatHasChanged == breakpointParameterEditor->levelSlider )
  {
    breakpointEditor->setSelectedBreakpointLevel(breakpointParameterEditor->levelSlider->getValue(), false);
    //breakpointEditor->updateModulationCurve();
  }
  else if( sliderThatHasChanged == breakpointParameterEditor->shapeSlider )
  {
    if( breakpointParameterEditor->shapeToAllButton->getToggleState() == false )
      breakpointEditor->setSelectedBreakpointShapeAmount(breakpointParameterEditor->shapeSlider->getValue());
    else
      breakpointEditor->setAllBreakpointShapeAmounts(breakpointParameterEditor->shapeSlider->getValue());
      //breakpointEditor->updateModulationCurve();
  }

  moduleToEdit->markStateAsDirty();
}
*/

void BreakpointModulatorEditor::paint(Graphics &g)
{
  Editor::paint(g);

  // draw rectangles for the parameter-groups:
  int x = breakpointGroupRectangle.getX();
  int y = breakpointGroupRectangle.getY();
  int w = breakpointGroupRectangle.getWidth();
  int h = breakpointGroupRectangle.getHeight();
  fillRectWithBilinearGradient(g, x, y, jmax(1,w), jmax(1,h),
    editorColourScheme.topLeft, editorColourScheme.topRight, 
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  g.setColour(editorColourScheme.outline);
  //g.drawRect(x, y, w, h);

  x = timeAndDepthGroupRectangle.getX();
  y = timeAndDepthGroupRectangle.getY();
  w = timeAndDepthGroupRectangle.getWidth();
  h = timeAndDepthGroupRectangle.getHeight();
  fillRectWithBilinearGradient(g, x, y, jmax(1,w), jmax(1,h),
    editorColourScheme.topLeft, editorColourScheme.topRight, 
    editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
  g.setColour(editorColourScheme.outline);
  //g.drawRect(x, y, w, h);

  g.setColour(editorColourScheme.outline);
  g.drawRect(0, 0, getWidth(), getHeight());

  //Editor::drawHeadline(g);
}

void BreakpointModulatorEditor::resized()
{
  AudioModuleEditor::resized();

  closeButton->setBounds(getWidth()-16, 0, 16, 16); // for what is this?

  int bottomSectionHeight = 44;
  int rightSectionWidth   = 100;
  int lineWidth           = 2;   // we should infer this from the skin

  int x = getWidth()-rightSectionWidth;
  int y = getHeadlineBottom()+4;
  int w = getWidth()-x;
  int h = getHeight()-y-bottomSectionHeight;

  breakpointGroupRectangle.setBounds(x, y, w, h+lineWidth);
  timeAndDepthGroupRectangle.setBounds(0, getHeight()-bottomSectionHeight, 
    getWidth()-rightSectionWidth+lineWidth, bottomSectionHeight);

  //x = 4; // old
  x = 0;
  y = getHeadlineBottom()+4;
  w = breakpointGroupRectangle.getX()-x   + 4;  // why do we need the +4 here?
  h = timeAndDepthGroupRectangle.getY()-y + 4;  // ..and here?
  breakpointEditor->setBounds(x, y, 
    w-breakpointZoomer->getZoomerSize(), 
    h-breakpointZoomer->getZoomerSize());
  breakpointZoomer->alignWidgetsToCoordinateSystem();

  // the right section:
  breakpointParameterEditor->setBounds(breakpointGroupRectangle);

  // the bottom section:
  globalEditor->setBounds(timeAndDepthGroupRectangle);

  // the snap-stuff
  x = globalEditor->getRight();
  w = getWidth()-x;
  y = globalEditor->getY()+2;
  snapXButton->setBounds(x+4, y+4,    32, 16);
  snapYButton->setBounds(x+4, y+14+4, 32, 16);
  x = snapXButton->getRight()-2;
  w = getWidth()-x;
  snapXComboBox->setBounds(x, y+4,    w-6, 16);
  snapYComboBox->setBounds(x, y+14+4, w-6, 16);
}

void BreakpointModulatorEditor::deSelectBreakpoint()
{
  breakpointEditor->setSelectedBreakpointIndex(-1);
  breakpointParameterEditor->deSelectBreakpoint();
  /*
  // make the numeric breakpoint edit widgets invisible - todo: move this into class BreakpointParameterEditor:
  breakpointParameterEditor->indexValueLabel->setText(juce::String(T("none selected")), false);
  breakpointParameterEditor->timeSlider->setVisible(false);
  breakpointParameterEditor->levelSlider->setVisible(false);
  breakpointParameterEditor->shapeLabel->setVisible(false);
  breakpointParameterEditor->shapeComboBox->setVisible(false);
  breakpointParameterEditor->shapeSlider->setVisible(false);
  breakpointParameterEditor->shapeToAllButton->setVisible(false);
  */
}

void BreakpointModulatorEditor::updateWidgetsAccordingToState(bool deSelectBreakpointBefore)
{
  if( modulatorToEdit == NULL )
    return;

  if( deSelectBreakpointBefore == true )
    deSelectBreakpoint();

  globalEditor->updateWidgetsAccordingToState();
  breakpointParameterEditor->updateWidgetsAccordingToState();

  snapXButton->setToggleState(breakpointEditor->isVerticalFineGridVisible(), false);
  snapXComboBox->selectItemByIndex(
    indexFromGridInterval(breakpointEditor->getVerticalFineGridInterval())-1, false);
  snapYButton->setToggleState(breakpointEditor->isHorizontalFineGridVisible(), false);
  snapYComboBox->selectItemByIndex(
    indexFromGridInterval(breakpointEditor->getHorizontalFineGridInterval())-1, false);

  // update the plot:
  breakpointEditor->updateMaximumRange(true);
  breakpointEditor->updatePlotCurveData();

  stateWidgetSet->updateStateNameField();
}

void BreakpointModulatorEditor::updateWidgetsAccordingToState()
{
  updateWidgetsAccordingToState(true);
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void BreakpointModulatorEditor::autoAdjustPlotRangeX()
{
  double maxTime = modulatorToEdit->getEndTime();
  breakpointEditor->setMaximumRangeMaxX(maxTime);
  breakpointEditor->setCurrentRangeMaxX(maxTime);
  breakpointZoomer->zoomToAllX();
}

void BreakpointModulatorEditor::autoAdjustPlotRangeY()
{
  double minLevel = modulatorToEdit->getMinLevel();
  double maxLevel = modulatorToEdit->getMaxLevel();
  breakpointEditor->setMaximumRangeY(minLevel, maxLevel);
  breakpointEditor->setCurrentRangeY(minLevel, maxLevel);
  breakpointZoomer->zoomToAllY();
}

double BreakpointModulatorEditor::gridIntervalFromIndex(int index)
{
  switch(index)
  {
  case  0: return 1.0/2.0;
  case  1: return 1.0/4.0;
  case  2: return 1.0/8.0;
  case  3: return 0.1;
  case  4: return 1.0/16.0;
  case  5: return 1.0/32.0;
  case  6: return 1.0/64.0;
  case  7: return 0.01;
  case  8: return 1.0/128.0;
  default: return 0.1;
  }
}

int BreakpointModulatorEditor::indexFromGridInterval(double interval)
{
  if     ( interval == 0.5       )  return 1;
  else if( interval == 0.25      )  return 2;
  else if( interval == 1.0/8.0   )  return 3;
  else if( interval == 0.1       )  return 4;
  else if( interval == 1.0/16.0  )  return 5;
  else if( interval == 1.0/32.0  )  return 6;
  else if( interval == 1.0/64.0  )  return 7;
  else if( interval == 0.01      )  return 8;
  else if( interval == 1.0/128.0 )  return 9;
  else                              return 4;
}

double BreakpointModulatorEditor::timeIntervalFromIndex(int index)
{
  switch(index)
  {
  case  2: return 1.0/16.0;
  case  3: return 1.0/12.0;
  case  4: return 1.0/8.0;
  case  5: return 1.0/6.0;
  case  6: return 1.0/4.0;
  case  7: return 1.0/3.0;
  case  8: return 1.0/2.0;
  case  9: return 1.0;
  case 10: return 1.5;
  case 11: return 2.0;
  case 12: return 3.0;
  case 13: return 4.0;
  case 14: return 6.0;
  case 15: return 8.0;
  case 16: return 12.0;
  case 17: return 16.0;
  default: return 1.0;
  }
}

int BreakpointModulatorEditor::indexFromTimeInterval(double interval)
{
  if     ( interval == 1.0/16.0  )  return  2;
  else if( interval == 1.0/12.0  )  return  3;
  else if( interval == 1.0/8.0   )  return  4;
  else if( interval == 1.0/6.0   )  return  5;
  else if( interval == 1.0/4.0   )  return  6;
  else if( interval == 1.0/3.0   )  return  7;
  else if( interval == 1.0/2.0   )  return  8;
  else if( interval == 1.0       )  return  9;
  else if( interval == 1.5       )  return 10;
  else if( interval == 2.0       )  return 11;
  else if( interval == 3.0       )  return 12;
  else if( interval == 4.0       )  return 13;
  else if( interval == 6.0       )  return 14;
  else if( interval == 8.0       )  return 15;
  else if( interval == 12.0      )  return 16;
  else if( interval == 16.0      )  return 17;
  else                              return  1; // no predefined interval
}

//=================================================================================================
// class BreakpointModulatorEditorCompact:

// construction/destruction:

BreakpointModulatorEditorCompact::BreakpointModulatorEditorCompact(CriticalSection *newPlugInLock,                                                                  
  BreakpointModulatorAudioModule* newModulatorToEdit) 
: AudioModuleEditor(newModulatorToEdit)
//: AudioModuleEditor(newPlugInLock, newModulatorToEdit)
{
  setLinkPosition(AudioModuleEditor::INVISIBLE);

  jassert(newModulatorToEdit != NULL ); // you must pass a valid module here
  modulatorModuleToEdit = newModulatorToEdit;
  modulatorToEdit       = newModulatorToEdit->wrappedBreakpointModulator;

  addWidget( editButton = new RButton("Edit") );
  editButton->addRButtonListener(this);
  editButton->setDescription("Open/close context menu with more options");
  editButton->setClickingTogglesState(true);

  numSamplesInPlot = 0;
  xValues          = NULL;
  yValues          = NULL;
  plot = new CurveFamilyPlotOld("Plot");
  plot->setDescription("Envelope");
  plot->setAxisLabels("", "");
  plot->setVerticalCoarseGrid(1.0, false);
  plot->setHorizontalCoarseGrid(1.0, false);
  plot->setAxisValuesPositionX(CoordinateSystemOld::INVISIBLE);
  plot->setAxisValuesPositionY(CoordinateSystemOld::INVISIBLE);
  addPlot(plot);

  popUpEditor = new BreakpointModulatorEditor(lock, modulatorModuleToEdit);
  popUpEditor->addChangeListener(this);
  popUpEditor->setAlwaysOnTop(true);
  popUpEditor->setOpaque(true);
  popUpEditor->closeButton->addRButtonListener(this);
  popUpEditor->closeButton->setVisible(true);
  addChildColourSchemeComponent(popUpEditor, false, false);
  popUpEditorX = -600;
  popUpEditorW =  600+48;
  popUpEditorY =  16;
  popUpEditorH =  300;

  //popUpEditor->setSize(200, 200);

  layout = STANDARD;
  layout = COMPACT;   // for test only
    
  updateWidgetsAccordingToState();
}

BreakpointModulatorEditorCompact::~BreakpointModulatorEditorCompact()
{
  delete popUpEditor; // this is not a child component
  delete xValues;
  delete yValues;
}

//-------------------------------------------------------------------------------------------------
// setup:

void BreakpointModulatorEditorCompact::setLayout(int newLayout)
{
  layout = newLayout;
  if( layout == COMPACT )
    setHeadlineStyle(NO_HEADLINE);
  else
    setHeadlineStyle(SUB_HEADLINE);
}

void BreakpointModulatorEditorCompact::setPopUpEditorBounds(int x, int y, int w, int h)
{
  popUpEditorX = x;
  popUpEditorW = w;
  popUpEditorY = y;
  popUpEditorH = h;
}

void BreakpointModulatorEditorCompact::setHeadlineText(const juce::String& newHeadlineText)
{
  AudioModuleEditor::setHeadlineText(newHeadlineText);
  popUpEditor->setHeadlineText(newHeadlineText);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void BreakpointModulatorEditorCompact::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( modulatorToEdit == NULL )
    return;

  if( buttonThatWasClicked == editButton )
  {
    if( editButton->getToggleState() == true )
    {
      int x = editButton->getScreenX();
      int y = editButton->getScreenY();
      popUpEditor->setBounds(x+popUpEditorX, y+popUpEditorY, popUpEditorW, popUpEditorH);
      popUpEditor->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
      popUpEditor->setVisible(true);
      popUpEditor->toFront(true);
    }
    else
      popUpEditor->setVisible(false);
  }
  else if( buttonThatWasClicked == popUpEditor->closeButton )
    editButton->setToggleState(false, true);
}

void BreakpointModulatorEditorCompact::changeListenerCallback(
  ChangeBroadcaster *objectThatHasChanged)
{
  /*
  updateWidgetsAccordingToState();
  sendChangeMessage();
  */
  if( modulatorToEdit == NULL )
    return;
  if( objectThatHasChanged == popUpEditor )
  {
    moduleToEdit->markStateAsDirty();  
    updatePlot();
  }
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void BreakpointModulatorEditorCompact::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  popUpEditor->updateWidgetsAccordingToState();
  updatePlot();
}

void BreakpointModulatorEditorCompact::updatePlot()
{
  if( modulatorModuleToEdit == NULL || modulatorToEdit == NULL )
    return;
  RAPT::rsBreakpointModulator<double> tmpModulator;
  tmpModulator.copyDataFrom(*modulatorToEdit);
  tmpModulator.fillBufferWithEnvelope(yValues, numSamplesInPlot, false);
  double xMin    = tmpModulator.getStartTime();
  double xMax    = tmpModulator.getEndTime();
  double yMin    = tmpModulator.getMinLevel();
  double yMax    = tmpModulator.getMaxLevel();
  double xRange  = xMax-xMin;
  double yRange  = yMax-yMin;
  double xMargin = 0.025;
  double yMargin = 0.1;
  plot->setMaximumRange(xMin-xMargin*xRange, xMax+xMargin*xRange, 
    yMin-yMargin*yRange, yMax+yMargin*yRange);
  plot->setCurrentRange(plot->getMaximumRange());
  for(int n=0; n<numSamplesInPlot; n++)
    xValues[n] = xMin + n*(xMax-xMin)/numSamplesInPlot;
  plot->setCurveValues(numSamplesInPlot, xValues, yValues);
}

void BreakpointModulatorEditorCompact::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  if( layout == COMPACT )
  {
    editButton->setBounds(w-32-4, y+4, 32, 16);
    y = editButton->getBottom()+4;
    plot->setBounds(x, y, w, h-y);
  }
  else
  {
    // something to do here....
  }

  numSamplesInPlot = plot->getWidth();
  if( xValues != NULL ) { delete[] xValues; xValues = NULL; }
  if( yValues != NULL ) { delete[] yValues; yValues = NULL; }
  xValues = new double[numSamplesInPlot];
  yValues = new double[numSamplesInPlot];
  fillWithIndex(xValues, numSamplesInPlot);
  fillWithZeros(yValues, numSamplesInPlot);
  plot->setMaximumRange(0.0, numSamplesInPlot, -1.1, 1.1);
  plot->setCurrentRange(plot->getMaximumRange());
  updatePlot();
}
