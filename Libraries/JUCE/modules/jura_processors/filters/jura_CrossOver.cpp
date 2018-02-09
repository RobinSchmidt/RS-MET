
// construction/destruction:

CrossOverAudioModule::CrossOverAudioModule(CriticalSection* newPlugInLock, 
  rosic::rsCrossOver4Way* crossOverToWrap) : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*lock);
  //jassert(crossOverToWrap != NULL); // you must pass a valid rosic-object to the constructor

  if(crossOverToWrap != nullptr)
    wrappedCrossOver   = crossOverToWrap;
  else
  {
    wrappedCrossOver = new rosic::rsCrossOver4Way;
    wrappedCrossOverIsOwned = true;
  }
  wantsTempoSyncInfo = false;  // mmmhh...maybe better set it to false in the baseclass and to true in subclasses that need it
  setModuleTypeName("CrossOver");
  createStaticParameters();
}

CrossOverAudioModule::~CrossOverAudioModule()
{
  if(wrappedCrossOverIsOwned)
    delete wrappedCrossOver;
}

AudioModuleEditor* CrossOverAudioModule::createEditor()
{
  return new CrossOverModuleEditor(lock, this);
}

// internal functions:

void CrossOverAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  AutomatableParameter* p;

  p = new AutomatableParameter(lock, "Mono", 0.0, 1.0,  1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedCrossOver, &rsCrossOver4Way::setMonoMode);

  p = new AutomatableParameter(lock, "OnOff_1_1", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setBandActive_0_0);

  p = new AutomatableParameter(lock, "Frequency_1_1", 20.0, 20000.0,  0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setCrossoverFrequency_0_0);

  p = new AutomatableParameter(lock, "Slope_1_1", 12.0, 96.0, 12.0, 24.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(this, &CrossOverAudioModule::setSlope_0_0);


  p = new AutomatableParameter(lock, "OnOff_2_1", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setBandActive_1_0);

  p = new AutomatableParameter(lock, "Frequency_2_1", 20.0, 20000.0,  0.0, 250.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setCrossoverFrequency_1_0);

  p = new AutomatableParameter(lock, "Slope_2_1", 12.0, 96.0, 12.0, 24.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(this, &CrossOverAudioModule::setSlope_1_0);


  p = new AutomatableParameter(lock, "OnOff_2_2", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setBandActive_1_1);

  p = new AutomatableParameter(lock, "Frequency_2_2", 20.0, 20000.0,  0.0, 6000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback(this, &CrossOverAudioModule::setCrossoverFrequency_1_1);

  p = new AutomatableParameter(lock, "Slope_2_2", 12.0, 96.0, 12.0, 24.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(this, &CrossOverAudioModule::setSlope_1_1);

  for(int i=0; i < (int) parameters.size(); i++ )
    parameters[i]->resetToDefaultValue(true, true);
}

//=================================================================================================

CrossOverPlotEditor::CrossOverPlotEditor(CriticalSection *newPlugInLock, CrossOverAudioModule* newCrossOverModuleToEdit)
{
  plugInLock = newPlugInLock;
  ScopedLock scopedLock(*plugInLock);


  setDescription("Vertical lines: adjust frequency, triangles: turn on/off");

  ParameterObserver::setIsGuiElement(true);

  crossOverModuleToEdit = newCrossOverModuleToEdit;  

  //crossOverModuleToEdit = NULL;

  // set up appearance:
  setAutoReRendering(false);
  setMaximumRange(15.625, 32000.0, -90.0, 18.0);
  setCurrentRange(15.625, 32000.0, -90.0, 18.0);
  setHorizontalCoarseGrid(12.0, true);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( true);
  setVerticalFineGridVisible(   false);
  rsPlot::setAxisValuesPositionX(rsPlotSettings::ABOVE_AXIS);
  rsPlot::setAxisValuesPositionY(rsPlotSettings::RIGHT_TO_AXIS);
  showPositionAsDescription = true;
  setAutoReRendering(true);

  //colourScheme.setCurveColouringStrategy(PlotColourScheme::UNIFORM); ...is managed in outlying editor

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  // initialize the pointers to Parameter objects with NULL:
  onOff21Parameter = NULL; 
  onOff22Parameter = NULL;
  freq11Parameter  = NULL; 
  freq21Parameter  = NULL; 
  freq22Parameter  = NULL; 
  slope11Parameter = NULL; 
  slope21Parameter = NULL; 
  slope22Parameter = NULL;

  // this stuff will be (re-) assigned in resized():
  numBins          = 0;
  frequencies      = NULL;
  magnitudes1      = NULL;
  magnitudes2      = NULL;
  magnitudes3      = NULL;
  magnitudes4      = NULL;
  allMagnitudes    = new double*[4];
  allMagnitudes[0] = magnitudes1;
  allMagnitudes[1] = magnitudes2;
  allMagnitudes[2] = magnitudes3;
  allMagnitudes[3] = magnitudes4;

  currentlyDraggedHandle = NONE;
  selectedLevel          = -1;
  selectedIndex          = -1;

  // activate automation for this ParameterObserver:
  ParameterObserver::setLocalAutomationSwitch(true);
}

CrossOverPlotEditor::~CrossOverPlotEditor(void)
{
  ScopedLock scopedLock(*plugInLock);

  // remove ourselves as listeners from the Parameter objects, such that they do not try to notify a nonexistent listener:
  ParameterObserver::setLocalAutomationSwitch(false);

  if( onOff21Parameter != NULL ) onOff21Parameter->deRegisterParameterObserver(this);
  if( onOff22Parameter != NULL ) onOff22Parameter->deRegisterParameterObserver(this);
  if( freq11Parameter  != NULL ) freq11Parameter->deRegisterParameterObserver( this);
  if( freq21Parameter  != NULL ) freq21Parameter->deRegisterParameterObserver( this);
  if( freq22Parameter  != NULL ) freq22Parameter->deRegisterParameterObserver( this);
  if( slope11Parameter != NULL ) slope11Parameter->deRegisterParameterObserver(this);
  if( slope21Parameter != NULL ) slope21Parameter->deRegisterParameterObserver(this);
  if( slope22Parameter != NULL ) slope22Parameter->deRegisterParameterObserver(this);

  deleteAndZero(frequencies);
  deleteAndZero(magnitudes1);
  deleteAndZero(magnitudes2);
  deleteAndZero(magnitudes3);
  deleteAndZero(magnitudes4);
  deleteAndZero(allMagnitudes);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

/*
void CrossOverPlotEditor::setCrossOverToEdit(rosic::CrossOver4Way* newCrossOverToEdit)
{
ScopedLock scopedLock(*plugInLock);
crossOverModuleToEdit = newCrossOverToEdit;
}
*/

void CrossOverPlotEditor::assignParameterOnOff(int treeLevel, int indexInLevel, Parameter* parameterToAssign)
{
  ScopedLock scopedLock(*plugInLock);
  if( treeLevel == 1 && indexInLevel == 0 )
  {
    onOff21Parameter = parameterToAssign;
    onOff21Parameter->registerParameterObserver(this);
  }
  else if( treeLevel == 1 && indexInLevel == 1 )
  {
    onOff22Parameter = parameterToAssign;
    onOff22Parameter->registerParameterObserver(this);
  }
}

void CrossOverPlotEditor::assignParameterFreq(int treeLevel, int indexInLevel, Parameter* parameterToAssign)
{
  ScopedLock scopedLock(*plugInLock);
  if( treeLevel == 0 && indexInLevel == 0 )
  {
    freq11Parameter = parameterToAssign;
    freq11Parameter->registerParameterObserver(this);
  }
  else if( treeLevel == 1 && indexInLevel == 0 )
  {
    freq21Parameter = parameterToAssign;
    freq21Parameter->registerParameterObserver(this);
  }
  else if( treeLevel == 1 && indexInLevel == 1 )
  {
    freq22Parameter = parameterToAssign;
    freq22Parameter->registerParameterObserver(this);
  }
}

void CrossOverPlotEditor::assignParameterSlope(int treeLevel, int indexInLevel, Parameter* parameterToAssign)
{
  ScopedLock scopedLock(*plugInLock);
  if( treeLevel == 0 && indexInLevel == 0 )
  {
    slope11Parameter = parameterToAssign;
    slope11Parameter->registerParameterObserver(this);
  }
  else if( treeLevel == 1 && indexInLevel == 0 )
  {
    slope21Parameter = parameterToAssign;
    slope21Parameter->registerParameterObserver(this);
  }
  else if( treeLevel == 1 && indexInLevel == 1 )
  {
    slope22Parameter = parameterToAssign;
    slope22Parameter->registerParameterObserver(this);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

void CrossOverPlotEditor::getSelectedTreeLevelAndIndex(int &treeLevel, int &indexInLevel)
{
  ScopedLock scopedLock(*plugInLock);
  treeLevel    = selectedLevel;
  indexInLevel = selectedIndex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void CrossOverPlotEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);

  updatePlot();

  // re-plot?
  sendChangeMessage();
}

void CrossOverPlotEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  ScopedLock scopedLock(*plugInLock);
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

void CrossOverPlotEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);

  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  setLocalAutomationSwitch(false);

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);

  // switch the wantsAutomationNotification flag on again:  
  setLocalAutomationSwitch(true);
}

void CrossOverPlotEditor::mouseMove(const MouseEvent &e)
{
  ScopedLock scopedLock(*plugInLock);

  int handleIndex = getDragHandleAt(e.x, e.y);
  if( handleIndex == FREQUENCY_1_1  )
    setMouseCursor(MouseCursor::LeftRightResizeCursor);
  else if( handleIndex != NONE )
  {
    if( e.y >= 16 )
      setMouseCursor(MouseCursor::LeftRightResizeCursor);
    else
      setMouseCursor(MouseCursor::PointingHandCursor);
  }
  else
    setMouseCursor(MouseCursor::NormalCursor);
}

void CrossOverPlotEditor::mouseDown(const MouseEvent &e)
{
  ScopedLock scopedLock(*plugInLock);

  if( crossOverModuleToEdit == NULL )
    return;

  currentlyDraggedHandle = getDragHandleAt(e.x, e.y);

  if( currentlyDraggedHandle == FREQUENCY_1_1  )
  {
    selectedLevel = 0;
    selectedIndex = 0;
  }
  else if( currentlyDraggedHandle == FREQUENCY_2_1  )
  {
    selectedLevel = 1;
    selectedIndex = 0;
    if( e.y < 16 )
    {
      onOff21Parameter->setValue( onOff21Parameter->getValue() == 0.0, true, true );
      currentlyDraggedHandle = NONE;
    }
  }
  else if( currentlyDraggedHandle == FREQUENCY_2_2  )
  {
    selectedLevel = 1;
    selectedIndex = 1;
    if( e.y < 16 )
    {
      onOff22Parameter->setValue( onOff22Parameter->getValue() == 0.0, true, true );
      currentlyDraggedHandle = NONE;
    }
  }
  else
  {
    selectedLevel = -1;
    selectedIndex = -1;
  }

  sendChangeMessage();
}

void CrossOverPlotEditor::mouseDrag(const juce::MouseEvent &e)
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit == NULL )
    return;

  // get the position of the event in components coordinates:
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();

  setupFilterAccordingToMousePosition(x, y);
  sendChangeMessage();
}

void CrossOverPlotEditor::mouseUp(const juce::MouseEvent &e)
{
  ScopedLock scopedLock(*plugInLock);
  currentlyDraggedHandle = NONE;
}

int CrossOverPlotEditor::getDragHandleAt(int x, int y)
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit == NULL )
    return NONE;

  // get the x-positions of the vertical lines (in pixels):
  double x11    = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(0, 0);
  double yDummy = 0.0;
  transformToComponentsCoordinates(x11, yDummy);

  double x21 = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 0);
  yDummy = 0.0;
  transformToComponentsCoordinates(x21, yDummy);

  double x22 = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 1);
  yDummy = 0.0;
  transformToComponentsCoordinates(x22, yDummy);

  // ckeck, if the passed x,y coordinates are near one of these vertical lines:
  double margin = 4.0;
  if(      fabs( (double)x - x11 ) < margin )  return FREQUENCY_1_1;
  else if( fabs( (double)x - x21 ) < margin )  return FREQUENCY_2_1;
  else if( fabs( (double)x - x22 ) < margin )  return FREQUENCY_2_2;
  else                                         return NONE;
}

void CrossOverPlotEditor::setupFilterAccordingToMousePosition(double mouseX, double mouseY)
{
  ScopedLock scopedLock(*plugInLock);

  double freq  = restrictDragHandleX(mouseX, currentlyDraggedHandle);  
  double dummy = mouseY;

  transformFromComponentsCoordinates(freq, dummy);

  if( currentlyDraggedHandle == FREQUENCY_1_1 && freq11Parameter != NULL ) 
    freq11Parameter->setValue(freq, true, true);
  if( currentlyDraggedHandle == FREQUENCY_2_1 && freq21Parameter != NULL ) 
    freq21Parameter->setValue(freq, true, true);
  if( currentlyDraggedHandle == FREQUENCY_2_2 && freq22Parameter != NULL ) 
    freq22Parameter->setValue(freq, true, true);
}

void CrossOverPlotEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  SpectrumDisplayOld::resized();

  // (re) allocate and fill the arrays for the magnitude plot
  numBins = getWidth();
  delete[] frequencies;
  delete[] magnitudes1;
  delete[] magnitudes2;
  delete[] magnitudes3;
  delete[] magnitudes4;
  frequencies = new double[numBins];
  magnitudes1 = new double[numBins];
  magnitudes2 = new double[numBins];
  magnitudes3 = new double[numBins];
  magnitudes4 = new double[numBins];
  allMagnitudes[0] = magnitudes1;
  allMagnitudes[1] = magnitudes2;
  allMagnitudes[2] = magnitudes3;
  allMagnitudes[3] = magnitudes4;

  getDisplayedFrequencies(frequencies, numBins);
  if( crossOverModuleToEdit != NULL )
  {
    crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes1, numBins, 0, true);
    crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes2, numBins, 1, true);
    crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes3, numBins, 2, true);
    crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes4, numBins, 3, true);
  }
  else
  {
    int k;
    for(k=0; k<numBins; k++)
    {
      frequencies[k] = 15.625;
      magnitudes1[k] = 0.0;
      magnitudes2[k] = 0.0;
      magnitudes3[k] = 0.0;
      magnitudes4[k] = 0.0;
    }
  }
  setSpectra(numBins, 4, frequencies, allMagnitudes);
}

void CrossOverPlotEditor::updatePlot()
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit == NULL || numBins == 0 )
    return;

  crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes1, numBins, 0, true);
  crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes2, numBins, 1, true);
  crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes3, numBins, 2, true);
  crossOverModuleToEdit->wrappedCrossOver->getMagnitudeResponse(frequencies, magnitudes4, numBins, 3, true);
  setSpectra(numBins, 4, frequencies, allMagnitudes);
}

void CrossOverPlotEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage, XmlElement *targetSVG)
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit == NULL )
    return;

  rsDataPlot::plotCurveFamily(g, targetImage, targetSVG);
  colourizeBackground(g, targetImage);
  g.setColour(plotColourScheme.getCurveColour(0));
  drawSelectionIndicator(g, targetImage);

  drawVerticalLineAtFrequency(g, targetImage, crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(0, 0), 0, 
    crossOverModuleToEdit->wrappedCrossOver->isBandActive(0, 0));
  drawVerticalLineAtFrequency(g, targetImage, crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 0), 1, 
    crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 0));
  drawVerticalLineAtFrequency(g, targetImage, crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 1), 1, 
    crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 1));

  //drawTriangleSwitchAtFrequency(g, targetImage, crossOverModuleToEdit->getCrossoverFrequency(0, 0), 0, crossOverModuleToEdit->isBandActive(0, 0));
  drawTriangleSwitchAtFrequency(g, targetImage, crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 0), 1, 
    crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 0));
  drawTriangleSwitchAtFrequency(g, targetImage, crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 1), 1, 
    crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 1));
}

void CrossOverPlotEditor::colourizeBackground(Graphics &g, juce::Image *targetImage)
{
  ScopedLock scopedLock(*plugInLock);

  int xL = 0;
  int xR = (int) getDragHandleX(FREQUENCY_2_1);
  float alpha = 0.125;

  juce::Rectangle<int> r(xL, 0, xR, getHeight());
  g.setColour(getLowBandColour().withMultipliedAlpha(alpha));
  g.fillRect(r);

  xL = xR;
  xR = (int) getDragHandleX(FREQUENCY_1_1);

  r = juce::Rectangle<int>(xL, 0, xR-xL, getHeight());
  g.setColour(getLowMidBandColour().withMultipliedAlpha(alpha));
  g.fillRect(r);

  xL = xR;
  xR = (int) getDragHandleX(FREQUENCY_2_2);

  r = juce::Rectangle<int>(xL, 0, xR-xL, getHeight());
  g.setColour(getHighMidBandColour().withMultipliedAlpha(alpha));
  g.fillRect(r);

  xL = xR;
  xR = getWidth();

  r = juce::Rectangle<int>(xL, 0, xR-xL, getHeight());
  g.setColour(getHighBandColour().withMultipliedAlpha(alpha));
  g.fillRect(r);

}

void CrossOverPlotEditor::drawSelectionIndicator(Graphics &g, juce::Image *targetImage)
{
  ScopedLock scopedLock(*plugInLock);

  if( selectedLevel == -1 || selectedIndex == -1 )
    return;

  double x = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(selectedLevel, selectedIndex);
  double y = 0.0;
  transformToComponentsCoordinates(x, y);
  y = 2.f;

  //Colour c = g.getCurrentColour();
  Colour c = plotColourScheme.getCurveColour(0);
  g.setColour(c.withMultipliedAlpha(0.5f));
  g.drawLine((float) x, (float) y, (float) x, (float) getHeight(), 6.f);
  g.setColour(c);
}

void CrossOverPlotEditor::drawVerticalLineAtFrequency(Graphics &g, juce::Image* targetImage, double frequency, int treeLevel, bool active)
{
  ScopedLock scopedLock(*plugInLock);

  double x = frequency;
  double y = 0.0; 
  transformToComponentsCoordinates(x, y);

  if( treeLevel == 0 )
    y = 0.0;
  else if( treeLevel == 1 )
    y = 16.0;

  if( active == false )
  {
    const float dashLengths[2] = {4.f, 6.f};
    g.drawDashedLine(Line<float>( (float) x, (float) y, (float) x, (float) getHeight()), 
      dashLengths, 2, 2.f); 
    /*g.drawDashedLine((float) x, (float) y, (float) x, (float) getHeight(), dashLengths, 2, 2.f); */
  }
  else
    g.drawLine(      (float) x, (float) (y-2), (float) x, (float) getHeight(), 2.f);
}

void CrossOverPlotEditor::drawTriangleSwitchAtFrequency(Graphics &g, juce::Image* targetImage, double frequency, int treeLevel, bool active)
{
  ScopedLock scopedLock(*plugInLock);

  double x = frequency;
  double y = 0.0; 
  transformToComponentsCoordinates(x, y);

  if( treeLevel == 0 )
    y = 0.0;
  else if( treeLevel == 1 )
    y = 2.0;

  drawTriangle(g, (float) (x-5), (float) y, (float) (x+5), (float) y, (float) x, 16.f, active);
}

double CrossOverPlotEditor::restrictDragHandleX(double x, int dragHandleIndex)
{
  ScopedLock scopedLock(*plugInLock);

  double xMin   = 0.0;
  double xMax   = getWidth();
  double margin = 5.0;

  if( dragHandleIndex == FREQUENCY_1_1 ) 
  {
    xMin = getDragHandleX(FREQUENCY_2_1) + margin;
    xMax = getDragHandleX(FREQUENCY_2_2) - margin;
  }
  else if( dragHandleIndex == FREQUENCY_2_1 ) 
    xMax = getDragHandleX(FREQUENCY_1_1) - margin;
  else if( dragHandleIndex == FREQUENCY_2_2 ) 
    xMin = getDragHandleX(FREQUENCY_1_1) + margin;
  else
    return x;

  return jlimit(xMin, xMax, x);
}

double CrossOverPlotEditor::getDragHandleX(int dragHandleIndex)
{
  ScopedLock scopedLock(*plugInLock);

  double y = 0.0; // dummy
  double x;
  switch( dragHandleIndex )
  {
  case FREQUENCY_1_1: x = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(0, 0); break;
  case FREQUENCY_2_1: x = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 0); break;
  case FREQUENCY_2_2: x = crossOverModuleToEdit->wrappedCrossOver->getCrossoverFrequency(1, 1); break;
  }
  transformToComponentsCoordinates(x, y);
  return x;
}

Colour CrossOverPlotEditor::getLowBandColour()
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 0) )
    return Colours::red;
  else
    //return Colour(255, 100, 0);
    return Colours::orange;
}

Colour CrossOverPlotEditor::getLowMidBandColour()
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 0) )
    return Colours::yellow;
  else
    return Colours::orange;
}

Colour CrossOverPlotEditor::getMidBandColour()
{
  return Colours::yellowgreen;
}

Colour CrossOverPlotEditor::getHighMidBandColour()
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 1) )
    return Colours::green;
  else
    return Colours::cyan.darker(1.0);
}

Colour CrossOverPlotEditor::getHighBandColour()
{
  ScopedLock scopedLock(*plugInLock);
  if( crossOverModuleToEdit->wrappedCrossOver->isBandActive(1, 1) )
    return Colours::blue;
  else
    return Colours::cyan.darker(1.0);
}

//=================================================================================================

CrossOverModuleEditor::CrossOverModuleEditor(CriticalSection *newPlugInLock, CrossOverAudioModule* newCrossOverAudioModule) 
  : AudioModuleEditor(newCrossOverAudioModule)
{
  ScopedLock scopedLock(*lock);

  //isTopLevelEditor = true;  // ?

  // assign the pointer to the rosic::CrossOver object to be used as aduio engine:
  jassert(newCrossOverAudioModule != NULL ); // you must pass a valid module here
  crossOverModuleToEdit = newCrossOverAudioModule;

  addWidget( frequency11Slider = new RSlider("Frequency11Slider") );
  frequency11Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Frequency_1_1") );
  frequency11Slider->setSliderName(juce::String("Frequency"));
  frequency11Slider->setDescription(juce::String("Crossover frequency"));
  frequency11Slider->setDescriptionField(infoField);
  frequency11Slider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( slope11Slider = new RSlider("Slope11Slider") );
  slope11Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Slope_1_1") );
  slope11Slider->setSliderName(juce::String("Slope"));
  slope11Slider->setDescription(juce::String("Slope for the separation filters"));
  slope11Slider->setDescriptionField(infoField);
  slope11Slider->setStringConversionFunction(&decibelsPerOctaveToString);


  addWidget( frequency21Slider = new RSlider("Frequency21Slider") );
  frequency21Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Frequency_2_1") );
  frequency21Slider->setSliderName(juce::String("Frequency"));
  frequency21Slider->setDescription(juce::String("Crossover frequency"));
  frequency21Slider->setDescriptionField(infoField);
  frequency21Slider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( slope21Slider = new RSlider("Slope21Slider") );
  slope21Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Slope_2_1") );
  slope21Slider->setSliderName(juce::String("Slope"));
  slope21Slider->setDescription(juce::String("Slope for the separation filters"));
  slope21Slider->setDescriptionField(infoField);
  slope21Slider->setStringConversionFunction(&decibelsPerOctaveToString);


  addWidget( frequency22Slider = new RSlider("Frequency22Slider") );
  frequency22Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Frequency_2_2") );
  frequency22Slider->setSliderName(juce::String("Frequency"));
  frequency22Slider->setDescription(juce::String("Crossover frequency"));
  frequency22Slider->setDescriptionField(infoField);
  frequency22Slider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( slope22Slider = new RSlider("Slope22Slider") );
  slope22Slider->assignParameter( crossOverModuleToEdit->getParameterByName("Slope_2_2") );
  slope22Slider->setSliderName(juce::String("Slope"));
  slope22Slider->setDescription(juce::String("Slope for the separation filters"));
  slope22Slider->setDescriptionField(infoField);
  slope22Slider->setStringConversionFunction(&decibelsPerOctaveToString);


  addWidget( monoButton = new RButton(juce::String("Mono")) );
  monoButton->assignParameter( crossOverModuleToEdit->getParameterByName("Mono") );
  monoButton->setDescription(juce::String("Switch into mono-mode (saves CPU)"));
  monoButton->setDescriptionField(infoField);
  monoButton->setClickingTogglesState(true);

  plotColourScheme.setCurveColouringStrategy(PlotColourScheme::UNIFORM);
  //plotEditor = new CrossOverPlotEditor(juce::String(T("SpectrumEditor")));
  plotEditor = new CrossOverPlotEditor(lock, crossOverModuleToEdit);
  //plotEditor->setDescription(juce::String(T("Drag vertical line to adjust crossover frequency")));
  plotEditor->setDescriptionField(infoField);
  plotEditor->addChangeListener(this);
  plotEditor->assignParameterOnOff( 1, 0, moduleToEdit->getParameterByName("OnOff_2_1"));
  plotEditor->assignParameterOnOff( 1, 1, moduleToEdit->getParameterByName("OnOff_2_2"));
  plotEditor->assignParameterFreq(  0, 0, moduleToEdit->getParameterByName("Frequency_1_1"));
  plotEditor->assignParameterFreq(  1, 0, moduleToEdit->getParameterByName("Frequency_2_1"));
  plotEditor->assignParameterFreq(  1, 1, moduleToEdit->getParameterByName("Frequency_2_2"));
  plotEditor->assignParameterSlope( 0, 0, moduleToEdit->getParameterByName("Slope_1_1"));
  plotEditor->assignParameterSlope( 1, 0, moduleToEdit->getParameterByName("Slope_2_1"));
  plotEditor->assignParameterSlope( 1, 1, moduleToEdit->getParameterByName("Slope_2_2"));
  //plotEditor->setCrossOverToEdit(      newCrossOverAudioModule->wrappedCrossOver);
  addPlot( plotEditor );

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(480, 240);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void CrossOverModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  ScopedLock scopedLock(*lock);
  if( objectThatHasChanged == plotEditor )
    updateWidgetVisibility();
  else
    AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
}

void CrossOverModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::updateWidgetsAccordingToState();
  updateWidgetVisibility();
}

void CrossOverModuleEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 4;
  int y = infoField->getY()-20;;
  int w = 160;
  int h = 16;

  frequency11Slider->setBounds(x, y, w, h); 
  frequency21Slider->setBounds(x, y, w, h); 
  frequency22Slider->setBounds(x, y, w, h); 

  x = frequency11Slider->getRight()+8;
  slope11Slider->setBounds(x, y, w, h); 
  slope21Slider->setBounds(x, y, w, h); 
  slope22Slider->setBounds(x, y, w, h); 

  w = 40;
  x = getWidth()-w-4;
  monoButton->setBounds(x, y, w, h);

  x = 0;
  y = getPresetSectionBottom();
  w = getWidth();
  h = frequency11Slider->getY() - y;
  plotEditor->setBounds(x+4, y+4, w-8, h-8);
}

void CrossOverModuleEditor::updateWidgetVisibility()
{
  ScopedLock scopedLock(*lock);

  if( crossOverModuleToEdit == NULL )
    return;
  if( crossOverModuleToEdit->wrappedCrossOver == NULL )
    return;

  frequency11Slider->setVisible(false);
  frequency21Slider->setVisible(false);
  frequency22Slider->setVisible(false);
  slope11Slider->setVisible(false);
  slope21Slider->setVisible(false);
  slope22Slider->setVisible(false);

  int treeLevel, index;
  plotEditor->getSelectedTreeLevelAndIndex(treeLevel, index);
  if( treeLevel == 0 && index == 0 )
  {
    frequency11Slider->setVisible(true);
    slope11Slider->setVisible(true);
  }
  else if( treeLevel == 1 && index == 0 )
  {
    frequency21Slider->setVisible(true);
    slope21Slider->setVisible(true);
  }
  else if( treeLevel == 1 && index == 1 )
  {
    frequency22Slider->setVisible(true);
    slope22Slider->setVisible(true);
  }

  plotEditor->updatePlot();
}