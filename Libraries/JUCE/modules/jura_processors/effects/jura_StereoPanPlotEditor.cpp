

StereoPanPlotEditor::StereoPanPlotEditor(const juce::String& name) 
: CurveFamilyPlotOld(name)
{
  setDescription("Shows the curves for the pan-law - red: L->L, blue: R->R, green: R->L, brown: L->R");


  ParameterObserver::setIsGuiElement(true);
  stereoPanToEdit = NULL;

  // set up the plot range:
  setMaximumRange(-1.2, 1.2, -0.2, 2.2);
  setCurrentRange(-1.2, 1.2, -0.2, 2.2);
  setHorizontalCoarseGrid(1.0, true);
  setVerticalCoarseGrid(  1.0, true);
  //setHorizontalFineGrid(  0.1, false);
  setAxisValuesPositionX(ABOVE_AXIS);
  setAxisLabels(juce::String(""), juce::String(""));

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  // initialize the two pointers to Parameter objects with NULL:
  panParameter     = NULL;
  panLawParameter = NULL;

  // this stuff will be (re-) assigned in resized():
  numValues   = 0;
  p           = NULL;
  gLL         = NULL;
  gRL         = NULL;
  gLR         = NULL;
  gRR         = NULL;
  allGains    = new double*[4];
  allGains[0] = gLL;
  allGains[1] = gRR;
  allGains[2] = gRL;
  allGains[3] = gLR;

  // activate automation for this ParameterObserver:
  ParameterObserver::setLocalAutomationSwitch(true);
}

StereoPanPlotEditor::~StereoPanPlotEditor(void)
{
  // remove ourselves as listeners from the Parameter objects, such that they do 
  // not try to notify a nonexistent listener:
  ParameterObserver::setLocalAutomationSwitch(false);
  if( panParameter != NULL )
    panParameter->deRegisterParameterObserver(this);
  if( panLawParameter != NULL )
    panLawParameter->deRegisterParameterObserver(this);
  deleteAndZero(p);
  deleteAndZero(gLL);
  deleteAndZero(gRL);
  deleteAndZero(gLR);
  deleteAndZero(gRR);
  deleteAndZero(allGains);
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void StereoPanPlotEditor::setStereoPanToEdit(rosic::StereoPan* newStereoPanToEdit)
{
  stereoPanToEdit = newStereoPanToEdit;
}

void StereoPanPlotEditor::assignParameterPan(Parameter *parameterToAssign)
{
  panParameter = parameterToAssign;
  if( panParameter != NULL )
    panParameter->registerParameterObserver(this);
}

void StereoPanPlotEditor::assignParameterPanLaw(Parameter *parameterToAssign)
{
  panLawParameter = parameterToAssign;
  if( panLawParameter != NULL )
    panLawParameter->registerParameterObserver(this);
}

void StereoPanPlotEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  updatePlot();
  sendChangeMessage();
}

void StereoPanPlotEditor::parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted)
{

}

void StereoPanPlotEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void StereoPanPlotEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  setLocalAutomationSwitch(false);

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);

  // switch the wantsAutomationNotification flag on again:  
  setLocalAutomationSwitch(true);
}

void StereoPanPlotEditor::mouseDown(const MouseEvent &e)
{
  setupPanAccordingToMousePosition(e.x);
}

void StereoPanPlotEditor::mouseDrag(const juce::MouseEvent &e)
{
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  setupPanAccordingToMousePosition(x);
  sendChangeMessage();
}

void StereoPanPlotEditor::setupPanAccordingToMousePosition(double mouseX)
{
  double p     = mouseX;
  double dummy = 0.0;
  transformFromComponentsCoordinates(p, dummy);
  panParameter->setValue(rosic::clip(p, -1.0, 1.0), true, true);
}

//-------------------------------------------------------------------------------------------------
// drawing:

void StereoPanPlotEditor::resized()
{
  CurveFamilyPlotOld::resized();

  // (re) allocate and fill the arrays for the magnitude plot
  numValues = getWidth();
  deleteAndZero(p);
  deleteAndZero(gLL);
  deleteAndZero(gRL);
  deleteAndZero(gLR);
  deleteAndZero(gRR);
  p   = new double[numValues];
  gLL = new double[numValues];
  gRL = new double[numValues];
  gLR = new double[numValues];
  gRR = new double[numValues];
  allGains[0] = gLL;
  allGains[1] = gRR;
  allGains[2] = gRL;
  allGains[3] = gLR;

  if( stereoPanToEdit != NULL )
    updatePlot();
}

void StereoPanPlotEditor::updatePlot()
{
  if( stereoPanToEdit == NULL )
    return;
  rosic::StereoPan spTmp = *stereoPanToEdit;
  spTmp.setGain(0.0);
  for(int n=0; n<numValues; n++)
  {
    double tmp = (double) n / (double) (numValues-1);
    tmp        = linToLin(tmp, 0.0, 1.0, -1.0, 1.0);
    p[n]       = tmp;
    spTmp.setPanoramaPosition(p[n]);
    gLL[n]     = spTmp.getLeftToLeftGain();
    gRL[n]     = spTmp.getRightToLeftGain();
    gLR[n]     = spTmp.getLeftToRightGain();
    gRR[n]     = spTmp.getRightToRightGain();
  }
  if( stereoPanToEdit->doesPanLawApplyCrossMix() )
  {
    setFunctionFamilyValues(numValues, 4, p, allGains);  
    setDescription("Shows the curves for the pan-law - red: gain left, blue: gain right, green: right to left, brown: left to right");
  }
  else
  {
    setFunctionFamilyValues(numValues, 2, p, allGains);
    setDescription("Shows the curves for the pan-law - red: gain left, blue: gain right");
  }

  int law = stereoPanToEdit->getPanLaw();
  if( law == rosic::StereoPan::LINEAR )
    setCurrentRange(-1.1, 1.1, -0.1, 2.2);
  else if( law == rosic::StereoPan::LINEAR_CLIPPED 
    || law == rosic::StereoPan::SQRT_CROSSMIX_SQUARE_NORMALIZED )
    setCurrentRange(-1.1, 1.1, -0.1, 1.2);
  else if( law == rosic::StereoPan::SINCOS || law == rosic::StereoPan::SQUARE_ROOT 
    || law == rosic::StereoPan::LINEAR_CLIPPED_SQUARE_NORMALIZED 
    || law == rosic::StereoPan::LINEAR_SQUARE_NORMALIZED )
    setCurrentRange(-1.1, 1.1, -0.1, 1.5);
  else
    setCurrentRange(-1.1, 1.1, -0.1, 1.2);
}

void StereoPanPlotEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage, 
                                                  XmlElement *targetSVG)
{
  if( stereoPanToEdit == NULL )
    return;

  CurveFamilyPlotOld::plotCurveFamily(g, targetImage, targetSVG);

  //if( colourScheme.plotColours.size() < 4 )
  //  return;

  //Colour graphColour = colourScheme.curves;
  Colour graphColour = plotColourScheme.getCurveColour(0);  

  rosic::StereoPan spTmp = *stereoPanToEdit;
  spTmp.setGain(0.0);

  // draw the draggable vertical line
  g.setColour(Colour(0xffa0a0a0)); // a somewhat bright grey
  //g.setColour(Colour(0xff505050));
  double x = spTmp.getPanoramaPosition();
  double y = 0.0;
  transformToImageCoordinates(x, y, targetImage);
  g.drawLine((float) x, 0.f, (float) x, (float) getHeight(), 2.0);

  if( spTmp.doesPanLawApplyCrossMix() )
  {
    x = spTmp.getPanoramaPosition();
    y = spTmp.getRightToLeftGain();
    transformToImageCoordinates(x, y, targetImage);
    //g.setColour(colourScheme.plotColours[2]); 
    g.fillRect(float(x-4), float(y-4), 8.f, 8.f);

    x = spTmp.getPanoramaPosition();
    y = spTmp.getLeftToRightGain();
    transformToImageCoordinates(x, y, targetImage);
    //g.setColour(colourScheme.plotColours[3]);
    g.fillRect(float(x-4), float(y-4), 8.f, 8.f);

  }
  x = spTmp.getPanoramaPosition();
  y = spTmp.getLeftToLeftGain();
  transformToImageCoordinates(x, y, targetImage);
  //g.setColour(colourScheme.plotColours[0]);
  g.fillEllipse(float(x-3), float(y-3), 6.f, 6.f);

  x = spTmp.getPanoramaPosition();
  y = spTmp.getRightToRightGain();
  transformToImageCoordinates(x, y, targetImage);
  //g.setColour(colourScheme.plotColours[1]);
  g.fillEllipse(float(x-3), float(y-3), 6.f, 6.f);
}



