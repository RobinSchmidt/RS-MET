
//-------------------------------------------------------------------------------------------------
// construction/destruction:

CoordinateSystem::CoordinateSystem(const String &newDescription) : RWidget(newDescription)
{
  //autoReRenderImage  = false;
  caption.setText(String());
  caption.setFont(Font(16), true);


  scaleX                        =  1.0;
  scaleY                        =  1.0;
  pixelsPerIntervalX            =  100.0;
  pixelsPerIntervalY            =  100.0;
  captionPosition               =  NO_CAPTION;
  axisPositionX                 =  ZERO;
  axisPositionY                 =  ZERO;
  axisLabelPositionX            =  ABOVE_AXIS;
  axisLabelPositionY            =  RIGHT_TO_AXIS;
  axisValuesPositionX           =  BELOW_AXIS;
  axisValuesPositionY           =  LEFT_TO_AXIS;

  axisLabelX                    =  String("x");
  axisLabelY                    =  String("y");

  horizontalCoarseGridIsVisible =  false;
  horizontalFineGridIsVisible	  =  false;
  verticalCoarseGridIsVisible   =  false;
  verticalFineGridIsVisible	    =  false;
  radialCoarseGridIsVisible     =  false;
  radialFineGridIsVisible       =  false;
  angularCoarseGridIsVisible    =  false;
  angularFineGridIsVisible      =  false;

  horizontalCoarseGridInterval  =  1.0;
  horizontalFineGridInterval    =  0.1;
  verticalCoarseGridInterval    =  1.0;
  verticalFineGridInterval      =  0.1;
  radialCoarseGridInterval      =  1.0;
  radialFineGridInterval        =  0.1;
  angularCoarseGridInterval     =  15.0;  // 15 degrees
  angularFineGridInterval       =  5.0;   // 5 degrees

  angleIsInDegrees              =  true;

  logScaledX	                  =  false;
  logBaseX	                    =  2.0;
  logScaledY	                  =  false;
  logBaseY	                    =  2.0;
  logScaledRadius               =  false;
  logBaseRadius                 =  2.0;

  maximumRange.setRangeX(-2.2, 2.2);
  maximumRange.setRangeY(-2.2, 2.2);
  currentRange.setRangeX(-2.2, 2.2);
  currentRange.setRangeY(-2.2, 2.2);

  // initialize the function-pointers for value->string conversion
  stringConversionFunctionX = &valueToString0;
  stringConversionFunctionY = &valueToString0;

  valuePopup    = true;
  useBitmapFont = true;

  caption.setColour(colourScheme.axes);

  /*
  // initialize the component-size and the image-size to 1x1 pixels, without
  // such initializations, a JUCE-breakpoint will be triggered or other screws
  // happen:
  backgroundImage = NULL;
  backgroundImage = new Image(Image::RGB, 1, 1, true);
  setBounds(0, 0, 1, 1);
  updateBackgroundImage();
  autoReRenderImage  = true;
  */

  // use a crosshair-cursor to precisely point to some coordinate for inspection:
  //currentMouseCursor = MouseCursor(MouseCursor::CrosshairCursor);
  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  // the Label for inspecting the current frequency ane level:
  inspectionField = new TextEditor( String("inspectionField") );
  inspectionField->setBounds(4, 4, 120, 40);
  inspectionField->setColour(TextEditor::backgroundColourId, Colours::white.withAlpha(0.7f) );
  inspectionField->setColour(TextEditor::outlineColourId,    Colours::black );
  inspectionField->setCaretVisible(false);
  inspectionField->setScrollbarsShown(false);
  inspectionField->setReadOnly(true);
  inspectionField->setMultiLine(true);
  addChildComponent( inspectionField );
}

CoordinateSystem::~CoordinateSystem()
{
  deleteAllChildren();
  //if( backgroundImage != NULL )
  //  delete backgroundImage;
}

//-------------------------------------------------------------------------------------------------
// component-overrides:

void CoordinateSystem::mouseDown(const MouseEvent &e)
{
  if( e.mods.isRightButtonDown() )
    openRightClickPopupMenu();
}

void CoordinateSystem::mouseEnter(const MouseEvent &e)
{
  RWidget::mouseEnter(e);
  if( valuePopup == true )
  {
    inspectionField->setVisible(true);
    inspectionField->setTopLeftPosition(e.x, e.y);
  }
}

void CoordinateSystem::mouseMove(const MouseEvent &e)
{
  if( valuePopup == true )
  {
    inspectionField->setVisible(true);
    inspectionField->setTopLeftPosition(e.x+8, e.y-16);

    // get frequency and amplitude-level at the current mouse-position:
    double x = e.x;
    double y = e.y;
    fromPixelCoordinates(x, y);

    // pass the strings to the inspectionLabel:
    inspectionField->setText(String("x=") + String(x) + String::charToString('\n') +
      String("y=") + String(y), false);
  }
}

void CoordinateSystem::mouseExit(const MouseEvent &e)
{
  RWidget::mouseExit(e);
  inspectionField->setVisible(false);
}

void CoordinateSystem::mouseDrag(const MouseEvent &e)
{
  Component::mouseDrag(e);
}

void CoordinateSystem::mouseUp(const MouseEvent &e)
{
  Component::mouseUp(e);
}

void CoordinateSystem::mouseDoubleClick(const MouseEvent &e)
{
  Component::mouseDoubleClick(e);
}

//void CoordinateSystem::mouseWheelMove(const MouseEvent &e, float wheelIncrementX, float wheelIncrementY)
void CoordinateSystem::mouseWheelMove(const MouseEvent& ev, const MouseWheelDetails& wheel)
{
  Component::mouseWheelMove(ev, wheel);
  //Component::mouseWheelMove(e, wheelIncrementX, wheelIncrementY);
}


void CoordinateSystem::resized()
{
  ThreadedDrawingPanel::resized();
  /*
  updateScaleFactors();
  if(autoReRenderImage == true)
    updateBackgroundImage();
  */
}


void CoordinateSystem::paint(juce::Graphics &g)
{
  ThreadedDrawingPanel::paint(g);
  /*
  if( backgroundImage != NULL )
  {
    g.drawImage(backgroundImage, 0, 
      0, 
      getWidth(), 
      getHeight(), 
      0, 
      0, 
      backgroundImage->getWidth(), 
      backgroundImage->getHeight(), 
      false);
  }
  else
    g.fillAll(Colours::red);
  */
}

/*
void CoordinateSystem::setDirty(bool shouldSetToDirty)
{
  ThreadedDrawingPanel::setDirty(shouldSetToDirty);
  // maybe more to do?
  //....
}
*/

//-------------------------------------------------------------------------------------------------
// range management:


//-------------------------------------------------------------------------------------------------
// appearance:


void CoordinateSystem::setColourScheme(const PlotColourScheme& newColourScheme)
{
  colourScheme = newColourScheme;
}

void CoordinateSystem::setColourSchemeFromXml(const XmlElement &xml)
{
  //colourScheme.setColourSchemeFromXml(xml);
}

bool CoordinateSystem::changeGraphColour(int index, Colour newColour)
{
  /*
  if( index < 0 || index >= colourScheme.plotColours.size() )
    return false;
  else
  {
    colourScheme.plotColours[index] = newColour;
    return true;
  }
  */

  return false;
}

/*
void CoordinateSystem::setAutoReRendering(bool shouldAutomaticallyReRender)
{
  autoReRenderImage = shouldAutomaticallyReRender;
}
*/

void CoordinateSystem::setCaption(const String &newCaption, int newPosition)
{
  captionPosition = newPosition;
  captionString   = newCaption;
  caption.setText(newCaption);
  setDirty();
}

void CoordinateSystem::setAxisLabels(const String &newLabelX, 
                                     const String &newLabelY, 
                                     int newLabelPositionX,
                                     int newLabelPositionY)
{
  axisLabelX         = newLabelX;
  axisLabelY         = newLabelY;
  axisLabelPositionX = newLabelPositionX;
  axisLabelPositionY = newLabelPositionY;
  setDirty();
}

void CoordinateSystem::setAxisLabelX(const String& newLabelX, 
                                     int newLabelPositionX)
{
  axisLabelX         = newLabelX;
  axisLabelPositionX = newLabelPositionX;
  setDirty();
}

void CoordinateSystem::setAxisLabelY(const String& newLabelY, 
                                     int newLabelPositionY)
{
  axisLabelY             = newLabelY;
  axisLabelPositionY = newLabelPositionY;
  setDirty();
}

void CoordinateSystem::setAxisValuesPositionX(int newValuesPositionX)
{
  if( newValuesPositionX == NO_ANNOTATION ||
    newValuesPositionX == BELOW_AXIS    ||
    newValuesPositionX == ABOVE_AXIS      )
  {
    axisValuesPositionX = newValuesPositionX;
  }
  setDirty();
}

void CoordinateSystem::setAxisValuesPositionY(int newValuesPositionY)
{
  if( newValuesPositionY == NO_ANNOTATION ||
    newValuesPositionY == LEFT_TO_AXIS  ||
    newValuesPositionY == RIGHT_TO_AXIS   )
  {
    axisValuesPositionY = newValuesPositionY;
  }
  setDirty();
}

void CoordinateSystem::setStringConversionFunctionX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  stringConversionFunctionX = newConversionFunctionX;
  setDirty();
}

void CoordinateSystem::setStringConversionFunctionY(
  String (*newConversionFunctionY) (double valueToBeConverted) )
{
  stringConversionFunctionY = newConversionFunctionY;
  setDirty();
}

void CoordinateSystem::setHorizontalCoarseGridVisible(bool shouldBeVisible)
{
  setHorizontalCoarseGrid(horizontalCoarseGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setHorizontalCoarseGridInterval(double newGridInterval)
{
  setHorizontalCoarseGrid(newGridInterval, horizontalCoarseGridIsVisible);
  setDirty();
}

void CoordinateSystem::setHorizontalCoarseGrid(double newGridInterval, 
                                               bool   shouldBeVisible)
{
  // for logarithmic scaling of an axis, we need the grid-intervals to be 
  // strictly greater than unity because the coordinate of a grid-line results
  // from the coordinate of the previous grid-line via multiplication - we 
  // would end up drawing an infinite number of grid-lines at the same 
  // coordinate with a unity-factor and denser and denser grid-lines when 
  // approaching zero with a factor lower than unity.
  if( logScaledY )
  {
    jassert(newGridInterval > 1.00001);
    if( newGridInterval <= 1.00001 )
    {
      horizontalCoarseGridInterval = 2.0;
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  horizontalCoarseGridIsVisible = shouldBeVisible;
  horizontalCoarseGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setHorizontalFineGridVisible(bool shouldBeVisible)
{
  setHorizontalFineGrid(horizontalFineGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setHorizontalFineGridInterval(double newGridInterval)
{
  setHorizontalFineGrid(newGridInterval, horizontalFineGridIsVisible);
  setDirty();
}

void CoordinateSystem::setHorizontalFineGrid(double newGridInterval, 
                                             bool   shouldBeVisible)
{
  if( logScaledY )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      horizontalFineGridInterval = pow(2.0, 1.0/3.0);
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  horizontalFineGridIsVisible = shouldBeVisible;
  horizontalFineGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setVerticalCoarseGridVisible(bool shouldBeVisible)
{
  setVerticalCoarseGrid(verticalCoarseGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setVerticalCoarseGridInterval(double newGridInterval)
{
  setVerticalCoarseGrid(newGridInterval, verticalCoarseGridIsVisible);
  setDirty();
}

void CoordinateSystem::setVerticalCoarseGrid(double newGridInterval, 
                                             bool   shouldBeVisible)
{
  if( logScaledX )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      verticalCoarseGridInterval = 2.0;
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  verticalCoarseGridIsVisible = shouldBeVisible;
  verticalCoarseGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setVerticalFineGridVisible(bool shouldBeVisible)
{
  setVerticalFineGrid(verticalFineGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setVerticalFineGridInterval(double newGridInterval)
{
  setVerticalFineGrid(newGridInterval, verticalFineGridIsVisible);
  setDirty();
}

void CoordinateSystem::setVerticalFineGrid(double newGridInterval, 
                                           bool   shouldBeVisible)
{
  if( logScaledX )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      verticalFineGridInterval = pow(2.0, 1.0/3.0);
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  verticalFineGridIsVisible = shouldBeVisible;
  verticalFineGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setRadialCoarseGridVisible(bool shouldBeVisible)
{
  setRadialCoarseGrid(radialCoarseGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setRadialCoarseGridInterval(double newGridInterval)
{
  setRadialCoarseGrid(newGridInterval, radialCoarseGridIsVisible);
  setDirty();
}

void CoordinateSystem::setRadialCoarseGrid(double newGridInterval, 
                                           bool   shouldBeVisible)
{
  if( logScaledRadius )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      radialCoarseGridInterval = 2.0;
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  radialCoarseGridIsVisible = shouldBeVisible;
  radialCoarseGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setRadialFineGridVisible(bool shouldBeVisible)
{
  setRadialFineGrid(radialFineGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setRadialFineGridInterval(double newGridInterval)
{
  setRadialFineGrid(newGridInterval, radialFineGridIsVisible);
  setDirty();
}

void CoordinateSystem::setRadialFineGrid(double newGridInterval, 
                                         bool   shouldBeVisible)
{
  if( logScaledRadius )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      radialFineGridInterval = pow(2.0, 1.0/3.0);
      return;
    }
  }
  else
  {
    jassert(newGridInterval > 0.000001); 
    // grid-intervals must be > 0
    if( newGridInterval <= 0.000001 )
      return;
  }

  radialFineGridIsVisible     = shouldBeVisible;
  radialFineGridInterval = newGridInterval;

  setDirty();
}

void CoordinateSystem::setAngularCoarseGridVisible(bool shouldBeVisible)
{
  setAngularCoarseGrid(angularCoarseGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setAngularCoarseGridInterval(double newGridInterval)
{
  setAngularCoarseGrid(newGridInterval, angularCoarseGridIsVisible);
  setDirty();
}

void CoordinateSystem::setAngularCoarseGrid(double newGridInterval, 
                                            bool   shouldBeVisible)
{
  jassert(newGridInterval > 0.000001); 
  // grid-intervals must be > 0
  if( newGridInterval <= 0.000001 )
    return;

  angularCoarseGridIsVisible = shouldBeVisible;
  angularCoarseGridInterval  = newGridInterval;

  setDirty();
}

void CoordinateSystem::setAngularFineGridVisible(bool shouldBeVisible)
{
  setAngularFineGrid(angularFineGridInterval, shouldBeVisible);
  setDirty();
}

void CoordinateSystem::setAngularFineGridInterval(double newGridInterval)
{
  setAngularFineGrid(newGridInterval, angularFineGridIsVisible);
  setDirty();
}

void CoordinateSystem::setAngularFineGrid(double newGridInterval, 
                                          bool   shouldBeVisible)
{
  jassert(newGridInterval > 0.000001); 
  // grid-intervals must be > 0
  if( newGridInterval <= 0.000001 )
    return;

  angularFineGridIsVisible     = shouldBeVisible;
  angularFineGridInterval = newGridInterval;

  setDirty();
}

void CoordinateSystem::setAngleUnitToDegrees(bool shouldBeInDegrees)
{
  angleIsInDegrees = shouldBeInDegrees;
}

bool CoordinateSystem::isHorizontalCoarseGridVisible() 
{ 
  return horizontalCoarseGridIsVisible; 
}

bool CoordinateSystem::isHorizontalFineGridVisible() 
{ 
  return horizontalFineGridIsVisible; 
}

bool CoordinateSystem::isVerticalCoarseGridVisible() 
{ 
  return verticalCoarseGridIsVisible; 
}

bool CoordinateSystem::isVerticalFineGridVisible() 
{ 
  return verticalFineGridIsVisible; 
}

bool CoordinateSystem::isRadialCoarseGridVisible() 
{ 
  return radialCoarseGridIsVisible; 
}

bool CoordinateSystem::isRadialFineGridVisible() 
{ 
  return radialFineGridIsVisible; 
}

bool CoordinateSystem::isAngularCoarseGridVisible() 
{ 
  return angularCoarseGridIsVisible; 
}

bool CoordinateSystem::isAngularFineGridVisible() 
{ 
  return angularFineGridIsVisible; 
}

double CoordinateSystem::getHorizontalCoarseGridInterval()
{
  return horizontalCoarseGridInterval;
}

double CoordinateSystem::getHorizontalFineGridInterval()
{
  return horizontalFineGridInterval;
}

double CoordinateSystem::getVerticalCoarseGridInterval()
{
  return verticalCoarseGridInterval;
}

double CoordinateSystem::getVerticalFineGridInterval()
{
  return verticalFineGridInterval;
}

double CoordinateSystem::getRadialCoarseGridInterval()
{
  return radialCoarseGridInterval;
}

double CoordinateSystem::getRadialFineGridInterval()
{
  return radialFineGridInterval;
}

double CoordinateSystem::getAngularCoarseGridInterval()
{
  return angularCoarseGridInterval;
}

double CoordinateSystem::getAngularFineGridInterval()
{
  return angularFineGridInterval;
}

void CoordinateSystem::setAxisPositionX(int newAxisPositionX)
{
  if(  newAxisPositionX == INVISIBLE || newAxisPositionX == ZERO || newAxisPositionX == TOP 
    || newAxisPositionX == BOTTOM       )
  {
    axisPositionX = newAxisPositionX;
    setDirty();
  }
}

void CoordinateSystem::setAxisPositionY(int newAxisPositionY)
{
  if(  newAxisPositionY == INVISIBLE || newAxisPositionY == ZERO || newAxisPositionY == LEFT 
    || newAxisPositionY == RIGHT        )
  {
    axisPositionY = newAxisPositionY;
    setDirty();
  }
}

void CoordinateSystem::setupAxisX(double newMin, double newMax, bool shouldBeLogScaled, 
  double newLogBase, int newAxisPosition, double newCoarseGridInterval, double newFineGridInterval)
{
  // axis settings seem not to make sense
  jassert(newMin < newMax);
  jassert(newMin > 0.0 || shouldBeLogScaled == false);
  jassert(shouldBeLogScaled == false ||
    (newCoarseGridInterval > 1.000001 && newFineGridInterval > 1.000001));
  if( newMin >= newMax )
    return;
  if( newMin < 0.0 && shouldBeLogScaled == true )
    return;
  if( shouldBeLogScaled == true && 
     (newCoarseGridInterval <= 1.000001 || newFineGridInterval <= 1.000001) )
    return;

  maximumRange.setRangeX(newMin, newMax);
  currentRange.setRangeX(newMin, newMax);
  logScaledX = shouldBeLogScaled;
  logBaseX     = newLogBase;
  if(  newAxisPosition == INVISIBLE || newAxisPosition == ZERO || newAxisPosition == TOP 
    || newAxisPosition == BOTTOM       )
  {
    axisPositionX = newAxisPosition;
  }
  verticalCoarseGridInterval = newCoarseGridInterval;
  verticalFineGridInterval   = newFineGridInterval;

  updateScaleFactors();
  setDirty();
}

void CoordinateSystem::setupAxisY(double newMin, double newMax, bool shouldBeLogScaled, 
  double newLogBase, int newAxisPosition, double newCoarseGridInterval, 
  double newFineGridInterval)
{
  // axis settings seem not to make sense
  jassert(newMin < newMax);
  jassert(newMin > 0.0 || shouldBeLogScaled == false);
  jassert(shouldBeLogScaled == false || (newCoarseGridInterval > 1.000001 &&
      newFineGridInterval > 1.000001));
  if( newMin >= newMax )
    return;
  if( newMin < 0.0 && shouldBeLogScaled == true )
    return;
  if( shouldBeLogScaled == true && 
     (newCoarseGridInterval <= 1.000001 || newFineGridInterval <= 1.000001) )
    return;

  maximumRange.setRangeY(newMin, newMax);
  currentRange.setRangeY(newMin, newMax);
  logScaledY = shouldBeLogScaled;
  logBaseY     = newLogBase;
  if( newAxisPosition == INVISIBLE ||
      newAxisPosition == ZERO      ||
      newAxisPosition == LEFT      ||
      newAxisPosition == RIGHT       )
  {
    axisPositionY = newAxisPosition;
  }
  horizontalCoarseGridInterval = newCoarseGridInterval;
  horizontalFineGridInterval   = newFineGridInterval;

  updateScaleFactors();
  setDirty();
}

void CoordinateSystem::useLogarithmicScale(bool   shouldBeLogScaledX, 
                                           bool   shouldBeLogScaledY, 
                                           double newLogBaseX, 
                                           double newLogBaseY)
{
  logScaledX = shouldBeLogScaledX;
  logScaledY = shouldBeLogScaledY;
  logBaseX     = newLogBaseX;
  logBaseY     = newLogBaseY;

  updateScaleFactors();
  setDirty();
}

void CoordinateSystem::useLogarithmicScaleX(bool   shouldBeLogScaledX, 
                                            double newLogBaseX)
{
  logScaledX = shouldBeLogScaledX;
  logBaseX     = newLogBaseX;
  updateScaleFactors();
  setDirty();
}

bool CoordinateSystem::isLogScaledX()
{
  return logScaledX;
}

void CoordinateSystem::useLogarithmicScaleY(bool   shouldBeLogScaledY, 
                                            double newLogBaseY)
{
  logScaledY = shouldBeLogScaledY;
  logBaseY     = newLogBaseY;
  updateScaleFactors();
  setDirty();
}

bool CoordinateSystem::isLogScaledY()
{
  return logScaledY;
}

void CoordinateSystem::setValueFieldPopup(bool shouldPopUp)
{
  valuePopup = shouldPopUp;
  inspectionField->setVisible(false);
}

//-------------------------------------------------------------------------------------------------
// functions for drawing and/or exporting the shown content:

void CoordinateSystem::openRightClickPopupMenu()
{
  // create a context menu to allow for export:
  PopupMenu menu;
  menu.addItem(1, "Export Image");
  const int result = menu.show();

  if (result == 0)
  {
    // user dismissed the menu without picking anything
  }
  else if (result == 1)
  {
    // user picked the Export item - open the export dialog window:
    openExportDialog(getWidth(), getHeight(), String("png"), File());
  }
}

void CoordinateSystem::addLineToSvgDrawing(XmlElement* theSVG, float x1, float y1, float x2, 
  float y2, float thickness, Colour colour, bool withArrowHead)
{
  if( theSVG == NULL )
    return;

  XmlElement* line = new XmlElement(String("line"));
  line->setAttribute(String("x1"), x1);
  line->setAttribute(String("y1"), y1);
  if( withArrowHead == true && y1 == y2 )
    line->setAttribute(String("x2"), x2-8);
  else
    line->setAttribute(String("x2"), x2);

  if( withArrowHead == true && x1 == x2 )
    line->setAttribute(String("y2"), y2+8);
  else
    line->setAttribute(String("y2"), y2);

  line->setAttribute(String("style"), String("stroke-width: ") + String(thickness) + 
    String("; stroke: #") + colour.toString().substring(2) + String(";") );
  theSVG->addChildElement(line);

  if( withArrowHead == true )
  {
    XmlElement* triangle = new XmlElement(String("path"));

    if( y1 == y2 ) // this is a horizontal rightward arrow 
    {
      triangle->setAttribute(String("d"), 
        String("M ")   + String(x2-8) + String(" ") + String(y2-4) + 
        String(", L ") + String(x2-8) + String(" ") + String(y2+4) + 
        String(", L ") + String(x2)   + String(" ") + String(y2)   + 
        String(", Z") );
      triangle->setAttribute(String("style"), String("stroke: none, fill: #") 
        + colour.toString().substring(2) + String(";") );
    }
    else if( x1 == x2 ) // this is an upward verzical rightward arrow 
    {
      triangle->setAttribute(String("d"), 
        String("M ")   + String(x2-4) + String(" ") + String(y2+8) + 
        String(", L ") + String(x2+4) + String(" ") + String(y2+8) + 
        String(", L ") + String(x2)   + String(" ") + String(y2)   + 
        String(", Z") );
      triangle->setAttribute(String("style"), String("stroke: none, fill: #") 
        + colour.toString().substring(2) + String(";") );
    }

    theSVG->addChildElement(triangle);
  }
}

void CoordinateSystem::addTextToSvgDrawing(XmlElement* theSVG, String theText, float x, float y, 
                                           Justification justification)
{
  if( theSVG == NULL )
    return;

  XmlElement* textContainer = new XmlElement(String("text"));
  auto text          = XmlElement::createTextElement(theText);

  String jString = String();
  if( justification.getFlags() == Justification::centredLeft )
    jString = String("start");
  else if( justification.getFlags() == Justification::centred )
    jString = String("middle");
  else if( justification.getFlags() == Justification::centredRight )
    jString = String("end");

  textContainer->setAttribute(String("x"), x);
  textContainer->setAttribute(String("y"), y);
  textContainer->setAttribute(String("style"), String("font-family: sans-serif;") +  
    String(" font-size: 12px;") + String(" stroke: none;") + String(" fill: black;") +
    String(" text-anchor: ") + jString + String(";") );
  textContainer->addChildElement(text);
  theSVG->addChildElement(textContainer);
}

void CoordinateSystem::drawBitmapText(Graphics &g, const String &text, double x, double y, 
  double w, double h, BitmapFont const* font, Justification justification)
{
  if( useBitmapFont == true )
  {
    if( font == NULL )
    {
      jassertfalse;
      return;
    }

    int hFlags = justification.getOnlyHorizontalFlags();
    if( justification.testFlags(hFlags & Justification::horizontallyCentred) )
      x = (x+x+w)/2.0; // average between x and x+w
    else if( justification.testFlags(hFlags & Justification::right) )
      x = x+w;

    int vFlags = justification.getOnlyVerticalFlags();
    if( justification.testFlags(vFlags & Justification::verticallyCentred) )
      y = (y+y+h)/2.0;
    else if( justification.testFlags(vFlags & Justification::bottom) )
      y = y+h;

    int  xInt = roundToInt(x);
    int  yInt = roundToInt(y);
    drawBitmapFontText(g, xInt, yInt, text, font, getTextColour(), -1, justification);
  }
  else
  {
    g.drawText(text, (int)x, (int)y, (int)w, (int)h, justification, false);
  }
}

Image* CoordinateSystem::getPlotAsImage(int width, int height)
{
  jassert(width  >= 1);
  jassert(height >= 1);    
  if( width < 1 || height < 1)
    return NULL;

  Image* thePlotImage = new Image(Image::RGB, width, height, true);

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*thePlotImage);

  drawCoordinateSystem(g, thePlotImage);

  return thePlotImage;
}

XmlElement* CoordinateSystem::getPlotAsSVG(int width, int height)
{
  jassert(width  >= 1);
  jassert(height >= 1);    
  if( width < 1 || height < 1)
    return NULL;

  // create an image object as dummy:
  Image* thePlotImage = new Image(Image::RGB, width, height, true);

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*thePlotImage);

  // create an XmlElement to be used for the SVG drawing:
  XmlElement* theSVG = new XmlElement(String("svg"));
  theSVG->setAttribute(String("width"),  width);
  theSVG->setAttribute(String("height"), height);

  // draw on the SVG:
  drawCoordinateSystem(g, thePlotImage, theSVG);

  // delete the dummy image:
  delete thePlotImage; 

  return theSVG;
}

void CoordinateSystem::openExportDialog(int defaultWidth, int defaultHeight, 
                                        const String &defaultFormat,
                                        const File& defaultTargetFile)
{
  /*
  ImageSavingDialog dialog(this, defaultWidth, defaultHeight, defaultFormat, defaultTargetFile);
  DialogWindow exportWindow(String(T("Export to Image or SVG Drawing")), Colours::white, true, true);
  exportWindow.showModalDialog(String(T("Export to Image or SVG Drawing")), &dialog, this, 
    Colours::white, true, false, false);
  */
}


void CoordinateSystem::drawComponent(Image* imageToDrawOnto)
{
  Graphics g(*imageToDrawOnto);
  drawCoordinateSystem(g);
  //repaint(); // necesarry?
}

/*
void CoordinateSystem::updateBackgroundImage()
{
  if( getWidth() < 1 || getHeight() < 1 )
    return;

  // allocate memory for the first time:
  if( backgroundImage == NULL )
  {
    backgroundImage = new Image(Image::RGB, getWidth(), getHeight(), true);

    if( backgroundImage == NULL )
      return; // memory allocation failed
  }

  // reallocate memory, if necessary (i.e. the size of the component differs
  // from the size of the current image):
  if( backgroundImage->getWidth()  != getWidth()  ||
    backgroundImage->getHeight() != getHeight()    )
  {
    // delete the old and create a new Image-object:
    if( backgroundImage != NULL )
    {
      delete backgroundImage;
      backgroundImage = NULL;
    }
    backgroundImage = new Image(Image::RGB, getWidth(), getHeight(), true);

    if( backgroundImage == NULL )
      return; // memory allocation failed
  }

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*backgroundImage);

  drawCoordinateSystem(g);

  // trigger a repaint:
  repaint();
}
*/

//-------------------------------------------------------------------------------------------------
// coordinate transformations:

inline double logB(double x, double b)
{
  // temporary - remove and use the function from RAPT
  return log(x) / log(b);
}

void CoordinateSystem::transformToImageCoordinates(
  double &x, double &y, const Image* theImage) const
{
  if( theImage == NULL )
  {
    toPixelCoordinates(x, y);
    return;
  }

  // transform the x,y values to coordinates inside this component:
  if( logScaledX )
  {
    jassert( x > 0.0 && currentRange.getMinX() > 0.0 );
    // caught a logarithm of a non-positive number
    if( x <= 0.0 || currentRange.getMinX() <= 0.0 )
      return;

    double imagePixelsPerIntervalX = 
      theImage->getWidth()/logB((currentRange.getMaxX()/currentRange.getMinX()), logBaseX);

    x = imagePixelsPerIntervalX * logB((x/currentRange.getMinX()), logBaseX);
  }
  else
  {
    double imageScaleX = theImage->getWidth()  / (currentRange.getMaxX()-currentRange.getMinX()); 

    x -= currentRange.getMinX();	   // shift origin left/right
    x *= imageScaleX;	       // scale to fit width
  }

  if( logScaledY )
  {
    jassert( y > 0.0 && currentRange.getMinY() > 0.0 );
    // caught a logarithm of a non-positive number
    if( y <= 0.0 || currentRange.getMinY() <= 0.0 )
      return;

    double imagePixelsPerIntervalY = 
      theImage->getHeight()/logB((currentRange.getMaxY()/currentRange.getMinY()), logBaseY);

    y = imagePixelsPerIntervalY * logB((y/currentRange.getMinY()), logBaseY);
  }
  else
  {
    double imageScaleY = theImage->getHeight() / (currentRange.getMaxY()-currentRange.getMinY());

    y -= currentRange.getMinY();	// shift origin up/down
    y *= imageScaleY;	         // scale to fit height
  }

  y  = theImage->getHeight()-y;	   // invert (pixels begin at top-left)
}

void CoordinateSystem::transformFromImageCoordinates(
  double &x, double &y, const Image *theImage) const
{
  if( theImage == NULL )
  {
    fromPixelCoordinates(x, y);
    return;
  }

  if( logScaledX )
  {
    double imagePixelsPerIntervalX = 
      theImage->getWidth()/logB((currentRange.getMaxX()/currentRange.getMinX()), logBaseX);
    x = currentRange.getMinX() * pow(logBaseX, (x/imagePixelsPerIntervalX));
  }
  else
  {
    double imageScaleX = theImage->getWidth()  / (currentRange.getMaxX()-currentRange.getMinX());
    x /= imageScaleX;              // scale to fit width
    x += currentRange.getMinX();    // shift origin left/right
  }
  if( logScaledY )
  {
    double imagePixelsPerIntervalY = 
      theImage->getHeight()/logB((currentRange.getMaxY()/currentRange.getMinY()), logBaseY);
    y = currentRange.getMinY() * pow(logBaseY, (y/imagePixelsPerIntervalY));
  }
  else
  {
    double imageScaleY = theImage->getHeight() / (currentRange.getMaxY()-currentRange.getMinY());
    y  = getHeight()-y; 
    y /= imageScaleY;             // scale to fit height
    y += currentRange.getMinY();    // shift origin up/down
  }
}

void CoordinateSystem::toPixelCoordinates(double &x, double &y) const
{
  // transform the x,y values to coordinates inside this component:
  if( logScaledX )
  {
    jassert( x > 0.0 && currentRange.getMinX() > 0.0 );
    // caught a logarithm of a non-positive number
    if( x <= 0.0 || currentRange.getMinX() <= 0.0 )
      return;

    x = pixelsPerIntervalX * logB((x/currentRange.getMinX()), logBaseX);
  }
  else
  {
    x -= currentRange.getMinX();	// shift origin left/right
    x *= scaleX;	         // scale to fit width
  }

  if( logScaledY )
  {
    jassert( y > 0.0 && currentRange.getMinY() > 0.0 );
    // caught a logarithm of a non-positive number
    if( y <= 0.0 || currentRange.getMinY() <= 0.0 )
      return;

    y = pixelsPerIntervalY * logB((y/currentRange.getMinY()), logBaseY);
  }
  else
  {
    y -= currentRange.getMinY();	// shift origin up/down
    y *= scaleY;	         // scale to fit height
  }

  y  = getHeight()-y;	   // invert (pixels begin at top-left)
}

void CoordinateSystem::toPixelCoordinates(float &x, float &y) const
{
  double xd = (double) x;
  double yd = (double) y;
  toPixelCoordinates(xd, yd);
  x = (float) xd;
  y = (float) yd;
}

void CoordinateSystem::fromPixelCoordinates(double &x, double &y) const
{
  if( logScaledX )
  {
    x = currentRange.getMinX() * pow(logBaseX, (x/pixelsPerIntervalX));
  }
  else
  {
    x /= scaleX;             // scale to fit width
    x += currentRange.getMinX();    // shift origin left/right
  }

  if( logScaledY )
  {
    y = currentRange.getMinY() * pow(logBaseY, (y/pixelsPerIntervalY));
  }
  else
  {
    y  = getHeight()-y; 
    y /= scaleY;             // scale to fit height
    y += currentRange.getMinY();    // shift origin up/down
  }
}

void CoordinateSystem::fromPixelCoordinates(float &x, float &y) const
{
  double xd = (double) x;
  double yd = (double) y;
  fromPixelCoordinates(xd, yd);
  x = (float) xd;
  y = (float) yd;
}

//-------------------------------------------------------------------------------------------------
// drawing functions

void CoordinateSystem::drawCoordinateSystem(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  g.setFont(Font(14));

  // draw the background:
  //g.fillAll(colourScheme.backgroundColour);
  if( targetImage != NULL )
  {
    Graphics g2(*targetImage);
    fillRectWithBilinearGradient(g2, 0, 0, targetImage->getWidth(), targetImage->getHeight(),
      colourScheme.topLeft, colourScheme.topRight, colourScheme.bottomLeft, colourScheme.bottomRight);
  }
  else
  {
    fillRectWithBilinearGradient(g, 0, 0, getWidth(), getHeight(),
      colourScheme.topLeft, colourScheme.topRight, colourScheme.bottomLeft, colourScheme.bottomRight);
  }

  // draw the grids, if desired:
  if( horizontalFineGridIsVisible )
    drawHorizontalGrid(g, horizontalFineGridInterval, logScaledY, colourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( verticalFineGridIsVisible )
    drawVerticalGrid(g, verticalFineGridInterval, logScaledX, colourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( radialFineGridIsVisible )
    drawRadialGrid(g, radialFineGridInterval, logScaledRadius, colourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( angularFineGridIsVisible )
    drawAngularGrid(g, angularFineGridInterval, colourScheme.fineGrid, 1.0f, targetImage, targetSVG);

  if( horizontalCoarseGridIsVisible )
    drawHorizontalGrid(g, horizontalCoarseGridInterval, logScaledY, colourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( verticalCoarseGridIsVisible )
    drawVerticalGrid(g, verticalCoarseGridInterval, logScaledX, colourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( radialCoarseGridIsVisible )
    drawRadialGrid(g, radialCoarseGridInterval, logScaledRadius, colourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( angularCoarseGridIsVisible )
    drawAngularGrid(g, angularCoarseGridInterval, colourScheme.coarseGrid, 1.0f, targetImage, targetSVG);

  // draw the coordinate system:
  if( axisPositionX != INVISIBLE )
    drawAxisX(g, targetImage, targetSVG);
  if( axisPositionY != INVISIBLE )
    drawAxisY(g, targetImage, targetSVG);

  // draw the labels on the axes:
  if( axisPositionX != INVISIBLE && axisLabelPositionX != NO_ANNOTATION )
    drawAxisLabelX(g, targetImage, targetSVG);
  if( axisPositionY != INVISIBLE && axisLabelPositionY != NO_ANNOTATION )
    drawAxisLabelY(g, targetImage, targetSVG);

  // draw the values on the axes:
  if( axisPositionX != INVISIBLE && axisValuesPositionX != NO_ANNOTATION )
    drawAxisValuesX(g, targetImage, targetSVG);
  if( axisPositionY != INVISIBLE && axisValuesPositionY != NO_ANNOTATION )
    drawAxisValuesY(g, targetImage, targetSVG);

  // draw the caption:
  drawCaption(g, targetImage, targetSVG);

  // draw an outlining rectangle:
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), 1);
}

void CoordinateSystem::drawCaption(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  //float x, y, w, h;
  //caption.getBounds(x, y, w, h);

  Rectangle<int> bounds = caption.getBounds();
  float x = (float) bounds.getX();
  float y = (float) bounds.getY();
  float w = (float) bounds.getWidth();

  g.setColour(colourScheme.axes);
  const BitmapFont* font = &BitmapFontRoundedBoldA10D0::instance;

  switch( captionPosition )
  {
  case NO_CAPTION: return;
  case TOP_CENTER:
    {
      //caption.drawAt(g, 0.5f*getWidth()-0.5f*w, 16);  
      drawBitmapText(g, captionString, 0.5f*getWidth()-0.5f*w, 16, getWidth(), 16,         
        font, Justification::centred);

      if( targetSVG != NULL )
      {
        float centerX;
        if( targetImage != NULL )
          centerX = 0.5f*targetImage->getWidth()-0.5f*w;
        else
          centerX = 0.5f*getWidth()-0.5f*w;
        addTextToSvgDrawing(targetSVG, captionString, centerX, 16, Justification::centred);
      }
    }
    break;
  case CENTER:
    {
      //caption.drawAt(g, 0.5f*getWidth()-0.5f*w, 0.5f*getHeight()-0.5f*h+10.f );

      x  = (float)getWidth()/2.f;
      y  = (float)getHeight()/2.f;
      x -= font->getTextPixelWidth(captionString, font->getDefaultKerning())/2;
      y -= font->getFontAscent()/2;
      drawBitmapText(g, captionString, x, y, getWidth(), 16, font, Justification::topLeft);

      //drawText(g, captionString, 0.5f*getWidth()-0.5f*w, 0.5f*getHeight()-0.5f*h+8.f, 
      //  getWidth(), 16, Justification::centred);
    }
    break;
  }
}


void CoordinateSystem::drawHorizontalGrid(Graphics &g, double interval, 
                                          bool exponentialSpacing, 
                                          Colour gridColour, 
                                          float lineThickness,
                                          Image* targetImage, XmlElement *targetSVG)
{
  if( exponentialSpacing == true )
  {
    jassert( interval >= 1.00001 );
    // grid spacing must be > 1 for exponentially spaced grid-lines
    if( interval < 1.00001 )
      return;
  }
  else
  {
    jassert( interval >= 0.000001 );
    // grid spacing must be > 0
    if( interval < 0.000001 )
      return;
  }

  g.setColour(gridColour);

  long	  i;
  double	startX, endX, startY, endY;
  double accumulator;

  String gridPathDataString;

  if( exponentialSpacing ) // draw grid with exponentially spaced lines
  {
    accumulator = interval*maximumRange.getMinY();
    while( accumulator < maximumRange.getMaxY() )
    {
      startX = currentRange.getMinX();
      endX   = currentRange.getMaxX();
      startY = accumulator;
      endY   = accumulator;

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing:
      if( targetSVG != NULL )
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      accumulator *= interval;
    }
  }
  else // draw grid with linearly spaced lines
  {
    i = 0;
    while( i*interval < currentRange.getMaxY() )
    {
      startX = currentRange.getMinX();
      endX   = currentRange.getMaxX();
      startY = i*interval;
      endY   = i*interval;

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing member:
      if( targetSVG != NULL )
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      i++;
    }
    i = 1;
    while( -i*interval > currentRange.getMinY() )
    {
      startX = currentRange.getMinX();
      endX   = currentRange.getMaxX();
      startY = -i*interval;
      endY   = -i*interval;

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing member:
      if( targetSVG != NULL )
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      i++;
    } // end while
  } // end else


  if( targetSVG != NULL )
  {
    XmlElement* gridPath = new XmlElement(String("path"));
    gridPath->setAttribute(String("d"), gridPathDataString);
    gridPath->setAttribute(String("style"), String("stroke-width: ") + String(lineThickness) + 
      String("; stroke: #") + gridColour.toString().substring(2) + String(";") );
    targetSVG->addChildElement(gridPath);
  }
}

void CoordinateSystem::drawVerticalGrid(Graphics &g, double interval, 
                                        bool exponentialSpacing, 
                                        Colour gridColour, 
                                        float lineThickness,
                                        Image* targetImage, XmlElement *targetSVG)
{
  if( exponentialSpacing == true )
  {
    jassert( interval >= 1.00001 );
    // grid spacing must be > 1 for exponentially spaced grid-lines
    if( interval < 1.00001 )
      return;
  }
  else
  {
    jassert( interval >= 0.000001 );
    // grid spacing must be > 0
    if( interval < 0.000001 )
      return;
  }

  g.setColour(gridColour);

  int    	i; 
  double	 startX, endX, startY, endY;
  double  accumulator;

  String gridPathDataString;

  if( exponentialSpacing ) // draw grid with exponentially spaced lines
  {
    accumulator = interval*maximumRange.getMinX();
    while( accumulator < maximumRange.getMaxX() )
    {
      startX = accumulator;
      endX   = accumulator;
      startY = currentRange.getMinY();
      endY   = currentRange.getMaxY();

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing:
      if( targetSVG != NULL )
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      accumulator *= interval;
    }
  }
  else // draw grid with linearly spaced lines
  {
    // draw vertical lines:
    i = 0;
    while( i*interval < currentRange.getMaxX() )
    {
      startX = i*interval;
      endX   = i*interval;
      startY = currentRange.getMinY();
      endY   = currentRange.getMaxY();

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing:
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      i++;
    }
    i = 1;
    while( -i*interval > currentRange.getMinX() )
    {
      startX = -i*interval;
      endX   = -i*interval;
      startY = currentRange.getMinY();
      endY   = currentRange.getMaxY();

      // transform:
      if( targetImage == NULL )
      {
        toPixelCoordinates(startX, startY);
        toPixelCoordinates(endX, endY);
      }
      else
      {
        transformToImageCoordinates(startX, startY, targetImage);
        transformToImageCoordinates(endX, endY, targetImage);
      }

      // draw:
      g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
        lineThickness);

      // add the line to the path which will be added to the SVG drawing:
      if( targetSVG != NULL )
      {
        gridPathDataString += String("M ") + String(startX) + String(" ") 
          + String(startY) + String(" ");
        gridPathDataString += String("L ") + String(endX) + String(" ")   
          + String(endY) + String(" ");
      }

      i++;
    } // end while
  } // end else

  if( targetSVG != NULL )
  {
    XmlElement* gridPath = new XmlElement(String("path"));
    gridPath->setAttribute(String("d"), gridPathDataString);
    gridPath->setAttribute(String("style"), String("stroke-width: ") + String(lineThickness) 
      + String("; stroke: #") + gridColour.toString().substring(2) + String(";") );
    targetSVG->addChildElement(gridPath);
  }
}

void CoordinateSystem::drawRadialGrid(Graphics &g, double interval, 
                                      bool exponentialSpacing, 
                                      Colour gridColour, 
                                      float lineThickness,
                                      Image* targetImage, XmlElement *targetSVG)
{
  if( exponentialSpacing == true )
  {
    jassert( interval >= 1.00001 );
    // grid spacing must be > 1 for exponentially spaced grid-lines
    if( interval < 1.00001 )
      return;
  }
  else
  {
    jassert( interval >= 0.000001 );
    // grid spacing must be > 0
    if( interval < 0.000001 )
      return;
  }

  g.setColour(gridColour);

  // calculate the radius of the largest circle to be drawn:
  double xTmp = jmax(fabs(currentRange.getMinX()), fabs(currentRange.getMaxX()) );
  double yTmp = jmax(fabs(currentRange.getMinY()), fabs(currentRange.getMaxY()) );
  double maxRadius = sqrt(xTmp*xTmp + yTmp*yTmp);

  // calculate the center-coordinates of the circles in terms of components
  // coordinates:
  double centerX = 0.0;
  double centerY = 0.0;
  double xScaler, yScaler;
  if( targetImage == NULL )
  {
    xScaler = getWidth()  / (currentRange.getMaxX()-currentRange.getMinX());
    yScaler = getHeight() / (currentRange.getMaxY()-currentRange.getMinY());
    toPixelCoordinates(centerX, centerY);
  }
  else
  {
    xScaler = targetImage->getWidth()  / (currentRange.getMaxX()-currentRange.getMinX());
    yScaler = targetImage->getHeight() / (currentRange.getMaxY()-currentRange.getMinY());
    transformToImageCoordinates(centerX, centerY, targetImage);
  }

  // draw the circles:
  int    i       = 1;
  double radius  = interval;
  double xL, xR, yT, yB;
  while( radius <= maxRadius )
  {
    // draw the circle (may deform to an ellipse depending on the scaling of the 
    // axes):
    xL = centerX - xScaler*radius;
    xR = centerX + xScaler*radius;
    yT = centerY - yScaler*radius;
    yB = centerY + yScaler*radius;
    g.drawEllipse((float)xL, (float)yT, (float)(xR-xL),(float)(yB-yT), 
      lineThickness);

    // add the circle to the svg-drawing:
    if( targetSVG != NULL )
    {
      XmlElement* ellipse = new XmlElement(String("ellipse"));
      ellipse->setAttribute(String("cx"), centerX);
      ellipse->setAttribute(String("cy"), centerY);
      ellipse->setAttribute(String("rx"), xScaler*radius);
      ellipse->setAttribute(String("ry"), yScaler*radius);
      ellipse->setAttribute(String("style"), 
        String("stroke-width: ") + String(lineThickness) + 
        String("; stroke: #") + gridColour.toString().substring(2) + String(";") + 
        String("fill: none;") );
      targetSVG->addChildElement(ellipse);
    }

    // calculate the next radius (in system-coordinates)
    i++;
    radius = interval * (double) i;
  }
}

void CoordinateSystem::drawAngularGrid(Graphics &g, double interval,
                                       Colour gridColour, 
                                       float lineThickness,
                                       Image* targetImage, XmlElement *targetSVG)
{
  g.setColour(gridColour);

  double angleIntervalInRadiant;
  if( angleIsInDegrees )
    angleIntervalInRadiant = interval*(PI/180.0);
  else
    angleIntervalInRadiant = interval;

  String gridPathDataString;

  double angle = 0.0;
  double startX, endX, startY, endY;
  int    i     = 0;
  while( angle <= PI )
  {
    endX   = cos(angle);
    endY   = sin(angle);
    startX = -endX;
    startY = -endY;

    // prolong (or shorten) the line such that it fits into the currently visible rectangle:
    fitLineToRectangle(startX, startY, endX, endY, currentRange.getMinX(), currentRange.getMinY(), 
      currentRange.getMaxX(), currentRange.getMaxY() );

    // transform:
    if( targetImage == NULL )
    {
      toPixelCoordinates(startX, startY);
      toPixelCoordinates(endX, endY);
    }
    else
    {
      transformToImageCoordinates(startX, startY, targetImage);
      transformToImageCoordinates(endX, endY, targetImage);
    }

    // draw:
    g.drawLine((float)startX, (float)startY, (float)endX, (float)endY, 
      lineThickness);

    // add the line to the SVG drawing:
    if( targetSVG != NULL )
    {

      // add the line to the path which will be added to the SVG drawing member:
      gridPathDataString += String("M ") + String(startX) + String(" ") 
        + String(startY) + String(" ");
      gridPathDataString += String("L ") + String(endX) + String(" ")   
        + String(endY) + String(" ");
    }

    i++;
    angle = angleIntervalInRadiant * (double) i;
  }

  if( targetSVG != NULL )
  {
    XmlElement* gridPath = new XmlElement(String("path"));
    gridPath->setAttribute(String("d"), gridPathDataString);
    gridPath->setAttribute(String("style"), String("stroke-width: ") + String(lineThickness) 
      + String("; stroke: #") + gridColour.toString().substring(2) + String(";") );
    targetSVG->addChildElement(gridPath);
  }
}

void CoordinateSystem::drawAxisX(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( logScaledX == true )
  {
    jassert( verticalCoarseGridInterval >= 1.00001 );
    if( verticalCoarseGridInterval < 1.00001 )
      return;
  }
  else
  {
    jassert( verticalCoarseGridInterval >= 0.000001 );
    // grid spacing must be > 0
    if( verticalCoarseGridInterval < 0.000001 )
      return;
  }

  if( axisPositionX == INVISIBLE ) 
    return;

  g.setColour(colourScheme.axes);

  double startX, endX, startY = 0, endY;

  startX = currentRange.getMinX();
  endX	  = currentRange.getMaxX();
  if( logScaledY )
    startY = currentRange.getMinY();
  else if( axisPositionX == ZERO )
    startY = 0.0;
  else if( axisPositionX == TOP ) 
    startY = currentRange.getMaxY();
  else if( axisPositionX == BOTTOM ) 
    startY = currentRange.getMinY();
  endY = startY;

  // transform:
  if( targetImage == NULL )
  {
    toPixelCoordinates(startX, startY);
    toPixelCoordinates(endX, endY);
  }
  else
  {
    transformToImageCoordinates(startX, startY, targetImage);
    transformToImageCoordinates(endX, endY, targetImage);
  }

  // include some margin for axes at the top and bottom:
  if( axisPositionX == TOP )
  {
    startY += 8;
    endY   += 8;
  }
  else if( axisPositionX == BOTTOM )
  {
    startY -= 8;
    endY   -= 8;
  }

  // draw:
  //g.drawArrow((float)startX, (float)startY, (float)endX, (float)endY, 
  //  2.0, 8.0, 8.0);
  g.drawArrow(Line<float>((float)startX, (float)startY, (float)endX, (float)endY), 1.0, 6.0, 6.0);

  if( targetSVG != NULL )
    addLineToSvgDrawing(targetSVG, (float)startX, (float)startY, (float) endX, (float)endY, 2.0, colourScheme.axes, true);
}

void CoordinateSystem::drawAxisY(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( logScaledY == true )
  {
    jassert( horizontalCoarseGridInterval >= 1.00001 );
    if( horizontalCoarseGridInterval < 1.00001 )
      return;
  }
  else
  {
    jassert( horizontalCoarseGridInterval >= 0.000001 );
    // grid spacing must be > 0
    if( horizontalCoarseGridInterval < 0.000001 )
      return;
  }

  if( axisPositionY == INVISIBLE ) 
    return;

  g.setColour(colourScheme.axes);

  double startX = 0, endX, startY, endY;

  startY = currentRange.getMinY();
  endY	  = currentRange.getMaxY();

  if( logScaledX )
    startX = currentRange.getMinX();
  else if( axisPositionY == ZERO )
    startX = 0.0;
  else if( axisPositionY == LEFT ) 
    startX = currentRange.getMinX();
  else if( axisPositionY == RIGHT ) 
    startX = currentRange.getMaxX();
  endX = startX;

  // transform:
  if( targetImage == NULL )
  {
    toPixelCoordinates(startX, startY);
    toPixelCoordinates(endX, endY);
  }
  else
  {
    transformToImageCoordinates(startX, startY, targetImage);
    transformToImageCoordinates(endX, endY, targetImage);
  }

  // include some margin for axes at the left and right:
  if( axisPositionY == LEFT )
  {
    startX += 8;
    endX   += 8;
  }
  else if( axisPositionY == RIGHT )
  {
    startX -= 8;
    endX   -= 8;
  }

  // draw:
  //g.drawArrow((float)startX, (float)startY, (float)endX, (float)endY, 
  //  2.0, 8.0, 8.0);
  g.drawArrow(Line<float>((float)startX, (float)startY, (float)endX, (float)endY), 1.0, 6.0, 6.0);

  if( targetSVG != NULL )
    addLineToSvgDrawing(targetSVG, (float)startX, (float)startY, (float) endX, (float)endY, 2.0, colourScheme.axes, true);
}

void CoordinateSystem::drawAxisLabelX(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisLabelPositionX == NO_ANNOTATION ) 
    return;

  g.setColour(colourScheme.axes);

  double posX, posY;

  // position for the label on the x-axis:
  posX = currentRange.getMaxX();
  if( logScaledY )
    posY = currentRange.getMinY();
  else if( axisPositionX == ZERO )
    posY = 0;
  else if( axisPositionX == TOP )
    posY = currentRange.getMaxY();
  else if( axisPositionX == BOTTOM )
    posY = currentRange.getMinY();

  // transform coordinates:
  if( targetImage == NULL )
    toPixelCoordinates(posX, posY);
  else
    transformToImageCoordinates(posX, posY, targetImage);

  // include some margin for axes at the top and bottom:
  posY += 2;
  if( axisPositionX == TOP )
    posY += 8;
  else if( axisPositionX == BOTTOM || axisLabelPositionX == ABOVE_AXIS )
    posY -= 28;

  BitmapFont font = BitmapFontRoundedBoldA10D0::instance;
  drawBitmapText(g, axisLabelX, (int)posX-100, (int)posY, 96, 16, &font, 
    Justification::centredRight);

  //drawTextAt(g, axisLabelX
  //posX -= valueFont->getTextPixelWidth(axisLabelX, valueFont->getDefaultKerning()) + 8;
  //drawBitmapFontText(g, (int)posX-4, (int)posY-12, axisLabelX, *valueFont, colourScheme.axesColour, -1, 
  //  Justification::centredRight);

  if( targetSVG != NULL )
    addTextToSvgDrawing(targetSVG, axisLabelX, 
    (float) (posX-4), (float) (posY+16), Justification::centredRight);
}

void CoordinateSystem::drawAxisLabelY(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisLabelPositionY == NO_ANNOTATION )
    return;

  g.setColour(colourScheme.axes);

  double posX, posY;

  if( logScaledX )
    posX = currentRange.getMinX();
  else if( axisPositionY == ZERO )
    posX = 0.0;
  else if( axisPositionY == LEFT ) 
    posX = currentRange.getMinX();
  else if( axisPositionY == RIGHT ) 
    posX = currentRange.getMaxX();

  posY  = currentRange.getMaxY();

  // transform coordinates:
  if( targetImage == NULL )
    toPixelCoordinates(posX, posY);
  else
    transformToImageCoordinates(posX, posY, targetImage);

  BitmapFont font = BitmapFontRoundedBoldA10D0::instance;

  // include some margin for axes at the left and right:
  if( axisPositionY == LEFT || axisLabelPositionY == RIGHT_TO_AXIS )
  {
    posX += 8;
    posY += 4;
    drawBitmapText(g, axisLabelY, (int) posX, (int) posY, 512, 16, &font,
      Justification::topLeft);

    //g.drawText(axisLabelY, (int) posX, (int) posY, 512, 16, 
    //  Justification::centredLeft, false);
    //drawBitmapFontText(g, (int)posX, (int)posY+8, axisLabelY, *valueFont, colourScheme.axesColour, -1,  
    //  Justification::centredLeft);

    if( targetSVG != NULL )
      addTextToSvgDrawing(targetSVG, axisLabelY, 
      (float) posX, (float) posY+12, Justification::centredLeft);
  }
  else
  {
    if( axisPositionY == ZERO )
      posX -=8;
    else
      posX -= 16;
    posX -= 512;

    //g.drawText(axisLabelY, (int) posX, (int) posY, 512, 16, 
    //  Justification::centredRight, false);
    drawBitmapFontText(g, (int)posX, (int)posY+8, axisLabelY, &font, colourScheme.axes);
    posX += 512;
    if( targetSVG != NULL )
      addTextToSvgDrawing(targetSVG, axisLabelY, 
      (float) posX, (float) posY+12, Justification::centredRight);
  }
}

void CoordinateSystem::drawAxisValuesX(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( logScaledX == true )
  {
    jassert( verticalCoarseGridInterval >= 1.00001 );
    if( verticalCoarseGridInterval < 1.00001 )
      return;
  }
  else
  {
    jassert( verticalCoarseGridInterval >= 0.000001 );
    // grid spacing must be > 0
    if( verticalCoarseGridInterval < 0.000001 )
      return;
  }

  if( axisValuesPositionX == NO_ANNOTATION )
    return;

  g.setColour(colourScheme.axes);

  long	  i;
  double	posX, posY, value;
  double accumulator;
  String numberString;

  // draw values on x-axis:
  if(logScaledX)
  {
    posX	      = verticalCoarseGridInterval*maximumRange.getMinX();
    accumulator = verticalCoarseGridInterval*maximumRange.getMinX();
    while( accumulator < maximumRange.getMaxX() )
    {
      posX = value = accumulator;
      if( logScaledY )
        posY = currentRange.getMinY();
      else
      {
        if( axisPositionX == ZERO )
          posY = 0;
        else if( axisPositionX == TOP )
          posY = currentRange.getMaxY();
        else if( axisPositionX == BOTTOM )
          posY = currentRange.getMinY();
      }

      // transform:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);

      // include some margin for axes at the top and bottom:
      if( axisPositionX == TOP )
        posY += 8;
      if( axisPositionX == BOTTOM )
        posY -= 8;

      // draw a small line:
      g.drawLine((float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 1.0, colourScheme.axes, false);

      // draw number:
      numberString = stringConversionFunctionX(value);
      if( axisValuesPositionX == ABOVE_AXIS || axisPositionX == BOTTOM )
      {
        drawBitmapText(g, numberString, posX-32, posY-20, 64, 16, &normalFont7px, 
          Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, (float)(posX-32), (float)(posY-8), 
          Justification::centred);
      }
      else
      {
        drawBitmapText(g, numberString, posX-32, posY+4, 64, 16, &normalFont7px, 
          Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, (float)(posX-32), (float)(posY+16), 
          Justification::centred);
      }

      accumulator *= verticalCoarseGridInterval;
    }
  }
  else	 // x is linerarly scaled
  {
    if( axisPositionY == LEFT  || 
      axisPositionY == RIGHT || 
      axisPositionY == INVISIBLE)
      i = 0;
    else
      i = 1;
    while( i*verticalCoarseGridInterval < currentRange.getMaxX() )
    {
      posX = value = i*verticalCoarseGridInterval;
      // "value" will not be transformed
      if( logScaledY )
        posY = currentRange.getMinY();
      else
      {
        if( axisPositionX == ZERO )
          posY = 0;
        else if( axisPositionX == TOP )
          posY = currentRange.getMaxY();
        else if( axisPositionX == BOTTOM )
          posY = currentRange.getMinY();
      }

      // transform coordinates:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);


      // include some margin for axes at the top and bottom:
      if( axisPositionX == TOP )
        posY += 8;
      if( axisPositionX == BOTTOM )
        posY -= 8;

      // draw a small line:
      g.drawLine((float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 1.0, colourScheme.axes, false);

      // draw the number:
      numberString = stringConversionFunctionX(value);
      if( axisValuesPositionX == ABOVE_AXIS || axisPositionX == BOTTOM )
      {
        drawBitmapText(g, numberString, (int)posX-32, (int)posY-20, 64, 16, 
          &normalFont7px, Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)posX, (float)(posY-8), Justification::centred);
      }
      else
      {
        drawBitmapText(g, numberString, (int)posX-32, (int)posY+4, 64, 16, 
          &normalFont7px, Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX), (float)(posY+16), Justification::centred);
      }

      i++;
    }
    i = 1;
    while( -i*verticalCoarseGridInterval > currentRange.getMinX() )
    {
      posX = value = -i*verticalCoarseGridInterval; 
      // "value" will not be transformed
      if( logScaledY )
        posY = currentRange.getMinY();
      else
      {
        if( axisPositionX == ZERO )
          posY = 0;
        else if( axisPositionX == TOP )
          posY = currentRange.getMaxY();
        else if( axisPositionX == BOTTOM )
          posY = currentRange.getMinY();
      }

      // transform coordinates:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);

      // include some margin for axes at the top and bottom:
      if( axisPositionX == TOP )
        posY += 8;
      if( axisPositionX == BOTTOM )
        posY -= 8;

      // draw a small line:
      g.drawLine((float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)posX, (float)(posY-4.0), (float)posX, (float)(posY+4.0), 1.0, colourScheme.axes, false);

      // draw the number:
      numberString = stringConversionFunctionX(value);
      if( axisValuesPositionX == ABOVE_AXIS || axisPositionX == BOTTOM )
      {
        drawBitmapText(g, numberString, (int)posX-32, (int)posY-20, 64, 16, 
          &normalFont7px, Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)posX, (float)(posY-8), Justification::centred);
      }
      else
      {
        drawBitmapText(g, numberString, (int)posX-32, (int)posY+4, 64, 16, 
          &normalFont7px, Justification::centred);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)posX, (float)(posY+16), Justification::centred);
      }

      i++;
    }
  }
}

void CoordinateSystem::drawAxisValuesY(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( logScaledY == true )
  {
    jassert( horizontalCoarseGridInterval >= 1.00001 );
    if( horizontalCoarseGridInterval < 1.00001 )
      return;
  }
  else
  {
    jassert( horizontalCoarseGridInterval >= 0.000001 );
    // grid spacing must be > 0
    if( horizontalCoarseGridInterval < 0.000001 )
      return;
  }

  if( axisValuesPositionY == NO_ANNOTATION )
    return;

  g.setColour(colourScheme.axes);

  long	  i;
  double	posX, posY, value;
  double accumulator;
  String numberString;

  // draw values on y-axis:
  if(logScaledY)
  {
    posY	      = horizontalCoarseGridInterval*maximumRange.getMinY();
    accumulator = horizontalCoarseGridInterval*maximumRange.getMinY();
    while( accumulator < maximumRange.getMaxY() )
    {
      if( logScaledX )
        posX = currentRange.getMinX();
      else
      {
        if( axisPositionY == ZERO )
          posX = 0.0;
        else if( axisPositionY == LEFT ) 
          posX = currentRange.getMinX();
        else if( axisPositionY == RIGHT ) 
          posX = currentRange.getMaxX();
      }
      posY = value = accumulator;

      // transform:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);

      // include some margin for axes at the left and right:
      if( axisPositionY == LEFT )
        posX += 8;
      else if( axisPositionY == RIGHT )
        posX -= 8;

      // draw a small line:
      g.drawLine((float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 1.0, colourScheme.axes, false);

      // draw number:
      numberString = stringConversionFunctionY(value);
      if( axisValuesPositionY == LEFT_TO_AXIS )
      {
        drawBitmapText(g, numberString, (int)posX-25, (int)posY-10, 20, 20, 
          &normalFont7px, Justification::centredRight);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX+8), (float)(posY+4), Justification::centredLeft);
      }
      else
      {
        drawBitmapText(g, numberString, (int)posX+4, (int)posY-10, 20, 20, 
          &normalFont7px, Justification::centredLeft);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX-8), (float)(posY+4), Justification::centredRight);
      }

      accumulator *= horizontalCoarseGridInterval;
    }
  }
  else // y is linearly scaled
  {
    if( axisPositionX == TOP    || 
      axisPositionX == BOTTOM || 
      axisPositionX == INVISIBLE )
      i = 0;
    else
      i = 1;
    while( i*horizontalCoarseGridInterval < currentRange.getMaxY() )
    {
      if( logScaledX )
        posX = currentRange.getMinX();
      else
      {
        if( axisPositionY == ZERO )
          posX = 0.0;
        else if( axisPositionY == LEFT ) 
          posX = currentRange.getMinX();
        else if( axisPositionY == RIGHT ) 
          posX = currentRange.getMaxX();
      }
      posY = value = i*horizontalCoarseGridInterval; 
      // "value" will not be transformed

      // transform coordinates:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);

      // include some margin for axes at the left and right:
      if( axisPositionY == LEFT )
        posX += 8;
      else if( axisPositionY == RIGHT )
        posX -= 8;

      // draw a small line:
      g.drawLine((float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 1.0, colourScheme.axes, false);

      // draw number:
      numberString = stringConversionFunctionY(value);
      if( axisValuesPositionY == RIGHT_TO_AXIS || axisPositionY == LEFT )
      {
        drawBitmapText(g, numberString, (int)posX+8, (int)posY-10, 64, 20, 
          &normalFont7px, Justification::centredLeft);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX+8), (float)(posY+4), Justification::centredLeft);
      }
      else
      {
        drawBitmapText(g, numberString, (int)posX-8-64, (int)posY-10, 64, 20, 
          &normalFont7px, Justification::centredRight);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX-8), (float)(posY+4), Justification::centredRight);
      }

      i++;
    }
    i = 1;
    while( -i*horizontalCoarseGridInterval > currentRange.getMinY() )
    {
      if( logScaledX )
        posX = currentRange.getMinX();
      else
      {
        if( axisPositionY == ZERO )
          posX = 0.0;
        else if( axisPositionY == LEFT ) 
          posX = currentRange.getMinX();
        else if( axisPositionY == RIGHT ) 
          posX = currentRange.getMaxX();
      }
      posY = value = -i*horizontalCoarseGridInterval; // "value" will not be transformed

      // transform coordinates:
      if( targetImage == NULL )
        toPixelCoordinates(posX, posY);
      else
        transformToImageCoordinates(posX, posY, targetImage);

      // include some margin for axes at the left and right:
      if( axisPositionY == LEFT )
        posX += 8;
      else if( axisPositionY == RIGHT )
        posX -= 8;

      // draw a small line:
      g.drawLine((float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 
        1.0);
      if( targetSVG != NULL )
        addLineToSvgDrawing(targetSVG, (float)(posX-4.0), (float)posY, (float)(posX+4.0), (float)posY, 1.0, colourScheme.axes, false);

      // draw number:
      numberString = stringConversionFunctionY(value);
      if( axisValuesPositionY == RIGHT_TO_AXIS || axisPositionY == LEFT )
      {
        drawBitmapText(g, numberString, (int)posX+8, (int)posY-10, 64, 20, 
          &normalFont7px, Justification::centredLeft);
        if( targetSVG != NULL )        
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX+8), (float)(posY+4), Justification::centredLeft);
      }
      else
      {
        drawBitmapText(g, numberString, (int)posX-8-64, (int)posY-10, 64, 20, 
          &normalFont7px, Justification::centredRight);
        if( targetSVG != NULL )
          addTextToSvgDrawing(targetSVG, numberString, 
          (float)(posX-8), (float)(posY+4), Justification::centredRight);
      }

      i++;
    }  // end while
  }
}

void CoordinateSystem::updateScaleFactors()
{
  if( !logScaledX )
    scaleX = getWidth()  / (currentRange.getMaxX()-currentRange.getMinX()); 
  // scaling factor for linear plots
  else
  {
    jassert(((currentRange.getMaxX()/currentRange.getMinX()) > 0.0));
      // caught a logarithm of a non-positive number, make sure that the 
      // minimum and the maximum for the x-coordinate are both strictly positive
      // for logarithmic axis-scaling

      if( (currentRange.getMaxX()/currentRange.getMinX()) > 0.0 )
        pixelsPerIntervalX = getWidth()/logB((currentRange.getMaxX()/currentRange.getMinX()), 
        logBaseX);
    // the number of pixels per interval for logarithmic plots - for
    // logBase==2 this is the number of pixels per octave:
      else
        pixelsPerIntervalX = 50; // some arbitrary fallback-value
  }
  if( !logScaledY )
    scaleY = getHeight() / (currentRange.getMaxY()-currentRange.getMinY());
  else
  {
    jassert(((currentRange.getMaxY()/currentRange.getMinY()) > 0.0));
      // caught a logarithm of a non-positive number, make sure that the 
      // minimum and the maximum for the y-coordinate are both strictly positive
      // for logarithmic axis-scaling

      if( (currentRange.getMaxY()/currentRange.getMinY()) > 0.0 )
        pixelsPerIntervalY = getHeight()/logB((currentRange.getMaxY()/currentRange.getMinY()), 
        logBaseY);
      else
        pixelsPerIntervalY = 50; // some arbitrary fallback-value
  }
}

double CoordinateSystem::getPlotHeight(Image *targetImage)
{
  if( targetImage == NULL )
    return getHeight();
  else
    return targetImage->getHeight();
}

double CoordinateSystem::getPlotWidth(Image *targetImage)
{
  if( targetImage == NULL )
    return getWidth();
  else
    return targetImage->getWidth();
}

//-------------------------------------------------------------------------------------------------
// state-management (storing and recall), still incomplete:

XmlElement* CoordinateSystem::getStateAsXml(const String& stateName) const
{
  XmlElement* xmlState = new XmlElement(stateName); 
  // the XmlElement which stores all the releveant state-information

  xmlState->setAttribute(String("MinX"), currentRange.getMinX());
  xmlState->setAttribute(String("MaxX"), currentRange.getMaxX());
  xmlState->setAttribute(String("MinY"), currentRange.getMinY());
  xmlState->setAttribute(String("MaxY"), currentRange.getMaxY());

  xmlState->setAttribute(String("HorizontalCoarseGridIsVisible"),    
    horizontalCoarseGridIsVisible);
  xmlState->setAttribute(String("HorizontalCoarseGridInterval"),
    horizontalCoarseGridInterval);
  xmlState->setAttribute(String("HorizontalFineGridIsVisible"),      
    horizontalFineGridIsVisible);
  xmlState->setAttribute(String("HorizontalFineGridInterval"),  
    horizontalFineGridInterval); 
  xmlState->setAttribute(String("VerticalCoarseGridIsVisible"),    
    verticalCoarseGridIsVisible);
  xmlState->setAttribute(String("VerticalCoarseGridInterval"),
    verticalCoarseGridInterval);
  xmlState->setAttribute(String("VerticalFineGridIsVisible"),      
    verticalFineGridIsVisible);
  xmlState->setAttribute(String("VerticalFineGridInterval"),  
    verticalFineGridInterval);

  return xmlState;
}

bool CoordinateSystem::setStateFromXml(const XmlElement &xmlState)
{
  bool success = true; // should report about success, not used yet

  currentRange.setMinX( xmlState.getDoubleAttribute(String("MinX"),getCurrentRangeMinX()) );
  currentRange.setMaxX( xmlState.getDoubleAttribute(String("MaxX"),getCurrentRangeMaxX()) );
  currentRange.setMinY( xmlState.getDoubleAttribute(String("MinY"),getCurrentRangeMinY()) );
  currentRange.setMaxY( xmlState.getDoubleAttribute(String("MaxY"),getCurrentRangeMaxY()) );

  horizontalCoarseGridIsVisible = xmlState.getBoolAttribute(
    String("HorizontalCoarseGridIsVisible"), isHorizontalCoarseGridVisible());
  horizontalCoarseGridInterval = xmlState.getDoubleAttribute(
    String("HorizontalCoarseGridInterval"), getHorizontalCoarseGridInterval());
  horizontalFineGridIsVisible = xmlState.getBoolAttribute(
    String("HorizontalFineGridIsVisible"), isHorizontalFineGridVisible());
  horizontalFineGridInterval = xmlState.getDoubleAttribute(
    String("HorizontalFineGridInterval"), getHorizontalFineGridInterval());
  verticalCoarseGridIsVisible = xmlState.getBoolAttribute(
    String("VerticalCoarseGridIsVisible"), isVerticalCoarseGridVisible());
  verticalCoarseGridInterval = xmlState.getDoubleAttribute(
    String("VerticalCoarseGridInterval"), getVerticalCoarseGridInterval());
  verticalFineGridIsVisible = xmlState.getBoolAttribute(
    String("VerticalFineGridIsVisible"), isVerticalFineGridVisible());
  verticalFineGridInterval = xmlState.getDoubleAttribute(
    String("VerticalFineGridInterval"), getVerticalFineGridInterval());

  //updateBackgroundImage();
  return success;
}






