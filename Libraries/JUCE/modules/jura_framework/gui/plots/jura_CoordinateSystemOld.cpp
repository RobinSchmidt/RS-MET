CoordinateSystemRangeOld::CoordinateSystemRangeOld(double initMinX, double initMaxX,
  double initMinY, double initMaxY)
{
  minX = initMinX;
  maxX = initMaxX;
  minY = initMinY;
  maxY = initMaxY;
}

CoordinateSystemRangeOld::~CoordinateSystemRangeOld()
{

}

void CoordinateSystemRangeOld::setMinX(double newMinX) 
{ 
  jassert(newMinX < maxX);
    if( newMinX < maxX )
      minX = newMinX;
}

void CoordinateSystemRangeOld::setMaxX(double newMaxX) 
{ 
  jassert(newMaxX > minX);
    if( newMaxX > minX )
      maxX = newMaxX;
}

void CoordinateSystemRangeOld::setRangeX(double newMinX, double newMaxX)
{
  jassert(newMaxX > newMinX);
    if( newMaxX > newMinX )
    {
      minX = newMinX;
      maxX = newMaxX;
    }
}

void CoordinateSystemRangeOld::setMinY(double newMinY) 
{ 
  jassert(newMinY < maxY);
    if( newMinY < maxY )
      minY = newMinY;
}

void CoordinateSystemRangeOld::setMaxY(double newMaxY) 
{ 
  jassert(newMaxY > minY);
    if( newMaxY > minY )
      maxY = newMaxY;
}

void CoordinateSystemRangeOld::setRangeY(double newMinY, double newMaxY)
{
  jassert(newMaxY > newMinY);
    if( newMaxY > newMinY )
    {
      minY = newMinY;
      maxY = newMaxY;
    }
}

void CoordinateSystemRangeOld::clipRange(CoordinateSystemRangeOld rangeToClipTo)
{
  if( minX < rangeToClipTo.getMinX() )
    minX = rangeToClipTo.getMinX();
  if( maxX > rangeToClipTo.getMaxX() )
    maxX = rangeToClipTo.getMaxX();
  if( minY < rangeToClipTo.getMinY() )
    minY = rangeToClipTo.getMinY();
  if( maxY > rangeToClipTo.getMaxY() )
    maxY = rangeToClipTo.getMaxY();
}

//=================================================================================================

// construction/destruction:

CoordinateSystemOld::CoordinateSystemOld(const String &newDescription) 
  : DescribedComponent(newDescription) //RWidget(newDescription)
{
  autoReRenderImage             = false;
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
  stringConversionForAxisX     = &valueToString0;
  stringConversionForAxisY     = &valueToString0;
  stringConversionForInfoLineX = &valueToString0;
  stringConversionForInfoLineY = &valueToString0;

  useBitmapFont = true;

  // initialize the component-size and the image-size to 1x1 pixels, without
  // such initializations, a JUCE-breakpoint will be triggered or other screws
  // happen:
  backgroundImage = NULL;
  backgroundImage = new Image(Image::RGB, 1, 1, true);
  //backgroundImage->duplicateIfShared(); // nah - call this on the copies
  setBounds(0, 0, 1, 1);
  updateBackgroundImage();

  autoReRenderImage  = true;

  // use a crosshair-cursor to precisely point to some coordinate for inspection:
  //currentMouseCursor = MouseCursor(MouseCursor::CrosshairCursor);
  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);
  showPositionAsDescription = false;
  showPopUpOnRightClick     = false;

  /*
  // the Label for inspecting the current frequency ane level:
  inspectionField = new TextEditor( String(T("inspectionField")) );
  inspectionField->setBounds(4, 4, 120, 40);
  inspectionField->setColour(TextEditor::backgroundColourId, Colours::white.withAlpha(0.7f) );
  inspectionField->setColour(TextEditor::outlineColourId,    Colours::black );
  inspectionField->setCaretVisible(false);
  inspectionField->setScrollbarsShown(false);
  inspectionField->setReadOnly(true);
  inspectionField->setMultiLine(true);
  addChildComponent( inspectionField );
  */
}

CoordinateSystemOld::~CoordinateSystemOld()
{
  deleteAllChildren();
  if( backgroundImage != NULL )
    delete backgroundImage;
}

//-------------------------------------------------------------------------------------------------
// component-overrides:

void CoordinateSystemOld::mouseDown(const MouseEvent &e)
{
  if( e.mods.isRightButtonDown() && showPopUpOnRightClick == true )
    openRightClickPopupMenu();
}

void CoordinateSystemOld::mouseEnter(const MouseEvent &e)
{
  if( showPositionAsDescription == true )
    setDescription( getInfoLineForPixelPosition(e.x, e.y) );
  DescribedComponent::mouseEnter(e);
}

void CoordinateSystemOld::mouseMove(const MouseEvent &e)
{
  if( showPositionAsDescription == true )
    setDescription( getInfoLineForPixelPosition(e.x, e.y) );
  DescribedComponent::mouseMove(e);
}

void CoordinateSystemOld::resized()
{
  updateMapping();
}

void CoordinateSystemOld::paint(juce::Graphics &g)
{
  if( backgroundImage != NULL )
    g.drawImage(*backgroundImage, 0, 0, getWidth(), getHeight(), 
      0, 0, backgroundImage->getWidth(), backgroundImage->getHeight(), false);
  else
    g.fillAll(Colours::red);
}

//-------------------------------------------------------------------------------------------------
// range management:

void CoordinateSystemOld::setMaximumRange(double newMinX, double newMaxX, 
                                       double newMinY, double newMaxY)
{
  maximumRange.setRangeX(newMinX, newMaxX);
  maximumRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRange(CoordinateSystemRangeOld newMaximumRange)
{
  maximumRange = newMaximumRange;
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeX(double newMinX, double newMaxX)
{
  maximumRange.setRangeX(newMinX, newMaxX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeY(double newMinY, double newMaxY)
{
  maximumRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeMinX(double newMinX)
{
  maximumRange.setMinX(newMinX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeMaxX(double newMaxX)
{
  maximumRange.setMaxX(newMaxX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeMinY(double newMinY)
{
  maximumRange.setMinY(newMinY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setMaximumRangeMaxY(double newMaxY)
{
  maximumRange.setMaxY(newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRange(double newMinX, double newMaxX, 
                                       double newMinY, double newMaxY)
{
  currentRange.setRangeX(newMinX, newMaxX);
  currentRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRange(CoordinateSystemRangeOld newRange)
{
  currentRange = newRange;
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeX(double newMinX, double newMaxX)
{
  currentRange.setRangeX(newMinX, newMaxX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeY(double newMinY, double newMaxY)
{
  currentRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeMinX(double newMinX)
{
  currentRange.setMinX(newMinX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeMaxX(double newMaxX)
{
  currentRange.setMaxX(newMaxX);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeMinY(double newMinY)
{
  currentRange.setMinY(newMinY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

void CoordinateSystemOld::setCurrentRangeMaxY(double newMaxY)
{
  currentRange.setMaxY(newMaxY);
  currentRange.clipRange(maximumRange);
  updateMapping();
}

String CoordinateSystemOld::getInfoLineForPixelPosition(int x, int y)
{
  double xd = (double) x;
  double yd = (double) y;
  transformFromComponentsCoordinates(xd, yd);
  String xString = stringConversionForInfoLineX(xd);
  String yString = stringConversionForInfoLineY(yd);
  return xString + String(", ") + yString;
}

//-------------------------------------------------------------------------------------------------
// appearance:

void CoordinateSystemOld::setColourScheme(const PlotColourScheme& newColourScheme)
{
  plotColourScheme = newColourScheme;
  updateBackgroundImage();
}

void CoordinateSystemOld::setColourSchemeFromXml(const XmlElement &xml)
{
  //colourScheme.setColourSchemeFromXml(xml);
}

void CoordinateSystemOld::setAutoReRendering(bool shouldAutomaticallyReRender)
{
  autoReRenderImage = shouldAutomaticallyReRender;
}

void CoordinateSystemOld::setCaption(const String &newCaption, int newPosition)
{
  captionPosition = newPosition;
  captionString   = newCaption;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabels(const String &newLabelX, const String &newLabelY,              
  int newLabelPositionX, int newLabelPositionY)
{
  axisLabelX         = newLabelX;
  axisLabelY         = newLabelY;
  axisLabelPositionX = newLabelPositionX;
  axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabelX(const String& newLabelX, int newLabelPositionX)
{
  axisLabelX         = newLabelX;
  axisLabelPositionX = newLabelPositionX;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabelY(const String& newLabelY, int newLabelPositionY)
{
  axisLabelY             = newLabelY;
  axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisValuesPositionX(int newValuesPositionX)
{
  if( newValuesPositionX == NO_ANNOTATION ||
    newValuesPositionX == BELOW_AXIS    ||
    newValuesPositionX == ABOVE_AXIS      )
  {
    axisValuesPositionX = newValuesPositionX;
  }
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisValuesPositionY(int newValuesPositionY)
{
  if( newValuesPositionY == NO_ANNOTATION ||
    newValuesPositionY == LEFT_TO_AXIS  ||
    newValuesPositionY == RIGHT_TO_AXIS   )
  {
    axisValuesPositionY = newValuesPositionY;
  }
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setStringConversionForAxisX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  stringConversionForAxisX = newConversionFunctionX;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setStringConversionForInfoLineX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  stringConversionForInfoLineX = newConversionFunctionX;
}

void CoordinateSystemOld::setStringConversionForAxisY(
  String (*newConversionFunctionY) (double valueToBeConverted) )
{
  stringConversionForAxisY = newConversionFunctionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setStringConversionForInfoLineY(
  String (*newConversionFunctionY) (double valueToBeConverted) )
{
  stringConversionForInfoLineY = newConversionFunctionY;
}

void CoordinateSystemOld::setHorizontalCoarseGridVisible(bool shouldBeVisible)
{
  setHorizontalCoarseGrid(horizontalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalCoarseGridInterval(double newGridInterval)
{
  setHorizontalCoarseGrid(newGridInterval, horizontalCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalCoarseGrid(double newGridInterval, bool   shouldBeVisible)
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
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGridVisible(bool shouldBeVisible)
{
  setHorizontalFineGrid(horizontalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGridInterval(double newGridInterval)
{
  setHorizontalFineGrid(newGridInterval, horizontalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGrid(double newGridInterval, 
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
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGridVisible(bool shouldBeVisible)
{
  setVerticalCoarseGrid(verticalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGridInterval(double newGridInterval)
{
  setVerticalCoarseGrid(newGridInterval, verticalCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGrid(double newGridInterval, 
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
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGridVisible(bool shouldBeVisible)
{
  setVerticalFineGrid(verticalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGridInterval(double newGridInterval)
{
  setVerticalFineGrid(newGridInterval, verticalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGrid(double newGridInterval, 
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
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGridVisible(bool shouldBeVisible)
{
  setRadialCoarseGrid(radialCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGridInterval(double newGridInterval)
{
  setRadialCoarseGrid(newGridInterval, radialCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGrid(double newGridInterval, 
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

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGridVisible(bool shouldBeVisible)
{
  setRadialFineGrid(radialFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGridInterval(double newGridInterval)
{
  setRadialFineGrid(newGridInterval, radialFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGrid(double newGridInterval, 
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
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularCoarseGridVisible(bool shouldBeVisible)
{
  setAngularCoarseGrid(angularCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularCoarseGridInterval(double newGridInterval)
{
  setAngularCoarseGrid(newGridInterval, angularCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularCoarseGrid(double newGridInterval, 
                                            bool   shouldBeVisible)
{
  jassert(newGridInterval > 0.000001); 
  // grid-intervals must be > 0
  if( newGridInterval <= 0.000001 )
    return;

  angularCoarseGridIsVisible = shouldBeVisible;
  angularCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularFineGridVisible(bool shouldBeVisible)
{
  setAngularFineGrid(angularFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularFineGridInterval(double newGridInterval)
{
  setAngularFineGrid(newGridInterval, angularFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularFineGrid(double newGridInterval, 
                                          bool   shouldBeVisible)
{
  jassert(newGridInterval > 0.000001); 
  // grid-intervals must be > 0
  if( newGridInterval <= 0.000001 )
    return;

  angularFineGridIsVisible     = shouldBeVisible;
  angularFineGridInterval = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

bool CoordinateSystemOld::isHorizontalCoarseGridVisible() 
{ 
  return horizontalCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isHorizontalFineGridVisible() 
{ 
  return horizontalFineGridIsVisible; 
}

bool CoordinateSystemOld::isVerticalCoarseGridVisible() 
{ 
  return verticalCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isVerticalFineGridVisible() 
{ 
  return verticalFineGridIsVisible; 
}

bool CoordinateSystemOld::isRadialCoarseGridVisible() 
{ 
  return radialCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isRadialFineGridVisible() 
{ 
  return radialFineGridIsVisible; 
}

bool CoordinateSystemOld::isAngularCoarseGridVisible() 
{ 
  return angularCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isAngularFineGridVisible() 
{ 
  return angularFineGridIsVisible; 
}

double CoordinateSystemOld::getHorizontalCoarseGridInterval()
{
  return horizontalCoarseGridInterval;
}

double CoordinateSystemOld::getHorizontalFineGridInterval()
{
  return horizontalFineGridInterval;
}

double CoordinateSystemOld::getVerticalCoarseGridInterval()
{
  return verticalCoarseGridInterval;
}

double CoordinateSystemOld::getVerticalFineGridInterval()
{
  return verticalFineGridInterval;
}

double CoordinateSystemOld::getRadialCoarseGridInterval()
{
  return radialCoarseGridInterval;
}

double CoordinateSystemOld::getRadialFineGridInterval()
{
  return radialFineGridInterval;
}

double CoordinateSystemOld::getAngularCoarseGridInterval()
{
  return angularCoarseGridInterval;
}

double CoordinateSystemOld::getAngularFineGridInterval()
{
  return angularFineGridInterval;
}

void CoordinateSystemOld::setAxisPositionX(int newAxisPositionX)
{
  if( newAxisPositionX == INVISIBLE ||
    newAxisPositionX == ZERO      ||
    newAxisPositionX == TOP       ||
    newAxisPositionX == BOTTOM       )
  {
    axisPositionX = newAxisPositionX;
    if(autoReRenderImage == true)
      updateBackgroundImage();
  }
}

void CoordinateSystemOld::setAxisPositionY(int newAxisPositionY)
{
  if( newAxisPositionY == INVISIBLE ||
    newAxisPositionY == ZERO      ||
    newAxisPositionY == LEFT      ||
    newAxisPositionY == RIGHT        )
  {
    axisPositionY = newAxisPositionY;
    if(autoReRenderImage == true)
      updateBackgroundImage();
  }
}

void CoordinateSystemOld::setupAxisX(double newMin, double newMax, bool shouldBeLogScaled, 
  double newLogBase, int newAxisPosition, double newCoarseGridInterval, 
  double newFineGridInterval)
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
  if( newAxisPosition == INVISIBLE ||
      newAxisPosition == ZERO      ||
      newAxisPosition == TOP       ||
      newAxisPosition == BOTTOM       )
  {
    axisPositionX = newAxisPosition;
  }
  verticalCoarseGridInterval = newCoarseGridInterval;
  verticalFineGridInterval   = newFineGridInterval;
  updateMapping();
}

void CoordinateSystemOld::setupAxisY(double newMin, double newMax, bool shouldBeLogScaled, 
  double newLogBase, int newAxisPosition, double newCoarseGridInterval, 
  double newFineGridInterval)
{
  // axis settings seem not to make sense
  jassert(newMin < newMax);
  jassert(newMin > 0.0 || shouldBeLogScaled == false);
  jassert(shouldBeLogScaled == false ||
    (newCoarseGridInterval > 1.000001 &&
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
  updateMapping();
}

void CoordinateSystemOld::useLogarithmicScale(bool shouldBeLogScaledX, bool shouldBeLogScaledY,                      
  double newLogBaseX, double newLogBaseY)
{
  logScaledX = shouldBeLogScaledX;
  logScaledY = shouldBeLogScaledY;
  logBaseX     = newLogBaseX;
  logBaseY     = newLogBaseY;
  updateMapping();
}

void CoordinateSystemOld::useLogarithmicScaleX(bool   shouldBeLogScaledX, 
                                            double newLogBaseX)
{
  logScaledX = shouldBeLogScaledX;
  logBaseX     = newLogBaseX;
  updateMapping();
}

bool CoordinateSystemOld::isLogScaledX()
{
  return logScaledX;
}

void CoordinateSystemOld::useLogarithmicScaleY(bool   shouldBeLogScaledY, 
                                            double newLogBaseY)
{
  logScaledY = shouldBeLogScaledY;
  logBaseY     = newLogBaseY;
  updateMapping();
}

bool CoordinateSystemOld::isLogScaledY()
{
  return logScaledY;
}

/*
void CoordinateSystemOld::setValueFieldPopup(bool shouldPopUp)
{
  valuePopup = shouldPopUp;
  inspectionField->setVisible(false);
}
*/

//-------------------------------------------------------------------------------------------------
// functions for drawing and/or exporting the shown content:

void CoordinateSystemOld::openRightClickPopupMenu()
{  
  //DEBUG_BREAK;
  // the code below for opening the context menu is outdated - change it to deal with the new RPopUpMenu

  /*
  // create a context menu to allow for export:
  RPopUpMenuOld menu;
  menu.setColourScheme(popUpColourScheme);  

  menu.addItem(1, "Export Image");
  const int result = menu.show();

  if (result == 0)
  {
    // user dismissed the menu without picking anything
  }
  else if (result == 1)
  {
    // user picked the Export item - open the export dialog window:
    openExportDialog(getWidth(), getHeight(), String(T("png")), File::nonexistent);
  }
  */
}

Image* CoordinateSystemOld::getPlotAsImage(int width, int height)
{
  jassert(width  >= 1);
  jassert(height >= 1); 
  if( width < 1 || height < 1)  
    return nullptr;
  Image* thePlotImage = new Image(Image::RGB, width, height, true);
  Graphics g(*thePlotImage);
  drawCoordinateSystem(g, thePlotImage);
  return thePlotImage;
}

// this needs to be refactored, maybe into a class rsPlotDrawerSvg:
XmlElement* CoordinateSystemOld::getPlotAsSVG(int width, int height)
{
  jassert(width  >= 1);
  jassert(height >= 1);
  if( width < 1 || height < 1)
    return nullptr;

  // we need dummy Image and Graphics objects:
  Image* thePlotImage = new Image(Image::RGB, width, height, true);
  Graphics g(*thePlotImage);

  // create an XmlElement to be used for the SVG drawing:
  XmlElement* theSVG = new XmlElement(String("svg"));
  theSVG->setAttribute(String("width"), width);
  theSVG->setAttribute(String("height"), height);
  drawCoordinateSystem(g, thePlotImage, theSVG);  // draw on the SVG
  delete thePlotImage;                            // delete the dummy image
  return theSVG;
}

void CoordinateSystemOld::openExportDialog(int defaultWidth, int defaultHeight, 
                                        const String &defaultFormat,
                                        const File& defaultTargetFile)
{
  ImageSavingDialog dialog(this, defaultWidth, defaultHeight, defaultFormat, defaultTargetFile);
  DialogWindow exportWindow(String("Export to Image or SVG Drawing"), Colours::white, true, true);
  exportWindow.showModalDialog(String("Export to Image or SVG Drawing"), &dialog, this, 
    Colours::white, true, false, false);
}

void CoordinateSystemOld::updateBackgroundImage()
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
  //repaint();
  repaintOnMessageThread();
}

//-------------------------------------------------------------------------------------------------
// coordinate transformations:

inline double logB(double x, double b)
{
  // temporary - remove and use the function from RAPT
  return log(x) / log(b);
}

void CoordinateSystemOld::transformToImageCoordinates(double &x, double &y, const Image* theImage)
{
  if( theImage == NULL )
  {
    transformToComponentsCoordinates(x, y);
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

void CoordinateSystemOld::transformFromImageCoordinates(double &x, double &y, const Image *theImage)
{
  if( theImage == NULL )
  {
    transformFromComponentsCoordinates(x, y);
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

void CoordinateSystemOld::transformToComponentsCoordinates(double &x, double &y)
{
  x = coordinateMapper.mapX(x);
  y = coordinateMapper.mapY(y);
}

void CoordinateSystemOld::transformToComponentsCoordinates(float &x, float &y)
{
  x = (float) coordinateMapper.mapX(x);
  y = (float) coordinateMapper.mapY(y);
}

void CoordinateSystemOld::transformFromComponentsCoordinates(double &x, double &y)
{
  x = coordinateMapper.unmapX(x);
  y = coordinateMapper.unmapY(y);
}

void CoordinateSystemOld::transformFromComponentsCoordinates(float &x, float &y)
{
  x = (float) coordinateMapper.unmapX(x);
  y = (float) coordinateMapper.unmapY(y);
}

void CoordinateSystemOld::updateMapping()
{
  updateScaleFactors();
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::updateCoordinateMapperOutputRange(Image* image, XmlElement* svg)
{
  double w = 0, h = 0;
  if(image != nullptr) {
    w = image->getWidth();
    h = image->getHeight(); }
  else if(svg != nullptr) {
    w = svg->getDoubleAttribute("width", 0);
    h = svg->getDoubleAttribute("height", 0);
    jassert(w != 0 && h != 0); } // svg must have width and height attributes
  else {
    w = getWidth();
    h = getHeight(); }

  //coordinateMapper.setOutputRange(0, w-1, h-1, 0);
  coordinateMapper.setOutputRange(0.5, w-0.5, h-0.5, 0.5);
}

//-------------------------------------------------------------------------------------------------
// drawing functions

void CoordinateSystemOld::drawCoordinateSystem(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  g.setFont(Font(14));

  updateCoordinateMapperOutputRange(targetImage, targetSVG);

  // draw the background:
  //g.fillAll(colourScheme.backgroundColour);
  if( targetImage != NULL )
  {
    Graphics g2(*targetImage);
    fillRectWithBilinearGradient(g2, 0, 0, targetImage->getWidth(), targetImage->getHeight(),
      plotColourScheme.topLeft, plotColourScheme.topRight, plotColourScheme.bottomLeft, plotColourScheme.bottomRight);
  }
  else
  {
    fillRectWithBilinearGradient(g, 0, 0, getWidth(), getHeight(),
      plotColourScheme.topLeft, plotColourScheme.topRight, plotColourScheme.bottomLeft, plotColourScheme.bottomRight);
  }

  // draw the grids, if desired:
  if( horizontalFineGridIsVisible )
    drawHorizontalGrid(g, horizontalFineGridInterval, logScaledY, plotColourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( verticalFineGridIsVisible )
    drawVerticalGrid(g, verticalFineGridInterval, logScaledX, plotColourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( radialFineGridIsVisible )
    drawRadialGrid(g, radialFineGridInterval, logScaledRadius, plotColourScheme.fineGrid, 1.0f, targetImage, targetSVG);
  if( angularFineGridIsVisible )
    drawAngularGrid(g, angularFineGridInterval, plotColourScheme.fineGrid, 1.0f, targetImage, targetSVG);

  if( horizontalCoarseGridIsVisible )
    drawHorizontalGrid(g, horizontalCoarseGridInterval, logScaledY, plotColourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( verticalCoarseGridIsVisible )
    drawVerticalGrid(g, verticalCoarseGridInterval, logScaledX, plotColourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( radialCoarseGridIsVisible )
    drawRadialGrid(g, radialCoarseGridInterval, logScaledRadius, plotColourScheme.coarseGrid, 1.0f, targetImage, targetSVG);
  if( angularCoarseGridIsVisible )
    drawAngularGrid(g, angularCoarseGridInterval, plotColourScheme.coarseGrid, 1.0f, targetImage, targetSVG);

  // draw the coordinate system:
  if( axisPositionX != INVISIBLE )
    drawAxisX(g, targetImage, targetSVG);
  if( axisPositionY != INVISIBLE )
    drawAxisY(g, targetImage, targetSVG);

  // draw the values on the axes:
  if( axisPositionX != INVISIBLE && axisValuesPositionX != NO_ANNOTATION )
    drawAxisValuesX(g, targetImage, targetSVG);
  if( axisPositionY != INVISIBLE && axisValuesPositionY != NO_ANNOTATION )
    drawAxisValuesY(g, targetImage, targetSVG);

  // draw the caption:
  drawCaption(g, targetImage, targetSVG);

  // draw an outlining rectangle:
  g.setColour(plotColourScheme.outline);
  g.drawRect(0, 0, getWidth(), getHeight(), 2);
}

void CoordinateSystemOld::drawCaption(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  // needs test:
  static const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  float w = (float) font->getTextPixelWidth(captionString);
  float h = (float) font->getFontHeight();
  switch( captionPosition )
  {
  case NO_CAPTION: return;
  case TOP_CENTER:
    {
      drawBitmapText(g, captionString, 0.5f*getWidth()-0.5f*w, 16, getWidth(), 16, font, 
        Justification::centred, plotColourScheme.text);
    }
    break;
  case CENTER:
    {
      float x = (float) (getWidth()  - w) * .5f;
      float y = (float) (getHeight() - h) * .5f;
      drawBitmapText(g, captionString, x, y, getWidth(), 16, font, Justification::topLeft,
        plotColourScheme.text);
    }
    break;
  }
}

void CoordinateSystemOld::drawHorizontalGrid(Graphics &g, double interval, bool exponentialSpacing, 
  Colour gridColour, float lineThickness, Image* targetImage, XmlElement *targetSVG)
{
  g.setColour(gridColour);
  if(targetSVG != nullptr)
    jura::drawHorizontalGrid(targetSVG, coordinateMapper, interval, lineThickness, gridColour);
  else
    jura::drawHorizontalGrid(g, coordinateMapper, interval, lineThickness);
}

void CoordinateSystemOld::drawVerticalGrid(Graphics &g, double interval, bool exponentialSpacing, 
  Colour gridColour, float lineThickness, Image* targetImage, XmlElement *targetSVG)
{

  g.setColour(gridColour);
  if(targetSVG != nullptr)
    jura::drawVerticalGrid(targetSVG, coordinateMapper, interval, lineThickness, gridColour);
  else
    jura::drawVerticalGrid(g, coordinateMapper, interval, lineThickness);
}

void CoordinateSystemOld::drawRadialGrid(Graphics &g, double interval, bool exponentialSpacing,
  Colour gridColour, float lineThickness, Image* targetImage, XmlElement *targetSVG)
{
  g.setColour(gridColour);
  if(targetSVG != nullptr)
    jura::drawRadialGrid(targetSVG, coordinateMapper, interval, lineThickness, gridColour);
  else
    jura::drawRadialGrid(g, coordinateMapper, interval, lineThickness);
}

void CoordinateSystemOld::drawAngularGrid(Graphics &g, double interval, Colour gridColour,                  
  float lineThickness, Image* targetImage, XmlElement *targetSVG)
{
  g.setColour(gridColour);
  if(targetSVG != nullptr)
    jura::drawAngularGrid(targetSVG, coordinateMapper, interval, lineThickness, gridColour);
  else
    jura::drawAngularGrid(g, coordinateMapper, interval, lineThickness);
}

void CoordinateSystemOld::drawAxisX(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisPositionX == INVISIBLE ) 
    return;

  double yPos = 0.0;
  if( axisPositionX == BOTTOM )
    yPos = coordinateMapper.unmapY(getHeight()-8);
  else if( axisPositionX == TOP )
    yPos = coordinateMapper.unmapY(8);

  if(targetSVG != nullptr)
    jura::drawAxisX(targetSVG, coordinateMapper, yPos, axisLabelX, plotColourScheme.axes);
  else 
  {
    g.setColour(plotColourScheme.axes);
    jura::drawAxisX(g, coordinateMapper, yPos, axisLabelX, plotColourScheme.text);
  }
}

void CoordinateSystemOld::drawAxisY(juce::Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisPositionX == INVISIBLE ) 
    return;

  double xPos = 0.0;
  if( axisPositionY == LEFT )
    xPos = coordinateMapper.unmapX(8);
  else if( axisPositionY == RIGHT )
    xPos = coordinateMapper.unmapX(getWidth()-8);

  if(targetSVG != nullptr)
  {
    jura::drawAxisY(targetSVG, coordinateMapper, xPos, axisLabelY, plotColourScheme.axes);
  }
  else 
  {
    g.setColour(plotColourScheme.axes);
    jura::drawAxisY(g, coordinateMapper, xPos, axisLabelY, plotColourScheme.text);
  }
}

void CoordinateSystemOld::drawAxisValuesX(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisValuesPositionX == NO_ANNOTATION )
    return;

  double yPos = 0.0;
  if( axisPositionX == BOTTOM )
    yPos = coordinateMapper.unmapY(getHeight()-8);
  else if( axisPositionX == TOP )
    yPos = coordinateMapper.unmapY(8);
  double spacing = verticalCoarseGridInterval;

  if(targetSVG != nullptr)
    jura::drawAxisValuesX(targetSVG, coordinateMapper, spacing, yPos, stringConversionForAxisX,
      plotColourScheme.axes);
  else {
    g.setColour(plotColourScheme.axes);
    jura::drawAxisValuesX(g, coordinateMapper, spacing, yPos, stringConversionForAxisX,
      plotColourScheme.text);
  }
}

void CoordinateSystemOld::drawAxisValuesY(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( axisValuesPositionY == NO_ANNOTATION )
    return;

  double xPos = 0.0; // allow for left or right, too
  if( axisPositionY == LEFT )
    xPos = coordinateMapper.unmapX(8);
  else if( axisPositionY == RIGHT )
    xPos = coordinateMapper.unmapX(getWidth()-8);
  double spacing = horizontalCoarseGridInterval;

  if(targetSVG != nullptr){
    jura::drawAxisValuesY(targetSVG, coordinateMapper, spacing, xPos, stringConversionForAxisY,
      plotColourScheme.axes);
  }
  else {
    g.setColour(plotColourScheme.axes);
    jura::drawAxisValuesY(g, coordinateMapper, spacing, xPos, stringConversionForAxisY,
      plotColourScheme.text);
  }
}

void CoordinateSystemOld::updateScaleFactors()
{
  coordinateMapper.mapperX.setLogScaled(logScaledX);  // logScaledX/Y are redundant here now
  coordinateMapper.mapperY.setLogScaled(logScaledY);  // ...get rid of them
  coordinateMapper.setInputRange(currentRange.getMinX(), currentRange.getMaxX(),
    currentRange.getMinY(), currentRange.getMaxY());
    // the currentRange member is actually also redundant now
}

// get rid of these:
double CoordinateSystemOld::getPlotHeight(Image *targetImage)
{
  if( targetImage == NULL )
    return getHeight();
  else
    return targetImage->getHeight();
}

double CoordinateSystemOld::getPlotWidth(Image *targetImage)
{
  if( targetImage == NULL )
    return getWidth();
  else
    return targetImage->getWidth();
}

//-------------------------------------------------------------------------------------------------
// state-management (storing and recall), still incomplete:

XmlElement* CoordinateSystemOld::getStateAsXml(const String& stateName) const
{
  XmlElement* xml = new XmlElement(stateName); 
  // the XmlElement which stores all the releveant state-information

  xml->setAttribute(String("MinX"), currentRange.getMinX());
  xml->setAttribute(String("MaxX"), currentRange.getMaxX());
  xml->setAttribute(String("MinY"), currentRange.getMinY());
  xml->setAttribute(String("MaxY"), currentRange.getMaxY());

  xml->setAttribute("HorizontalCoarseGridIsVisible", horizontalCoarseGridIsVisible);
  xml->setAttribute("HorizontalCoarseGridInterval",  horizontalCoarseGridInterval);
  xml->setAttribute("HorizontalFineGridIsVisible",   horizontalFineGridIsVisible);
  xml->setAttribute("HorizontalFineGridInterval",    horizontalFineGridInterval); 
  xml->setAttribute("VerticalCoarseGridIsVisible",   verticalCoarseGridIsVisible);
  xml->setAttribute("VerticalCoarseGridInterval",    verticalCoarseGridInterval);
  xml->setAttribute("VerticalFineGridIsVisible",     verticalFineGridIsVisible);
  xml->setAttribute("VerticalFineGridInterval",      verticalFineGridInterval);

  return xml;
}

bool CoordinateSystemOld::setStateFromXml(const XmlElement &xml)
{
  bool success = true; // should report about success, not used yet

  currentRange.setMinX( xml.getDoubleAttribute("MinX", getCurrentRangeMinX()) );
  currentRange.setMaxX( xml.getDoubleAttribute("MaxX", getCurrentRangeMaxX()) );
  currentRange.setMinY( xml.getDoubleAttribute("MinY", getCurrentRangeMinY()) );
  currentRange.setMaxY( xml.getDoubleAttribute("MaxY", getCurrentRangeMaxY()) );

  horizontalCoarseGridIsVisible = xml.getBoolAttribute(  "HorizontalCoarseGridIsVisible", isHorizontalCoarseGridVisible());
  horizontalCoarseGridInterval  = xml.getDoubleAttribute("HorizontalCoarseGridInterval", getHorizontalCoarseGridInterval());
  horizontalFineGridIsVisible   = xml.getBoolAttribute(  "HorizontalFineGridIsVisible", isHorizontalFineGridVisible());
  horizontalFineGridInterval    = xml.getDoubleAttribute("HorizontalFineGridInterval", getHorizontalFineGridInterval());
  verticalCoarseGridIsVisible   = xml.getBoolAttribute(  "VerticalCoarseGridIsVisible", isVerticalCoarseGridVisible());
  verticalCoarseGridInterval    = xml.getDoubleAttribute("VerticalCoarseGridInterval", getVerticalCoarseGridInterval());
  verticalFineGridIsVisible     = xml.getBoolAttribute(  "VerticalFineGridIsVisible", isVerticalFineGridVisible());
  verticalFineGridInterval      = xml.getDoubleAttribute("VerticalFineGridInterval", getVerticalFineGridInterval());

  updateBackgroundImage();
  return success;
}
