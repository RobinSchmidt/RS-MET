// construction/destruction:

CoordinateSystemOld::CoordinateSystemOld(const String &newDescription) 
  : DescribedComponent(newDescription) //RWidget(newDescription)
{
  autoReRenderImage = false;

  // initialize the function-pointers for value->string conversion
  stringConversionForInfoLineX = &valueToString0;
  stringConversionForInfoLineY = &valueToString0;

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
  updateMapperOutputRange();
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
  plotSettings.maximumRange.setRangeX(newMinX, newMaxX);
  plotSettings.maximumRange.setRangeY(newMinY, newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRange(rsPlotRange newMaximumRange)
{
  plotSettings.maximumRange = newMaximumRange;
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeX(double newMinX, double newMaxX)
{
  plotSettings.maximumRange.setRangeX(newMinX, newMaxX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeY(double newMinY, double newMaxY)
{
  plotSettings.maximumRange.setRangeY(newMinY, newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeMinX(double newMinX)
{
  plotSettings.maximumRange.setMinX(newMinX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeMaxX(double newMaxX)
{
  plotSettings.maximumRange.setMaxX(newMaxX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeMinY(double newMinY)
{
  plotSettings.maximumRange.setMinY(newMinY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setMaximumRangeMaxY(double newMaxY)
{
  plotSettings.maximumRange.setMaxY(newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRange(double newMinX, double newMaxX, 
                                       double newMinY, double newMaxY)
{
  plotSettings.currentRange.setRangeX(newMinX, newMaxX);
  plotSettings.currentRange.setRangeY(newMinY, newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRange(rsPlotRange newRange)
{
  plotSettings.currentRange = newRange;
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeX(double newMinX, double newMaxX)
{
  plotSettings.currentRange.setRangeX(newMinX, newMaxX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeY(double newMinY, double newMaxY)
{
  plotSettings.currentRange.setRangeY(newMinY, newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeMinX(double newMinX)
{
  plotSettings.currentRange.setMinX(newMinX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeMaxX(double newMaxX)
{
  plotSettings.currentRange.setMaxX(newMaxX);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeMinY(double newMinY)
{
  plotSettings.currentRange.setMinY(newMinY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
}

void CoordinateSystemOld::setCurrentRangeMaxY(double newMaxY)
{
  plotSettings.currentRange.setMaxY(newMaxY);
  plotSettings.currentRange.clipRange(plotSettings.maximumRange);
  updateMapperInputRange();
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
  plotSettings.captionPosition = newPosition;
  plotSettings.captionString   = newCaption;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabels(const String &newLabelX, const String &newLabelY,              
  int newLabelPositionX, int newLabelPositionY)
{
  plotSettings.axisLabelX         = newLabelX;
  plotSettings.axisLabelY         = newLabelY;
  plotSettings.axisLabelPositionX = newLabelPositionX;
  plotSettings.axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabelX(const String& newLabelX, int newLabelPositionX)
{
  plotSettings.axisLabelX         = newLabelX;
  plotSettings.axisLabelPositionX = newLabelPositionX;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisLabelY(const String& newLabelY, int newLabelPositionY)
{
  plotSettings.axisLabelY             = newLabelY;
  plotSettings.axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisValuesPositionX(int newValuesPositionX)
{
  if( newValuesPositionX == rsPlotSettings::NO_ANNOTATION ||
    newValuesPositionX == rsPlotSettings::BELOW_AXIS    ||
    newValuesPositionX == rsPlotSettings::ABOVE_AXIS      )
  {
    plotSettings.axisValuesPositionX = newValuesPositionX;
  }
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAxisValuesPositionY(int newValuesPositionY)
{
  if( newValuesPositionY == rsPlotSettings::NO_ANNOTATION ||
    newValuesPositionY == rsPlotSettings::LEFT_TO_AXIS  ||
    newValuesPositionY == rsPlotSettings::RIGHT_TO_AXIS   )
  {
    plotSettings.axisValuesPositionY = newValuesPositionY;
  }
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setStringConversionForAxisX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  plotSettings.stringConversionForAxisX = newConversionFunctionX;
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
  plotSettings.stringConversionForAxisY = newConversionFunctionY;
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
  setHorizontalCoarseGrid(plotSettings.horizontalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalCoarseGridInterval(double newGridInterval)
{
  setHorizontalCoarseGrid(newGridInterval, plotSettings.horizontalCoarseGridIsVisible);
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
  if( plotSettings.logScaledY )
  {
    jassert(newGridInterval > 1.00001);
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.horizontalCoarseGridInterval = 2.0;
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

  plotSettings.horizontalCoarseGridIsVisible = shouldBeVisible;
  plotSettings.horizontalCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGridVisible(bool shouldBeVisible)
{
  setHorizontalFineGrid(plotSettings.horizontalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGridInterval(double newGridInterval)
{
  setHorizontalFineGrid(newGridInterval, plotSettings.horizontalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setHorizontalFineGrid(double newGridInterval, 
                                             bool   shouldBeVisible)
{
  if( plotSettings.logScaledY )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.horizontalFineGridInterval = pow(2.0, 1.0/3.0);
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

  plotSettings.horizontalFineGridIsVisible = shouldBeVisible;
  plotSettings.horizontalFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGridVisible(bool shouldBeVisible)
{
  setVerticalCoarseGrid(plotSettings.verticalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGridInterval(double newGridInterval)
{
  setVerticalCoarseGrid(newGridInterval, plotSettings.verticalCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalCoarseGrid(double newGridInterval, 
                                             bool   shouldBeVisible)
{
  if( plotSettings.logScaledX )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.verticalCoarseGridInterval = 2.0;
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

  plotSettings.verticalCoarseGridIsVisible = shouldBeVisible;
  plotSettings.verticalCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGridVisible(bool shouldBeVisible)
{
  setVerticalFineGrid(plotSettings.verticalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGridInterval(double newGridInterval)
{
  setVerticalFineGrid(newGridInterval, plotSettings.verticalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setVerticalFineGrid(double newGridInterval, 
                                           bool   shouldBeVisible)
{
  if( plotSettings.logScaledX )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.verticalFineGridInterval = pow(2.0, 1.0/3.0);
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

  plotSettings.verticalFineGridIsVisible = shouldBeVisible;
  plotSettings.verticalFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGridVisible(bool shouldBeVisible)
{
  setRadialCoarseGrid(plotSettings.radialCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGridInterval(double newGridInterval)
{
  setRadialCoarseGrid(newGridInterval, plotSettings.radialCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialCoarseGrid(double newGridInterval, 
                                           bool   shouldBeVisible)
{
  if( plotSettings.logScaledRadius )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.radialCoarseGridInterval = 2.0;
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

  plotSettings.radialCoarseGridIsVisible = shouldBeVisible;
  plotSettings.radialCoarseGridInterval  = newGridInterval;

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGridVisible(bool shouldBeVisible)
{
  setRadialFineGrid(plotSettings.radialFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGridInterval(double newGridInterval)
{
  setRadialFineGrid(newGridInterval, plotSettings.radialFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setRadialFineGrid(double newGridInterval, 
                                         bool   shouldBeVisible)
{
  if( plotSettings.logScaledRadius )
  {
    jassert(newGridInterval > 1.00001);
    // for logarithmic scaling, we need the grid-intervals to be > 1
    if( newGridInterval <= 1.00001 )
    {
      plotSettings.radialFineGridInterval = pow(2.0, 1.0/3.0);
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

  plotSettings.radialFineGridIsVisible     = shouldBeVisible;
  plotSettings.radialFineGridInterval = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularCoarseGridVisible(bool shouldBeVisible)
{
  setAngularCoarseGrid(plotSettings.angularCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularCoarseGridInterval(double newGridInterval)
{
  setAngularCoarseGrid(newGridInterval, plotSettings.angularCoarseGridIsVisible);
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

  plotSettings.angularCoarseGridIsVisible = shouldBeVisible;
  plotSettings.angularCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularFineGridVisible(bool shouldBeVisible)
{
  setAngularFineGrid(plotSettings.angularFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::setAngularFineGridInterval(double newGridInterval)
{
  setAngularFineGrid(newGridInterval, plotSettings.angularFineGridIsVisible);
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

  plotSettings.angularFineGridIsVisible = shouldBeVisible;
  plotSettings.angularFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

bool CoordinateSystemOld::isHorizontalCoarseGridVisible() 
{ 
  return plotSettings.horizontalCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isHorizontalFineGridVisible() 
{ 
  return plotSettings.horizontalFineGridIsVisible; 
}

bool CoordinateSystemOld::isVerticalCoarseGridVisible() 
{ 
  return plotSettings.verticalCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isVerticalFineGridVisible() 
{ 
  return plotSettings.verticalFineGridIsVisible; 
}

bool CoordinateSystemOld::isRadialCoarseGridVisible() 
{ 
  return plotSettings.radialCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isRadialFineGridVisible() 
{ 
  return plotSettings.radialFineGridIsVisible; 
}

bool CoordinateSystemOld::isAngularCoarseGridVisible() 
{ 
  return plotSettings.angularCoarseGridIsVisible; 
}

bool CoordinateSystemOld::isAngularFineGridVisible() 
{ 
  return plotSettings.angularFineGridIsVisible; 
}

double CoordinateSystemOld::getHorizontalCoarseGridInterval()
{
  return plotSettings.horizontalCoarseGridInterval;
}

double CoordinateSystemOld::getHorizontalFineGridInterval()
{
  return plotSettings.horizontalFineGridInterval;
}

double CoordinateSystemOld::getVerticalCoarseGridInterval()
{
  return plotSettings.verticalCoarseGridInterval;
}

double CoordinateSystemOld::getVerticalFineGridInterval()
{
  return plotSettings.verticalFineGridInterval;
}

double CoordinateSystemOld::getRadialCoarseGridInterval()
{
  return plotSettings.radialCoarseGridInterval;
}

double CoordinateSystemOld::getRadialFineGridInterval()
{
  return plotSettings.radialFineGridInterval;
}

double CoordinateSystemOld::getAngularCoarseGridInterval()
{
  return plotSettings.angularCoarseGridInterval;
}

double CoordinateSystemOld::getAngularFineGridInterval()
{
  return plotSettings.angularFineGridInterval;
}

void CoordinateSystemOld::setAxisPositionX(int newAxisPositionX)
{
  if( newAxisPositionX == rsPlotSettings::INVISIBLE ||
    newAxisPositionX == rsPlotSettings::ZERO      ||
    newAxisPositionX == rsPlotSettings::TOP       ||
    newAxisPositionX == rsPlotSettings::BOTTOM       )
  {
    plotSettings.axisPositionX = newAxisPositionX;
    if(autoReRenderImage == true)
      updateBackgroundImage();
  }
}

void CoordinateSystemOld::setAxisPositionY(int newAxisPositionY)
{
  if( newAxisPositionY == rsPlotSettings::INVISIBLE ||
    newAxisPositionY == rsPlotSettings::ZERO      ||
    newAxisPositionY == rsPlotSettings::LEFT      ||
    newAxisPositionY == rsPlotSettings::RIGHT        )
  {
    plotSettings.axisPositionY = newAxisPositionY;
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

  plotSettings.maximumRange.setRangeX(newMin, newMax);
  plotSettings.currentRange.setRangeX(newMin, newMax);
  plotSettings.logScaledX = shouldBeLogScaled;
  if( newAxisPosition == rsPlotSettings::INVISIBLE ||
      newAxisPosition == rsPlotSettings::ZERO      ||
      newAxisPosition == rsPlotSettings::TOP       ||
      newAxisPosition == rsPlotSettings::BOTTOM       )
  {
    plotSettings.axisPositionX = newAxisPosition;
  }
  plotSettings.verticalCoarseGridInterval = newCoarseGridInterval;
  plotSettings.verticalFineGridInterval   = newFineGridInterval;
  updateMapperInputRange();
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

  plotSettings.maximumRange.setRangeY(newMin, newMax);
  plotSettings.currentRange.setRangeY(newMin, newMax);
  plotSettings.logScaledY = shouldBeLogScaled;
  if( newAxisPosition == rsPlotSettings::INVISIBLE ||
      newAxisPosition == rsPlotSettings::ZERO      ||
      newAxisPosition == rsPlotSettings::LEFT      ||
      newAxisPosition == rsPlotSettings::RIGHT       )
  {
    plotSettings.axisPositionY = newAxisPosition;
  }
  plotSettings.horizontalCoarseGridInterval = newCoarseGridInterval;
  plotSettings.horizontalFineGridInterval   = newFineGridInterval;
  updateMapperInputRange();
}

void CoordinateSystemOld::useLogarithmicScale(bool shouldBeLogScaledX, bool shouldBeLogScaledY,                      
  double newLogBaseX, double newLogBaseY)
{
  plotSettings.logScaledX = shouldBeLogScaledX;
  plotSettings.logScaledY = shouldBeLogScaledY;
  updateMapperInputRange();
}

void CoordinateSystemOld::useLogarithmicScaleX(bool   shouldBeLogScaledX, 
                                            double newLogBaseX)
{
  plotSettings.logScaledX = shouldBeLogScaledX;
  updateMapperInputRange();
}

bool CoordinateSystemOld::isLogScaledX()
{
  return plotSettings.logScaledX;
}

void CoordinateSystemOld::useLogarithmicScaleY(bool   shouldBeLogScaledY, 
                                            double newLogBaseY)
{
  plotSettings.logScaledY = shouldBeLogScaledY;
  updateMapperInputRange();
}

bool CoordinateSystemOld::isLogScaledY()
{
  return plotSettings.logScaledY;
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
  jassert(width  >= 1 && height >= 1);
  if( width < 1 || height < 1) return nullptr;
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  Image* img = new Image(Image::RGB, width, height, true);
  Graphics g(*img);
  drawer.drawPlotBackground(g);
  return img;

  // old:
  /*
  jassert(width  >= 1);
  jassert(height >= 1); 
  if( width < 1 || height < 1)  
    return nullptr;
  Image* thePlotImage = new Image(Image::RGB, width, height, true);
  Graphics g(*thePlotImage);
  drawCoordinateSystem(g, thePlotImage);
  return thePlotImage;
  */

}

XmlElement* CoordinateSystemOld::getPlotAsSVG(int width, int height)
{
  jassert(width  >= 1 && height >= 1);
  if( width < 1 || height < 1) return nullptr;
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);

  XmlElement* svg = new XmlElement(String("svg"));
  svg->setAttribute("width", width);
  svg->setAttribute("height", height);
  drawer.drawPlotBackground(svg);
  return svg;

  // old:
  /*
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
  */
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

  // create graphics object associated with the image and the draw on it:
  Graphics g(*backgroundImage);
  drawCoordinateSystem(g);
  repaintOnMessageThread();
}

//-------------------------------------------------------------------------------------------------
// coordinate transformations:

void CoordinateSystemOld::transformToImageCoordinates(double &x, double &y, const Image* theImage)
{
  if( theImage == nullptr ) { transformToComponentsCoordinates(x, y); return; }
  setupCoordinateMapper(coordinateMapper, theImage);  // not very elegant to set it back and forth
  x = coordinateMapper.mapX(x);
  y = coordinateMapper.mapY(y);
  setupCoordinateMapper(coordinateMapper, this);
}

void CoordinateSystemOld::transformFromImageCoordinates(double &x, double &y, 
  const Image *theImage)
{
  if( theImage == NULL ) { transformFromComponentsCoordinates(x, y); return; }
  setupCoordinateMapper(coordinateMapper, theImage);
  x = coordinateMapper.unmapX(x);
  y = coordinateMapper.unmapY(y);
  setupCoordinateMapper(coordinateMapper, this);
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

//-------------------------------------------------------------------------------------------------
// drawing functions

void CoordinateSystemOld::drawCoordinateSystem(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  // not very elegant - try to refactor (it would be best to get rid of the image and svg 
  // parameters)

  // figure out desired width and height:
  double w, h;
  if(targetSVG != nullptr) {
    double w = targetSVG->getDoubleAttribute("width", 0);
    double h = targetSVG->getDoubleAttribute("height", 0);
    jassert(w != 0 && h != 0); } // svg must have width and height attributes 
  else if(targetImage != nullptr) {
    w = targetImage->getWidth();
    h = targetImage->getHeight(); }
  else {    
    w = getWidth();
    h = getHeight(); }

  // create drawer and draw:
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, w, h);
  if(targetSVG != nullptr)
    drawer.drawPlotBackground(targetSVG);
  else if(targetImage != nullptr) {
    Graphics g2(*targetImage);
    fillRectWithBilinearGradient(g2, 0, 0, targetImage->getWidth(), targetImage->getHeight(),
      plotColourScheme.topLeft, plotColourScheme.topRight, plotColourScheme.bottomLeft,
      plotColourScheme.bottomRight);
    drawer.drawPlotBackground(g2); }
  else {    
    fillRectWithBilinearGradient(g, 0, 0, getWidth(), getHeight(),
      plotColourScheme.topLeft, plotColourScheme.topRight, 
      plotColourScheme.bottomLeft, plotColourScheme.bottomRight);
    drawer.drawPlotBackground(g); }
}

void CoordinateSystemOld::updateMapperOutputRange(Image* image, XmlElement* svg)
{
  if(image != nullptr)    setupCoordinateMapper(coordinateMapper, image);
  else if(svg != nullptr) setupCoordinateMapper(coordinateMapper, svg);
  else                    setupCoordinateMapper(coordinateMapper, this);

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void CoordinateSystemOld::updateMapperInputRange()
{
  coordinateMapper.mapperX.setLogScaled(plotSettings.logScaledX);  // logScaledX/Y are redundant here now
  coordinateMapper.mapperY.setLogScaled(plotSettings.logScaledY);  // ...get rid of them
  coordinateMapper.setInputRange(
    plotSettings.currentRange.getMinX(), plotSettings.currentRange.getMaxX(),
    plotSettings.currentRange.getMinY(), plotSettings.currentRange.getMaxY());

  if(autoReRenderImage == true)
    updateBackgroundImage();
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
  return plotSettings.getStateAsXml();
}

bool CoordinateSystemOld::setStateFromXml(const XmlElement &xml)
{
  plotSettings.setStateFromXml(xml);
  updateBackgroundImage();
  return true; // make function void
}
