// construction/destruction:

rsPlot::rsPlot(const String &newDescription) 
  : DescribedComponent(newDescription) //RWidget(newDescription)
{
  autoReRenderImage = false;

  // initialize the function-pointers for value->string conversion
  stringConversionForInfoLineX = &valueToString0;
  stringConversionForInfoLineY = &valueToString0;

  // initialize the component-size and the image-size to 1x1 pixels, without
  // such initializations, a JUCE-breakpoint will be triggered or other screws
  // happen:
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
}

rsPlot::~rsPlot()
{
  deleteAllChildren(); // ?  
  delete backgroundImage;
}

//-------------------------------------------------------------------------------------------------
// component-overrides:

void rsPlot::mouseDown(const MouseEvent &e)
{
  if( e.mods.isRightButtonDown() && showPopUpOnRightClick == true )
    openRightClickPopupMenu();
}

void rsPlot::mouseEnter(const MouseEvent &e)
{
  if( showPositionAsDescription == true )
    setDescription( getInfoLineForPixelPosition(e.x, e.y) );
  DescribedComponent::mouseEnter(e);
}

void rsPlot::mouseMove(const MouseEvent &e)
{
  if( showPositionAsDescription == true )
    setDescription( getInfoLineForPixelPosition(e.x, e.y) );
  DescribedComponent::mouseMove(e);
}

void rsPlot::resized()
{
  updateMapperOutputRange();
}

void rsPlot::paint(juce::Graphics &g)
{
  if( backgroundImage != NULL )
    g.drawImage(*backgroundImage, 0, 0, getWidth(), getHeight(), 
      0, 0, backgroundImage->getWidth(), backgroundImage->getHeight(), false);
  else
    g.fillAll(Colours::red);
}

void rsPlot::rsPlotAppearanceChanged(rsPlotSettings* plotSettings)
{
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::rsPlotVisibleRangeChanged(rsPlotSettings* plotSettings)
{
  updateMapperInputRange();
}

//-------------------------------------------------------------------------------------------------
// range management:

void rsPlot::setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  plotSettings.setMaximumRange(newMinX, newMaxX, newMinY, newMaxY);
}

void rsPlot::setMaximumRange(rsPlotRange newMaximumRange)
{
  plotSettings.setMaximumRange(newMaximumRange);
}

void rsPlot::setMaximumRangeX(double newMinX, double newMaxX)
{
  plotSettings.setMaximumRangeX(newMinX, newMaxX);
}

void rsPlot::setMaximumRangeY(double newMinY, double newMaxY)
{
  plotSettings.setMaximumRangeY(newMinY, newMaxY);
}

void rsPlot::setMaximumRangeMinX(double newMinX)
{
  plotSettings.setMaximumRangeMinX(newMinX);
}

void rsPlot::setMaximumRangeMaxX(double newMaxX)
{
  plotSettings.setMaximumRangeMaxX(newMaxX);
}

void rsPlot::setMaximumRangeMinY(double newMinY)
{
  plotSettings.setMaximumRangeMinY(newMinY);
}

void rsPlot::setMaximumRangeMaxY(double newMaxY)
{
  plotSettings.setMaximumRangeMaxY(newMaxY);
}

void rsPlot::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  plotSettings.setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
}

void rsPlot::setCurrentRange(rsPlotRange newRange)
{
  plotSettings.setCurrentRange(newRange);
}

void rsPlot::setCurrentRangeX(double newMinX, double newMaxX)
{
  plotSettings.setCurrentRangeX(newMinX, newMaxX);
}

void rsPlot::setCurrentRangeY(double newMinY, double newMaxY)
{
  plotSettings.setCurrentRangeY(newMinY, newMaxY);
}

void rsPlot::setCurrentRangeMinX(double newMinX)
{
  plotSettings.setCurrentRangeMinX(newMinX);
}

void rsPlot::setCurrentRangeMaxX(double newMaxX)
{
  plotSettings.setCurrentRangeMaxX(newMaxX);
}

void rsPlot::setCurrentRangeMinY(double newMinY)
{
  plotSettings.setCurrentRangeMinY(newMinY);
}

void rsPlot::setCurrentRangeMaxY(double newMaxY)
{
  plotSettings.setCurrentRangeMaxY(newMaxY);
}

String rsPlot::getInfoLineForPixelPosition(int x, int y)
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

void rsPlot::setColourScheme(const PlotColourScheme& newColourScheme)
{
  plotColourScheme = newColourScheme;
  updateBackgroundImage();
}

void rsPlot::setColourSchemeFromXml(const XmlElement &xml)
{
  //colourScheme.setColourSchemeFromXml(xml);
}

void rsPlot::setAutoReRendering(bool shouldAutomaticallyReRender)
{
  autoReRenderImage = shouldAutomaticallyReRender;
}

void rsPlot::setCaption(const String &newCaption, int newPosition)
{
  plotSettings.captionPosition = newPosition;
  plotSettings.captionString   = newCaption;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAxisLabels(const String &newLabelX, const String &newLabelY,              
  int newLabelPositionX, int newLabelPositionY)
{
  plotSettings.axisLabelX         = newLabelX;
  plotSettings.axisLabelY         = newLabelY;
  plotSettings.axisLabelPositionX = newLabelPositionX;
  plotSettings.axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAxisLabelX(const String& newLabelX, int newLabelPositionX)
{
  plotSettings.axisLabelX         = newLabelX;
  plotSettings.axisLabelPositionX = newLabelPositionX;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAxisLabelY(const String& newLabelY, int newLabelPositionY)
{
  plotSettings.axisLabelY             = newLabelY;
  plotSettings.axisLabelPositionY = newLabelPositionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAxisValuesPositionX(int newValuesPositionX)
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

void rsPlot::setAxisValuesPositionY(int newValuesPositionY)
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

void rsPlot::setStringConversionForAxisX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  plotSettings.stringConversionForAxisX = newConversionFunctionX;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setStringConversionForInfoLineX(
  String (*newConversionFunctionX) (double valueToBeConverted) )
{
  stringConversionForInfoLineX = newConversionFunctionX;
}

void rsPlot::setStringConversionForAxisY(
  String (*newConversionFunctionY) (double valueToBeConverted) )
{
  plotSettings.stringConversionForAxisY = newConversionFunctionY;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setStringConversionForInfoLineY(
  String (*newConversionFunctionY) (double valueToBeConverted) )
{
  stringConversionForInfoLineY = newConversionFunctionY;
}

void rsPlot::setHorizontalCoarseGridVisible(bool shouldBeVisible)
{
  setHorizontalCoarseGrid(plotSettings.horizontalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setHorizontalCoarseGridInterval(double newGridInterval)
{
  setHorizontalCoarseGrid(newGridInterval, plotSettings.horizontalCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setHorizontalCoarseGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledY);
  plotSettings.horizontalCoarseGridIsVisible = shouldBeVisible;
  plotSettings.horizontalCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setHorizontalFineGridVisible(bool shouldBeVisible)
{
  setHorizontalFineGrid(plotSettings.horizontalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setHorizontalFineGridInterval(double newGridInterval)
{
  setHorizontalFineGrid(newGridInterval, plotSettings.horizontalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setHorizontalFineGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledY);
  plotSettings.horizontalFineGridIsVisible = shouldBeVisible;
  plotSettings.horizontalFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalCoarseGridVisible(bool shouldBeVisible)
{
  setVerticalCoarseGrid(plotSettings.verticalCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalCoarseGridInterval(double newGridInterval)
{
  setVerticalCoarseGrid(newGridInterval, plotSettings.verticalCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalCoarseGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledX);
  plotSettings.verticalCoarseGridIsVisible = shouldBeVisible;
  plotSettings.verticalCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalFineGridVisible(bool shouldBeVisible)
{
  setVerticalFineGrid(plotSettings.verticalFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalFineGridInterval(double newGridInterval)
{
  setVerticalFineGrid(newGridInterval, plotSettings.verticalFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setVerticalFineGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledX);
  plotSettings.verticalFineGridIsVisible = shouldBeVisible;
  plotSettings.verticalFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialCoarseGridVisible(bool shouldBeVisible)
{
  setRadialCoarseGrid(plotSettings.radialCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialCoarseGridInterval(double newGridInterval)
{
  setRadialCoarseGrid(newGridInterval, plotSettings.radialCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialCoarseGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledRadius);
  plotSettings.radialCoarseGridIsVisible = shouldBeVisible;
  plotSettings.radialCoarseGridInterval  = newGridInterval;

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialFineGridVisible(bool shouldBeVisible)
{
  setRadialFineGrid(plotSettings.radialFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialFineGridInterval(double newGridInterval)
{
  setRadialFineGrid(newGridInterval, plotSettings.radialFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setRadialFineGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, plotSettings.logScaledRadius);
  plotSettings.radialFineGridIsVisible = shouldBeVisible;
  plotSettings.radialFineGridInterval = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularCoarseGridVisible(bool shouldBeVisible)
{
  setAngularCoarseGrid(plotSettings.angularCoarseGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularCoarseGridInterval(double newGridInterval)
{
  setAngularCoarseGrid(newGridInterval, plotSettings.angularCoarseGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularCoarseGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, false);
  plotSettings.angularCoarseGridIsVisible = shouldBeVisible;
  plotSettings.angularCoarseGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularFineGridVisible(bool shouldBeVisible)
{
  setAngularFineGrid(plotSettings.angularFineGridInterval, shouldBeVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularFineGridInterval(double newGridInterval)
{
  setAngularFineGrid(newGridInterval, plotSettings.angularFineGridIsVisible);
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAngularFineGrid(double newGridInterval, bool shouldBeVisible)
{
  sanityCheckGridSpacing(newGridInterval, false);
  plotSettings.angularFineGridIsVisible = shouldBeVisible;
  plotSettings.angularFineGridInterval  = newGridInterval;
  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::setAxisPositionX(int newAxisPositionX)
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

void rsPlot::setAxisPositionY(int newAxisPositionY)
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

void rsPlot::setupAxisX(double newMin, double newMax, bool shouldBeLogScaled, 
  int newAxisPosition, double newCoarseGridInterval, double newFineGridInterval)
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

  plotSettings.setMaximumRangeX(newMin, newMax);
  plotSettings.setCurrentRangeX(newMin, newMax);
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

void rsPlot::setupAxisY(double newMin, double newMax, bool shouldBeLogScaled, 
  int newAxisPosition, double newCoarseGridInterval, double newFineGridInterval)
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

  plotSettings.setMaximumRangeY(newMin, newMax);
  plotSettings.setCurrentRangeY(newMin, newMax);
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

void rsPlot::useLogarithmicScale(bool shouldBeLogScaledX, bool shouldBeLogScaledY)
{
  plotSettings.logScaledX = shouldBeLogScaledX;
  plotSettings.logScaledY = shouldBeLogScaledY;
  updateMapperInputRange();
}

void rsPlot::useLogarithmicScaleX(bool shouldBeLogScaledX)
{
  plotSettings.logScaledX = shouldBeLogScaledX;
  updateMapperInputRange();
}

void rsPlot::useLogarithmicScaleY(bool shouldBeLogScaledY)
{
  plotSettings.logScaledY = shouldBeLogScaledY;
  updateMapperInputRange();
}

/*
void rsPlot::setValueFieldPopup(bool shouldPopUp)
{
  valuePopup = shouldPopUp;
  inspectionField->setVisible(false);
}
*/

void rsPlot::sanityCheckGridSpacing(double& spacing, bool logScaled)
{
  if( logScaled ) {
    jassert(spacing > 1.00001); // for logarithmic scaling, we need the grid-intervals to be > 1
    if( spacing <= 1.00001 ) 
      spacing = 2.0; }
  else {
    jassert(spacing > 0.000001); // grid-intervals must be > 0
    if( spacing <= 0.000001 )
      spacing = 1.0; }
}

//-------------------------------------------------------------------------------------------------
// functions for drawing and/or exporting the shown content:
void rsPlot::openRightClickPopupMenu()
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

Image* rsPlot::getPlotAsImage(int width, int height)
{
  jassert(width  >= 1 && height >= 1);
  if( width < 1 || height < 1) return nullptr;
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  Image* img = new Image(Image::RGB, width, height, true);
  Graphics g(*img);
  drawer.drawPlotBackground(g);
  return img;
}

XmlElement* rsPlot::getPlotAsSVG(int width, int height)
{
  jassert(width  >= 1 && height >= 1);
  if( width < 1 || height < 1) return nullptr;
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  XmlElement* svg = new XmlElement(String("svg"));
  svg->setAttribute("width", width);
  svg->setAttribute("height", height);
  drawer.drawPlotBackground(svg);
  return svg;
}

void rsPlot::openExportDialog(int defaultWidth, int defaultHeight, 
  const String &defaultFormat, const File& defaultTargetFile)
{
  ImageSavingDialog dialog(this, defaultWidth, defaultHeight, defaultFormat, defaultTargetFile);
  DialogWindow exportWindow(String("Export to Image or SVG Drawing"), Colours::white, true, true);
  exportWindow.showModalDialog(String("Export to Image or SVG Drawing"), &dialog, this, 
    Colours::white, true, false, false);
}

void rsPlot::updateBackgroundImage()
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
// drawing functions

//void rsPlot::drawCoordinateSystem(Graphics &g, Image *targetImage, XmlElement *targetSVG)
void rsPlot::drawCoordinateSystem(Graphics &g)
{
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, getWidth(), getHeight());
  drawer.drawPlotBackground(g);
}

void rsPlot::updateMapperOutputRange(Image* image, XmlElement* svg)
{
  if(image != nullptr)    setupCoordinateMapper(coordinateMapper, image);
  else if(svg != nullptr) setupCoordinateMapper(coordinateMapper, svg);
  else                    setupCoordinateMapper(coordinateMapper, this);

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

void rsPlot::updateMapperInputRange()
{
  coordinateMapper.mapperX.setLogScaled(plotSettings.logScaledX);
  coordinateMapper.mapperY.setLogScaled(plotSettings.logScaledY);
  coordinateMapper.setInputRange(
    getCurrentRangeMinX(), getCurrentRangeMaxX(),
    getCurrentRangeMinY(), getCurrentRangeMaxY());

  if(autoReRenderImage == true)
    updateBackgroundImage();
}

//-------------------------------------------------------------------------------------------------
// state-management (storing and recall), still incomplete:

XmlElement* rsPlot::getStateAsXml(const String& stateName) const
{
  return plotSettings.getStateAsXml();
}

bool rsPlot::setStateFromXml(const XmlElement &xml)
{
  plotSettings.setStateFromXml(xml);
  updateBackgroundImage();
  return true; // make function void
}
