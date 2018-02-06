
CurveFamilyPlotOld::CurveFamilyPlotOld(const String& name)
: CoordinateSystemOld(name)
{
  numCurves		    = 0;
  numCurvesToDraw	= 0;
  numValues		    = 0;
  familyValuesX 	= NULL;
  familyValuesY 	= NULL;
  valuesX1		= NULL;
  valuesY1		= NULL;
  fillAreaUnderFunction = false; // ...or perhaps we should fill it? ...maybe with a gradient even?
  isFunctionFamily	    = false;
  plotImage		          = NULL;
  plotImage          		= new Image(Image::RGB, 1, 1, true);
  highlightedCurve      = -1;

  /*
  colourScheme.plotColours.clear();
  colourScheme.plotColours.add(Colour(0.7f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.3f, 1.0f, 0.60f, (uint8) 255));
  colourScheme.plotColours.add(Colour(1.0f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.6f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.2f, 1.0f, 0.60f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.9f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.5f, 1.0f, 0.60f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.1f, 1.0f, 0.60f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.8f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.4f, 1.0f, 0.75f, (uint8) 255));
  colourScheme.plotColours.add(Colour(0.0f, 1.0f, 0.75f, (uint8) 255));
  */
}

CurveFamilyPlotOld::~CurveFamilyPlotOld()
{
  deleteAllChildren();
  if( plotImage != NULL )
    delete plotImage;
}

//-------------------------------------------------------------------------------------------------
// data input:

void CurveFamilyPlotOld::setCurveFamilyValues(int newNumValues, int newNumCurves, 
					   double** newFamilyValuesX, double** newFamilyValuesY)
{
  if( newNumValues != 0 && newNumCurves != 0 
    && newFamilyValuesX != NULL && newFamilyValuesY != NULL )
  {
    numValues	     = newNumValues;
    numCurves	     = newNumCurves;
    familyValuesX    = newFamilyValuesX;
    familyValuesY    = newFamilyValuesY;
    isFunctionFamily = false;
    updatePlotImage();
  }
  else
  {
    numValues	  = 0;
    numCurves	  = 0;
    valuesX1	  = NULL;
    valuesY1	  = NULL;
    familyValuesX = NULL;
    familyValuesY = NULL;
  }
}

void CurveFamilyPlotOld::setCurveValues(int newNumValues, double* newValuesX, double* newValuesY)
{
  if( newNumValues != 0 && newValuesX != NULL && newValuesY != NULL )
  {
    numValues	     = newNumValues;
    numCurves	     = 1;
    valuesX1	     = newValuesX;
    valuesY1	     = newValuesY;
    familyValuesX    = &valuesX1;
    familyValuesY    = &valuesY1;
    isFunctionFamily = false;
    updatePlotImage();
  }
  else
  {
    numValues	  = 0;
    numCurves	  = 0;
    valuesX1	  = NULL;
    valuesY1	  = NULL;
    familyValuesX = NULL;
    familyValuesY = NULL;
  }
}

void CurveFamilyPlotOld::setFunctionFamilyValues(int newNumValues, int newNumCurves, 
					      double* newValuesX, double** newFamilyValuesY)
{
  if( newNumValues != 0 && newNumCurves != 0 
    && newValuesX != NULL && newFamilyValuesY != NULL )
  {

    numValues	     = newNumValues;
    numCurves	     = newNumCurves;
    valuesX1	     = newValuesX;
    familyValuesX    = &valuesX1;
    familyValuesY    = newFamilyValuesY;
    isFunctionFamily = true;
    updatePlotImage();
  }
  else
  {
    numValues	  = 0;
    numCurves	  = 0;
    valuesX1	  = NULL;
    valuesY1	  = NULL;
    familyValuesX = NULL;
    familyValuesY = NULL;
  }
}

void CurveFamilyPlotOld::setHighlightedCurve(int newHighlightedCurve)
{
  if( newHighlightedCurve >= -1 )
    highlightedCurve = newHighlightedCurve;
}

//-------------------------------------------------------------------------------------------------
// drawing, painting:

void CurveFamilyPlotOld::paint(juce::Graphics &g)
{
  if( plotImage != NULL )
    g.drawImage(*plotImage, 0, 0, plotImage->getWidth(), plotImage->getHeight(), 
                            0, 0, plotImage->getWidth(), plotImage->getHeight(), false);
  else
    g.fillAll(Colours::red);

  /*
  // test to see, where juce deems the pixel coordinates to be:
  g.setColour(Colours::red);
  //g.drawLine(10, 0, getWidth()-10, 0, 1.f); 
  //g.drawLine(10, 0.5, getWidth()-10, 0.5, 1.f); 
  g.drawLine(10, 10.5, getWidth()-10, 10.5, 1.f); 
  g.drawLine(10, 100.0f, getWidth()-10, 100.0f, 1.f); 
  g.drawLine(10, 200.5f, getWidth()-10, 200.5f, 1.f); 
  // the lines at N+0.5 (for integer N) are 1 pixel wide, the lines at N are two pixels wide (and 
  // the color is scaled down by 1/2), so it seems we should use the 0.5...width-0.5 instead of
  // 0...width-1
  */
}

void CurveFamilyPlotOld::updatePlotImage(bool redrawCoordinateSystem)
{
  if( getWidth() < 1 || getHeight() < 1 )
    return;

  // redraw the underlying CoordinateSystem if so chosen:
  if( redrawCoordinateSystem == true )
    CoordinateSystemOld::updateBackgroundImage();

  // the coordinate-system itself has been drawn on our (inherited) backgroundImage member variable
  // now we plot the curve(s) on top of that:

  //copyImage(backgroundImage, plotImage);  // doesn't work anymore
  *plotImage = backgroundImage->createCopy(); // test - experimental - ok, seems to work
  Graphics g(*plotImage);
  plotCurveFamily(g); 

  repaintOnMessageThread();
  //repaint();  
}

void CurveFamilyPlotOld::updateBackgroundImage()
{
  CoordinateSystemOld::updateBackgroundImage();
  updatePlotImage(false);
}

XmlElement* CurveFamilyPlotOld::getPlotAsSVG(int width, int height)
{
  jassert( width  >= 1 );
  jassert( height >= 1);
    
  if( width < 1 || height < 1)  
    return NULL;

  Image* thePlot = new Image(Image::RGB, width, height, true);

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*thePlot);

  // create an XmlElement to be used for the SVG drawing:
  XmlElement* theSVG = new XmlElement(String("svg"));
  theSVG->setAttribute(String("width"), width);
  theSVG->setAttribute(String("height"), height);

  // fill the background:
  //g.fillAll(colourScheme.backgroundColour);

  // draw the grids, if desired:
  if( horizontalFineGridIsVisible )
    drawHorizontalGrid(g, horizontalFineGridInterval, logScaledY, plotColourScheme.fineGrid, 1.0f, 
      thePlot, theSVG);
  if( verticalFineGridIsVisible )
    drawVerticalGrid(g, verticalFineGridInterval, logScaledX, plotColourScheme.fineGrid, 1.0f, 
      thePlot, theSVG);
  if( radialFineGridIsVisible )
    drawRadialGrid(g, radialFineGridInterval, logScaledRadius, plotColourScheme.fineGrid, 1.0f, 
      thePlot, theSVG);
  if( angularFineGridIsVisible )
    drawAngularGrid(g, angularFineGridInterval, plotColourScheme.fineGrid, 1.0f, thePlot, theSVG);

  if( horizontalCoarseGridIsVisible )
    drawHorizontalGrid(g, horizontalCoarseGridInterval, logScaledY, plotColourScheme.coarseGrid, 
      1.0f, thePlot, theSVG);
  if( verticalCoarseGridIsVisible )
    drawVerticalGrid(g, verticalCoarseGridInterval, logScaledX, plotColourScheme.coarseGrid, 
      1.0f, thePlot, theSVG);
  if( radialCoarseGridIsVisible )
    drawRadialGrid(g, radialCoarseGridInterval, logScaledRadius, plotColourScheme.coarseGrid, 
      1.0f, thePlot, theSVG);
  if( angularCoarseGridIsVisible )
    drawAngularGrid(g, angularCoarseGridInterval, plotColourScheme.coarseGrid, 1.0f, thePlot, theSVG);


  // draw the function family values:
  plotCurveFamily(g, thePlot, theSVG); 


  // draw the coordinate system:
  if( axisPositionX != INVISIBLE )
    drawAxisX(g, thePlot, theSVG);
  if( axisPositionY != INVISIBLE )
    drawAxisY(g, thePlot, theSVG);

  // draw the labels on the axes:
  if( axisPositionX != INVISIBLE && axisLabelPositionX != NO_ANNOTATION )
    drawAxisLabelX(g, thePlot, theSVG);
  if( axisPositionY != INVISIBLE && axisLabelPositionY != NO_ANNOTATION )
    drawAxisLabelY(g, thePlot, theSVG);

  // draw the values on the axes:
  if( axisPositionX != INVISIBLE && axisValuesPositionX != NO_ANNOTATION )
    drawAxisValuesX(g, thePlot, theSVG);
  if( axisPositionY != INVISIBLE && axisValuesPositionY != NO_ANNOTATION )
    drawAxisValuesY(g, thePlot, theSVG);

  // draw the caption:
  drawCaption(g, thePlot, theSVG);

  return theSVG;
}

//  BAAHHH! code duplication smells bad!!! refactor this!
Image* CurveFamilyPlotOld::getPlotAsImage(int width, int height)
{
  jassert(width  >= 1);
  jassert(height >= 1);
  if( width < 1 || height < 1)  
    return NULL;

  Image* thePlot = new Image(Image::RGB, width, height, true);

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*thePlot);

  // fill the background:
  //g.fillAll(colourScheme.backgroundColour);

  CoordinateSystemOld::drawCoordinateSystem(g, thePlot);

  // draw the grids, if desired:
  if( horizontalFineGridIsVisible )
    drawHorizontalGrid(g, horizontalFineGridInterval, logScaledY, plotColourScheme.fineGrid, 
      1.0f, thePlot);
  if( verticalFineGridIsVisible )
    drawVerticalGrid(g, verticalFineGridInterval, logScaledX, plotColourScheme.fineGrid, 
      1.0f, thePlot);
  if( radialFineGridIsVisible )
    drawRadialGrid(g, radialFineGridInterval, logScaledRadius, plotColourScheme.fineGrid, 
      1.0f, thePlot);
  if( angularFineGridIsVisible )
    drawAngularGrid(g, angularFineGridInterval, plotColourScheme.fineGrid, 1.0f, thePlot);

  if( horizontalCoarseGridIsVisible )
    drawHorizontalGrid(g, horizontalCoarseGridInterval, logScaledY, plotColourScheme.coarseGrid, 
      1.0f, thePlot);
  if( verticalCoarseGridIsVisible )
    drawVerticalGrid(g, verticalCoarseGridInterval, logScaledX, plotColourScheme.coarseGrid, 
      1.0f, thePlot);
  if( radialCoarseGridIsVisible )
    drawRadialGrid(g, radialCoarseGridInterval, logScaledRadius, plotColourScheme.coarseGrid, 
      1.0f, thePlot);
  if( angularCoarseGridIsVisible )
    drawAngularGrid(g, angularCoarseGridInterval, plotColourScheme.coarseGrid, 1.0f, thePlot);

  // draw the function family values:
  plotCurveFamily(g, thePlot, NULL); 

  // draw the coordinate system:
  if( axisPositionX != INVISIBLE )
    drawAxisX(g, thePlot);
  if( axisPositionY != INVISIBLE )
    drawAxisY(g, thePlot);

  // draw the labels on the axes:
  if( axisPositionX != INVISIBLE && axisLabelPositionX != NO_ANNOTATION )
    drawAxisLabelX(g, thePlot);
  if( axisPositionY != INVISIBLE && axisLabelPositionY != NO_ANNOTATION )
    drawAxisLabelY(g, thePlot);

  // draw the values on the axes:
  if( axisPositionX != INVISIBLE && axisValuesPositionX != NO_ANNOTATION )
    drawAxisValuesX(g, thePlot);
  if( axisPositionY != INVISIBLE && axisValuesPositionY != NO_ANNOTATION )
    drawAxisValuesY(g, thePlot);

  // draw the caption:
  drawCaption(g, thePlot, NULL);

  // draw an outlining rectangle:
  g.setColour(plotColourScheme.outline);
  g.drawRect(0, 0, thePlot->getWidth(), thePlot->getHeight(), 1);

  return thePlot;
}

void CurveFamilyPlotOld::plotCurveFamily(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  Colour graphColour = getCurveColour(0);
  g.setColour(graphColour); 

  if( highlightedCurve == -1 )
  {
    for(int k=0; k<numCurves; k++)
    {
      graphColour = getCurveColour(k);
      g.setColour(graphColour); 
      plotCurve(g, targetImage, targetSVG, k);
    }
  }
  else
  {
    for(int k=0; k<numCurves; k++)
    {
      graphColour = getCurveColour(k);
      if( k != highlightedCurve )
        graphColour = graphColour.withMultipliedAlpha(0.5f);  
      g.setColour(graphColour); 
      plotCurve(g, targetImage, targetSVG, k);
    }
  }

}

void CurveFamilyPlotOld::plotCurve(Graphics &g, Image* targetImage, XmlElement *targetSVG, int index)
{
  int	   i;	// indices for the loops through the curves and through the values
  double x, y;	// the current x- and y-value (in the current curve)

  // make sure that the arrays are valid:
  if( familyValuesX==NULL || familyValuesY==NULL )
    return;

  // create a mask to multiply the first index for the function-family plots:
  int mask = 1;
  if( isFunctionFamily )
    mask = 0;

  String curvePathDataString = String::empty;

  // start values for a line segment:
  x = familyValuesX[index*mask][0]; 
  y = familyValuesY[index][0];

  if( targetImage == NULL )
    transformToComponentsCoordinates(x, y);
  else
    transformToImageCoordinates(x, y, targetImage);

  // create a Path object for the actual plot:
  Path	 funcPath;

  // start the path at the first (transformed) data-point:
  funcPath.startNewSubPath((float) x, (float) y);

  // move the 'pen' for the svg-drawing to the start-position:
  if( targetSVG != NULL )
    curvePathDataString += String(x) + String(" ") + String(y) + String(", ");

  for(i=1; i<numValues; i++)
  {
    // read out the tables:
    x = familyValuesX[index*mask][i];
    y = familyValuesY[index][i];

    // transform:
    if( targetImage == NULL )
      transformToComponentsCoordinates(x, y);
    else
      transformToImageCoordinates(x, y, targetImage);

    // add a line-segment to the transformed (x,y):
    funcPath.lineTo((float) x, (float) y);

    // add the line-segment the svg-drawing also:
    if( targetSVG != NULL )
      curvePathDataString += String(x) + String(" ") + String(y) + String(", ");
  }

  // set the colour for the current curve:
  //Colour graphColour = getCurveColour(index);
  //g.setColour(graphColour); 

  // draw the path:
  if(fillAreaUnderFunction)
  {
    // close the path:
    x = familyValuesX[index*mask][numValues-1];
    y = 0;
    if( targetImage == NULL )
      transformToComponentsCoordinates(x, y);
    else
      transformToImageCoordinates(x, y, targetImage);
    funcPath.lineTo((float) x, (float) y);

    x = familyValuesX[index*mask][0];
    y = 0;
    if( targetImage == NULL )
      transformToComponentsCoordinates(x, y);
    else
      transformToImageCoordinates(x, y, targetImage);
    funcPath.lineTo((float) x, (float) y);

    x = familyValuesX[index*mask][0];
    y = familyValuesY[index][0];
    if( targetImage == NULL )
      transformToComponentsCoordinates(x, y);
    else
      transformToImageCoordinates(x, y, targetImage);
    funcPath.lineTo((float) x, (float) y);

    funcPath.closeSubPath();

    //g.setColour(graphColour);
    g.fillPath(funcPath);
  }
  else
  {
    g.strokePath(funcPath, PathStrokeType(2.f));

    if( targetSVG != NULL )
    {
      // create a XmlElement for the path, setup its attributes and add it ot the svg-drawing:
      Colour colour = Colours::black; // preliminary
      XmlElement* curvePath = new XmlElement(String("polyline"));
      curvePath->setAttribute(String("points"), curvePathDataString);
      curvePath->setAttribute(String("style"), String("stroke-width: ") + String(1.0) + 
        String("; stroke: #") + colour.toString().substring(2) + 
        String("; fill: none;") );
      targetSVG->addChildElement(curvePath);
    }
  }
}

void CurveFamilyPlotOld::plotFamilyValuesAsDots(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  /*
  // make sure that the arrays are valid:
  if( valuesX==NULL || valuesY==NULL )
  return;

  g.setColour(graphColor);

  float x0 = 0.0;
  float y0 = 0.0;
  transformToComponentsCoordinates(x0, y0);
  // transformed zero values for x and y, will be used to draw the lines to
  // the x- and y-axes if this option is chosen

  float x,y;	   // current x and y value
  long	i;	   // counter through the values
  for(i=0; i<numValues; i++)
  {
  //read out the tables:
  x = valuesX[i];
  y = valuesY[i];

  //transform:
  transformToComponentsCoordinates(x, y);

  //add a dot at postion x,y:
  g.fillEllipse(x-1, y-1, 3, 3);

  // draw lines to x- and y-axis if the option is selected:
  if( drawLineToXAxis )
  g.drawLine(x, y, x, y0);
  if( drawLineToYAxis )
  g.drawLine(x0, y, x, y);
  }
  */
}

//-------------------------------------------------------------------------------------------------
// others:

void CurveFamilyPlotOld::resized()
{
  CoordinateSystemOld::resized();
  CoordinateSystemOld::updateBackgroundImage();

  delete plotImage;
  plotImage = new Image(backgroundImage->getFormat(), backgroundImage->getWidth(), 
    backgroundImage->getHeight(), false); 

  if(autoReRenderImage == true)
    updatePlotImage(false);
}

void CurveFamilyPlotOld::invalidatePointers()
{
  numCurves	= 0;
  familyValuesX = NULL;
  familyValuesY = NULL;
  valuesX1	= NULL;
  valuesY1	= NULL;
}
