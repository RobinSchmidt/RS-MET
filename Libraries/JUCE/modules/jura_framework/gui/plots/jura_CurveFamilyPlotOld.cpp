
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
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  XmlElement* svg = CoordinateSystemOld::getPlotAsSVG(width, height);
  for(int k = 0; k < numCurves; k++)
  {
    //drawer.drawWithLines(svg, numValues, familyValuesX[k], familyValuesY[k]);
    // what about color?
  }
  //drawer.drawForeground(svg);
  return svg;


  /*
  jassert( width  >= 1 );
  jassert( height >= 1);
    
  if( width < 1 || height < 1)  
    return NULL;

  Image* thePlot = new Image(Image::RGB, width, height, true); // obsolete?

  // create a graphics object which is associated with the image to perform
  // the drawing-operations
  Graphics g(*thePlot);

  // create an XmlElement to be used for the SVG drawing:
  XmlElement* theSVG = new XmlElement(String("svg"));
  theSVG->setAttribute(String("width"), width);
  theSVG->setAttribute(String("height"), height);

  updateMapperOutputRange(nullptr, theSVG);

  //CoordinateSystemOld::drawCoordinateSystem...
  //// draw the function family values:
  //plotCurveFamily(g, thePlot, theSVG); 

  return theSVG;
  */
}

Image* CurveFamilyPlotOld::getPlotAsImage(int width, int height)
{
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  Image* img = CoordinateSystemOld::getPlotAsImage(width, height);
  Graphics g(*img);
  for(int k = 0; k < numCurves; k++)
  {
    drawer.drawWithLines(g, numValues, familyValuesX[k], familyValuesY[k]);
    // preliminary, what about color? currently, the curve is black (on black background)
  }
  //drawer.drawForeground(g);
  return img;
}

void CurveFamilyPlotOld::plotCurveFamily(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( familyValuesX == nullptr || familyValuesY == nullptr )
    return;

  // new (does not yet support image or svg drawing)
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, getWidth(), getHeight());
  for(int k = 0; k < numCurves; k++) // or numCurvesToDraw?
  {
    Colour graphColour = getCurveColour(k);
    if( k != highlightedCurve || highlightedCurve == -1 )
      graphColour = graphColour.withMultipliedAlpha(0.5f);  
    g.setColour(graphColour); 
    //drawer.fillFunction( g, numValues, familyValuesX[k], familyValuesY[k]); // test
    drawer.drawWithLines(g, numValues, familyValuesX[k], familyValuesY[k]);
    //drawer.drawAsDots(   g, numValues, familyValuesX[k], familyValuesY[k], 5.f, true, false);
  }

  // old:
  /*
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
  */

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

/*
void CurveFamilyPlotOld::plotFamilyValuesAsDots(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{

}
*/

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
