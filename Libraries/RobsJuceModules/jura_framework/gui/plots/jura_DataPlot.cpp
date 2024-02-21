rsDataPlot::rsDataPlot(const String& name)
: rsPlot(name)
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

rsDataPlot::~rsDataPlot()
{
  deleteAllChildren();
  if( plotImage != NULL )
    delete plotImage;
}

//-------------------------------------------------------------------------------------------------
// data input:

void rsDataPlot::setCurveFamilyValues(int newNumValues, int newNumCurves, 
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

void rsDataPlot::setCurveValues(int newNumValues, double* newValuesX, double* newValuesY)
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

void rsDataPlot::setFunctionFamilyValues(int newNumValues, int newNumCurves, 
  double* newValuesX, double** newFamilyValuesY)
{
  if( newNumValues != 0 && newNumCurves != 0 
    && newValuesX != NULL && newFamilyValuesY != NULL )
  {

    numValues        = newNumValues;
    numCurves        = newNumCurves;
    valuesX1         = newValuesX;
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

void rsDataPlot::setHighlightedCurve(int newHighlightedCurve)
{
  if( newHighlightedCurve >= -1 )
    highlightedCurve = newHighlightedCurve;
}

//-------------------------------------------------------------------------------------------------
// drawing, painting:

void rsDataPlot::paint(juce::Graphics &g)
{
  if( plotImage != NULL )
    g.drawImage(*plotImage, 0, 0, plotImage->getWidth(), plotImage->getHeight(), 
                            0, 0, plotImage->getWidth(), plotImage->getHeight(), false);
  else
    g.fillAll(Colours::red);
}

void rsDataPlot::updatePlotImage(bool redrawCoordinateSystem)
{
  if( getWidth() < 1 || getHeight() < 1 )
    return;

  // redraw the underlying CoordinateSystem if so chosen:
  if( redrawCoordinateSystem == true )
    rsPlot::updateBackgroundImage();

  // the coordinate-system itself has been drawn on our (inherited) backgroundImage member variable
  // now we plot the curve(s) on top of that:

  //copyImage(backgroundImage, plotImage);  // doesn't work anymore
  *plotImage = backgroundImage->createCopy(); // test - experimental - ok, seems to work
  Graphics g(*plotImage);                     // but we should really copy only the data, not create
  plotCurveFamily(g);                         // a new object

  repaintOnMessageThread();
}

void rsDataPlot::updateBackgroundImage()
{
  rsPlot::updateBackgroundImage();
  updatePlotImage(false);
}

XmlElement* rsDataPlot::getPlotAsSVG(int width, int height)
{
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  XmlElement* svg = rsPlot::getPlotAsSVG(width, height);
  for(int k = 0; k < numCurves; k++)
  {
    drawer.drawWithLines(svg, numValues, familyValuesX[k], familyValuesY[k]);
    // what about color?
  }
  //drawer.drawForeground(svg);
  return svg;
}

Image* rsDataPlot::getPlotAsImage(int width, int height)
{
  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, width, height);
  Image* img = rsPlot::getPlotAsImage(width, height);
  Graphics g(*img);
  for(int k = 0; k < numCurves; k++)
  {
    drawer.drawWithLines(g, numValues, familyValuesX[k], familyValuesY[k]);
    // preliminary, what about color? currently, the curve is black (on black background)
  }
  //drawer.drawForeground(g);
  return img;
}

void rsDataPlot::plotCurveFamily(Graphics &g, Image* targetImage, XmlElement *targetSVG)
{
  if( familyValuesX == nullptr || familyValuesY == nullptr )
    return;

  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, getWidth(), getHeight());
  for(int k = 0; k < numCurves; k++) // or numCurvesToDraw?
  {
    Colour graphColour = getCurveColour(k);
    if( k != highlightedCurve && highlightedCurve != -1 )
      graphColour = graphColour.withMultipliedAlpha(0.5f);  // maybe highlight by filling area?
    g.setColour(graphColour); 
    //drawer.fillFunction( g, numValues, familyValuesX[k], familyValuesY[k]); // test
    drawer.drawWithLines(g, numValues, familyValuesX[k], familyValuesY[k]);
    //drawer.drawAsDots(   g, numValues, familyValuesX[k], familyValuesY[k], 5.f, true, false);
  }

}

//-------------------------------------------------------------------------------------------------
// others:

void rsDataPlot::resized()
{
  rsPlot::resized();
  rsPlot::updateBackgroundImage();

  delete plotImage;
  plotImage = new Image(backgroundImage->getFormat(), backgroundImage->getWidth(), 
    backgroundImage->getHeight(), false); 

  if(autoReRenderImage == true)
    updatePlotImage(false);
}

void rsDataPlot::invalidatePointers()
{
  numCurves	= 0;
  familyValuesX = NULL;
  familyValuesY = NULL;
  valuesX1	= NULL;
  valuesY1	= NULL;
}
