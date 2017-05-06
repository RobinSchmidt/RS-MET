

//=================================================================================================
// class OscilloscopeAudioModule:

OscilloscopeAudioModule::OscilloscopeAudioModule(CriticalSection *newPlugInLock, 
  rosic::SyncedWaveformDisplayBuffer *displayBufferToUse)
 : AudioModule(newPlugInLock)
{
  jassert(displayBufferToUse != NULL); // you must pass a valid rosic-object to the constructor
  waveformDisplayBuffer = displayBufferToUse;
  moduleName = juce::String("Oscilloscope");
  setActiveDirectory(getApplicationDirectory() + juce::String("/OscilloscopePresets") );
  initializeAutomatableParameters();
  displayIsFrozen = false;
}

void OscilloscopeAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( waveformDisplayBuffer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: waveformDisplayBuffer->setSyncMode(         (int) value  );  break;
  //case 1: waveformDisplayBuffer->setMidSideMode(      value != 0.0 );  break;
  case 2: waveformDisplayBuffer->setTimeWindowLength( value        );  break;
  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load which 
                        // renders it meaningless - fix this!
}

void OscilloscopeAudioModule::initializeAutomatableParameters()
{
  Parameter *p = new Parameter(lock, "SyncMode", 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String("Free Running"));
  p->addStringValue(juce::String("Lowpass Zeros"));
  //p->setValue(1.0);
  addObservedParameter(p);
  addObservedParameter( new Parameter(lock, "MidSideMode",       0.0,   1.0, 1.0,    0.0,  Parameter::BOOLEAN) );
  addObservedParameter( new Parameter(lock, "TimeWindowLength",  0.001, 1.5, 0.001,  0.1,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(lock, "MinAmplitude",     -2.0,  +2.0, 0.1,   -1.5,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(lock, "MaxAmplitude",     -2.0,  +2.0, 0.1,    1.5,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(lock, "FrameRate",        10.0,  50.0, 1.0,   15.0,  Parameter::LINEAR)  );
  addObservedParameter( new Parameter(lock, "Freeze",            0.0,   1.0, 1.0,    1.0,  Parameter::BOOLEAN) );

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class SpectrumAnalyzerAudioModule:

SpectrumAnalyzerAudioModule::SpectrumAnalyzerAudioModule(CriticalSection *newPlugInLock, rosic::SpectrumAnalyzer *spectrumAnalyzerToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(spectrumAnalyzerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedSpectrumAnalyzer = spectrumAnalyzerToWrap;
  moduleName = juce::String("SpectrumAnalyzer");
  setActiveDirectory(getApplicationDirectory() + juce::String("/SpectrumAnalyzerPresets") );
  initializeAutomatableParameters();
  displayIsFrozen  =  false;
}

void SpectrumAnalyzerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedSpectrumAnalyzer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: wrappedSpectrumAnalyzer->setBlockSize(   (int) pow(2.0, value+8)     );  break;
  case 1: wrappedSpectrumAnalyzer->setMidSideMode(                value != 0.0 );  break;

  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load which 
                        // renders it meaningless - fix this!
}

void SpectrumAnalyzerAudioModule::initializeAutomatableParameters()
{
  addObservedParameter( new ParameterPowersOfTwo(lock, "FFTSize", 256.0, 32768.0, 0.0, 1024.0) );

  addObservedParameter( new Parameter(lock, "MidSideMode",        0.0,       1.0, 1.0,     0.0,   Parameter::BOOLEAN)      );
  addObservedParameter( new Parameter(lock, "LinearFrequency",    0.0,       1.0, 1.0,     0.0,   Parameter::BOOLEAN)      );
  addObservedParameter( new Parameter(lock, "MinFrequency",      15.625, 32000.0, 0.0,    15.625, Parameter::EXPONENTIAL)  );
  addObservedParameter( new Parameter(lock, "MaxFrequency",      15.625, 32000.0, 0.0, 32000.0,   Parameter::EXPONENTIAL)  );
  addObservedParameter( new Parameter(lock, "MinLevel",        -100.0,      10.0, 0.0,  -100.0,   Parameter::LINEAR)       );
  addObservedParameter( new Parameter(lock, "MaxLevel",        -100.0,      10.0, 0.0,    10.0,   Parameter::LINEAR)       );
  addObservedParameter( new Parameter(lock, "FrameRate",        10.0,       50.0, 1.0,   15.0,    Parameter::LINEAR)       );
  addObservedParameter( new Parameter(lock, "Freeze",            0.0,        1.0, 1.0,    1.0,    Parameter::BOOLEAN)      );

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class MultiAnalyzerAudioModule:

MultiAnalyzerAudioModule::MultiAnalyzerAudioModule(CriticalSection *newPlugInLock, 
                                                   OscilloscopeAudioModule *newOscilloscopeModule, 
                                                   SpectrumAnalyzerAudioModule *newSpectrumAnalyzerModule)
                                                    : AudioModule(newPlugInLock)
{
  oscilloscopeModule     = newOscilloscopeModule;
  spectrumAnalyzerModule = newSpectrumAnalyzerModule;

  jassert( oscilloscopeModule     != NULL );
  jassert( spectrumAnalyzerModule != NULL );

  addChildAudioModule(oscilloscopeModule);
  addChildAudioModule(spectrumAnalyzerModule);

  moduleName = juce::String("MultiAnalyzer");
  setActiveDirectory(getApplicationDirectory() + juce::String("/MultiAnalyzerPresets") );

  initializeAutomatableParameters();
}

void MultiAnalyzerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: mode = (int) value;  break;
  } 

  markStateAsDirty(); // this feature is de-activated because the state will be marked as dirty immediately after preset-load due to some
  // async scrollBarMove callback - this renders the feature meaningless - fix this!
}

void MultiAnalyzerAudioModule::initializeAutomatableParameters()
{
  Parameter *p = new Parameter(lock, "Mode", 0.0, 1.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String("Oscilloscope"));
  p->addStringValue(juce::String("SpectrumAnalyzer"));
  p->setValue(1.0, false, false);
  addObservedParameter(p);

  //addObservedParameter( new Parameter("FrameRate",  10.0,  50.0, 1.0,  25.0, Parameter::LINEAR)  );

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

SpectrumAnalyzerDisplay::SpectrumAnalyzerDisplay(const juce::String& name) : SpectrumDisplayOld(name)
{
  setMaximumRange(15.625, 32000.0, -100.0, 10.0);
  setCurrentRange(15.625, 32000.0, -100.0, 10.0);
  showPositionAsDescription = true;
  // we must set this here to make the MultiAnalyzer intialize its range correctly because the inherited (lower) values would trigger a
  // callback to SpectrumAnalyzerModuleEditor::coordinateSystemChanged in the constructor 
  // (due to spectrumDisplay->setMaximumRange(minX, maxX, minY, maxY)) which will will there be misinterpreted as a result of
  // spectrumDisplay->setCurrentRange(minX, maxX, minY, maxY)
  // mmm...this smells bad...i need to clean this up
}

SpectrumAnalyzerDisplay::~SpectrumAnalyzerDisplay(void)
{
  deleteAllChildren();
}

void SpectrumAnalyzerDisplay::useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase)
{
  SpectrumDisplayOld::useLogarithmicScaleX(shouldBeLogScaledX, newLogBase);
  MessengingCoordinateSystemOld::useLogarithmicScaleX(shouldBeLogScaledX);
}

void SpectrumAnalyzerDisplay::paint(juce::Graphics &g)
{
  CoordinateSystemOld::paint(g);
  plotCurveFamily(g);
}

void SpectrumAnalyzerDisplay::updateBackgroundImage()
{
  CoordinateSystemOld::updateBackgroundImage();
}

void SpectrumAnalyzerDisplay::plotCurveFamily(Graphics &g, juce::Image* targetImage, XmlElement *targetSVG)
{
  // make sure that the arrays are valid:
  if( familyValuesX == NULL || familyValuesY == NULL )
    return;

  int	   k, i;	           // indices for the loops through the curves and through the values
  double x1, y1, x2, y2;   // for temporary storage of the frequency and amplitude
  double xOld, yOld;
  double lowFreq, highFreq;
  int    minBin, maxBin;


  //Colour graphColour = colourScheme.curves; 
  Colour graphColour = plotColourScheme.getCurveColour(0);  

  for(k=0; k<numCurves; k++)
  {
    juce::String curvePathDataString = juce::String::empty;

    graphColour = plotColourScheme.getCurveColour(k);  
    g.setColour(graphColour); 

    // just for dbug:
    //int y = (k+1)*100;
    //g.drawLine(0, y, getWidth(), y, 16);

    i = 0;
    xOld = 0;
    yOld = getHeight();
    int width;
    if( targetImage == 0 )
      width = getWidth();
    else
      width = targetImage->getWidth();
    while( i < width )
    {
      // get the lower frequency represented by this pixel column:
      x1 = (double) i;
      y1 = 0.0;
      if( targetImage == NULL )
        transformFromComponentsCoordinates(x1, y1);
      else
        transformFromImageCoordinates(x1, y1, targetImage);
      lowFreq = x1;

      // find the upper frequency:
      i++;
      x2 = (double) i;
      y2 = 0.0;
      if( targetImage == NULL )
        transformFromComponentsCoordinates(x2, y2);
      else
        transformFromImageCoordinates(x2, y2, targetImage);
      highFreq = x2;

      while( !getRepresentingBins(lowFreq, highFreq, k, minBin, maxBin)  &&  i < getWidth()-1  )
      {
        x1      = x2;
        y1      = y2;
        lowFreq = x1;

        i++;
        x2 = (double) i;
        y2 = 0.0;
        if( targetImage == NULL )
          transformFromComponentsCoordinates(x2, y2);
        else
          transformFromImageCoordinates(x2, y2, targetImage);

        highFreq = x2;
      }

      // get the x,y- values of the bin representing the range between the lower and upper 
      // frequency and transfrom this to components coordinates
      getRepresentingBins(lowFreq, highFreq, k, minBin, maxBin);

      if( minBin < 0 || minBin > numValues-1 || maxBin < 0 || maxBin > numValues-1 )
        continue;

      x1 = jmax(1.0, familyValuesX[0][minBin]);
      y1 = familyValuesY[k][minBin];
      y1 = amp2dBWithCheck(y1, 0.00001);
      if( targetImage == NULL )
        transformToComponentsCoordinates(x1, y1);
      else
        transformToImageCoordinates(x1, y1, targetImage);

      x2 = jmax(1.0, familyValuesX[0][maxBin]);
      y2 = familyValuesY[k][maxBin];
      y2 = amp2dBWithCheck(y2, 0.00001);
      if( targetImage == NULL )
        transformToComponentsCoordinates(x2, y2);
      else
        transformToImageCoordinates(x2, y2, targetImage);

      if( minBin > maxBin )
      {
        std::swap(x1, x2);
        std::swap(y1, y2);
      }

      // draw a line form the old coordinates (from the previous iteration) and the new ones:
      g.drawLine((float) xOld, (float) yOld, (float) x1, (float) y1);
      g.drawLine((float) x1,   (float) y1,   (float) x2, (float) y2);

      // add the line-segment the svg-drawing also:
      if( targetSVG != NULL )
      {
        curvePathDataString += juce::String(x1) + juce::String(" ") + juce::String(y1) + juce::String(", ");
        curvePathDataString += juce::String(x2) + juce::String(" ") + juce::String(y2) + juce::String(", ");
      }

      // remember the old coordinates for the next iteration:
      xOld = x2;
      yOld = y2;
    }

    // add the generated juce::String (which represents the curve) to the svg-drawing, if the respective
    // flag is true:
    if( targetSVG != NULL )
    {
      // create a XmlElement for the path, setup its attributes and add it ot the svg-drawing:
      XmlElement* curvePath = new XmlElement(juce::String("polyline"));
      curvePath->setAttribute(juce::String("points"), curvePathDataString);
      //Colour currentGraphColour = graphColour;
      //if( curveColours.size() != 0 )
      //  currentGraphColour = *curveColours[k % curveColours.size()];
      curvePath->setAttribute(juce::String("style"), juce::String("stroke-width: ") + juce::String(1.0) + 
        juce::String("; stroke: #") + graphColour.toString().substring(2) + 
        juce::String("; fill: none;") );
      targetSVG->addChildElement(curvePath);
    }
  }
}

bool SpectrumAnalyzerDisplay::getRepresentingBins(double lowFreq, double highFreq, int k, 
  int &minBin, int &maxBin)
{
  int    lowBin  = numValues/2;
  int    highBin = numValues/2;
  int    offset  = lowBin/2;
  double lowBinFreq, highBinFreq;
  while( offset >= 1 )
  {
    if( familyValuesX[0][lowBin] < lowFreq  )
      lowBin += offset;
    else if( familyValuesX[0][lowBin] > lowFreq  ) 
      lowBin -= offset;
    lowBinFreq = familyValuesX[0][lowBin];

    if( familyValuesX[0][highBin] < highFreq  )
      highBin += offset;
    else if( familyValuesX[0][highBin] > highFreq  ) 
      highBin -= offset;
    highBinFreq = familyValuesX[0][highBin];

    offset /= 2;
  }

  while( familyValuesX[0][lowBin] < lowFreq && lowBin <= numValues-1 )
    lowBin++;  
  while( familyValuesX[0][highBin] >= highFreq && highBin >= 0)
    highBin--;

  lowBinFreq  = familyValuesX[0][lowBin];
  highBinFreq = familyValuesX[0][highBin];

  // determine the bin with the maximum magnitude in between (and including) the range 
  // lowBin...highBin and return the index of this bin:
  if( highBin < lowBin )
  {
    minBin = -1;
    maxBin = -1;
    return false;
  }
  else if( highBin == lowBin )
  {
    minBin = maxBin = highBin;

    if( familyValuesX[0][lowBin] < lowFreq || familyValuesX[0][highBin] >= highFreq )
      return false;
    else
      return true;
  }
  /*
  else if( highBin-lowBin <= 2  )
  {
  }
  */
  else
  {
    int length = highBin-lowBin+1;
    offset = rosic::arrayMinIndex(&(familyValuesY[k][lowBin]), length);
    minBin = lowBin + offset;
    offset = rosic::arrayMaxIndex(&(familyValuesY[k][lowBin]), length);
    maxBin = lowBin + offset;
    return true;
  }
}

//=================================================================================================


OscilloscopeDisplay::OscilloscopeDisplay(const juce::String& name) 
  : CurveFamilyPlotOld(name)
{
  setAutoReRendering(false); 
  // we want to manually trigger the re-rendering of the background-image inside this class to avoid unnecesarry calls because we need to 
  // perform some extra actions before rendering here

  // initialize the appearance settings:
  setAxisPositionY(CoordinateSystemOld::ZERO);

  //setMaximumRange(-0.0375, 1.5, -1.5, 1.5);
  //setCurrentRange(-0.0375, 1.5, -1.5, 1.5);
  setMaximumRange(0.0, 1.5, -1.5, 1.5);
  setCurrentRange(0.0, 1.5, -1.5, 1.5);

  setVerticalCoarseGrid(1.0, true);
  setVerticalFineGrid(0.1, true);
  setHorizontalCoarseGrid(1.0, true);
  setHorizontalFineGrid(0.1, true);
  setAxisLabelX(juce::String::empty);
  setAxisLabelY(juce::String::empty);
  setAxisValuesPositionY(CoordinateSystemOld::RIGHT_TO_AXIS);
  //setAxisLabelX(juce::String(T("Time [seconds]")));
  //setAxisLabelY(juce::String(T("Amplitude")));
  setStringConversionForInfoLineX(&secondsToStringWithUnit4);
  setStringConversionForInfoLineY(&amplitudeRawAndInDecibels);
  showPositionAsDescription = true;

  updateBackgroundImage();

  timeAxis = NULL;
  peakData = NULL;
}

void OscilloscopeDisplay::setCurrentRangeX(double newMinX, double newMaxX)
{
  setCurrentRange(newMinX, newMaxX, currentRange.getMinY(), currentRange.getMaxY());
}

void OscilloscopeDisplay::setCurrentRangeY(double newMinY, double newMaxY)
{
  setCurrentRange(currentRange.getMinX(), currentRange.getMaxX(), newMinY, newMaxY);
}

void OscilloscopeDisplay::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  newMaxX = jlimit(0.0001, 1.5, newMaxX);
  newMinX = 0.0 - 0.025*newMaxX;
  MessengingCoordinateSystemOld::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
  adjustGrid();
  updateBackgroundImage();
  //repaint();
}

void OscilloscopeDisplay::setWaveformData(int newNumSamples, int newNumChannels, float** newData, double* newTimeAxis)
{
  numValues = newNumSamples;
  numCurves = newNumChannels;
  peakData  = newData;
  timeAxis  = newTimeAxis;
  updatePlotImage();
}

/*
void OscilloscopeDisplay::resized()
{
CurveFamilyPlotOld::resized();
updatePlotImage(true);
}
*/

void OscilloscopeDisplay::plotCurveFamily(Graphics &g, juce::Image* targetImage, XmlElement *targetSVG)
{
  /*
  // for debug/optimize:
  static int callCount = 0;
  callCount++;
  return;
  // OK - function is called as often as expected
  */



  /*
  if( peakData == NULL ) 
  return;

  double x1, y1, x2, y2;  
  for(int c=0; c<numCurves; c++)
  {
  g.setColour(plotColourScheme.getCurveColour(c)); 
  for(int n=0; n<numValues-1; n++)  
  {
  x1 = timeAxis[n];
  y1 = peakData[c][n];
  x2 = timeAxis[n+1];
  y2 = peakData[c][n+1];
  transformToImageCoordinates(x1, y1, targetImage);
  transformToImageCoordinates(x2, y2, targetImage);
  g.drawLine((float) x1, (float) y1, (float) x2, (float) y2, 1.f); 
  // uses a lot of CPU - try using a juce::Path - nah, Path uses even more
  }
  }
  */



  /*
  double x, y;  
  juce::Path path;
  for(int c = 0; c < numCurves; c++)
  {
  path.clear();
  x = timeAxis[0];
  y = peakData[c][0];  
  transformToImageCoordinates(x, y, targetImage); 
  path.startNewSubPath((float) x, (float) y);
  for(int n = 1; n < numValues; n++)  
  {
  x = timeAxis[n];
  y = peakData[c][n];
  transformToImageCoordinates(x, y, targetImage); 
  path.lineTo((float) x, (float) y);
  }
  g.setColour(plotColourScheme.getCurveColour(c)); 
  g.strokePath(path, PathStrokeType(1.f));
  }
  */
}

void OscilloscopeDisplay::adjustGrid()
{
  // ahhemm...this code is pretty naive and it certainly can be done more clever
  if( currentRange.getMaxX() > 1.1 )
  {
    setVerticalCoarseGridInterval(1.0);
    setVerticalFineGridInterval(0.1);
    setStringConversionForAxisX(&valueToString0);
  }
  else if( currentRange.getMaxX() > 0.11 )
  {
    setVerticalCoarseGridInterval(0.1);
    setVerticalFineGridInterval(0.01);
    setStringConversionForAxisX(&valueToString1);
  }
  else if( currentRange.getMaxX() > 0.011 )
  {
    setVerticalCoarseGridInterval(0.01);
    setVerticalFineGridInterval(0.001);
    setStringConversionForAxisX(&valueToString2);
  }
  else if( currentRange.getMaxX() > 0.0011 )
  {
    setVerticalCoarseGridInterval(0.001);
    setVerticalFineGridInterval(0.0001);
    setStringConversionForAxisX(&valueToString3);
  }
  else if( currentRange.getMaxX() > 0.00011 )
  {
    setVerticalCoarseGridInterval(0.0001);
    setVerticalFineGridInterval(0.00001);
    setStringConversionForAxisX(&valueToString4);
  }

  double yRange = currentRange.getMaxY() - currentRange.getMinY();
  if( yRange > 2.0 )
  {
    setHorizontalCoarseGrid(1.0, true);
    setHorizontalFineGrid(0.25, true);
    setStringConversionForAxisY(&valueToString0);
  }
  else if( yRange > 1.0 )
  {
    setHorizontalCoarseGrid(0.5, true);
    setHorizontalFineGrid(0.125, true);
    setStringConversionForAxisY(&valueToString1);
  }
  else if( yRange > 0.5 )
  {
    setHorizontalCoarseGrid(0.25, true);
    setHorizontalFineGrid(0.0625, true);
    setStringConversionForAxisY(&valueToString2);
  }
  else if( yRange > 0.25 )
  {
    setHorizontalCoarseGrid(0.125, true);
    setHorizontalFineGrid(0.03125, true);
    setStringConversionForAxisY(&valueToString3);
  }
  else if( yRange > 0.125 )
  {
    setHorizontalCoarseGrid(0.0625, true);
    setHorizontalFineGrid(0.015625, true);
    setStringConversionForAxisY(&valueToString4);
  }
}


