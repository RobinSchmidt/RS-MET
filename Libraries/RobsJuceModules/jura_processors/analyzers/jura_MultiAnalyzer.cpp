

//=================================================================================================
// class OscilloscopeAudioModule:

OscilloscopeAudioModule::OscilloscopeAudioModule(CriticalSection *newPlugInLock,
  rosic::OscilloscopeBufferOld* displayBufferToUse)
 : AudioModule(newPlugInLock)
{
  jassert(displayBufferToUse != NULL); // you must pass a valid rosic-object to the constructor
  waveformBuffer = displayBufferToUse;
  setModuleTypeName("Oscilloscope"); // use WaveDisplay
  initializeAutomatableParameters();
  displayIsFrozen = false;
}

void OscilloscopeAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( waveformBuffer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case 0: waveformBuffer->setSyncMode(         (int) value  );  break;
  case 1: waveformBuffer->setMidSideMode(      value != 0.0 );  break;
  case 2: waveformBuffer->setTimeWindowLength( value        );  break;
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
  addObservedParameter( new Parameter(lock, "Freeze",            0.0,   1.0, 1.0,    0.0,  Parameter::BOOLEAN) );

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
  setModuleTypeName("SpectrumAnalyzer");
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
  addObservedParameter( new Parameter(lock, "Freeze",            0.0,        1.0, 1.0,    0.0,    Parameter::BOOLEAN)      );

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================
// class MultiAnalyzerAudioModule:

MultiAnalyzerAudioModule::MultiAnalyzerAudioModule(CriticalSection *newPlugInLock
  //, OscilloscopeAudioModule *newOscilloscopeModule
  //, SpectrumAnalyzerAudioModule *newSpectrumAnalyzerModule
  )
  : AudioModule(newPlugInLock)
{
  //oscilloscopeModule     = newOscilloscopeModule;
  //spectrumAnalyzerModule = newSpectrumAnalyzerModule;

  //jassert( oscilloscopeModule     != NULL );
  //jassert( spectrumAnalyzerModule != NULL );

  //waveformBuffer = new rosic::SyncedWaveformDisplayBuffer;
  waveformBuffer          = new rosic::OscilloscopeBufferOld(400); // displayWidth must be passed - check value
  wrappedSpectrumAnalyzer = new rosic::SpectrumAnalyzer;

  oscilloscopeModule      = new OscilloscopeAudioModule(lock, waveformBuffer);
  spectrumAnalyzerModule  = new SpectrumAnalyzerAudioModule(lock, wrappedSpectrumAnalyzer);

  addChildAudioModule(oscilloscopeModule);
  addChildAudioModule(spectrumAnalyzerModule);

  setModuleTypeName("MultiAnalyzer");
  initializeAutomatableParameters();
}

MultiAnalyzerAudioModule::~MultiAnalyzerAudioModule()
{
  removeChildAudioModule(oscilloscopeModule, true);
  removeChildAudioModule(spectrumAnalyzerModule, true);
  //delete oscilloscopeModule;
  //delete spectrumAnalyzerModule;
  delete waveformBuffer;
  delete wrappedSpectrumAnalyzer;
}

AudioModuleEditor* MultiAnalyzerAudioModule::createEditor(int type)
{
  return new MultiAnalyzerModuleEditor(lock, this);
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
  //p->addStringValue(juce::String("Oscilloscope"));
  p->addStringValue("WaveformDisplay");
  p->addStringValue("SpectrumAnalyzer");
  p->setValue(1.0, false, false);
  addObservedParameter(p);

  //addObservedParameter( new Parameter("FrameRate",  10.0,  50.0, 1.0,  25.0, Parameter::LINEAR)  );

  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

SpectrumAnalyzerDisplay::SpectrumAnalyzerDisplay(const juce::String& name) : rsSpectrumPlot(name)
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
  rsSpectrumPlot::useLogarithmicScaleX(shouldBeLogScaledX, newLogBase);
  rsObservablePlot::useLogarithmicScaleX(shouldBeLogScaledX);
}

void SpectrumAnalyzerDisplay::paint(juce::Graphics &g)
{
  rsPlot::paint(g);
  plotCurveFamily(g);
}

void SpectrumAnalyzerDisplay::updateBackgroundImage()
{
  rsPlot::updateBackgroundImage();
}

void SpectrumAnalyzerDisplay::plotCurveFamily(Graphics &g, juce::Image* targetImage, 
  XmlElement *targetSVG)
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
    juce::String curvePathDataString = juce::String();

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
      fromPixelCoordinates(x1, y1);
      lowFreq = x1;

      // find the upper frequency:
      i++;
      x2 = (double) i;
      y2 = 0.0;
      fromPixelCoordinates(x2, y2);
      highFreq = x2;

      while( !getRepresentingBins(lowFreq, highFreq, k, minBin, maxBin)  &&  i < getWidth()-1  )
      {
        x1      = x2;
        y1      = y2;
        lowFreq = x1;

        i++;
        x2 = (double) i;
        y2 = 0.0;
        fromPixelCoordinates(x2, y2);
        highFreq = x2;
      }

      // get the x,y- values of the bin representing the range between the lower and upper
      // frequency and transfrom this to components coordinates
      getRepresentingBins(lowFreq, highFreq, k, minBin, maxBin);

      if( minBin < 0 || minBin > numValues-1 || maxBin < 0 || maxBin > numValues-1 )
        continue;

      x1 = jmax(1.0, familyValuesX[0][minBin]);
      y1 = familyValuesY[k][minBin];
      y1 = RAPT::rsAmpToDbWithCheck(y1, 0.00001);
      toPixelCoordinates(x1, y1);

      x2 = jmax(1.0, familyValuesX[0][maxBin]);
      y2 = familyValuesY[k][maxBin];
      y2 = RAPT::rsAmpToDbWithCheck(y2, 0.00001);
      toPixelCoordinates(x2, y2);

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

    /*
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
    */
  }
}

bool SpectrumAnalyzerDisplay::getRepresentingBins(double lowFreq, double highFreq, int k,
  int &minBin, int &maxBin)
{
  int    lowBin  = numValues/2;
  int    highBin = numValues/2;
  int    offset  = lowBin/2;
  //double lowBinFreq, highBinFreq;
  while( offset >= 1 )
  {
    if( familyValuesX[0][lowBin] < lowFreq  )
      lowBin += offset;
    else if( familyValuesX[0][lowBin] > lowFreq  )
      lowBin -= offset;
    //lowBinFreq = familyValuesX[0][lowBin];

    if( familyValuesX[0][highBin] < highFreq  )
      highBin += offset;
    else if( familyValuesX[0][highBin] > highFreq  )
      highBin -= offset;
    //highBinFreq = familyValuesX[0][highBin];

    offset /= 2;
  }

  while( familyValuesX[0][lowBin] < lowFreq && lowBin <= numValues-1 )
    lowBin++;
  while( familyValuesX[0][highBin] >= highFreq && highBin >= 0)
    highBin--;

  //lowBinFreq  = familyValuesX[0][lowBin];
  //highBinFreq = familyValuesX[0][highBin];

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
    offset = RAPT::rsArrayTools::minIndex(&(familyValuesY[k][lowBin]), length);
    minBin = lowBin + offset;
    offset = RAPT::rsArrayTools::maxIndex(&(familyValuesY[k][lowBin]), length);
    maxBin = lowBin + offset;
    return true;
  }
}

//=================================================================================================


OscilloscopeDisplay::OscilloscopeDisplay(const juce::String& name)
  : rsDataPlot(name)
{
  setAutoReRendering(false);
  // we want to manually trigger the re-rendering of the background-image inside this class to avoid unnecesarry calls because we need to
  // perform some extra actions before rendering here

  // initialize the appearance settings:
  setAxisPositionY(rsPlotSettings::ZERO);

  //setMaximumRange(-0.0375, 1.5, -1.5, 1.5);
  //setCurrentRange(-0.0375, 1.5, -1.5, 1.5);
  setMaximumRange(0.0, 1.5, -1.5, 1.5);
  setCurrentRange(0.0, 1.5, -1.5, 1.5);

  setVerticalCoarseGrid(1.0, true);
  setVerticalFineGrid(0.1, true);
  setHorizontalCoarseGrid(1.0, true);
  setHorizontalFineGrid(0.1, true);
  setAxisLabelX(juce::String());
  setAxisLabelY(juce::String());
  setAxisValuesPositionY(rsPlotSettings::RIGHT_TO_AXIS);
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
  setCurrentRange(newMinX, newMaxX, plotSettings.getCurrentRangeMinY(), 
    plotSettings.getCurrentRangeMaxY());
}

void OscilloscopeDisplay::setCurrentRangeY(double newMinY, double newMaxY)
{
  setCurrentRange(plotSettings.getCurrentRangeMinX(), plotSettings.getCurrentRangeMaxX(), 
    newMinY, newMaxY);
}

void OscilloscopeDisplay::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  newMaxX = jlimit(0.0001, 1.5, newMaxX);
  newMinX = 0.0 - 0.025*newMaxX;
  rsObservablePlot::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
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
rsDataPlot::resized();
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

  if(peakData == NULL)
    return;

  double x1, y1, x2, y2;
  for(int c = 0; c < numCurves; c++)
  {
    g.setColour(plotColourScheme.getCurveColour(c));
    for(int n = 0; n < numValues-1; n++)
    {
      x1 = timeAxis[n];
      y1 = peakData[c][n];
      x2 = timeAxis[n+1];
      y2 = peakData[c][n+1];
      toPixelCoordinates(x1, y1);
      toPixelCoordinates(x2, y2);
      g.drawLine((float)x1, (float)y1, (float)x2, (float)y2, 1.f);
      // uses a lot of CPU - try using a juce::Path - nah, Path uses even more

      // y stays always fixed and x seems to be way too small (after transformation)
    }
  }

  /*
  double x, y;
  juce::Path path;
  for(int c = 0; c < numCurves; c++)
  {
    path.clear();
    x = timeAxis[0];
    y = peakData[c][0];
    transformToImageCoordinates(x, y, targetImage);
    path.startNewSubPath((float)x, (float)y);
    for(int n = 1; n < numValues; n++)
    {
      x = timeAxis[n];
      y = peakData[c][n];
      transformToImageCoordinates(x, y, targetImage);
      path.lineTo((float)x, (float)y);
    }
    g.setColour(plotColourScheme.getCurveColour(c));
    g.strokePath(path, PathStrokeType(1.f));
  }
  */
}

void OscilloscopeDisplay::adjustGrid()
{
  // ahhemm...this code is pretty naive and it certainly can be done more clever
  if( plotSettings.getCurrentRangeMaxX() > 1.1 )
  {
    setVerticalCoarseGridInterval(1.0);
    setVerticalFineGridInterval(0.1);
    setStringConversionForAxisX(&valueToString0);
  }
  else if( plotSettings.getCurrentRangeMaxX() > 0.11 )
  {
    setVerticalCoarseGridInterval(0.1);
    setVerticalFineGridInterval(0.01);
    setStringConversionForAxisX(&valueToString1);
  }
  else if( plotSettings.getCurrentRangeMaxX() > 0.011 )
  {
    setVerticalCoarseGridInterval(0.01);
    setVerticalFineGridInterval(0.001);
    setStringConversionForAxisX(&valueToString2);
  }
  else if( plotSettings.getCurrentRangeMaxX() > 0.0011 )
  {
    setVerticalCoarseGridInterval(0.001);
    setVerticalFineGridInterval(0.0001);
    setStringConversionForAxisX(&valueToString3);
  }
  else if( plotSettings.getCurrentRangeMaxX() > 0.00011 )
  {
    setVerticalCoarseGridInterval(0.0001);
    setVerticalFineGridInterval(0.00001);
    setStringConversionForAxisX(&valueToString4);
  }

  double yRange = plotSettings.getCurrentRangeMaxY() - plotSettings.getCurrentRangeMinY();
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


//=========================================================================================================================================
// class AudioModuleEditorAnimated:

AudioModuleEditorAnimated::AudioModuleEditorAnimated(CriticalSection *newPlugInLock, AudioModule* newAudioModuleToEdit)
  : AudioModuleEditor(newAudioModuleToEdit)
{
  addWidget( frameRateSlider = new RSlider("FrameRateSlider") );
  frameRateSlider->assignParameter( moduleToEdit->getParameterByName("FrameRate") );

  frameRateSlider->setSliderName(juce::String("FPS"));
  frameRateSlider->setDescription(juce::String("Number of frames per second"));
  frameRateSlider->setDescriptionField(infoField);
  frameRateSlider->addListener(this);
  frameRateSlider->setStringConversionFunction(valueToString0); // \todo: fpsToString

  addWidget( freezeButton = new RButton(juce::String("Freeze")) );
  freezeButton->assignParameter( moduleToEdit->getParameterByName("Freeze") );
  freezeButton->setDescription(juce::String("Freeze the display"));
  freezeButton->setDescriptionField(infoField);
  freezeButton->setClickingTogglesState(true);
  freezeButton->addRButtonListener(this);
  freezeButton->setToggleState(false, false);

  updateWidgetsAccordingToState();
}

void AudioModuleEditorAnimated::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == freezeButton )
    setFreeze(freezeButton->getToggleState());
  else
    AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
}

void AudioModuleEditorAnimated::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  if( sliderThatHasChanged == frameRateSlider )
    setFrameRate(frameRateSlider->getValue());
}

void AudioModuleEditorAnimated::setFrameRate(double newFrameRate)
{
  frameRateSlider->setValue(newFrameRate, false, false);
  updateTimerSettings();
}

void AudioModuleEditorAnimated::setFreeze(bool shouldBeFrozen)
{
  freezeButton->setToggleState(shouldBeFrozen, false);
  updateTimerSettings();
}

void AudioModuleEditorAnimated::timerCallback()
{
  if( freezeButton->getToggleState() != true )
    updateEditorContent();
}

void AudioModuleEditorAnimated::setVisible(bool shouldBeVisible)
{
  AudioModuleEditor::setVisible(shouldBeVisible);
  updateTimerSettings();
}

void AudioModuleEditorAnimated::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  updateTimerSettings();
}

void AudioModuleEditorAnimated::updateTimerSettings()
{
  if( isVisible() && freezeButton->getToggleState() == false )
    startTimer( roundToInt(1000.0 / frameRateSlider->getValue()) );
  else
    stopTimer();
}

void AudioModuleEditorAnimated::resized()
{
  freezeButton->setBounds(getWidth()-52-4, getHeight()-20, 52, 16);
  frameRateSlider->setBounds(freezeButton->getX()-80-4, getHeight()-20, 80, 16);
}

//=========================================================================================================================================
// class OscilloscopeModuleEditor:

OscilloscopeModuleEditor::OscilloscopeModuleEditor(CriticalSection *newPlugInLock, OscilloscopeAudioModule* newOscilloscopeAudioModule)
  : AudioModuleEditorAnimated(newPlugInLock, newOscilloscopeAudioModule)
{
  setHeadlineStyle(Editor::NO_HEADLINE);

  jassert(newOscilloscopeAudioModule != NULL ); // you must pass a valid module here
  oscilloscopeAudioModule = newOscilloscopeAudioModule;

  addWidget( midSideButton = new RButton(juce::String("Mid/Side")) );
  midSideButton->setDescription(juce::String("Switch to Mid/Side mode"));
  midSideButton->assignParameter( moduleToEdit->getParameterByName("MidSideMode") );
  midSideButton->setClickingTogglesState(true);
  midSideButton->setToggleState(false, false);
  midSideButton->addRButtonListener(this);

  addWidget( syncModeComboBox = new RNamedComboBox(juce::String("SyncModeComboBox"), juce::String("Sync:")) );
  syncModeComboBox->setDescription(juce::String("Selects the syncronization mode"));
  syncModeComboBox->assignParameter( moduleToEdit->getParameterByName("SyncMode") );
  syncModeComboBox->setDescriptionField(infoField);
  //syncModeComboBox->addListener(this);

  addPlot( oscilloscopeDisplay = new OscilloscopeDisplay(juce::String("OscilloscopeDisplay")) );
  oscilloscopeDisplay->addCoordinateSystemOldObserver(this);

  addChildColourSchemeComponent( oscilloscopeZoomer = new rsPlotZoomer() );
  oscilloscopeZoomer->setZoomerSize(16);
  oscilloscopeZoomer->setCoordinateSystem(oscilloscopeDisplay);
  oscilloscopeZoomer->hideScrollBarX(true);
  oscilloscopeZoomer->setVerticalMouseWheelMode(rsPlotZoomer::horizontalZoomViaVerticalMouseWheel);

  /*
  timeAxis = NULL;
  xL       = NULL;
  xR       = NULL;
  px       = new float*[2];
  px[0]    = xL;
  px[1]    = xR;
  */

  updateWidgetsAccordingToState();
}

OscilloscopeModuleEditor::~OscilloscopeModuleEditor()
{
  /*
  delete[] timeAxis;
  delete[] xL;
  delete[] xR;
  delete[] px;
  */
}

void OscilloscopeModuleEditor::coordinateSystemChanged(rsObservablePlot *coordinateSystemThatHasChanged)
{
  if( oscilloscopeAudioModule == NULL )
    return;

  if( coordinateSystemThatHasChanged == oscilloscopeDisplay )
  {
    moduleToEdit->getParameterByName("TimeWindowLength")->setValue( oscilloscopeDisplay->getCurrentRangeMaxX(), true, true );
    moduleToEdit->getParameterByName("MinAmplitude")->setValue(     oscilloscopeDisplay->getCurrentRangeMinY(), true, true );
    moduleToEdit->getParameterByName("MaxAmplitude")->setValue(     oscilloscopeDisplay->getCurrentRangeMaxY(), true, true );


    //oscilloscopeAudioModule->wrappedOscilloscope->setTimeWindowLength(oscilloscopeDisplay->getCurrentRangeMaxX());
  }
}

void OscilloscopeModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  double maxX = moduleToEdit->getParameterByName("TimeWindowLength")->getValue();
  double minY = moduleToEdit->getParameterByName("MinAmplitude")->getValue();
  double maxY = moduleToEdit->getParameterByName("MaxAmplitude")->getValue();

  oscilloscopeAudioModule->setIgnoreDirtification(true);       // prevents the subsequent call from spawning a state dirtification
  oscilloscopeDisplay->setCurrentRange(0.0, maxX, minY, maxY); // 1st argument irrelevant
  oscilloscopeZoomer->updateScrollbars();
  oscilloscopeAudioModule->setIgnoreDirtification(false);
}

void OscilloscopeModuleEditor::resized()
{
  AudioModuleEditorAnimated::resized();

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = h - 24;

  midSideButton->setBounds(4, y, 60, 20);
  syncModeComboBox->setBounds(midSideButton->getRight()+4, y, 140, 20);

  y = getHeadlineBottom();
  h = midSideButton->getY()-y;
  oscilloscopeDisplay->setBounds(x, y, w-16, h-16);
  oscilloscopeDisplay->setWaveformData(0, 0, NULL, NULL);
  oscilloscopeZoomer->alignWidgetsToCoordinateSystem();

  /*
  // later - take stereo into account - maybe wrap this stuff into a function
  int N = oscilloscopeAudioModule->waveformBuffer->getDisplayBufferLength();
  delete[] timeAxis;
  delete[] xL;
  delete[] xR;
  timeAxis = new double[N];
  xL       = new float[N];
  xR       = new float[N];
  px[0]    = xL;
  px[1]    = xR;
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getTimeAxis(),      timeAxis, N);
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getDisplayBuffer(), xL,       N);
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getDisplayBuffer(), xR,       N);
  oscilloscopeDisplay->setWaveformData(N, 1, px, timeAxis); // 1 only for debug/optimize - later: 2
  */

  updateEditorContent();
  //int numSamples   = oscilloscopeAudioModule->waveformBuffer->getViewBufferLength();
  //int numChannels  = oscilloscopeAudioModule->waveformBuffer->getNumChannels();
  //double* timeAxis = oscilloscopeAudioModule->waveformBuffer->getTimeAxis();
  //float** values   = oscilloscopeAudioModule->waveformBuffer->getCurrentDisplayBuffers();
  //oscilloscopeDisplay->setWaveformData(numSamples, numChannels, values, timeAxis);
  ////jassertfalse; // we need to adapt the commented code above to work with class OscilloscopeBufferOld
  //// i think, we can just call updateEditorContent()
}

void OscilloscopeModuleEditor::updateEditorContent()
{
  /*
  // todo: let the display buffer directly store float arrays - avoid the copy/conversion
  int N = oscilloscopeAudioModule->waveformBuffer->getDisplayBufferLength();
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getTimeAxis(),      timeAxis, N);
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getDisplayBuffer(), xL,       N);
  convertBuffer(oscilloscopeAudioModule->waveformBuffer->getDisplayBuffer(), xR,       N);
  oscilloscopeDisplay->setWaveformData(N, 2, px, timeAxis);
  */

  oscilloscopeAudioModule->waveformBuffer->updateDisplayBuffers();
  int numSamples   = oscilloscopeAudioModule->waveformBuffer->getViewBufferLength();
  int numChannels  = oscilloscopeAudioModule->waveformBuffer->getNumChannels();
  double* timeAxis = oscilloscopeAudioModule->waveformBuffer->getTimeAxis();
  float** values   = oscilloscopeAudioModule->waveformBuffer->getCurrentDisplayBuffers();
  oscilloscopeDisplay->setWaveformData(numSamples, numChannels, values, timeAxis);
}


//=========================================================================================================================================
// class SpectrumAnalyzerModuleEditor:

SpectrumAnalyzerModuleEditor::SpectrumAnalyzerModuleEditor(CriticalSection *newPlugInLock,
  SpectrumAnalyzerAudioModule* newSpectrumAnalyzerAudioModule)
  : AudioModuleEditorAnimated(newPlugInLock, newSpectrumAnalyzerAudioModule)
{
  setHeadlineStyle(Editor::NO_HEADLINE);

  jassert(newSpectrumAnalyzerAudioModule != NULL ); // you must pass a valid module here
  spectrumAnalyzerAudioModule = newSpectrumAnalyzerAudioModule;

  addWidget( midSideButton = new RButton(juce::String("Mid/Side")) );
  midSideButton->assignParameter( moduleToEdit->getParameterByName("MidSideMode") );
  midSideButton->setDescription(juce::String("Switch to Mid/Side mode"));
  midSideButton->setClickingTogglesState(true);
  midSideButton->setToggleState(false, false);
  midSideButton->addRButtonListener(this);

  addWidget( fftSizeComboBox = new RNamedComboBox(juce::String("fftSizeComboBox"), juce::String("Size:")) );
  fftSizeComboBox->assignParameter( moduleToEdit->getParameterByName("FFTSize") );
  fftSizeComboBox->setDescription(juce::String("Selects the FFT-Size."));
  fftSizeComboBox->setDescriptionField(infoField);

  addWidget( linearFrequencyAxisButton = new RButton(juce::String("Linear")) );
  linearFrequencyAxisButton->assignParameter( moduleToEdit->getParameterByName("LinearFrequency") );
  linearFrequencyAxisButton->setDescription(juce::String("Toggle linear scaling of the frequency axis"));
  linearFrequencyAxisButton->setClickingTogglesState(true);
  linearFrequencyAxisButton->setToggleState(true, false);
  linearFrequencyAxisButton->addRButtonListener(this);

  addPlot( spectrumDisplay = new SpectrumAnalyzerDisplay(juce::String("SpectrumDisplay")) );
  spectrumDisplay->addCoordinateSystemOldObserver(this);
  double minX = moduleToEdit->getParameterByName("MinFrequency")->getMinValue();
  double maxX = moduleToEdit->getParameterByName("MaxFrequency")->getMaxValue();
  double minY = moduleToEdit->getParameterByName("MinLevel")->getMinValue();
  double maxY = moduleToEdit->getParameterByName("MaxLevel")->getMaxValue();
  spectrumDisplay->setMaximumRange(minX, maxX, minY, maxY);
  minX = moduleToEdit->getParameterByName("MinFrequency")->getValue();
  maxX = moduleToEdit->getParameterByName("MaxFrequency")->getValue();
  minY = moduleToEdit->getParameterByName("MinLevel")->getValue();
  maxY = moduleToEdit->getParameterByName("MaxLevel")->getValue();
  spectrumDisplay->setCurrentRange(minX, maxX, minY, maxY);

  addChildColourSchemeComponent( spectrumZoomer = new rsPlotZoomer() );
  spectrumZoomer->setZoomerSize(16);
  spectrumZoomer->setCoordinateSystem(spectrumDisplay);
  spectrumZoomer->setVerticalMouseWheelMode(rsPlotZoomer::horizontalZoomViaVerticalMouseWheel);

  updateWidgetsAccordingToState();
}

void SpectrumAnalyzerModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == linearFrequencyAxisButton )
  {
    spectrumDisplay->useLogarithmicScaleX( !linearFrequencyAxisButton->getToggleState());
    spectrumZoomer->zoomToAllXY();
  }
  else
    AudioModuleEditorAnimated::rButtonClicked(buttonThatWasClicked);
}

void SpectrumAnalyzerModuleEditor::coordinateSystemChanged(rsObservablePlot *coordinateSystemThatHasChanged)
{
  if( spectrumAnalyzerAudioModule == NULL )
    return;

  if( coordinateSystemThatHasChanged == spectrumDisplay )
  {
    // here seems to be a bug - the pointer has a differentv address than spectrumDisplay when the message was send
    // from the spectrumDisplay - maybe something about the address of the inherited MessengingCoordinateSystem?
    moduleToEdit->getParameterByName("MinFrequency")->setValue( spectrumDisplay->getCurrentRangeMinX(), true, true );
    moduleToEdit->getParameterByName("MaxFrequency")->setValue( spectrumDisplay->getCurrentRangeMaxX(), true, true );
    moduleToEdit->getParameterByName("MinLevel")->setValue(     spectrumDisplay->getCurrentRangeMinY(), true, true );
    moduleToEdit->getParameterByName("MaxLevel")->setValue(     spectrumDisplay->getCurrentRangeMaxY(), true, true );
  }
}

void SpectrumAnalyzerModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  double minX = moduleToEdit->getParameterByName("MinFrequency")->getValue();
  double maxX = moduleToEdit->getParameterByName("MaxFrequency")->getValue();
  double minY = moduleToEdit->getParameterByName("MinLevel")->getValue();
  double maxY = moduleToEdit->getParameterByName("MaxLevel")->getValue();

  spectrumAnalyzerAudioModule->setIgnoreDirtification(true);       // prevents the subsequent call from spawning a state dirtification
  spectrumDisplay->setCurrentRange(minX, maxX, minY, maxY);
  spectrumZoomer->updateScrollbars();
  spectrumAnalyzerAudioModule->setIgnoreDirtification(false);
}

void SpectrumAnalyzerModuleEditor::resized()
{
  AudioModuleEditorAnimated::resized();

  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = h - 24;

  midSideButton->setBounds(4, y, 60, 20);
  fftSizeComboBox->setBounds(midSideButton->getRight()+4, y, 80, 20);
  // methodComboBox...
  linearFrequencyAxisButton->setBounds(fftSizeComboBox->getRight()+4, y, 48, 20);

  y = getHeadlineBottom();
  h = midSideButton->getY()-y;
  spectrumDisplay->setBounds(x, y, w-16+RWidget::outlineThickness, h-16); // 16: zoomerSize
  spectrumZoomer->alignWidgetsToCoordinateSystem();
}


void SpectrumAnalyzerModuleEditor::updateEditorContent()
{
  if( spectrumAnalyzerAudioModule == NULL )
    return;
  if( spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer == NULL )
    return;

  // trigger an update of the display-buffers:
  if( !spectrumAnalyzerAudioModule->displayIsFrozen )
    spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer->updateDisplayBuffers();

  // retrieve the data from the analyzer and pass it to the display:
  spectrumDisplay->setSpectra(
    spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer->getNumNonRedundantBins(),
    spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer->getNumChannels(),
    spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer->getBinFrequencies(),
    spectrumAnalyzerAudioModule->wrappedSpectrumAnalyzer->getCurrentSpectra() );
}

//=========================================================================================================================================
// class MultiAnalyzerModuleEditor:

MultiAnalyzerModuleEditor::MultiAnalyzerModuleEditor(CriticalSection *newPlugInLock, MultiAnalyzerAudioModule* newMultiAnalyzerAudioModule)
  : AudioModuleEditor(newMultiAnalyzerAudioModule)
{
  jassert(newMultiAnalyzerAudioModule != NULL ); // you must pass a valid module here
  multiAnalyzerAudioModule = newMultiAnalyzerAudioModule;

  /*
  addWidget( frameRateSlider = new RSlider (T("FrameRateSlider")) );
  frameRateSlider->assignParameter( moduleToEdit->getParameterByName(T("FrameRate")) );
  frameRateSlider->setSliderName(juce::String(T("FPS")));
  frameRateSlider->setDescription(juce::String(T("Number of frames per second")));
  frameRateSlider->setDescriptionField(infoField);
  frameRateSlider->setStringConversionFunction(&rojue::valueToString0); // \todo: fpsToString

  // pull Freeze button into this class also...

  addWidget( freezeButton = new RButton(juce::String(T("Freeze"))) );
  freezeButton->setDescription(juce::String(T("Freeze the display")));
  freezeButton->setDescriptionField(infoField);
  freezeButton->setClickingTogglesState(true);
  freezeButton->addRButtonListener(this);
  freezeButton->setToggleState(false, false);
  //spectrumAnalyzerAudioModule->displayIsFrozen = false;
  */

  // embedded sub-editors:
  oscilloscopeEditor = new OscilloscopeModuleEditor(lock, multiAnalyzerAudioModule->oscilloscopeModule);
  oscilloscopeEditor->setDescriptionField(infoField, true);
  addChildEditor(oscilloscopeEditor, true, false);

  spectrumAnalyzerEditor = new SpectrumAnalyzerModuleEditor(lock, multiAnalyzerAudioModule->spectrumAnalyzerModule);
  spectrumAnalyzerEditor->setDescriptionField(infoField, true);
  addChildEditor(spectrumAnalyzerEditor, true, false);

  // buttons for tabbing between the subeditors:
  addWidget( oscilloscopeButton = new RButton("Wave") );
  oscilloscopeButton->setDescription("Switch to waveform display mode");
  oscilloscopeButton->setDescriptionField(infoField);
  oscilloscopeButton->addRButtonListener(this);

  addWidget( spectrumAnalyzerButton = new RButton("Spectrum") );
  spectrumAnalyzerButton->setDescription("Switch to spectrum analyzer mode");
  spectrumAnalyzerButton->setDescriptionField(infoField);
  spectrumAnalyzerButton->addRButtonListener(this);

  // todo: add meters, spectrogram, pitch

  loadPreferencesFromFile();
  updateWidgetsAccordingToState();

  setSize(600, 300);
}

void MultiAnalyzerModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( multiAnalyzerAudioModule == NULL )
    return;

  if( buttonThatWasClicked == oscilloscopeButton )
    multiAnalyzerAudioModule->setMode(MultiAnalyzerAudioModule::WAVEFORM);
  else if( buttonThatWasClicked == spectrumAnalyzerButton )
    multiAnalyzerAudioModule->setMode(MultiAnalyzerAudioModule::SPECTRUM);
  else
    AudioModuleEditor::rButtonClicked(buttonThatWasClicked);

  updateSubEditorVisibilitiesAndTabButtonStates();
}

void MultiAnalyzerModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int w, h;
  int x = 0;
  int y = getPresetSectionBottom();

  spectrumAnalyzerButton->setBounds(x+4, y+4, 64, 20);
  x = spectrumAnalyzerButton->getRight();
  oscilloscopeButton->setBounds(x+4, y+4, 48, 20);

  /*
  x = getWidth() - 52 - 4;
  freezeButton->setBounds(x, y+4, 52, 16);
  x = freezeButton->getX() - 80 - 4;
  frameRateSlider->setBounds(x, y+4, 80, 16);
  */

  x = 0;
  y = spectrumAnalyzerButton->getBottom() - RWidget::outlineThickness;
  w = getWidth();
  h = infoField->getY()-y;
  oscilloscopeEditor->setBounds(x, y, w, h);
  spectrumAnalyzerEditor->setBounds(x, y, w, h);
}

void MultiAnalyzerModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();

  multiAnalyzerAudioModule->setIgnoreDirtification(true);
  oscilloscopeEditor->updateWidgetsAccordingToState();
  spectrumAnalyzerEditor->updateWidgetsAccordingToState();
  updateSubEditorVisibilitiesAndTabButtonStates();
  multiAnalyzerAudioModule->setIgnoreDirtification(false);
}

void MultiAnalyzerModuleEditor::updateSubEditorVisibilitiesAndTabButtonStates()
{
  oscilloscopeButton->setToggleState(false, false);
  spectrumAnalyzerButton->setToggleState(false, false);

  oscilloscopeEditor->setVisible(false);
  spectrumAnalyzerEditor->setVisible(false);

  int mode = multiAnalyzerAudioModule->getMode();
  switch( mode )
  {
  case MultiAnalyzerAudioModule::WAVEFORM:
  {
    oscilloscopeButton->setToggleState(true, false);
    oscilloscopeEditor->setVisible(true);
  }
  break;
  case MultiAnalyzerAudioModule::SPECTRUM:
  {
    spectrumAnalyzerButton->setToggleState(true, false);
    spectrumAnalyzerEditor->setVisible(true);
  }
  break;
  }
}

