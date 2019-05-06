rsWaveformPlot::rsWaveformPlot(const String& name)
: rsDataPlot(name)
{
  sampleRate      = 44100.0;
  timeFormat      = HOUR_MINUTE_SECOND;
  //timeFormat      = SAMPLES;
  numChannels     = 0;
  numSampleFrames = 0;
  peakData        = NULL;

  // initialize the range - time-window: 1 ms, amplitude: -1.1...+1.1
  setMaximumRange(0.0, 0.001, -1.2, +1.2);
  setCurrentRange(0.0, 0.001, -1.0, +1.0);

  // intialize the peakData with some dummy-data:
  double dummyWave[4];
  dummyWave[0] = 0.0;
  dummyWave[1] = 0.0;
  dummyWave[2] = 0.0;
  dummyWave[3] = 0.0;
  double*  p   = &dummyWave[0];
  double** pp  = &p;
  setWaveform(pp, 4, 1);
  // ...this is ugly - try to get rid of that

  //setValueFieldPopup(false);
  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);
}

rsWaveformPlot::~rsWaveformPlot(void)
{
  if( peakData != NULL )
    delete[] peakData;
  //deleteAllChildren();
}

rsPlotRange rsWaveformPlot::getMaximumMeaningfulRange(
  double relativeMarginLeft, double relativeMarginRight,
  double relativeMarginTop,  double relativeMarginBottom)
{
  // require at least 1.0% margin:
  jassert( relativeMarginLeft   >= 0.0 );
  jassert( relativeMarginRight  >= 0.0 );
  jassert( relativeMarginTop    >= 0.0 );
  jassert( relativeMarginBottom >= 0.0 );
  if( relativeMarginLeft   < 0.0 )
    relativeMarginLeft   = 0.0;
  if( relativeMarginRight  < 0.0 )
    relativeMarginRight  = 0.0;
  if( relativeMarginTop    < 0.0 )
    relativeMarginTop    = 0.0;
  if( relativeMarginBottom < 0.0 )
    relativeMarginBottom = 0.0;
  // use jmax

  rsPlotRange r = getMaximumRange();

  double minX   = 0.0;
  double maxX   = jmax(0.001, (double) (numSampleFrames-1) / sampleRate);
  // ToDo: encapsulate this in maximumRange, and support different time-axis units

  double width  = maxX-minX;
  double minY   = r.getMinY();
  double maxY   = r.getMaxY();
  double height = maxY-minY;

  r.setMinX(minX - 0.01*relativeMarginLeft   * width);
  r.setMaxX(maxX + 0.01*relativeMarginRight  * width);
  r.setMinY(minY - 0.01*relativeMarginBottom * height);
  r.setMaxY(maxY + 0.01*relativeMarginTop    * height);

  return r;

  // this function needs some sophistication....
}

void rsWaveformPlot::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  repaint();
}

bool rsWaveformPlot::setWaveform(double** newWaveformData, int newNumSampleFrames,
                                  int newNumChannels)
{
  // copy the pointers into member variables:
  numSampleFrames = newNumSampleFrames;
  numChannels     = newNumChannels;

  if( numSampleFrames <= 0 || numChannels <= 0 )
    return false;

  // create the decimated (peak) data, which is to be displayed when the display-width
  // in pixels is smaller than the number of data-points:
  int peakArraySize = 2*nextPowerOfTwo(newNumSampleFrames);
  if( peakData != NULL )
    delete[] peakData;
  peakData = new float[numChannels*peakArraySize];

  // catch memory allocation errors:
  if( peakData == NULL )
  {
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, String("Memory Allocation Failed"),
      String("Memory allocation failed in rsWaveformPlot::setWaveform"), String("OK") );
    numSampleFrames = 0;
    numChannels     = 0;
    return false;
  }

  // copy (and cast) the undecimated data inot the first sub-array:
  int c, n; // channel and sample-number
  for(c=0; c<newNumChannels; c++)
  {
    for(n=0; n<numSampleFrames; n++)
      peakData[c*peakArraySize+n] = (float) newWaveformData[c][n];
  }

  // for the rest of the array copy the last value in order to have a well defined value which does
  // not confuse or corrupt the min/max calculation (leaving these values undefined may result in
  // rendering of ridiciously long lines):
  for(c=0; c<newNumChannels; c++)
  {
    for(n=numSampleFrames; n<peakArraySize; n++)
      peakData[c*peakArraySize+n] = (float) newWaveformData[c][numSampleFrames-1];
  }

  // create the decimated versionds of the data and report sucess:
  createDecimatedData();
  return true;
}

bool rsWaveformPlot::setWaveform(float** newWaveformData, int newNumSampleFrames,
                                  int newNumChannels)
{
  // same procedure as the corresponding function for doubles:
  numSampleFrames = newNumSampleFrames;
  numChannels     = newNumChannels;
  if( numSampleFrames <= 0 || numChannels <= 0 )
    return false;
  int peakArraySize = 2*nextPowerOfTwo(newNumSampleFrames);
  if( peakData != NULL )
    delete[] peakData;
  peakData = new float[numChannels*peakArraySize];
  if( peakData == NULL )
  {
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, String("Memory Allocation Failed"),
      String("Memory allocation failed in rsWaveformPlot::setWaveform"), String("OK") );
    numSampleFrames = 0;
    numChannels     = 0;
    return false;
  }
  int c, n;
  for(c=0; c<newNumChannels; c++)
  {
    for(n=0; n<numSampleFrames; n++)
      peakData[c*peakArraySize+n] = (float) newWaveformData[c][n];
  }
  for(c=0; c<newNumChannels; c++)
  {
    for(n=numSampleFrames; n<peakArraySize; n++)
      peakData[c*peakArraySize+n] = newWaveformData[c][numSampleFrames-1];
  }
  createDecimatedData();
  return true;
}

bool rsWaveformPlot::setWaveform(const AudioSampleBuffer& newWaveformBuffer)
{
  // copy the pointers into member variables:
  numSampleFrames = newWaveformBuffer.getNumSamples();
  numChannels     = newWaveformBuffer.getNumChannels();

  if( numSampleFrames <= 0 || numChannels <= 0 )
    return false;

  // create the decimated (peak) data, which is to be displayed when the display-width
  // in pixels is smaller than the number of data-points:
  int peakArraySize = 2*nextPowerOfTwo(numSampleFrames);
  if( peakData != NULL )
    delete[] peakData;
  peakData = new float[numChannels*peakArraySize];

  // copy (and cast) the undecimated data inot the first sub-array:
  int c, n;           // channel and sample-number
  const float* readPointer;
  for(c=0; c<numChannels; c++)
  {
    readPointer = newWaveformBuffer.getReadPointer(c);

    for(n=0; n<numSampleFrames; n++)
      peakData[c*peakArraySize+n] = (float) readPointer[n];
  }

  // for the rest of the array copy the last value in order to have a well defined value which does
  // not confuse or corrupt the min/max calculation (leaving these values undefined may result in
  // rendering of ridiciously long lines):
  for(c=0; c<numChannels; c++)
  {
    readPointer = newWaveformBuffer.getReadPointer(c);
    for(n=numSampleFrames; n<peakArraySize; n++)
      peakData[c*peakArraySize+n] = readPointer[numSampleFrames-1];
  }

  createDecimatedData();

  return true;
}

void rsWaveformPlot::createDecimatedData()
{
  int peakArraySize     = 2*nextPowerOfTwo(numSampleFrames);

  // outer loop over the channels:
  for(int c=0; c<numChannels; c++)
  {
    int   readOffset       = c*peakArraySize;               // here starts the sub-array from which we read
    int   writeOffset      = readOffset + peakArraySize/2;  // here we start to write the next sub-array
    int   nextSubArraySize = peakArraySize/4;
    int   nRead            = readOffset;
    int   nWrite           = writeOffset;
    int   minIndex, maxIndex;
    float minValue, maxValue;

    // inner loop over the samples:
    while( writeOffset < readOffset + peakArraySize && nextSubArraySize >= 2 )
    {
      nRead  = readOffset;
      nWrite = writeOffset;
      while( nWrite < writeOffset+nextSubArraySize-2 )
      {
        minIndex = arrayMinIndex(&(peakData[nRead]), 4);
        maxIndex = arrayMaxIndex(&(peakData[nRead]), 4);
        minValue = peakData[nRead+minIndex];
        maxValue = peakData[nRead+maxIndex];

        if( minIndex < maxIndex )
        {
          peakData[nWrite]   = minValue;
          peakData[nWrite+1] = maxValue;
        }
        else
        {
          peakData[nWrite]   = maxValue;
          peakData[nWrite+1] = minValue;
        }

        nRead  += 4;
        nWrite += 2;
      }
      readOffset        = writeOffset; // we read where we just wrote
      writeOffset      += nextSubArraySize;
      nextSubArraySize /= 2;
    }
  }

  // it's not very elegant to have this here - but we want to start at fully zoomed out view, 
  // whenever we get new waveform data via one of the setWaveform functions and this function gets
  // called by all of them:
  double timeScaler = 1.0;
  if(timeFormat != SAMPLES) 
    timeScaler = 1.0/sampleRate;     // time axis is in seconds seconds
  setMaximumRangeX(0.0, timeScaler * (numSampleFrames+1) ); // why +1?
  setCurrentRangeX(0.0, timeScaler * (numSampleFrames+1) );
}

/*
void rsWaveformPlot::paint(juce::Graphics &g)
{
  CoordinateSystem::paint(g);
  plotWaveform(g);
}
*/

void rsWaveformPlot::plotCurveFamily(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  plotWaveform(g, targetImage, targetSVG);
}

void rsWaveformPlot::plotWaveform(Graphics &g, Image *targetImage, XmlElement *targetSVG)
{
  //return; // test
  // test:
  //g.fillAll(Colours::black); // has no effect in the SamplePlayer ..the background is gray
  //return;                    // even, if we just fill black and return - still gray - wtf?

  // make sure that the arrays are valid:
  if( peakData == nullptr || numSampleFrames <= 0 || numChannels <= 0 )
    return;


  int firstVisibleFrame = (int) floor(plotSettings.getCurrentRangeMinX() * sampleRate);
  int lastVisibleFrame  = (int) ceil(plotSettings.getCurrentRangeMaxX() * sampleRate);
  firstVisibleFrame     = jlimit(0, numSampleFrames-1, firstVisibleFrame);
  lastVisibleFrame      = jlimit(0, numSampleFrames-1, lastVisibleFrame);
  int numVisibleFrames  = lastVisibleFrame - firstVisibleFrame;
  int decimationFactor  = (int) floor((double) numVisibleFrames / (double) (4*getWidth()) );
  decimationFactor      = jmax(1, decimationFactor);
  int subArrayIndex     = (int) floor(log2(decimationFactor));
  decimationFactor      = (int) pow(2.0, (double) subArrayIndex);
  int peakArraySize     = 2*nextPowerOfTwo(numSampleFrames);
  //int subArraySize      = peakArraySize/(2*decimationFactor);
  double timeScaler = 1.0;
  if(timeFormat != SAMPLES) 
    timeScaler = 1.0/sampleRate;     // time axis is in seconds seconds

  //Colour graphColour    = Colours::blue;
  bool  drawDots        = numVisibleFrames <= 0.1 * getWidth();
  float dotRadius       = 3.f;

  // adjust the pointer at which we want to read out from the array, i.e. skip
  // forward to the appropriate sub-array:
  int readOffset  = 0;
  int offsetToAdd = peakArraySize/2;
  for(int i=0; i<subArrayIndex; i++)
  {
    readOffset  += offsetToAdd;
    offsetToAdd /= 2;
  }

  // outer loop over the channels:
  for(int c=0 ; c<numChannels; c++)
  {
    // set the colour for the current channel:
    Colour graphColour = plotColourScheme.getCurveColour(c);
    g.setColour(graphColour);

    // calculate the range for sample-index n in order to draw only the actually visible range:
    int nMin = (int) floor((double) firstVisibleFrame / (double) decimationFactor);
    int nMax = (int) ceil((double) lastVisibleFrame / (double) decimationFactor);
    nMin = jlimit(0, numSampleFrames-1, nMin);
    nMax = jlimit(0, numSampleFrames-1, nMax);

    // inner loop over the samples:
    for(int n=nMin; n<nMax; n++)
    {
      double x1, y1, x2, y2;

      // read out the tables:
      x1 = (float) ( timeScaler * (double) n*decimationFactor);
      y1 = peakData[c*peakArraySize + readOffset+n];
      x2 = (float) ( timeScaler * (double) ((n+1)*decimationFactor));
      y2 = peakData[c*peakArraySize + readOffset+n+1];

      // transform:
      toPixelCoordinates(x1, y1);
      toPixelCoordinates(x2, y2);

      // draw:
      g.drawLine((float)x1, (float)y1, (float)x2, (float)y2, 2.f);
      if( drawDots == true )
      {
         g.fillEllipse((float) (x1-dotRadius), (float) (y1-dotRadius),
           (float) (2*dotRadius), (float) (2*dotRadius) );
      }
    }
  }
}

void rsWaveformPlot::updateWaveformCurve()
{
  /*
  if( numSampleFrames == 0 || numChannels == 0 )
    return;

  CoordinateSystemRange r;
  r.setMinX(-0.1*numSampleFrames);
  r.setMaxX(numSampleFrames+0.1*numSampleFrames);
  r.setMinY(-1.2);
  r.setMaxY(+1.2);
  setRange(r);

  setCoarseGrid(true, numSampleFrames/8.0, true, 1.0/2.0);
  setFineGrid(true, numSampleFrames/32.0, true, 1.0/16.0);

  // something more to do here.....

  repaint();
  */
}
