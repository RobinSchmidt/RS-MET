
//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveformDisplay::WaveformDisplay(AudioFileBuffer *newBuffer)
: InteractiveCoordinateSystem("WaveformDisplay"), AudioFileBufferUser(newBuffer)
{
  ScopedLock pointerLock(audioFileBufferPointerLock);

  Component::setName(String("WaveformDisplay"));
  firstChannelToPlot = 0;
  lastChannelToPlot  = 7; // restrict plot to 8 channels by default
  //setValueFieldPopup(false);

  if( bufferToUse != NULL )
  {
    double maxTime = jmax(bufferToUse->getLengthInSeconds(), 0.0001); // 1/10 ms minimum range
    setMaximumRangeX(0.0, maxTime);
    setCurrentRangeX(0.0, maxTime);
    minVisibleTime = 0.0;
    maxVisibleTime = maxTime;
    //numChannels    = bufferToUse->getNumChannels();
    //numSamples     = bufferToUse->getNumSamples();
    //sampleRate     = bufferToUse->getFileSampleRate(); // mmmhh. could be used directly in plot...
  }
  else
  {
    // use some fallback-values when we got a NULL pointer
    setMaximumRange(0.0, 0.001, -1.0, +1.0);
    setCurrentRange(0.0, 0.001, -1.0, +1.0);
    minVisibleTime = 0.0;
    maxVisibleTime = 0.001;
    //numChannels    = 0;
    //numSamples     = 0;
    //sampleRate     = 44100.0;
  }
}

//-------------------------------------------------------------------------------------------------
// setup:

void WaveformDisplay::assignAudioFileBuffer(AudioFileBuffer *newBuffer)
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( newBuffer != bufferToUse )
  {
    AudioFileBufferUser::assignAudioFileBuffer(newBuffer);
    setRangeToBufferLength();
    /*
    if( bufferToUse != NULL )
    {
      double maxTime = jmax(bufferToUse->getLengthInSeconds(), 0.0001); // 1/10 ms minimum range
      setMaximumRangeX(0.0, maxTime);
      setCurrentRangeX(0.0, maxTime);
      minVisibleTime = 0.0;
      maxVisibleTime = maxTime;
      //numChannels    = bufferToUse->getNumChannels();
      //numSamples     = bufferToUse->getNumSamples();
      //sampleRate     = bufferToUse->getFileSampleRate(); // mmmhh. could be used directly in plot...
    }
    */
  }
}

void WaveformDisplay::setRangeToBufferLength()
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
  {
    double maxTime = jmax(bufferToUse->getLengthInSeconds(), 0.0001); // 1/10 ms minimum range
    setMaximumRangeX(0.0, maxTime);
    setCurrentRangeX(0.0, maxTime);
    minVisibleTime = 0.0;
    maxVisibleTime = maxTime;
  }
}


/*
void WaveformDisplay::setSampleRate(double newSampleRate)
{
  jassert(newSampleRate > 0.0);  // zero or negative sample-rates are not supported ;-)
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  repaint();
}
*/

void WaveformDisplay::setFirstChannelToPlot(int newFirstChannelIndex)
{
  jassert(newFirstChannelIndex >= 0); // negative channel-index? no! counting starts at 0
  if( newFirstChannelIndex >= 0 )
    firstChannelToPlot = newFirstChannelIndex;
}

void WaveformDisplay::setLastChannelToPlot(int newLastChannelIndex)
{
  jassert(newLastChannelIndex >= 0); // negative channel-index? no! counting starts at 0
  if( newLastChannelIndex >= 0 )
    lastChannelToPlot = newLastChannelIndex;
}

void WaveformDisplay::plotOnlyOneChannel(int channelToPlotIndex)
{
  setFirstChannelToPlot(channelToPlotIndex);
  setLastChannelToPlot(channelToPlotIndex);
}

void WaveformDisplay::setVisibleTimeRange(double newMinTimeInSeconds, double newMaxTimeInSeconds)
{
  minVisibleTime = newMinTimeInSeconds;
  maxVisibleTime = newMaxTimeInSeconds;
}

//-------------------------------------------------------------------------------------------------
// others:

/*
void WaveformDisplay::resized()
{
  Panel::resized();
  ThreadedDrawingComponent::resized();
}

void WaveformDisplay::setDirty(bool shouldBeSetToDirty)
{
  Panel::setDirty(shouldBeSetToDirty);
  ThreadedDrawingComponent::setDirty(shouldBeSetToDirty);
}
*/

//-------------------------------------------------------------------------------------------------
// drawing:

void WaveformDisplay::drawComponent(Image* imageToDrawOnto)
{
  Graphics g(*imageToDrawOnto);

  drawCoordinateSystem(g, imageToDrawOnto);

  //g.fillAll(Colours::white); // preliminary - call CoordinateSystem background drawing stuff here

  plotWaveform(imageToDrawOnto);
  g.drawRect(0, 0, imageToDrawOnto->getWidth(), imageToDrawOnto->getHeight(), 1);
}

void WaveformDisplay::plotWaveform(Image *targetImage)
{
  if( targetImage == NULL )
    return;

  Graphics g(*targetImage);

  // make sure that the pointer to the bufferToUse is not modified during this function:
  bool pointerLockAcquired = audioFileBufferPointerLock.tryEnter();
  if( pointerLockAcquired == false )
    return;

  // check if we have a valid buffer - if not, return - we don't need to aquire the readlock for
  // the data here, because we use getMinMaxSamples to retrieve the data which in itself aquires
  // the read-lock:
  if( bufferToUse == NULL )
	{
		audioFileBufferPointerLock.exit();
		return;
	}

  //ScopedReadLock dataLock(bufferToUse->audioDataReadWriteLock);
	bufferToUse->acquireReadLock();
  // we have our read-lock - do the work:

  int    numChannels = bufferToUse->getNumChannels();
  int    numSamples  = bufferToUse->getNumSamples();
  double sampleRate  = bufferToUse->getFileSampleRate();

  double pixelWidth = (double) getWidth();
  if( targetImage != NULL )
    pixelWidth = (double) targetImage->getWidth();

  double tSecMin  = currentRange.getMinX();           // min-time in seconds
  double tSecMax  = currentRange.getMaxX();           // max-time in seconds
  double tSecInc  = (tSecMax-tSecMin) / pixelWidth;   // increment per pixel
  double tSmpMin  = tSecMin * sampleRate;             // min-time in samples
  double tSmpMax  = tSecMax * sampleRate;             // max-time in samples
  double tSmpInc  = (tSmpMax-tSmpMin) / pixelWidth;   // increment per pixel
  int cMin        = jmax(0, firstChannelToPlot);
  int cMax        = jmin(numChannels-1, lastChannelToPlot);
  //int curveDrawn  = 0;
  //float dotRadius = 3.f;

  // outer loop over the channels:
  double tSec, tSmp;
  double x1, x2, y1, y2;
  int    nMin, nMax;
	g.setColour(Colours::blue);
  for(int c=cMin ; c<=cMax; c++)
  {
    // inner loop over pixels:
    tSec = tSecMin;
    tSmp = tSmpMin;
    while( tSmp <= tSmpMax )
    {
      nMin = jlimit(0, numSamples-1, (int) floor(tSmp)      );
      nMax = jlimit(0, numSamples-1, (int) ceil( tSmp+tSmpInc) );
      x1   = (float) tSec;
      x2   = (float) (tSec+tSecInc);

      //bufferToUse->getMinMaxSamples(c, nMin, nMax-nMin+1, y1, y2);
      bufferToUse->getMinMaxSamplesWithoutLock(c, nMin, nMax-nMin+1, y1, y2);

      toPixelCoordinates(x1, y1);
      toPixelCoordinates(x2, y2);
      g.drawLine((float)x1, (float)y1, (float)x2, (float)y2);

      tSec += tSecInc;
      tSmp += tSmpInc;
    }
  }

	bufferToUse->releaseReadLock();
  audioFileBufferPointerLock.exit();
}

void WaveformDisplay::restrictToVisibleSection(double &tMin, double &tMax,
                                               double marginInPercent) const
{
  double marginAbsolute = 0.01*marginInPercent*(maxVisibleTime-minVisibleTime);
  tMin = jmax(tMin, minVisibleTime-marginAbsolute);
  tMax = jmin(tMax, maxVisibleTime+marginAbsolute);
}
