//-------------------------------------------------------------------------------------------------
// construction/destruction:

DualWaveformDisplay::DualWaveformDisplay(AudioFileBuffer* newBuffer) 
: Component("DualWaveformDisplay"), AudioFileBufferUser(newBuffer)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  //lockUsedBufferPointer();

  addAndMakeVisible( waveDisplayL = new WaveformDisplay(newBuffer) );
  //waveDisplayL->setAxisPositionX(CoordinateSystem::INVISIBLE);
  //waveDisplayL->setAxisPositionY(CoordinateSystem::INVISIBLE);  
  waveDisplayL->plotOnlyOneChannel(0);

  addAndMakeVisible( waveDisplayR = new WaveformDisplay(newBuffer) );
  //waveDisplayR->setAxisPositionX(CoordinateSystem::INVISIBLE);
  //waveDisplayR->setAxisPositionY(CoordinateSystem::INVISIBLE);  
  waveDisplayR->plotOnlyOneChannel(1);

  //unlockUsedBufferPointer();
}

DualWaveformDisplay::~DualWaveformDisplay()
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  //lockUsedBufferPointer();
  deleteAllChildren();
  //unlockUsedBufferPointer();
}

//-------------------------------------------------------------------------------------------------
// setup:

void DualWaveformDisplay::assignAudioFileBuffer(AudioFileBuffer* newBuffer)
{
  lockUsedBufferPointer();
  AudioFileBufferUser::assignAudioFileBuffer(newBuffer);
  waveDisplayL->assignAudioFileBuffer(newBuffer);
  waveDisplayR->assignAudioFileBuffer(newBuffer);
  unlockUsedBufferPointer();

  resized(); // to update the arrangement of one or two displays
}

//-------------------------------------------------------------------------------------------------
// appearance stuff:

void DualWaveformDisplay::resized()
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  lockUsedBufferPointer();

  int x = 0;
  int w = getWidth();
  int h = getHeight();

  /*
  // handle the situation when no audio buffer is assigned:
  if( bufferToUse == NULL )
  {
    waveDisplayL->setBounds(0,   0, w, h  );
    waveDisplayR->setBounds(0, h/2, w, h/2);
    return;
  }
  */

  int numChannels = waveDisplayL->getNumChannels();

  // show either one waveform display over the full height (for mono clips) or two displays of
  // half height (for stereo clips)
  if( numChannels < 2 )
  {
    waveDisplayL->setBounds(0, 0, w, h);
    waveDisplayR->setBounds(0, 0, 0, 0);
  }
  else
  {
    waveDisplayL->setBounds(0,   0, w, h/2);
    waveDisplayR->setBounds(0, h/2, w, h/2);
  }

  unlockUsedBufferPointer();
}

//-------------------------------------------------------------------------------------------------
// the CoordinateSystem mimics:

double DualWaveformDisplay::getMaximumRangeMinX() const
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  return waveDisplayL->getMaximumRangeMinX();
}

double DualWaveformDisplay::getMaximumRangeMaxX() const
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  return waveDisplayL->getMaximumRangeMaxX();
}

double DualWaveformDisplay::getCurrentRangeMinX() const
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  return waveDisplayL->getCurrentRangeMinX();
}

double DualWaveformDisplay::getCurrentRangeMaxX() const
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  return waveDisplayL->getCurrentRangeMaxX();
}

void DualWaveformDisplay::setMaximumRange(double newMinX, double newMaxX, 
                                          double newMinY, double newMaxY)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setMaximumRange(newMinX, newMaxX, newMinY, newMaxY);
  waveDisplayR->setMaximumRange(newMinX, newMaxX, newMinY, newMaxY);
}

void DualWaveformDisplay::setMaximumRangeX(double newMinX, double newMaxX)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setMaximumRangeX(newMinX, newMaxX);
  waveDisplayR->setMaximumRangeX(newMinX, newMaxX);
}

void DualWaveformDisplay::setMaximumRangeY(double newMinY, double newMaxY)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setMaximumRangeY(newMinY, newMaxY);
  waveDisplayR->setMaximumRangeY(newMinY, newMaxY);
}

void DualWaveformDisplay::setCurrentRange(double newMinX, double newMaxX,                                               
                                          double newMinY, double newMaxY)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
  waveDisplayR->setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
}

void DualWaveformDisplay::setCurrentRangeX(double newMinX, double newMaxX)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setCurrentRangeX(newMinX, newMaxX);
  waveDisplayR->setCurrentRangeX(newMinX, newMaxX);
}

void DualWaveformDisplay::setCurrentRangeMinX(double newMinX)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setCurrentRangeMinX(newMinX);
  waveDisplayR->setCurrentRangeMinX(newMinX);
}

void DualWaveformDisplay::setCurrentRangeMaxX(double newMaxX)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setCurrentRangeMaxX(newMaxX);
  waveDisplayR->setCurrentRangeMaxX(newMaxX);
}

void DualWaveformDisplay::setCurrentRangeY(double newMinY, double newMaxY)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setCurrentRangeY(newMinY, newMaxY);
  waveDisplayR->setCurrentRangeY(newMinY, newMaxY);
}

void DualWaveformDisplay::setDrawingThread(TimeSliceThread* newDrawingThread)
{
  waveDisplayL->setDrawingThread(newDrawingThread);
  waveDisplayR->setDrawingThread(newDrawingThread);
}

void DualWaveformDisplay::setVisibleTimeRange(double newMinTimeInSeconds, 
                                              double newMaxTimeInSeconds)
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setVisibleTimeRange(newMinTimeInSeconds, newMaxTimeInSeconds);
  waveDisplayR->setVisibleTimeRange(newMinTimeInSeconds, newMaxTimeInSeconds);
}

void DualWaveformDisplay::updatePlotImage()
{
  //ScopedLock pointerLock(audioFileBufferPointerLock);
  waveDisplayL->setDirty();
  waveDisplayR->setDirty();
}

void DualWaveformDisplay::lockUsedBufferPointer()
{
  AudioFileBufferUser::lockUsedBufferPointer();
  waveDisplayL->lockUsedBufferPointer();
  waveDisplayR->lockUsedBufferPointer();
}

void DualWaveformDisplay::unlockUsedBufferPointer()
{
  waveDisplayR->unlockUsedBufferPointer();
  waveDisplayL->unlockUsedBufferPointer();
  AudioFileBufferUser::unlockUsedBufferPointer();
}

void DualWaveformDisplay::acquireAudioFileBufferWriteLock()
{
  AudioFileBufferUser::acquireAudioFileBufferWriteLock();
  waveDisplayL->acquireAudioFileBufferWriteLock();
  waveDisplayR->acquireAudioFileBufferWriteLock();
}

void DualWaveformDisplay::releaseAudioFileBufferWriteLock()
{
  waveDisplayR->releaseAudioFileBufferWriteLock();
  waveDisplayL->releaseAudioFileBufferWriteLock();
  AudioFileBufferUser::releaseAudioFileBufferWriteLock();
}