//#include "rosic_WaveformDisplayBuffer.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

WaveformDisplayBuffer::WaveformDisplayBuffer()
{
  sampleRate          = 44100.0;
  timeWindowLength    = 0.1;
  displayWidth        = 1;
  displayBufferLength = 2*displayWidth;   // times two because of min/max storage

  timeAxisValues = NULL;
  inputBuffer    = NULL;
  displayBuffer  = NULL;
  allocateBuffers();
  clearBuffers();
  setSampleRate(sampleRate);

  xMin =  INF;
  xMax = -INF;
  xOld = 0.0;
  dw   = 0.0;
  nMin = 0;
  nMax = 0;
  n    = 0;
  w    = 0;
  inc  = 0.0;
  wd   = 0.0;
}

WaveformDisplayBuffer::~WaveformDisplayBuffer()
{
  freeBuffers();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void WaveformDisplayBuffer::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  updateTimeVariables();
}

void WaveformDisplayBuffer::setTimeWindowLength(double newTimeWindowLength)
{
  if( newTimeWindowLength > 0.0 )
    timeWindowLength = newTimeWindowLength;
  updateTimeVariables();
}

void WaveformDisplayBuffer::setDisplayWidth(int newDisplayWidth)
{
  if( newDisplayWidth == displayWidth )
    return;
  displayWidth        = RAPT::rsMax(newDisplayWidth, 1);
  displayBufferLength = 2*displayWidth;
  allocateBuffers();
  updateTimeVariables();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

double* WaveformDisplayBuffer::getDisplayBuffer()
{
  updateDisplayBuffer();
  return displayBuffer;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void WaveformDisplayBuffer::updateDisplayBuffer()
{
  RAPT::rsArrayTools::copy(inputBuffer, displayBuffer, displayBufferLength);
  RAPT::rsArrayTools::circularShift(displayBuffer, displayBufferLength, -w);
    // shifts the most recently added min/max pair to the end of the linear buffer
}

void WaveformDisplayBuffer::updateTimeVariables()
{
  numSamplesShown = timeWindowLength*sampleRate;
  inc = displayBufferLength / numSamplesShown;


  RAPT::rsArrayTools::fillWithRangeLinear(timeAxisValues, displayBufferLength, 0.0, (numSamplesShown)/(sampleRate));
  //fillWithRangeLinear(timeAxisValues, displayBufferLength, 0.0, (2*numSamplesShown-1)/(2*sampleRate));
  //fillWithRangeLinear(timeAxisValues, displayBufferLength, 0.0, (numSamplesShown-1)/(sampleRate));


  // \todo: check, if this is exact... might also be (numSamplesShown-1)/(sampleRate) or something
  // maybe use convenient values for the sample-rate (like 100 or something)

  //double test = timeAxisValues[displayBufferLength-1];

  clearBuffers();
}

void WaveformDisplayBuffer::clearBuffers()
{
  RAPT::rsArrayTools::fillWithZeros(inputBuffer,   displayBufferLength);
  RAPT::rsArrayTools::fillWithZeros(displayBuffer, displayBufferLength);
}

void WaveformDisplayBuffer::allocateBuffers()
{
  freeBuffers();
  timeAxisValues = new double[displayBufferLength];
  inputBuffer    = new double[displayBufferLength];
  displayBuffer  = new double[displayBufferLength];
}

void WaveformDisplayBuffer::freeBuffers()
{
  delete[] timeAxisValues;
  delete[] inputBuffer;
  delete[] displayBuffer;
}

void WaveformDisplayBuffer::feedInputBuffer(double *inBuffer, int length)
{
  for(int n = 0; n < length; n++)
    feedInputSample(inBuffer[n]);
}

void WaveformDisplayBuffer::feedInputBuffer(float* inBuffer, int length)
{
  for(int n = 0; n < length; n++)
    feedInputSample((double) inBuffer[n]);
}

void WaveformDisplayBuffer::feedInputSample(double x)
{
  w %= displayBufferLength; // wraparound
  if( numSamplesShown > (double) displayWidth )
  {
    // on-the-fly min/max decimation:
    if( x < xMin )
    {
      xMin = x;
      nMin = n;
    }
    if( x > xMax )
    {
      xMax = x;
      nMax = n;
    }
    dw += inc;
    n  += 1;
    if( dw >= 2.0 )
    {
      if( nMin < nMax )
      {
        inputBuffer[w]   = xMin;
        inputBuffer[w+1] = xMax;
      }
      else
      {
        inputBuffer[w]   = xMax;
        inputBuffer[w+1] = xMin;
      }
      w    += 2;
      dw   -= 2.0;
      n     = 0;
      xMin  =  INF;
      xMax  = -INF;
    }
  }
  else if( numSamplesShown < (double) displayWidth )
  {
    // on-the-fly linear interpolation:
    int start = (int) floor(wd);
    wd        = wd+inc;
    int end   = (int) floor(wd);
    double d  = 0.0;
    double c  = 1.0 / (end-start);
    for(int i = start; i <= end; i++)
    {
      inputBuffer[i % displayBufferLength] = (1.0-d)*xOld + d*x;
      d += c;
    }
    w  = (int) floor(wd) + 1;
    wd = fmod(wd, (double) displayBufferLength);
  }
  else
  {
    // exceptional trivial case - store incoming data as is, using the arithmetic mean in between the samples:
    inputBuffer[w]   = 0.5 * (x + xOld);
    inputBuffer[w+1] = x;
    w += 2;
  }
  xOld = x;
}


//=========================================================================================================================================
// class SyncedWaveformDisplayBuffer:

SyncedWaveformDisplayBuffer::SyncedWaveformDisplayBuffer()
{
  syncMode      = FREE_RUNNING;
  syncThreshold = 1.0;
  xOld2         = 0.0;
  yOld          = 0.0;
  bufferFull    = false;
  syncBackward  = false;

  // initialize the filter for the zero-crossing detection:
  lowpass.setMode(FourPoleFilterParameters::BANDPASS_RBJ);
  lowpass.useTwoStages(true);
  lowpass.setFrequency(30.0);
  lowpass.setQ(1.0/sqrt(2.0));
}

SyncedWaveformDisplayBuffer::~SyncedWaveformDisplayBuffer()
{

}

void SyncedWaveformDisplayBuffer::setSyncMode(int newSyncMode)
{
  syncMode = newSyncMode;
}

void SyncedWaveformDisplayBuffer::setSampleRate(double newSampleRate)
{
  WaveformDisplayBuffer::setSampleRate(newSampleRate);
  lowpass.setSampleRate(sampleRate);
}

double* SyncedWaveformDisplayBuffer::getDisplayBuffer()
{
  if( syncMode == FREE_RUNNING )
    return WaveformDisplayBuffer::getDisplayBuffer();
  else
    return displayBuffer;
}

void SyncedWaveformDisplayBuffer::feedInputSample(double x)
{
  if( syncMode == FREE_RUNNING || (timeWindowLength > syncThreshold) )
    WaveformDisplayBuffer::feedInputSample(x);
  else if( syncMode == ZEROS )
    feedInputSampleWithSync(x, x);
  else if( syncMode == LOWPASS_ZEROS )
    feedInputSampleWithSync(x, lowpass.getSample(x));
}

void SyncedWaveformDisplayBuffer::feedInputSampleWithSync(double x, double y)
{
  // treat case where we sync the last pixel in the display:
  if( syncBackward == true )
  {
    WaveformDisplayBuffer::feedInputSample(x);
    if( yOld <= 0.0 && y > 0.0 )
      updateDisplayBuffer();
    yOld = y;
    return;
  }

  // treat case where we sync the first pixel:
  if( yOld <= 0.0 && y > 0.0 )
  {
    if( bufferFull == true )
    {
      w    = 0;
      wd   = 0.0;
      xOld = xOld2; // dunno why, but it is needed to make it work
    }
    else
    {
      // maybe do a circular shift? or not?
    }
    bufferFull = false;
  }
  xOld2 = x;
  yOld  = y;

  if( bufferFull == false && !isNextIndexBehindBufferEnd() )
    WaveformDisplayBuffer::feedInputSample(y);
  else
  {
    if( bufferFull == false ) // seems to mess up 1st sample in N=400 cases (interpolating)
    {
      double b0 = inputBuffer[0];
      double b1 = inputBuffer[1];

      WaveformDisplayBuffer::feedInputSample(y);
      wd = w = displayBufferLength; // counteract wraparound in feedInputSample
      RAPT::rsArrayTools::copy(inputBuffer, displayBuffer, displayBufferLength);

      displayBuffer[0] = b0;
      displayBuffer[1] = b1;
      bufferFull = true;
    }
  }
}

bool SyncedWaveformDisplayBuffer::isNextIndexBehindBufferEnd()
{
  if( numSamplesShown < (double) displayWidth )
    return (wd+inc) > (double) (displayBufferLength-1);
  else
    return (w+2) > (displayBufferLength-1);
}
