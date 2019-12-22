//#include "rosic_OscilloscopeBufferOld.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OscilloscopeBufferOld::OscilloscopeBufferOld(int displayWidthInPixels)
{
  sampleRate          = 44100.0;
  timeWindowStart     = 0.0;
  timeWindowLength    = 0.1;
  syncThreshold       = 1.0;  
  displayWidth        = displayWidthInPixels;
  decimationFactor    = 1;
  syncMode            = FREE_RUNNING;
  numChannels         = 2;
  viewBufferLength    = 2*displayWidth;   // times two because of min/max storage
  sampleCounter       = 0;
  midSideMode         = false;

  // create the buffers:
  int c;
  circularBufferFlat = new float[maxNumChannels*maxBufferSize];
  circularBuffer     = new float*[maxNumChannels];
  for(c=0; c<maxNumChannels; c++)
    circularBuffer[c] = &(circularBufferFlat[c*maxBufferSize]);

  linearBufferFlat = new float[maxNumChannels*maxBufferSize];
  linearBuffer     = new float*[maxNumChannels];
  for(c=0; c<maxNumChannels; c++)
    linearBuffer[c] = &(linearBufferFlat[c*maxBufferSize]);

  viewBufferFlat = new float[maxNumChannels*viewBufferLength];
  viewBuffer     = new float*[maxNumChannels];
  for(c=0; c<maxNumChannels; c++)
    viewBuffer[c] = &(viewBufferFlat[c*viewBufferLength]);

  filteredCircularBuffer = new float[maxBufferSize];
  timeAxisValues         = new double[viewBufferLength];

  clearBuffers();

  // initialize the filter for the zero-crossing detection:
  lowpass.setMode(FourPoleFilterParameters::BANDPASS_RBJ);
  lowpass.useTwoStages(true);
  lowpass.setFrequency(30.0);
  lowpass.setQ(1.0/sqrt(2.0));

  setSampleRate(sampleRate);

  calculateTimeAxisValues();
}

OscilloscopeBufferOld::~OscilloscopeBufferOld()
{
  delete[] circularBufferFlat;
  delete[] circularBuffer;
  delete[] linearBufferFlat;
  delete[] linearBuffer;
  delete[] viewBufferFlat;
  delete[] viewBuffer;
  delete[] filteredCircularBuffer;
  delete[] timeAxisValues;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void OscilloscopeBufferOld::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  lowpass.setSampleRate(sampleRate);
}

void OscilloscopeBufferOld::setTimeWindowLength(double newTimeWindowLength)
{
  if( newTimeWindowLength > 0.0 )
    timeWindowLength = newTimeWindowLength;

  if( (int) ceil(timeWindowLength*sampleRate) > maxBufferSize )
    timeWindowLength = floor((double)maxBufferSize/sampleRate);

  //calculateTimeAxisValues();
}

void OscilloscopeBufferOld::setTimeWindowStart(double newTimeWindowStart)
{
  if( newTimeWindowStart > 0.0 )
    timeWindowStart = newTimeWindowStart;

  //calculateTimeAxisValues();
}

void OscilloscopeBufferOld::setDisplayWidth(int newDisplayWidth)
{
  if( newDisplayWidth == displayWidth )
    return;


  // to actually use this, we need a mutex

  /*
  displayWidth     = newDisplayWidth;
  viewBufferLength = 2*displayWidth;

  delete[] viewBufferFlat;
  delete[] viewBuffer;
  viewBufferFlat = new float[maxNumChannels*viewBufferLength];
  viewBuffer     = new float*[maxNumChannels];
  for(c=0; c<maxNumChannels; c++)
    viewBuffer[c] = &(viewBufferFlat[c*viewBufferLength]);
  */
}

void OscilloscopeBufferOld::setSyncMode(int newSyncMode)
{
  syncMode = newSyncMode;
}

//-------------------------------------------------------------------------------------------------
// others:

void OscilloscopeBufferOld::clearBuffers()
{
  for(int c=0; c<numChannels; c++)
  {
    int n;
    for(n=0; n<maxBufferSize; n++)
    {
      circularBuffer[c][n]      = 0.0;
      linearBuffer[c][n]        = 0.0;
      filteredCircularBuffer[n] = 0.0;
    }
    for(n=0; n<viewBufferLength; n++)
    {
      viewBuffer[c][n]  = 0.0;
      timeAxisValues[n] = 0.0;
    }
  }
  lowpass.reset();
}

void OscilloscopeBufferOld::updateDisplayBuffers()
{
  int c, nr, nw; // channel, read and write positions

  // get the number of samples, we must look into the past (this function call will also set up
  // the member-variable 'decimationFactor'):
  int hindsight = getHindsight();

  // copy some chunk of the circular buffer into the linear buffer for peak-extraction:
  int readStart = sampleCounter - hindsight;
  for(c=0; c<numChannels; c++)
  {
    for(nw=0; nw<(hindsight+1); nw++)
    {
      nr = wrapAround(readStart + nw);  
      linearBuffer[c][nw] = circularBuffer[c][nr];  
    }
  }

  // copy the peaks (min and max values) of the linear buffer into the peak array:
  if( decimationFactor == 1 ) 
  {
    // interpolate the linear buffer:
    for(c=0; c<numChannels; c++)
    {
      for(nw=0; nw<displayWidth; nw++)
      {
        viewBuffer[c][2*nw]   = linearBuffer[c][nw];
        viewBuffer[c][2*nw+1] = 0.5f * (linearBuffer[c][nw] + linearBuffer[c][nw+1]);
      }
    }
  }
  else if( decimationFactor == 2 ) 
  {
    // just copy from the linear buffer:
    for(c=0; c<numChannels; c++)
    {
      for(nw=0; nw<viewBufferLength; nw++)
        viewBuffer[c][nw] = linearBuffer[c][nw];
    }
  }
  else 
  {
    // put the min- and max- values into the peak-buffer:
    int   minIndex, maxIndex;
    float minValue, maxValue;
    for(c=0; c<numChannels; c++)
    {
      for(nw=0; nw<displayWidth; nw++)
      {
        nr      = decimationFactor*nw;
        minIndex = RAPT::rsArrayTools::minIndex(&(linearBuffer[c][nr]), decimationFactor);
        maxIndex = RAPT::rsArrayTools::maxIndex(&(linearBuffer[c][nr]), decimationFactor);
        minValue = linearBuffer[c][nr+minIndex];
        maxValue = linearBuffer[c][nr+maxIndex];

        if( minIndex < maxIndex )
        {
          viewBuffer[c][2*nw]   = minValue;
          viewBuffer[c][2*nw+1] = maxValue;
        }
        else
        {
          viewBuffer[c][2*nw]   = maxValue;
          viewBuffer[c][2*nw+1] = minValue;
        }
      }
    }
  }
}

void OscilloscopeBufferOld::calculateTimeAxisValues()
{
  double timeIncrement = (double) decimationFactor / sampleRate;
  for(int n=0; n<displayWidth; n++)
  {
    timeAxisValues[2*n]   = timeWindowStart + (double) n * timeIncrement;
    timeAxisValues[2*n+1] = timeAxisValues[2*n] + 0.5*timeIncrement;
  }
}

int OscilloscopeBufferOld::getHindsight()
{
  int hindsight = (int) ceil(sampleRate*timeWindowLength);
  hindsight     = RAPT::rsClip(hindsight, displayWidth, maxBufferSize-1);

  // update the decimation factor (integer division is approriate here, hindsight is the length 
  // of the time-slice in samples):
  decimationFactor = (hindsight / displayWidth) + 1; // why +1 ?

  // refine the hindsight to the next zero crossing of a lowpassed version of the signal, if sync
  // is active:
  if( syncMode == LOWPASS_ZEROS )
  {
    bool hindsightHasBeenRefined = false;
    int  offset    = 0;
    int  maxOffset = (int) floor(0.2*sampleRate); // restricts the zero-crossing search
                                                  // inside some reasonable time-window
    while( hindsightHasBeenRefined == false )
    {
      if( filteredCircularBuffer[wrapAround(sampleCounter-hindsight-offset)]   >= 0.0 && 
          filteredCircularBuffer[wrapAround(sampleCounter-hindsight-offset-1)] <  0.0 )
      {
        hindsight              += offset;
        hindsightHasBeenRefined = true;
      }
      else if( offset >= maxOffset )
        hindsightHasBeenRefined = true; // fallback if we don't find a zero crossing inside a resonable time window
      else
        offset++;
    }
  }

  // update the time-axis (why here? - refactor -  maybe this causes the gum-effect?):
  calculateTimeAxisValues();

  return hindsight;
}

/*
\todo use a better algorithm with on-the-fly decimation/interpolation
W: display-width in pixels (integer)
N: number of samples in time-window (real)
d = N/W: decimation-factor (samples per pixel)
r = W/N: resolution (pixels per sample)
x: current incoming sample
-maintain a circular buffer of length W
-use w as buffer index and n as sample index
-we have to consider 3 cases:

1: N = W -> d = r = 1, easiest case but rare
buf[w] = x;
w = wrapAround(w+1, W);
n = wrapAround(n+1, N);  // not needed?

2: N > W -> d > 1, r < 1, we must min/max decimate for display
if( x < xMin )
 xMin = x;
 nMin = n;
endif
if( x > xMax )
 xMax = x;
 nMax = n;
endif
dw = dw+r;      // increment for w
if( dw >= 2.0 )
  buf[w]   = xMin;
  buf[w+1] = xMax; // hmm...may write beyond end of buffer - we'll see
  // ... needs to be refined to take into account min before max or max before min

  w  = wrapAround(floor(w+dw), W);
  dw = dw - 2.0;
endif
n = wrapAround(n+1, N); // not needed?


3: N < W -> d < 1, r > 1, we must interpolate for display



*/
