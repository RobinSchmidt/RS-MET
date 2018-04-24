using namespace RSLib;

// Construction/Destruction:

rsBasicDelayLine::rsBasicDelayLine()
{
  maxDelay  = 3;
  delayLine = new double[maxDelay+1];
  tapIn     = 0;
  tapOut    = 0;
  reset();
}

rsBasicDelayLine::~rsBasicDelayLine()
{
  if( delayLine != NULL )
    delete[] delayLine;
}

// Setup:

void rsBasicDelayLine::setMaximumDelayInSamples(int newMaxDelay)
{
  if( newMaxDelay > maxDelay )
  {
    delete[] delayLine;
    maxDelay  = rsNextPowerOfTwo(newMaxDelay + 1) - 1;
    delayLine = new double[maxDelay+1];
  }
}

void rsBasicDelayLine::setDelayInSamples(int newDelay)
{
  int delay = rsMax(0, newDelay);
  if( delay > maxDelay )
    setMaximumDelayInSamples(delay);

  // adjust tapOut-pointer:
  tapOut = tapIn - delay;
  if( tapOut < 0 )
    tapOut += maxDelay+1;
}

// Misc:

void rsBasicDelayLine::reset()
{
  for(int i = 0; i < maxDelay+1; i++)  
    delayLine[i] = 0.0;
}

//-------------------------------------------------------------------------------------------------

// Construction/Destruction:

rsDelayLine::rsDelayLine() 
{
  delayInSeconds = 0.001;
  sampleRate     = 44100.0; 
  setDelayInSeconds(delayInSeconds);
}

rsDelayLine::~rsDelayLine()
{

}

// Setup:

void rsDelayLine::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  setDelayInSeconds(delayInSeconds);
}

void rsDelayLine::setDelayInSamples(int newDelayInSamples)
{
  rsBasicDelayLine::setDelayInSamples(newDelayInSamples);
  delayInSeconds = (double) newDelayInSamples / sampleRate;
}

void rsDelayLine::setDelayInSeconds(double newDelayInSeconds)
{
  delayInSeconds = rsMax(0.0, newDelayInSeconds);
  setDelayInSamples( rsRoundToInt(sampleRate*delayInSeconds) );
}

void rsDelayLine::setDelayInMilliseconds(double newDelayInMilliseconds)
{
  setDelayInSeconds(0.001*newDelayInMilliseconds);
}

//-------------------------------------------------------------------------------------------------

// construction/destruction:

rsFractionalDelayLine::rsFractionalDelayLine(int maximumDelayInSamples)
{
  length      = maximumDelayInSamples + 1;
  delayBuffer = new double[length+interpolatorMargin];

  tapIn      = 0;
  tapOut     = 0;
  delayTime  = 0.25;
  sampleRate = 44100.0;
  bpm        = 120.0;
  tempoSync  = false;

  interpolator.setInterpolationMethod(rsInterpolator::WARPED_ALLPASS);

  interpolator.setInterpolationMethod(rsInterpolator::LINEAR);  // for debug

  setDelayTime(delayTime);
  clearDelayBuffer();
}

rsFractionalDelayLine::~rsFractionalDelayLine()
{
  if( delayBuffer != NULL )
  {
    delete[] delayBuffer;
    delayBuffer = NULL;
  }
}

// parameter settings (set-functions):

void rsFractionalDelayLine::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
  }
}

void rsFractionalDelayLine::setDelayTime(double newDelayTime)
{
  //delayTime = rsClip(newDelayTime, 0.0, 4.25);  
  delayTime = newDelayTime;
  setupDelayInSamples();
}

void rsFractionalDelayLine::setSyncMode(bool shouldTempoSync)
{
  tempoSync = shouldTempoSync;
  setupDelayInSamples();
}

void rsFractionalDelayLine::setTempoInBPM(double newTempoInBPM)
{
  if(newTempoInBPM >= 0.0)
  {
    bpm = newTempoInBPM;
    setupDelayInSamples();
  }
  else
    rsError("Tempo < 0");
}

void rsFractionalDelayLine::setInterpolationMethod(int newMethod)
{
  if(newMethod <= rsInterpolator::WARPED_ALLPASS)
    interpolator.setInterpolationMethod(newMethod);
  else
    rsError("Unknown interpolation method");
}

// others:

void rsFractionalDelayLine::clearDelayBuffer()
{
  for(int i=0; i<length+interpolatorMargin; i++)
    delayBuffer[i] = 0.0;
  interpolator.reset();
}

void rsFractionalDelayLine::setupDelayInSamples()
{
  double delayInSeconds;
  if( tempoSync )
    delayInSeconds = rsBeatsToSeconds(delayTime, bpm);
  else
    delayInSeconds = delayTime;
    // get rid of the temposync stuff

  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = rsLimitToRange(delayInSamples, (double)(interpolatorMargin-1), 
    (double) (length-1-interpolatorMargin));
  //delayInSamples = clip(delayInSamples, (double)(interpolatorMargin-1), 
  //  (double) (length-1-interpolatorMargin)); // old


  // update member delayTime to the clipped value:
  delayInSeconds = delayInSamples / sampleRate;   
  if( tempoSync )
    delayTime = rsSecondsToBeats(delayInSeconds, bpm);
  else
    delayTime = delayInSeconds;

  // calculate the integer and fractional parts of the delay:
  double tmp   = floor(delayInSamples);
  int    dInt  = (int) tmp;
  double dFrac = delayInSamples - tmp;
  frac         = 1.0 - dFrac; // because we look backwards

                              // adjust tapOut-pointer:
  tapOut = tapIn - dInt - 1;
  if( frac >= 1.0 )
  {
    frac    = 0.0;
    tapOut += 1;
  }
  tapOut = wrapAround(tapOut);
}
