//#include "rosic_FractionalDelayLine.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FractionalDelayLine::FractionalDelayLine(int maximumDelayInSamples)
{
  length      = maximumDelayInSamples + 1;
  delayBuffer = new double[length+interpolatorMargin];


  tapIn      = 0;
  tapOut     = 0;
  delayTime  = 0.25;
  sampleRate = 44100.0;
  bpm        = 120.0;
  tempoSync  = false;

  interpolator.setInterpolationMethod(Interpolator::WARPED_ALLPASS);

  interpolator.setInterpolationMethod(Interpolator::LINEAR);  // for debug

  setDelayTime(delayTime);

  clearDelayBuffer();
}

FractionalDelayLine::~FractionalDelayLine()
{
  if( delayBuffer != NULL )
  {
    delete[] delayBuffer;
    delayBuffer = NULL;
  }
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void FractionalDelayLine::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
  }
}

void FractionalDelayLine::setDelayTime(double newDelayTime)
{
  delayTime = clip(newDelayTime, 0.0, 4.25);  
  setupDelayInSamples();
}

void FractionalDelayLine::setSyncMode(bool shouldTempoSync)
{
  tempoSync = shouldTempoSync;
  setupDelayInSamples();
}

void FractionalDelayLine::setTempoInBPM(double newTempoInBPM)
{
  if( newTempoInBPM >= 0.0 )
  {
    bpm = newTempoInBPM;
    setupDelayInSamples();
  }
  else
    DEBUG_BREAK;
}

void FractionalDelayLine::setInterpolationMethod(int newMethod)
{
  if( newMethod <= Interpolator::WARPED_ALLPASS )
    interpolator.setInterpolationMethod(newMethod);
  else
    DEBUG_BREAK;
}

//-------------------------------------------------------------------------------------------------
// others:

void FractionalDelayLine::clearDelayBuffer()
{
  for(int i=0; i<length+interpolatorMargin; i++)
    delayBuffer[i] = 0.0;
  interpolator.reset();
}

void FractionalDelayLine::setupDelayInSamples()
{
  double delayInSeconds;
  if( tempoSync )
    delayInSeconds = RAPT::rsBeatsToSeconds(delayTime, bpm);
  else
    delayInSeconds = delayTime;

  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = clip(delayInSamples, (double)(interpolatorMargin-1), 
    (double) (length-1-interpolatorMargin));

  // update member delayTime to the clipped value:
  delayInSeconds = delayInSamples / sampleRate;   
  if( tempoSync )
    delayTime = secondsToBeats(delayInSeconds, bpm);
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