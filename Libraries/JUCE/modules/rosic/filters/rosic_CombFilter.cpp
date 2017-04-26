#include "rosic_CombFilter.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombFilter::CombFilter(int bufferLengthToAllocate)
{
  if( bufferLengthToAllocate < 512 )
    DEBUG_BREAK;  // you need allocate some reasonable amount of memory

  length         = bufferLengthToAllocate;  // usable/nominal length  
  buffer         = new(std::nothrow) double[length+interpolatorMargin+1];
  tapIn          = 0;
  tapOut         = 0;
  frac           = 0.0;
  yOld           = 0.0;
  sampleRate     = 44100.f;
  dryWetRatio    = 1.0;
  frequency      = 1000.0;  
  feedbackFactor = 0.0;
  wetPolarity    = 1.f;
  delayInSamples = 0.0;
 
  interpolator.setInterpolationMethod(Interpolator::CUBIC_ZERO_DERIVATIVE);

  setupDelayInSamples();
  clearBuffer();
}

CombFilter::~CombFilter()
{
  delete[] buffer;
}

//-------------------------------------------------------------------------------------------------
// setup:

void CombFilter::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
  }
}

void CombFilter::setFrequency(double newFrequency)
{
  if(newFrequency > 0.f)
  {
    frequency = newFrequency;
    setupDelayInSamples();
  }
}

void CombFilter::setInterpolationMethod(int newMethod)
{
  if( newMethod <= Interpolator::WARPED_ALLPASS )
    interpolator.setInterpolationMethod(newMethod);
  else
    DEBUG_BREAK;
}

//-------------------------------------------------------------------------------------------------
// others:

void CombFilter::clearBuffer()
{
  yOld = 0.0;
  interpolator.reset();
  for(int i=0; i<length+interpolatorMargin; i++)
    buffer[i] = 0.0;
}

void CombFilter::setupDelayInSamples()
{
  double delayInSeconds = 0.5 / frequency;
  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = clip(delayInSamples, (double)(interpolatorMargin-1), 
    (double) (length-1-interpolatorMargin));

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
  tapOut = wrapAround(tapOut, length);
}