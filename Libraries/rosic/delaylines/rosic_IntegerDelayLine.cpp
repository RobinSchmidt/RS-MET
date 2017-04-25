#include "rosic_IntegerDelayLine.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

IntegerDelayLine::IntegerDelayLine(int maximumDelayInSamples) 
: BasicIntegerDelayLine(maximumDelayInSamples)
{
  delayInSeconds        = 0.001;
  sampleRate            = 44100.0; 
  setDelayInSeconds(delayInSeconds);
}

IntegerDelayLine::~IntegerDelayLine()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void IntegerDelayLine::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
    sampleRate = newSampleRate;
  setDelayInSeconds(delayInSeconds);
}

int IntegerDelayLine::setDelayInSamples(int newDelayInSamples)
{
  int delayInSamples = BasicIntegerDelayLine::setDelayInSamples(newDelayInSamples);
  delayInSeconds     = (double) delayInSamples / sampleRate;  // update the delay in seconds
  return delayInSamples;
}

void IntegerDelayLine::setDelayInSeconds(double newDelayInSeconds)
{
  delayInSeconds = rmax(0.0, newDelayInSeconds);
  setDelayInSamples( roundToInt(sampleRate*delayInSeconds) );
}

void IntegerDelayLine::setDelayInMilliseconds(double newDelayInMilliseconds)
{
  setDelayInSeconds(0.001*newDelayInMilliseconds);
}

