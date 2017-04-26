#include "rosic_BasicIntegerDelayLine.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BasicIntegerDelayLine::BasicIntegerDelayLine(int maximumDelayInSamples)
{
  //maximumDelayInSamples += 1;
  maximumDelayInSamples  = nextPowerOfTwo(maximumDelayInSamples);
  bitMask                = maximumDelayInSamples-1;
  delayLine              = new double[(bitMask+1)];
  tapIn                  = 0;
  tapOut                 = 0;
  clearDelayBuffer();
}

BasicIntegerDelayLine::~BasicIntegerDelayLine()
{
  if( delayLine != NULL )
    delete[] delayLine;
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

int BasicIntegerDelayLine::setDelayInSamples(int newDelayInSamples)
{
  int delayInSamples = rmax(0, newDelayInSamples);

  // cap the delay in samples to the maximum allowed delay-time:
  if( delayInSamples > bitMask )
    delayInSamples = bitMask;

  // adjust tapOut-pointer:
  tapOut = tapIn - delayInSamples;
  if( tapOut < 0 )
    tapOut += (bitMask+1);

  return delayInSamples;
}

//-------------------------------------------------------------------------------------------------
// others:

void BasicIntegerDelayLine::clearDelayBuffer()
{
  for(int i=0; i<(bitMask+1); i++)
    delayLine[i] = 0.0;
}