//#include "rosic_Noisifier.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Noisifier::Noisifier()
{
  envelopeAmount = 1.0;
  passGain       = 1.0;
  noiseGain      = 1.0;
  reset();
}

//-------------------------------------------------------------------------------------------------
// setup

void Noisifier::setSampleRate(double newSampleRate)
{
  envFollower.setSampleRate(newSampleRate);
}

//-------------------------------------------------------------------------------------------------
// others:

void Noisifier::reset()
{
  envFollower.reset();
  noiseGenerator.trigger();
}

