#include "rosic_Flanger.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Flanger::Flanger()
{
  pitch  = freqToPitch(1000.0);
  reset();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Flanger::setSampleRate(double newSampleRate)
{
  ModulationEffect::setSampleRate(newSampleRate);
  combL.setSampleRate(newSampleRate);
  combR.setSampleRate(newSampleRate);
}

//-------------------------------------------------------------------------------------------------
// others:

void Flanger::reset()
{
  combL.clearBuffer();
  combR.clearBuffer();
}

