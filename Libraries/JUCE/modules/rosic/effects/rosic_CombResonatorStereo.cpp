#include "rosic_CombResonatorStereo.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombResonatorStereo::CombResonatorStereo()
{
  pitch        = freqToPitch(1000.0);
  pitchOffsetL = 0.0;
  pitchOffsetR = 0.0;
  pan1         = 0.0;
  pan2         = 1.0;
  gain         = 1.0;
  reset();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void CombResonatorStereo::setSampleRate(double newSampleRate)
{
  combL.setSampleRate(newSampleRate);
  combR.setSampleRate(newSampleRate);
}

//-------------------------------------------------------------------------------------------------
// others:

void CombResonatorStereo::reset()
{
  combL.clearBuffer();
  combR.clearBuffer();
}

