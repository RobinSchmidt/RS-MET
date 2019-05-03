//#include "rosic_AlgoVerb.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AlgoVerb::AlgoVerb()
{
  impulseAmplitude = 0.f;

  dryVol   = 0.0;
  wetVol   = 1.0;
  earlyVol = 1.0;  
  lateVol  = 1.0;;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AlgoVerb::setSampleRate(float newSampleRate)
{
  fdn.setSampleRate(newSampleRate);
}
  
//-------------------------------------------------------------------------------------------------
// others:

void AlgoVerb::reset()
{
  fdn.reset();

}