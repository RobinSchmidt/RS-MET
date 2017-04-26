#include "rosic_BesselFilterForGainSignal.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BesselFilterForGainSignal::BesselFilterForGainSignal()
{
  reset();
}

//-------------------------------------------------------------------------------------------------
// others:

void BesselFilterForGainSignal::reset()
{
  for(int i=0; i<8; i++)
    w[i] = 0.0;
}

