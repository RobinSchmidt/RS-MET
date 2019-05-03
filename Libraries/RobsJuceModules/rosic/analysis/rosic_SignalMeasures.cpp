#include "rosic_SignalMeasures.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SignalMeasures::SignalMeasures()
{
  leftLevel = rightLevel = midLevel = sideLevel = crossCorrelation = 0.0;
}

