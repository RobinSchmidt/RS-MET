#include "rosic_InstantaneousEnvelopeDetector.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

InstantaneousEnvelopeDetector::InstantaneousEnvelopeDetector()
{
  quadratureNetwork.setApproximationMethod(PrototypeDesigner::ELLIPTIC);
  quadratureNetwork.setOrder(15);
  quadratureNetwork.setRipple(0.1);
  quadratureNetwork.setStopbandRejection(40.0);
  reset();
}

InstantaneousEnvelopeDetector::~InstantaneousEnvelopeDetector()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

