// construction/destruction:

InstantaneousEnvelopeDetector::InstantaneousEnvelopeDetector()
{
  quadratureNetwork.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  quadratureNetwork.setOrder(15);
  quadratureNetwork.setRipple(0.1);
  quadratureNetwork.setStopbandRejection(40.0);
  reset();
}

InstantaneousEnvelopeDetector::~InstantaneousEnvelopeDetector()
{

}
