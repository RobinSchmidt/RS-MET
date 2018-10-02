//-------------------------------------------------------------------------------------------------
// construction/destruction:

QuadratureNetwork::QuadratureNetwork() 
{
  order = 12;
  gain  = 1.0;
  designer.setApproximationMethod(rsPrototypeDesigner::ELLIPTIC);  
  designer.setMode(rsPrototypeDesigner::LOWPASS_PROTOTYPE);
  designer.setPrototypeOrder(order);
  designer.setSampleRate(44100.0);
  designer.setFrequency(11025.0);
  designer.setRipple(0.1);
  designer.setStopbandRejection(40.0);

  updateCoefficients();
  reset();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void QuadratureNetwork::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients();
}

void QuadratureNetwork::setOrder(int newOrder)
{
  if( newOrder <= maxOrder )
  {
    order = newOrder;
    designer.setPrototypeOrder(newOrder);
    updateCoefficients();
  }
}

void QuadratureNetwork::setRipple(double newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

void QuadratureNetwork::setStopbandRejection(double newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

//-------------------------------------------------------------------------------------------------
// others:

void QuadratureNetwork::reset()
{
  for(int i=0; i<maxOrder; i++)
  {
    x[i] = 0.0;
    y[i] = 0.0;
  }
}

void QuadratureNetwork::updateCoefficients()
{
  // design the prototype halfband lowpass filter:
  designer.getPolesAndZeros(poles, zeros);

  // evaluate lowpass transfer function at z=1 and obtain the gain factor:
  double desiredGain = 1.0;
  if( RAPT::rsIsEven(order) )
    desiredGain = - RAPT::rsDbToAmp(designer.getPassbandRipple());
  Complex num = 1.0;
  Complex den = 1.0;
  for(int i=0; i<order; i++)
  {
    num *= (1.0-zeros[i]);
    den *= (1.0-poles[i]);
  }
  double actualGain = (num/den).getRadius();
  gain = 2.0*desiredGain/actualGain; // factor 2 due to cos(x) = (exp(j*x)+exp(-j*x)) / 2

  // rotate pole/zero pattern by 90 degrees counterclockwise:
  Complex j(0.0, 1.0); 
  for(int i=0; i<order; i++)
  {
    poles[i] *= j;
    zeros[i] *= j;
  }
}
