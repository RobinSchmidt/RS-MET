template<class TSig, class TPar>
QuadratureNetwork<TSig, TPar>::QuadratureNetwork() 
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

// setup:

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients();
}

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::setOrder(int newOrder)
{
  if( newOrder <= maxOrder )
  {
    order = newOrder;
    designer.setPrototypeOrder(newOrder);
    updateCoefficients();
  }
}

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::setRipple(TPar newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::setStopbandRejection(TPar newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

// others:

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::reset()
{
  for(int i = 0; i < maxOrder; i++)
  {
    x[i] = 0.0;
    y[i] = 0.0;
  }
}

template<class TSig, class TPar>
void QuadratureNetwork<TSig, TPar>::updateCoefficients()
{
  // design the prototype halfband lowpass filter:
  designer.getPolesAndZeros(poles, zeros);

  // evaluate lowpass transfer function at z=1 and obtain the gain factor:
  TPar desiredGain = 1.0;
  if( isEven(order) )
    desiredGain = - dB2amp(designer.getPassbandRipple());
  ComplexPar num = 1.0;
  ComplexPar den = 1.0;
  for(int i = 0; i < order; i++)
  {
    num *= (1.0-zeros[i]);
    den *= (1.0-poles[i]);
  }
  double actualGain = (num/den).getRadius();
  gain = 2.0*desiredGain/actualGain; // factor 2 due to cos(x) = (exp(j*x)+exp(-j*x)) / 2

  // rotate pole/zero pattern by 90 degrees counterclockwise:
  ComplexPar j(0.0, 1.0); 
  for(int i = 0; i < order; i++)
  {
    poles[i] *= j;
    zeros[i] *= j;
  }
}
