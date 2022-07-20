template<class TSig, class TPar>
rsQuadratureNetwork<TSig, TPar>::rsQuadratureNetwork() 
{
  order = 12;
  gain  = 1.0;
  designer.setApproximationMethod(rsPrototypeDesigner<TPar>::ELLIPTIC);  
  designer.setMode(rsPrototypeDesigner<TPar>::LOWPASS_PROTOTYPE);
  designer.setPrototypeOrder(order);
  designer.setSampleRate(TPar(44100));
  designer.setFrequency(TPar(11025));
  designer.setRipple(TPar(0.1));
  designer.setStopbandRejection(40.0);

  updateCoefficients();
  reset();
}

// setup:

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::setOrder(int newOrder)
{
  if( newOrder <= maxOrder )
  {
    order = newOrder;
    designer.setPrototypeOrder(newOrder);
    updateCoefficients();
  }
}

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::setRipple(TPar newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::setStopbandRejection(TPar newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

// others:

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::reset()
{
  for(int i = 0; i < maxOrder; i++)
  {
    x[i] = 0.0;
    y[i] = 0.0;
  }
}

template<class TSig, class TPar>
void rsQuadratureNetwork<TSig, TPar>::updateCoefficients()
{
  // design the prototype halfband lowpass filter:
  designer.getPolesAndZeros(poles, zeros);

  // evaluate lowpass transfer function at z=1 and obtain the gain factor:
  TPar desiredGain = TPar(1);
  if( rsIsEven(order) )
    desiredGain = - rsDbToAmp(designer.getPassbandRipple());
  ComplexPar num = TPar(1);
  ComplexPar den = TPar(1);
  for(int i = 0; i < order; i++)
  {
    num *= (TPar(1)-zeros[i]);
    den *= (TPar(1)-poles[i]);
  }
  TPar actualGain = abs(num/den);
  gain = TPar(2)*desiredGain/actualGain; // factor 2 due to cos(x) = (exp(j*x)+exp(-j*x)) / 2

  // rotate pole/zero pattern by 90 degrees counterclockwise:
  ComplexPar j(0.0, 1.0); 
  for(int i = 0; i < order; i++)
  {
    poles[i] *= j;
    zeros[i] *= j;
  }
}


/*

ToDo:
-Plot phase responses of real and imag output
-Try to form linear combinations of re and im outputs - can we obtain arbitrary phse-shifts with 
 that?
-What happens, if we inject the input not only into the real part but with weights into the real
 and imaginary inputs of the filter (using a sin/cos rule such that sum-of-squares of weights is
 always 1)? Maybe that's equivalen to forming linear combinations of the outputs?

Ideas:
-Maybe also implement a bona-fide Hilbert transform wilter, i.e. a filter which takes the original
 signal as is as real part and produces a fitting imaginary part. The impuse response af a digital
 ideal Hilbert filter is h[k] = (1-cos(pi*k)) / (pi*k) = 2/(pi*k) for odd k, 0 for even k. The 
 analog Hilbert filter has h(t) = 1 / (pi*t). Applying the Hilbert transform 3 times yields the 
 inverse transform, so applying it 4 times yields the identity: sin -> cos -> -sin -> -cos -> sin.
 https://en.wikipedia.org/wiki/Hilbert_transform
 https://de.wikipedia.org/wiki/Hilbert-Transformation


*/