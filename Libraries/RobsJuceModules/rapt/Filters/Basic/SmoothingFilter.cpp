template<class TSig, class TPar>
rsSmoothingFilter<TSig, TPar>::rsSmoothingFilter()
{
  y1.resize(1);
  coeffs.resize(1);
  reset();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate)
{
  decay = sampleRate * timeConstant;
  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setOrder(int newOrder)
{
  order = rsMax(1, newOrder);

  y1.resize(order);
  coeffs.resize(order);

  reset();
  // todo: if newOrder > oldOrder, init only the vector values in y1 above oldOrder-1 to 0
  // not all of them

  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setNumSamplesToReachHalf(TPar numSamples)
{
  decay = numSamples * TPar(LN2_INV);
  updateCoeffs();
}

//template<class TSig, class TPar>
//void rsSmoothingFilter<TSig, TPar>::setShape(int newShape)
//{
//  shape = newShape;
//  updateCoeffs();
//}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setShapeParameter(TPar newParam)
{
  shapeParam = newParam;
  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::reset()
{
  setStates(0);
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setStates(TSig value)
{
  for(int i = 0; i < order; i++)
    y1[i] = value;
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::updateCoeffs()
{
  TPar tmp;
  if(shapeParam != 0)
  {
    for(int i = 0; i < order; i++)
    {
      tmp  = decay / (TPar) pow(i+1, shapeParam);  // tau[n] = tau[0] / n^p // p == shapeParam
      //if(i > 0)
      //  tmp *= (shapeParam+1);
      coeffs[i] = exp(-order/tmp); 
    }
    // maybe try, if it responds different to modulations of the time-constants are in reverse
    // order (from short to long instead of long to short)
  }
  else
  {
    // all filter stages use the same time-constant
    tmp = exp(-order/decay); // amounts to divide the time-constant by the order
    for(int i = 0; i < order; i++)
      coeffs[i] = tmp;
  }
}


/*

ToDo:

- Make an envelope follower based on this filter for use in dynamics processors

- The env-follower should be availbale as modulation source in chainer

- Look inot this  https://cytomic.com/files/dsp/DynamicSmoothing.pdf  for an idea of a dynamic
  smoothing filter. The results are indeed smoother that for a normal filter with fixed cutoff. The
  idea is that in a multimode filter the (absolute value of) the bandpass output can be seen as a 
  measure for how much the signal is changing (why?) and that is used to modulate the cutoff 
  frequency of the filter. If the signal is changing a lot, the opens the filter. Or put another 
  way, when not much is happening, the filter can use a lower cutoff and therby smooth more 
  aggressively. It is tested there on a noisy signal with steps in it like a typical control 
  signal.

- In Notes/SmoothingFilter.txt are some more ideas.

*/