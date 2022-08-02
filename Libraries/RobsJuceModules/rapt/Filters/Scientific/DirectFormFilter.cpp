// construction/destruction:

template<class TSig, class TCoef>
rsDirectFormFilter<TSig, TCoef>::rsDirectFormFilter(int maximumOrder)
{
  maxOrder = maximumOrder;
  order    = maxOrder;
  w        = new TSig[order+1];
  a        = new TCoef[order+1];
  b        = new TCoef[order+1];
  initializeCoefficients();
  reset();  
}

template<class TSig, class TCoef>
rsDirectFormFilter<TSig, TCoef>::~rsDirectFormFilter()
{
  delete[] w;
  delete[] a;
  delete[] b;
}

// parameter settings:

template<class TSig, class TCoef>
void rsDirectFormFilter<TSig, TCoef>::setCoefficients(TCoef* newCoeffsA, TCoef* newCoeffsB, 
  int newOrder)
{
  initializeCoefficients();  // why that?
  order = newOrder;
  for(int k = 0; k <= newOrder; k++)
  {
    a[k] = newCoeffsA[k];
    b[k] = newCoeffsB[k];
  }
}

template<class TSig, class TCoef>
void rsDirectFormFilter<TSig, TCoef>::setGlobalGainFactor(TCoef newFactor)
{
  for(int k = 0; k <= order; k++)
    b[k] *= newFactor;
}

template<class TSig, class TCoef>
TCoef rsDirectFormFilter<TSig, TCoef>::getMagnitudeResponseAt(TCoef omega)
{
  // todo: move computation into rsFilterAnalyzer

  TCoef ca = 0.0;
  TCoef sa = 0.0;
  TCoef cb = 0.0; 
  TCoef sb = 0.0;
  TCoef sk, ck;
  for(int k = 0; k <= order; k++)
  {
    sk  = sin(k*omega);  // \todo optimize by trigonometric recursion
    ck  = cos(k*omega);  
    ca += ck * a[k];
    sa += sk * a[k];
    cb += ck * b[k];
    sb += sk * b[k];
  }
  return sqrt( (cb*cb + sb*sb) / (ca*ca + sa*sa) );
}

template<class TSig, class TCoef>
void rsDirectFormFilter<TSig, TCoef>::getMagnitudeResponse(TCoef* frequencies, TCoef* magnitudes, 
  int numBins, TCoef sampleRate, bool inDecibels, bool accumulate)
{
  TCoef H;
  for(int k = 0; k < numBins; k++)
  {
    H = getMagnitudeResponseAt(TCoef(2*PI)*frequencies[k] / sampleRate);

    // \todo move this decision function somewhere where it can be shared - other filter use that, 
    // too
    if( !inDecibels && !accumulate )
      magnitudes[k]  = H;
    else if( !inDecibels && accumulate )
      magnitudes[k] *= H;
    else if( inDecibels && !accumulate )
      magnitudes[k]  = rsAmpToDbWithCheck(H, TCoef(0.0000001)); // clip at -140 dB
    else if( inDecibels && accumulate )
      magnitudes[k] += rsAmpToDbWithCheck(H, TCoef(0.0000001));
  }
}

// others:

template<class TSig, class TCoef>
void rsDirectFormFilter<TSig, TCoef>::reset()
{
  for(int i = 0; i <= maxOrder; i++)
    w[i] = 0.0;
}

template<class TSig, class TCoef>
void rsDirectFormFilter<TSig, TCoef>::initializeCoefficients()
{
  b[0] = 1.0;
  a[0] = 1.0;
  for(int i = 1; i <= maxOrder; i++)
  {
    b[i] = 0.0;
    a[i] = 0.0;
  }
}


/*

Ideas:
-Allow to use it in TDF-1 form, too. I think, it can be done by re-interpreting the state-buffer. 
 We just need another function getSampleTDF1 (and maybe rename getSample to getSampleDF2, but keep
 getSample as alias).
-Maybe make a subclass that allows for non-delay-canonical computation, i.e. DF-1 and TDF-2. This 
 needs a 2nd state buffer, that's why it should go into a subclass.
-Try to replace the unit delays with multi-sample delays, i.e. z^-1 is replaced by z^-M for some 
 natural number M. I think, we'll get spectral replications of the frequency responses, i.e. a 
 lowpass turns into a sort of comb-filter? Maybe that can be useful to create higher order 
 comb-filters with flatter passbands? If we get useful results, maybe try to do the same thing
 with more modulatable filter structures such as biquad cascades where biquads may be implemented
 as ZDF-SVF.

*/