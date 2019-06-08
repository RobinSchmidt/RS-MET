
template<class T>
void rsNonUniformOnePole<T>::setOmega(T newOmega)
{
  w = newOmega;
  T dummy;       // will be set to 0
  //RAPT::rsFirstOrderFilterBase<T, T>::coeffsLowpassIIT(w, &a, &b, &dummy);

  RAPT::rsFirstOrderFilterBase<T, T>::coeffsLowpassIIT(w, &a, &dummy, &b);


  //coeffsLowpassIIT(T w, T* b0, T* b1, T* a1)

  // update s? or maybe in reset?

  dummy = 0;
}

template<class T>
T rsNonUniformOnePole<T>::getSample(T x, T dt)
{
  T bdt = pow(b, dt);  // b^dt
  // can later be optimized as bdt = exp(log(b) * dt) where log(b) can be precomputed

  
  // update scaler:
  s = a + bdt*s;       //  Eq. 16, with w = 0

  
  // update state:
  //y = (a*x + bdt*y) / rsAbs(s);   //  Eq. 13
  //return y;

  // verify, if the scaler really has to be applied when updating the state or if it should be 
  // applied only to the returned output

  //// alternative version:
  y = a*x + bdt*y;
  return y / rsAbs(s);
  // with dt == 1, there's no difference - s is always 1 - but probably it's correct to apply it
  // to the output because the paper says, with several parallel filtersm the normalization should
  // *not* be done in the individual complex one-poles but once for the whole filter - which makes 
  // sense only if it's applied to the output
}

template<class T>
T rsNonUniformOnePole<T>::getSample2(T x, T dt)
{
  T bdt = pow(b, dt); // b^dt - express exp(log(b) * dt), precompute log(b)

  // some coeffs r0,r1 that can be precomputed:
  T bm1 = b-1;
  T r0  = bm1*bm1 / (a*b);
  T r1  = a / bm1;

  // compute additional compensation term:
  T R   = (bdt-1)/(r0*dt);
  T Phi = (R-r1*b)*x - (R-r1*bdt)*x1;

  y = a*x + bdt*y + Phi;
  return y;
}
// needs testing and debugging


template<class T>
void rsNonUniformOnePole<T>::reset()
{
  x1 = T(0);
  y  = T(0);
  //s = a;
  s = a / (T(1) - b);
  // the general formula is s = a / (1 - b * exp(j*wr)) where wr is the reference frequency where 
  // we  want unit gain, which is zero in this (lowpass) case
  // there's an alternative formula: s = a ...figure out, which one should be used
  // ...comparing results with a uniform 1-pole filter, it seems like s = a/(1-b) is the correct 
  // one - in this case, s comes out as 1 - but probably only in the special case of a lowpass
}

template class rsNonUniformOnePole<double>;




/*



*/