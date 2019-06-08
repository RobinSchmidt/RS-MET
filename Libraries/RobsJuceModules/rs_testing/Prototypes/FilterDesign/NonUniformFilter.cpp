
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
  T bdt = pow(b, dt);
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
// this function implements the "spatial invariant scaling" normalization method (section 3.2.4)
// todo: implement the "piecewise resampling method" from section 3.2.3, too - maybe we need to do
// something else in reset() then?

template<class T>
void rsNonUniformOnePole<T>::reset()
{
  y = T(0);
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