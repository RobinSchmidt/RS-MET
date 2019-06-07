
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

  
  // update scaler:
  s = a + bdt*s;       //  Eq. 16, with w = 0

  
  // update state:
  y = (a*x + bdt*y) / abs(s);   //  Eq. 13
  // verify, if the scaler really has to be applied when updating the state or if it should be 
  // applied only to the returned output

  return y;
}

template<class T>
void rsNonUniformOnePole<T>::reset()
{
  y = T(0);
  s = a; 
  //s = a / (T(1) - b); 
  // the general formula is s = a / (1 - b * exp(j*wr)) where wr is the reference frequency where 
  // we  want unit gain, which is zero in this (lowpass) case
  // there's an alternative formula: s = a ...figure out, which one should be used
}

template class rsNonUniformOnePole<double>;




/*



*/