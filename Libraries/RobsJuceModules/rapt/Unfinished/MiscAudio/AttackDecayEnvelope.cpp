template<class T>
rsAttackDecayFilter<T>::rsAttackDecayFilter()
{
  reset();
  updateCoeffs();
}

template<class T>
void rsAttackDecayFilter<T>::updateCoeffs()
{
  T tauAttack;
  expDiffScalerAndTau2(decaySamples, attackSamples+T(1), &tauAttack, &s);
  ca = exp(-1.0/tauAttack);       // = exp(-alpha), pole radius
  cd = exp(-1.0/decaySamples);
  coeffsDirty = false;
}

/*
Notes:

The formula for rsAttackDecayFilter::getGainAtDC was derived as follows:

  gd = 1 / (1 - cd)     DC gain of decay filter
  ga = 1 / (1 - ca)     DC gain of attack filter
  g  = s * (gd - ga)    DC gain of the whole filter

and then algebraically simplified.

ToDo:
When setting up a nonzero sustain, the location and height of the peak shifts - with increasing 
sustain, the peak comes later and gets higher. Try to figure out, how we have to adjust the filter 
coeffs ca,cd,s to compensate for this effect. I think, we need to consider the output signal of the 
filter, when the input is a weighted sum of a unit-impulse and unit step, which is for an exp-decay
filter just the same weighted sum of exp(-t/tau) and 1-exp(-t/tau). We need to set up the equation 
and churn through the same analysis and algebra (take derivative, set zero, etc.) as was originally
done to derive the formulas in expDiffScalerAndTau2.


*/