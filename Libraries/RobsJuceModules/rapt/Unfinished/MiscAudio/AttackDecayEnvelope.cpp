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

The formula rsAttackDecayFilter::getGainAtDC was derived as follows:

  gd = 1 / (1 - cd)     DC gain of decay filter
  ga = 1 / (1 - ca)     DC gain of attack filter
  g  = s * (gd - ga)    DC gain of the whole filter

and then algebraically simplified.


*/