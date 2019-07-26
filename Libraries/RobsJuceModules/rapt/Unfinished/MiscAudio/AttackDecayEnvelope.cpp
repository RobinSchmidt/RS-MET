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