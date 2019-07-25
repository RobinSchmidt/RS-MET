
template<class T>
rsAttackDecayFilter<T>::rsAttackDecayFilter()
{
  reset();
  updateCoeffs();
}

template<class T>
void rsAttackDecayFilter<T>::updateCoeffs()
{
  // preliminary - todo: 
  // -compute the two required time-constants (formulas are somewhere in the modal filter)
  // -use IIT to compute the coeffs
  //ca = 0.99; 
  //cd = 0.999;


  T tauAttack;
  expDiffScalerAndTau2(decaySamples, attackSamples, &tauAttack, &s);

  // verify these:
  ca = exp(-1.0/tauAttack);              // = exp(-alpha), pole radius
  cd = exp(-1.0/decaySamples);




  coeffsDirty = false;
}