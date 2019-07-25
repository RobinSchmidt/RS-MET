

template<class T>
void rsAttackDecayFilter<T>::updateCoeffs()
{
  // preliminary - todo: 
  // -compute the two required time-constants (formulas are somewhere in the modal filter)
  // -use IIT to compute the coeffs
  ca = 0.99; 
  cd = 0.999;

  coeffsDirty = false;
}