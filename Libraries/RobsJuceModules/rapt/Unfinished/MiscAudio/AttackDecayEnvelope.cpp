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

  T attackSamples2 = rsMin(T(0.99) * decaySamples, attackSamples + T(1));
  // Why the 2nd +1? to avoid numerical problems when is goes down to zero? Then maybe 
  // using max would be better...or is there some offset of 1 sample that is being compesated?


  expDiffScalerAndTau2(decaySamples, attackSamples2, &tauAttack, &s);
  ca = exp(-1.0/tauAttack);       // = exp(-alpha), pole radius
  cd = exp(-1.0/decaySamples);
  coeffsDirty = false;
}


template<class T>
T rsAttackDecayEnvelope<T>::getExactAccuCompensation()
{
  if( this->yd == 0 )
    return 1;


  // some precomputations:
  T lcd = log(this->cd);
  T lca = log(this->ca);
  T R   = this->cd / this->ca;
  T rlR = T(1)/log(R);

  // objective function of which we want to find a root:
  auto f = [=](T x)->T
  {
    //return getPeakForInputImpulse(x) - T(1); // naive, not optimized

    T a0 = x + this->ya * this->ca;
    T d0 = x + this->yd * this->cd;
    T D  = d0*lcd;
    T A  = a0*lca;

    // we need to catch the case D==0
    //if( D == 0 )
    //  return 1;


    T np = log(A/D) * rlR;    // = rsLogB(A/D, R), peak location in samples
    T dp = d0 * exp(lcd*np);  // = d0*pow(cd, np), decay output at peak
    T ap = a0 * exp(lca*np);  // = a0*pow(ca, np), attack output at peak
    T ep = this->s * (dp-ap); // env output at peak
    return ep - T(1);         // subtract target value of 1
  };

  return rsRootFinder<T>::bisection(f, T(0), T(1), T(0));
  //return rsRootFinder<T>::falsePosition(f, T(0), T(1), T(0));
  // false position seems to work quite well but maybe try Newton, Brent, Ridders, ... and
  // maybe use a higher tolerance, try better initial interval - maybe (0,1-yd) - then test
  // how many iterations are typically taken
  // oh - damn - it is not guaranteed to converge - this is not yet usable
}


template<class T>
T rsAttackDecayEnvelope<T>::getAccuCompensatedImpulse()
{
  switch(accuFormula)
  {
  case AccuFormula::one_minus_yd: return T(1) - this->yd;
  case AccuFormula::exact:        return getExactAccuCompensation();
  default: return T(1);
  }
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

ToDo:
-maybe to implement a release time different from the decay time, just switch the decay time 
 constant to a different value in release phase or let the constant input not go down to zero 
 immediately but via another RC filter...but that would work well only if release >= decay

*/
