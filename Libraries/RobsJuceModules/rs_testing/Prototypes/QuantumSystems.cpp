

//=================================================================================================

// setup

template<class T>
void rsQuantumSpin<T>::randomizeState(rsNoiseGenerator<T>* prng) 
{
  T Pu = prng->getSample();             // probability of "up"
  T Pd = T(1) - Pu;                     // probability of "down"
  T r  = T(1) / sqrt(Pu*Pu + Pd*Pd);    // normalizer
  Pu  *= r;
  Pd  *= r;
  T pu = T(2.0*PI) * prng->getSample(); // phase of au
  T pd = T(2.0*PI) * prng->getSample(); // phase of ad

  au = std::polar(Pu, pu);
  ad = std::polar(Pd, pd);

  //normalize(); // should already be normalized thanks to *= r
}


// measurements

template<class T>
T rsQuantumSpin<T>::measureObservable(const rsSpinOperator<T>& M, rsNoiseGenerator<T>* prng)
{
  rsQuantumSpin<T> E1 = M.getEigenvector1();
  T P1 = getStateProbability(*this, E1);
  T rnd = prng->getSample();
  if(rnd <= P1) {
    //T P2 = getStateProbability(*this, M.getEigenvector2()); // should be 1-P1
    setState(E1);
    return M.getEigenvalue1().real(); // is real if M is Hermitian
  }
  else {
    rsQuantumSpin<T> E2 = M.getEigenvector2();
    setState(E2);
    return M.getEigenvalue2().real(); // is real if M is Hermitian
  }
}
// in general, we'll have an NxN matrix and the probability to be in state k is given by
//  (v * Ek) * (v * Ek) where v is our N dimensional complex state vector and Ek is the k-th
// eigenvector of M. To collapse into one of the N states, we'll have to look at into which 
// interval the random variable falls. For example, if P1 = 0.2, P2 = 0.5, P3 = 0.3 for a
// 3D state, we'll fall into E1 if rnd in 0..0.2, into E2 if rnd in 0.2...0.7 and into E3 if
// rnd in 0.7..1.0 (and return eigenvalue e1, e1 or e3 respectively) - it's probably a good idea
// to pass precomputed eigenvalues and -vectors into such a N-dimensional measurement function


template<class T>
T rsQuantumSpin<T>::measureSpinZ(rsNoiseGenerator<T>* prng)
{
  //T Pu  = getUpProbability(*this); // optimize this!
  T Pu = getStateProbability(*this, up());       // (1) Eq 2.2
  //T Pu = getSquaredNorm(au); // same result as Pu = getUpProbability(*this) but more efficient
  T rnd = prng->getSample();
  if(rnd <= Pu) {       // should it be <= or < ?
    prepareUpState();
    return +1; }
  else {
    prepareDownState();
    return -1; }
}

template<class T>
T rsQuantumSpin<T>::measureSpinX(rsNoiseGenerator<T>* prng)
{
  T Pr = getRightProbability(*this); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Pr) 
  {
    prepareRightState();
    return +1; 
  }
  else 
  {
    prepareLeftState();
    return -1; 
  }
}

template<class T>
T rsQuantumSpin<T>::measureSpinY(rsNoiseGenerator<T>* prng)
{
  T Pi = getInProbability(*this); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Pi) {
    prepareOutState();
    return -1; }
  else {
    prepareInState();
    return +1; }
}


/*
template<class T>
T rsQuantumSpin<T>::measureRightComponent(rsNoiseGenerator<T>* prng)
{
  T Pr = getRightProbability(*this);
  T rnd = prng->getSample();
  if(rnd <= Pr) {
    prepareRightState();
    return +1; }
  else {
    prepareLeftState();
    return -1; }
}

template<class T>
T rsQuantumSpin<T>::measureInComponent(rsNoiseGenerator<T>* prng)
{
  T Pi = getInProbability(*this);
  T rnd = prng->getSample();
  if(rnd <= Pi) {
    prepareInState();
    return +1; }
  else {
    prepareOutState();
    return -1; }
}
*/



// todo: make a function that computes the two angles that describe the state as point on
// the Bloch sphere https://en.wikipedia.org/wiki/Bloch_sphere