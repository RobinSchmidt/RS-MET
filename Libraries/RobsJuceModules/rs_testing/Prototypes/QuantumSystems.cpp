

//=================================================================================================

// setup

template<class T>
void rsQuantumSpin<T>::randomizeState() 
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
T rsQuantumSpin<T>::measureObservable(const rsSpinOperator<T>& M)
{
  rsQuantumSpin<T> E1 = M.getEigenvector1();
  T P1 = getStateProbability(*this, E1);       // (1) Eq 2.2
  T rnd = prng->getSample();
  if(rnd <= P1) {
    rsQuantumSpin<T> E2 = M.getEigenvector2();
    //T P2 = getStateProbability(*this, E2); // should be 1-P1 - move to tests
    setState(E2);
    return M.getEigenvalue2().real(); // should be real anyway, if M is Hermitian
  }
  else {
    setState(E1);
    return M.getEigenvalue1().real(); // should be real anyway, if M is Hermitian
  }
  // my intuition says, the branches should be theo other way around but that doesn't fit with 
  // measureUpComponent - why is this so?
}

template<class T>
T rsQuantumSpin<T>::measureUpComponent()
{
  //T Pu  = getUpProbability(*this); // optimize this!
  //T Pu = getStateProbability(*this, up());       // (1) Eq 2.2
  T Pu = getSquaredNorm(au); // same result as Pu = getUpProbability(*this) but more efficient
  T rnd = prng->getSample();
  if(rnd <= Pu) {       // should it be <= or < ?
    prepareUpState();
    return +1; }
  else {
    prepareDownState();
    return -1; }
}

template<class T>
T rsQuantumSpin<T>::measureRightComponent()
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
T rsQuantumSpin<T>::measureInComponent()
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



// todo: make a function that computes the two angles that describe the state as point on
// the Bloch sphere https://en.wikipedia.org/wiki/Bloch_sphere