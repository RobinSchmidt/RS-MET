template<class T>
void rsQuantumSpinFunctions<T>::randomizeState(Vec& v, PRNG* prng)
{
  T Pu = prng->getSample();             // probability of "up"
  T Pd = T(1) - Pu;                     // probability of "down"
  T r  = T(1) / sqrt(Pu*Pu + Pd*Pd);    // normalizer
  Pu  *= r;
  Pd  *= r;
  T pu = T(2.0*PI) * prng->getSample(); // phase of v.y
  T pd = T(2.0*PI) * prng->getSample(); // phase of v.x

  v.y = std::polar(Pu, pu);
  v.x = std::polar(Pd, pd);
}

template<class T>
T rsQuantumSpinFunctions<T>::getExpectedMeasurement(const Mat& M, const Vec& A)
{
  // computation via (1) Eq 3.26 - this is just naively applying the general formula for 
  // an expectation value:
  rsQuantumSpin<T> E1 = M.eigenvector1();
  rsQuantumSpin<T> E2 = M.eigenvector2();
  std::complex<T>  e1 = M.eigenvalue1();
  std::complex<T>  e2 = M.eigenvalue2();
  T P1 = rsQuantumSpin<T>::getStateProbability(A, E1);
  T P2 = rsQuantumSpin<T>::getStateProbability(A, E2);
  T E  = (e1*P1 + e2*P2).real();
  return E;


  // ...but the same value can be computed more efficiently by (1) Eq 4.14:
  //T E = (A * (M * A)).real();
  // use a function E = bracket(A, M, A) that takes two kets and a matrix, converts the 1st into 
  // a bra and computes the inner product

  return E;
}







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
  T pu = T(2.0*PI) * prng->getSample(); // phase of v.y
  T pd = T(2.0*PI) * prng->getSample(); // phase of v.x

  y = std::polar(Pu, pu);
  x = std::polar(Pd, pd);

  //normalize(); // should already be normalized thanks to *= r
}


// measurements

template<class T>
T rsQuantumSpin<T>::measureObservable(const rsSpinOperator<T>& M, rsNoiseGenerator<T>* prng)
{
  rsQuantumSpin<T> E1 = M.eigenvector1();
  T P1 = getStateProbability(*this, E1);
  T rnd = prng->getSample();
  //if(rnd >= P1) { // new
  if(rnd <= P1) { // should it be <= or < ?  old
    //T P2 = getStateProbability(*this, M.getEigenvector2()); // should be 1-P1
    setState(E1);
    return M.eigenvalue1().real(); // is real if M is Hermitian
  }
  else {
    rsQuantumSpin<T> E2 = M.eigenvector2();
    setState(E2);
    return M.eigenvalue2().real(); // is real if M is Hermitian
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
T rsQuantumSpin<T>::measureSpin(T nx, T ny, T nz, rsNoiseGenerator<T>* prng)
{
  rsSpinOperator<T> M;
  M.setValues(nz, nx-i*ny, nx+i*ny, -nz); // (1) Eq 3.23
  return measureObservable(M, prng);
}

template<class T>
T rsQuantumSpin<T>::measureSpinZ(rsNoiseGenerator<T>* prng)
{
  //T Pd  = getStateProbability(*this, down());       // (1) Eq 2.2
  //T Pd  = getDownProbability(*this);
  T Pd  = getSquaredNorm(x);
  T rnd = prng->getSample();
  if(rnd <= Pd) {       // should it be <= or < ?
    prepareDownState();
    return -1; }
  else {
    prepareUpState();
    return +1; }
}

template<class T>
T rsQuantumSpin<T>::measureSpinX(rsNoiseGenerator<T>* prng)
{
  T Pl  = getLeftProbability(*this); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Pl) {
    prepareLeftState();
    return -1; }
  else {
    prepareRightState();
    return +1; }
}

template<class T>
T rsQuantumSpin<T>::measureSpinY(rsNoiseGenerator<T>* prng)
{
  T Po = getOutProbability(*this); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Po) {
    prepareOutState();
    return -1; }
  else {
    prepareInState();
    return +1; }
}

//=================================================================================================

template<class T>
T rsSpinOperator<T>::getExpectedMeasurement(
  const rsSpinOperator<T>& M, const rsQuantumSpin<T>& A)
{
  // computation via (1) Eq 3.26 - this is just naively applying the general formula for 
  // an expectation value:
  rsQuantumSpin<T> E1 = M.eigenvector1();
  rsQuantumSpin<T> E2 = M.eigenvector2();
  std::complex<T>  e1 = M.eigenvalue1();
  std::complex<T>  e2 = M.eigenvalue2();
  T P1 = rsQuantumSpin<T>::getStateProbability(A, E1);
  T P2 = rsQuantumSpin<T>::getStateProbability(A, E2);
  T E  = (e1*P1 + e2*P2).real();
  return E;
  

  // ...but the same value can be computed more efficiently by (1) Eq 4.14:
  //T E = (A * (M * A)).real();
  // use a function E = bracket(A, M, A) that takes two kets and a matrix, converts the 1st into 
  // a bra and computes the inner product

  return E;
}


// todo: make a function that computes the two angles that describe the state as point on
// the Bloch sphere https://en.wikipedia.org/wiki/Bloch_sphere





//=================================================================================================
/*                            Background Information



Notes:
I think, the relationship to what is called the "wavefunction" in quantum mechanics is as 
follows: The wavefunction is in general some function from a set S into the complex numbers. 
For position and momentum of a particle, that set S is the set of real numbers R (for 1D) or 
R^2 or R^3 for 2D and 3D space. In our case, the set S is just S = { up, down } - which can be
renamed to S = { 0, 1 } or S = { |0>, |1> } to get the common qubit noation. The "collapse of 
the wavefunction" occurs when we assign the pure up/down states in the measurement operation. In
a general state, both of these pure states have a complex number associated with them - the 
"probability amplitude" - and the square of its magnitude gives the actual probability. When the
wavefunctions is collapsed due to a measurement one the values becomes 1 and the other 0 - it 
becomes a delta distribution - although its spikey nature is not really obvious in this simple 
case where we have only two possible input values into the function.


*/