
// Operator facory functions


template<class T>
rsVector3D<rsMatrix2x2<std::complex<T>>> rsQuantumSpin<T>::pauliVector() 
{
  rsVector3D<rsMatrix2x2<Complex>> sigma;
  setToPauliX(sigma.x);
  setToPauliY(sigma.y);
  setToPauliZ(sigma.z);
  return sigma;
}

// state setup:

template<class T>
void rsQuantumSpin<T>::normalizeState(Vec& A)
{
  T r = sqrt(T(1) / getTotalProbability(A));  
  // or 1/sqrt(t) instead of sqrt(1/t)? - which one is better numerically?
  A.x *= r;
  A.y *= r;
}

template<class T>
void rsQuantumSpin<T>::randomizeState(Vec& v, PRNG* prng)
{
  T Pu = prng->getSample();             // probability of "up"
  T Pd = T(1) - Pu;                     // probability of "down"
  T r  = T(1) / sqrt(Pu*Pu + Pd*Pd);    // normalizer
  Pu  *= r;
  Pd  *= r;
  T pu = T(2.0*PI) * prng->getSample(); // phase of v.x
  T pd = T(2.0*PI) * prng->getSample(); // phase of v.y

  v.x = std::polar(Pu, pu);
  v.y = std::polar(Pd, pd);
  // exchange x,y - have amount parameter
}

// inquiry:

template<class T>
T rsQuantumSpin<T>::getStateProbability(const Vec& A, const Vec& t)
{
  std::complex<T> r = bracket(A,t) * bracket(t,A); // (1), Eq 3.11 (with lambda_i replaced by t)
  return r.real();                                 // imag should be zero
}

template<class T>
T rsQuantumSpin<T>::getExpectedMeasurement(const Mat& M, const Vec& A)
{
  //// computation via (1) Eq 3.26 - this is just naively applying the general formula for 
  //// an expectation value:
  //Vec E1 = M.eigenvector1();
  //Vec E2 = M.eigenvector2();
  //std::complex<T> e1 = M.eigenvalue1();
  //std::complex<T> e2 = M.eigenvalue2();
  //T P1 = getStateProbability(A, E1);
  //T P2 = getStateProbability(A, E2);
  //T E  = (e1*P1 + e2*P2).real();
  //return E;

  // ...but the same value can be computed more efficiently by (1) Eq 4.14:
  return bracket(A, M*A).real();
}

// measurement:

template<class T>
T rsQuantumSpin<T>::measureObservable(Vec& A, const Mat& M, rsNoiseGenerator<T>* prng)
{
  Vec E1 = M.eigenvector1();
  T P1 = getStateProbability(A, E1);
  T rnd = prng->getSample();
  if(rnd <= P1) { // should it be <= or < ? 
    //T P2 = getStateProbability(A, M.getEigenvector2()); // should be 1-P1
    A = E1;
    return M.eigenvalue1().real(); // is real if M is Hermitian
  }
  else {
    Vec E2 = M.eigenvector2();
    A = E2;
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
T rsQuantumSpin<T>::measureSpin(Vec& A, T nx, T ny, T nz, rsNoiseGenerator<T>* prng)
{
  Mat M;
  M.setValues(nz, nx-i*ny, nx+i*ny, -nz); // (1) Eq 3.23
  return measureObservable(A, M, prng);
}

template<class T>
T rsQuantumSpin<T>::measureSpinZ(Vec& A, rsNoiseGenerator<T>* prng)
{
  //T Pd  = getStateProbability(A, down());       // (1) Eq 2.2
  //T Pd = getDownProbability(A);
  T Pd  = rsAbsSquared(A.y);   // more efficient
  T rnd = prng->getSample();
  if(rnd <= Pd) {       // should it be <= or < ?
    prepareDownState(A);
    return -1; }
  else {
    prepareUpState(A);
    return +1; }
}

template<class T>
T rsQuantumSpin<T>::measureSpinX(Vec& A, rsNoiseGenerator<T>* prng)
{
  T Pl  = getLeftProbability(A); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Pl) {
    prepareLeftState(A);
    return -1; }
  else {
    prepareRightState(A);
    return +1; }
}

template<class T>
T rsQuantumSpin<T>::measureSpinY(Vec& A, rsNoiseGenerator<T>* prng)
{
  T Po = getOutProbability(A); // optimizable?
  T rnd = prng->getSample();
  if(rnd <= Po) {
    prepareOutState(A);
    return -1; }
  else {
    prepareInState(A);
    return +1; }
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

misc:
what is used as the quantum "bit" is actually one of spin's 
components, such as the z-component (a pure "up" represents binary 1 and a pure "down" state
represents binary 0. Mixed states (with respect to the z-axis) represent a superposition of 0 and
1 (a pure left or right or in or out state is "mixed" with respect to the z axis)
done ..maybe make a simpler class rsQuantumBit that deals exclusively with the up/down direction
and doesn't consider the others at all - it should also not call the states "up" and "down" but
|1> and |0> respectively - it's more abstract and needs less "physical" features

*/