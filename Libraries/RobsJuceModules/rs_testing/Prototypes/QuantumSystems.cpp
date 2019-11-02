
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

template<class T>
rsMatrix2x2<std::complex<T>> rsQuantumSpin<T>::densityMatrix(
  std::vector<T> p, std::vector<rsVector2D<std::complex<T>>> s)
{
  rsAssert(p.size() == s.size(), "sizes must match");
  Mat r(T(0), T(0), T(0), T(0));
  for(size_t i = 0; i < p.size(); i++)
    r = r + Complex(p[i]) * projector(s[i]);
  return r;
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
  return sandwich(A, M, A).real();
}

template<class T>
T rsQuantumSpin<T>::getUncertaintySquared(const Mat& M, const Vec& A)
{
  T Me = getExpectedMeasurement(M, A);         // expectation value of observable M in state A
  Mat Mc = M - Complex(Me) * Mat::identity();  // centered version of M, has expectation zero
  return sandwich(A, Mc*Mc, A).real();         // (1) pg 141
}

template<class T>
T rsQuantumSpin<T>::getUncertaintyProduct(const Mat& M, const Mat& L, const Vec& A)
{
  T   Me, Le;     // expectation values of measurables associated with M,L
  Mat Mc, Lc, C;  // centered operators and their commutator
  Me = getExpectedMeasurement(M, A); Mc = M - Complex(Me) * Mat::identity();
  Le = getExpectedMeasurement(L, A); Lc = L - Complex(Le) * Mat::identity();
  C  = Mat::commutator(Mc, Lc);
  return T(0.5) * abs(sandwich(A, C, A));  // (1) Eq 5.13 (with centered matrices)
}

// measurement:

template<class T>
T rsQuantumSpin<T>::measureObservable(Vec& A, const Mat& M, rsNoiseGenerator<T>* prng)
{
  Vec E1 = M.getEigenvector1();
  T P1 = getStateProbability(A, E1);
  T rnd = prng->getSample();
  if(rnd <= P1) { // should it be <= or < ? 
    //T P2 = getStateProbability(A, M.getEigenvector2()); // should be 1-P1
    A = E1;                        // collapse the state into eigenstate 1
    return M.getEigenvalue1().real(); // return eigenvalue 1, (it's real if M is Hermitian)
  }
  else {
    //Vec E2 = M.eigenvector2();
    A = M.getEigenvector2();          // collapse into eigenstate 2
    return M.getEigenvalue2().real(); // return eigenvalue 2
  }
}
// in general, we'll have an NxN matrix and the probability to be in state k is given by
//  <v|Ek> * <v|Ek> where v is our N dimensional complex state vector and Ek is the k-th
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
//=================================================================================================

template<class T>
void rsQuantumParticle<T>::initializeWaveFunction(std::vector<Complex>& Psi_0)
{
  Psi = Psi_0;
  Psi_t.resize(Psi.size());
  Psi_xx.resize(Psi.size());
}

template<class T>
void rsQuantumParticle<T>::updateWaveFunction(T dt)
{
  int n;
  int Nx = (int)Psi.size();
  T dx   = (rsLast(x) - x[0]) / (x.size()-1);
  Complex i(0, 1); // imaginary unit

  // compute second spatial derivative of wave function Psi by central differences 
  // (treating the ends cyclically):
  for(n = 0; n < Nx; n++)
    Psi_xx[n] = (Psi[wrap(n-1,Nx)] + Psi[wrap(n+1,Nx)] - 2.*Psi[n])/(dx*dx);

  // compute time derivative of wave function via the Schroedinger equation:
  // Psi_t = (i*hBar)/(2*m)*Psi_xx - (i/hBar)*V*Psi:

  //for(n = 0; n < Nx; n++) 
  //  Psi_t[n] = ((i*hBar)/(2*m))    * Psi_xx[n] // term for free particle
  //            -((i/hBar)* V(x[n])) * Psi[n];   // term from the potential

  //for(n = 0; n < Nx; n++) 
  //  Psi_t[n] = ((i*hBar)/(2*m))    * Psi_xx[n] // term for free particle
  //            +((i/hBar)* V(x[n])) * Psi[n];   // term from the potential

  //for(n = 0; n < Nx; n++) 
  //  Psi_t[n] = ((-hBar/2.)*Psi_xx[n]  +  (1./(hBar))*V(x[n])*Psi[n])/i; // 10.13

  // 10.13:
  for(n = 0; n < Nx; n++)
  {
    Psi_t[n]  = -0.5 * Psi_xx[n];
    Psi_t[n] +=  V(x[n]) * Psi[n];
    Psi_t[n] /= i;
  }

  // update wave function using a forward Euler step:
  for(n = 0; n < Nx; n++) 
    Psi[n]  = Psi[n] + dt * Psi_t[n];
}
// how else (other than cyclically) could we treat the ends?
// optimize: get rid of the call to "wrap" by treating the ends separately outside the loop 
// susskind, pg 119
// somehow, the wvaefunctions shown here
// https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
// look different - is there something wrong with my potential function? and/or the implementation 
// of the schroedinger equation? something seems to be wrong - triple-check everything
// ...maybe get rid of hBar - set it to 1, also get rid of the mass here
// https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#One-dimensional_examples

//=================================================================================================


template<class T>
void rsQuantumComputer<T>::applyGate(const QGate& g, int t)
{
  std::vector<Vec> tmp = qbits;  
  // silly - make this tmp a member - hmm - but the code int the paper doesn't use a temporary
  // array - it only uses local temporaries during the update of the two amplitudes - i think, it 
  // makes no difference - each index appears exactly once in the loop, so there should be no
  // issues of overwriting and then reading a state again

  //int iMax = pow(2, numQBits-1);
  int iMax = rsPowInt(2, numQBits-1);

  //for(int i = 0; i < numStates; i++)
  for(int i = 0; i < iMax; i++) 
    // the paper says it goes to 2^(n-1) but then the last gives an access violation (n=4,t=3)
  {
    int a = nth_cleared(i, t);
    int b = a | (1 << t);
    // with numQBits = 4, numStates = 16, b gets larger than 15 - but that would mean an access
    // violation - something must be wrong here (the code is taken from (1)) - they call it with
    // global_id - yeas - that seems to correspond to teh loop index:
    // https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/get_global_id.html
    // could it have to do with the number representation? maybe we should use unsigned ints?
    // ...nope - that doesn't seem to make any difference
    // int in OpenCL is a 32-bit signed int in two's complement - that should be fine...
    // https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/scalarDataTypes.html
    // todo: maybe figure out a function for the n-th digit clear myself (maybe an ineeficient 
    // one)
    // ...maybe we must additionally mask by numStates? but the paper says no such thing - might 
    // this be done implieitly by OpenCL? a sort of auto-wrap around? or do the bitwise 
    // operators work differently?
    // https://www.khronos.org/registry/OpenCL//sdk/2.2/docs/man/html/bitwiseOperators.html
    // maybe check out the actual source code - maybe the code in the paper is buggy?

    // oh no - the text says: "The structure of the algorithm is a for loop through have the 
    // number of amplitudes." .....and that "have" probably means "half"?
    // ...but even that gives an access violation on the last index (if we use 
    // i <= pow(2, numQBits-1)) - perhaps it should be < instead of <=
    // ..ok - no access violation anymore - but with so much guesswork necessarry, who knowsm 
    // what else is wrong (i think, i probably should not use a tempoprary array)
    // ...how can we verify, if the algo does the right thing?


    qbits[a] = g.a * tmp[a] + g.b * tmp[b];
    qbits[b] = g.c * tmp[a] + g.d * tmp[b];
    int dummy = 0;
  }

  // see also https://quantum-journal.org/papers/q-2018-01-31-49/pdf/?
}