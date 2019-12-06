template<class T>
void rsHeatEquation1D<T>::setMaxCycleLength(int newLength) 
{ 
  rodArray1.resize(newLength);
  rodArray2.resize(newLength);
  rsArray::fillWithZeros(&rodArray1[0], newLength);
  rsArray::fillWithZeros(&rodArray2[0], newLength);
  reset(); // possibly re-adjust pointers
}

template<class T>
void rsHeatEquation1D<T>::setHeatDistribution(T* d, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  for(int i = 0; i < N; i++)
    rodIn[i] = rodOut[i] = d[i];
}

template<class T>
void rsHeatEquation1D<T>::setRandomHeatDistribution(int seed, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  // todo: allow different probability densities for the prng



  rsNoiseGenerator<T> prng;
  prng.setSeed(seed);
  for(int i = 0; i < N; i++)
    rodIn[i] = prng.getSample();
}


template<class T>
void rsHeatEquation1D<T>::setTwoValueDistribution(T highFraction, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  int nh = (int) round(highFraction*N);
  for(int i = 0; i < nh; i++) rodIn[i] = +1.0;
  for(int i = nh; i < N; i++) rodIn[i] = -1.0;
}


template<class T>
void rsHeatEquation1D<T>::normalizeHeatDistribution(T targetMean, T targetVariance)
{
  // set mean to desired target value (maybe factor out):
  //int N = (int) rodArray1.size();
  int N = (int) rodLength;
  typedef rsArray AR;
  T mean = AR::mean(rodIn, N);
  AR::add(rodIn, -mean, rodIn, N); 
  // ...hmm..this actually just set the mean to zero...which is the most reasonable target value
  // nayway (in case of sound generation)...but what if the user wants some other mean..we'll see

  // set variance to desired target variance:
  // ....
}

//=================================================================================================

template<class T>
T rsWaveEquation1D<T>::getCourantNumber(T timeStep) const
{
  int N = getNumGridPoints()-1; // see (1), section 5.2.8
  T c   = waveSpeed;            // what unit?
  T L   = T(1);                 // we use the unit interval [0,1], so the spatial length is 1
  T gam = c / L;                // "gamma" - (1), Eq. 6.5
  T k   = timeStep;             // temporal sampling interval
  T h   = L / N;                // spatial sampling interval
  return gam*k / h;             // "lambda", the Courant number
}
// see:
// https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

template<class T>
T rsWaveEquation1D<T>::getOmegaForWaveNumber(T waveNumber, T timeStep) const
{
  T beta   = waveNumber;
  T k      = timeStep;
  T h      = getGridSpacing();
  T lambda = getCourantNumber(timeStep);
  return (T(2)/k) * asin(lambda * sin(beta*h/2)); // (1), Eq. 6.43
  // does this depend on the scheme? ...probably - we may need different formulas for different 
  // schemes
}

template<class T>
void rsWaveEquation1D<T>::setInitialConditions(T* newPositions, T* newVelocities,
  int length, T timeStep)
{
  rsAssert(length == getNumGridPoints(), "array length should match number of grid points");
  T k = timeStep;
  for(int l = 0; l < length; l++) {
    u[l]  = newPositions[l];
    u1[l] = u[l] + k*newVelocities[l];  // (1), Eq 6.36
  }
}

template<class T>
void rsWaveEquation1D<T>::updateState(T timeStep)
{
  computeInteriorPoints(timeStep);
  //computeInteriorPointsSimple();      // simplified formula for Courant number == 1
  computeBoundaryPoints(timeStep);
  updateStateArrays();
}

template<class T>
void rsWaveEquation1D<T>::computeInteriorPoints(T timeStep)
{
  // We implement the scheme in (1), Eq. 6.34: u_tt = g^2 * u_xx where u_tt and u_xx are central 
  // difference approximations to the second temporal and spatial derivative respectively. g is a 
  // constant (gamma).

  // these intermediate variables are mostly for clarity and consistency with the mathematical 
  // notation in (1) - production code may get away without them ...or maybe the compiler 
  // optimizes them away anyway:
  int N = getNumGridPoints()-1; // see (1), section 5.2.8
  T c   = waveSpeed;     // what unit?
  T L   = T(1);          // we use the unit interval [0,1], so the spatial length is 1
  T gam = c / L;         // "gamma" - (1), Eq. 6.5
  T k   = timeStep;      // temporal sampling interval
  T h   = L / N;         // spatial sampling interval
  T lam = gam*k / h;     // "lambda", the Courant number
  rsAssert(lam <= T(1), "scheme unstable with these settings!"); // (1), Eq. 6.40
  // actually lamda == 1 is most desirable - in this special case, the numerical solution becomes 
  // exact
  // maybe factor out into getLamda() or getCourantNumber(timeStep)

  // compute updated solution at interior points (factor out to allow switching between schemes):
  T l2 = lam*lam;  // lambda-squared
  for(int l = 1; l <= N-1; l++)                           // (1), section 5.2.8
    tmp[l] = 2*(1-l2)*u[l] + l2*(u[l-1]+u[l+1]) - u1[l];  // (1), Eq 6.35
}
// todo: allow to switch between different schemes
// it's a bit surprising (in a good way), that a central difference in time leads to an explicit 
// rather than implicit recursion
// -> figure out, why -> derive the recursion from the operators

template<class T>
void rsWaveEquation1D<T>::computeInteriorPointsSimple()
{
  for(size_t l = 1; l < u.size()-1; l++) 
    tmp[l] = u[l+1] + u[l-1] - u1[l];    // (1), Eq. 6.54
}

template<class T>
void rsWaveEquation1D<T>::computeBoundaryPoints(T timeStep)
{
  int N  = getNumGridPoints()-1; // see (1), section 5.2.8
  tmp[0] = tmp[N] = T(0);        // endpoints fixed at zero - "Dirichlet" conditions
  // todo: allow to let client code choose from various boundary conditions (Dirichlet, Neumann, 
  // mixed, etc.)
}

template<class T>
void rsWaveEquation1D<T>::updateStateArrays()
{
  int N = getNumGridPoints()-1;             // see (1), section 5.2.8
  RAPT::rsArray::copy(&u[0],   &u1[0], N);  // u goes into u1
  RAPT::rsArray::copy(&tmp[0], &u[0],  N);  // tmp goes into u
}

// todo: implement 6.38, 6.45, 146: bottom, 149: u_l^{n+1} =..., 6.59, 
// 151: implicit scheme, 6.62, 6.66: recursion

// see:
// https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition