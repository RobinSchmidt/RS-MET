#include "PhysicsExperiments.h"

void doublePendulum()
{
  int N  = 150000;    // number of samples
  int fs = 44100;

  // create and set up the double pendulum object:
  rsDoublePendulumDD dp;
  dp.setLength1(1.0);
  dp.setLength2(0.5);
  dp.setMass1(1.0);
  dp.setMass2(0.5);
  dp.setStepSize(0.01);


  // let the pendulum swing and record the locations of both bobs over time:
  vector<double> a1(N), a2(N);                 // angles
  vector<double> m1(N), m2(N);                 // momenta
  vector<double> x1(N), y1(N), x2(N), y2(N);   // x- and y-coordinates
  for(int n = 0; n < N; n++)
  {
    dp.getAngles(&a1[n], &a2[n]);
    dp.getMomenta(&m1[n], &m2[n]);
    dp.getPendulumCoordinates(&x1[n], &y1[n], &x2[n], &y2[n]);
    dp.updateState();
  }

  GNUPlotter plt;
  //plt.addDataArrays(N, &a1[0]);
  //plt.addDataArrays(N, &a2[0]);

  plt.addDataArrays(N, &m1[0]);  // maybe momenta should be scaled by the mass
  plt.addDataArrays(N, &m2[0]);

  //plt.addDataArrays(N, &x1[0]);
  //plt.addDataArrays(N, &x2[0]);

  //plt.addDataArrays(N, &y1[0]);
  //plt.addDataArrays(N, &y2[0]);
  //plt.plot();


  // normalize observables and write into wavefiles:
  RAPT::rsArray::normalize(&m1[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumM1.wav", &m1[0], N, fs, 16);
  RAPT::rsArray::normalize(&m2[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumM2.wav", &m2[0], N, fs, 16);
  RAPT::rsArray::normalize(&x1[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumX1.wav", &x1[0], N, fs, 16);
  RAPT::rsArray::normalize(&x2[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumX2.wav", &x2[0], N, fs, 16);

  // Observations:
  // Notation: m1, m2: masses of the bobs; l1, l2: lengths of the arms
  // We start with the 1st arm horizontal and the 2nd arm upward
  // when  l2 << l1 and m2 << m1, like l1=1, l2=0.01, m1=1, m2=0.01, x2 looks like a distorted
  // sinewave and x2 almost the same.
}
// I think, in order to ensure reproducible results when using chaotic systems, it is necessarry to
// ensure that the rounding behavior is always the same (for example on different platforms, 
// compilers, etc.). There's a compiler switch for the floating point model which can be set to 
// "precise", but i'm not sure, if that's enough. The standard library functions may have different
// implementations. In any case, I think there shouldn't be any algebraically equivalent different
// formulas between versions of a software that uses such chaotic systems.
//
// Idea: 
// -start with a simple linear system with well understood and musically useful behavior, like the 
//  (driven) damped harmonic oscillator: 
//  m*a = -k*x - r*v + F 
//  where m: mass, a: acceleration, k: spring hardness, x: excursion, r: damping/resistence, 
//        v: velocity, F: driving force
// -add nonlinear terms one by one and study their influence on the system, for example a 
//  progressive spring hardening term: -k3*x^3 - this should distort the waveform and increase the 
//  frequency at higher amplitudes
// -generally, instead of using just the linear contributions (m*a, k*x, r*v), we could use 
//  polynomials in x, v, a
// -maybe we could also consider more derivatives like the increase/decrease of acceleration 
//  ("jerk"?)

void heatEquation1D()
{
  int fs = 44100;
  int N  = 2000;    // number of samples

  rsHeatEquation1D<double> hteq;
  hteq.setMaxCycleLength(2048);
  hteq.setDiffusionCoefficient(1.0);
  //hteq.setRandomHeatDistribution(2, 50);
  hteq.setTwoValueDistribution(0.45, 50); 
  hteq.normalizeHeatDistribution();

  std::vector<double> y(N);
  for(int n = 0; n < N; n++)
    y[n] = hteq.getSample();



  rsPlotVector(y);
  //rosic::writeToMonoWaveFile("HeatEquation1D.wav", &y[0], N, fs);

  // it's buzzy and there's a parasitic oscillation at the Nyquist freq.
  // -buzz is probably because of end-handling (try cyclic end-handling to get rid of the buzz)
  // -i think, the parasitic oscillation was due to choosing 1.0 as diffusion coeff - maybe it must
  //  be strictly less than 1
  //  hmmm...with 0.95, there's still s little bit of tha oscillation in the transient
  // -different seeds give wildly different sounds - maybe try a random phase-spectrum with defined
  //  magnitude sprectrum
  // -with cyclic end-handling, it converges to some non-zero DC

  
  //GNUPlotter plt;
}


// maybe to really challenge the blep/blamp class, try to hardsync a sinewave and try to anti-alias
// more higher order derivatives with bladratics, blubics, blartics, etc.

void particleForceDistanceLaw()
{
  // Plots the force vs the distance of the rsPartcielSystem class for various choices of the
  // parameters.

  rsParticleSystem<float> ps(1);
  ps.setForceLawExponent(2);
  //ps.setForceLawOffset(0.1);

  static const int N = 1000;
  static const int numExponents = 7;
  float exponents[numExponents] = { -3, -2, -1, 0, 1, 2, 3 };
  float size1 = 0.01f;
  float size2 = 0.01f;

  float test = ps.getForceByDistance(1, size1, size2);

  float dMin = 0;
  float dMax = 2;
  float d[N];
  float f[numExponents][N];
  RAPT::rsArray::fillWithRangeLinear(d, N, dMin, dMax);
  for(int p = 0; p < numExponents; p++)
  {
    ps.setForceLawExponent(exponents[p]);
    for(int n = 0; n < N; n++)
      f[p][n] = ps.getForceByDistance(d[n], size1, size2);
  }

  GNUPlotter plt;
  plt.addDataArrays(N, d, f[0], f[1], f[2], f[3], f[4], f[5], f[6]);
  plt.plot();

  // hmm...somehow the force-distance laws currently implemented are not really good. Ideally, we 
  // would like the law to satisfy:
  // f(d=0) = y0 (y0: user parameter)
  // f(d=1) = y1 (a1: another user parameter)
  // f(d->inf) -> f^p (p: power, also user parameter, asymptote should work also for p < 0)
  // for y0 = inf, y1 = 1, p = -2, the physically correct (for gravitaion, etc) inverse-square law
  // should come out

  // try: f(d) = (a+d)^p, 1/(a + d^-p), a / (b + d^-p)

  // or maybe the function f(d) = d / (c + d^-p) that we currently effectively use is not so bad
  // it may correspond to cloud-like particles that exert no force on each other when they are at
  // the same location

  // for wolfram: derivative of x / (c + x^p) with respect to x
  // finds also the location where the derivative is 0, evaluating f(x) at that location gives the
  // maximum, solving for c gives c in terms of a desired force-limit

  // maybe give each particle a size and then use f = d / (size1 + size2 + d^p) - that corresponds 
  // to a mental model where each particle is a spherically symmetric cloud of matter, the clouds
  // can peneterate and fall through each other
  // for the force of a uniform cloud, see section 5.7 here:
  // http://www.feynmanlectures.caltech.edu/II_05.html
  // the function above is qualitatively similar bur more smooth, i think (verify)

}

void getTwoParticleTrajectories(rsParticleSystem<float>& ps, int N, float* x1, float* y1, float* z1,
  float* x2, float* y2, float* z2, float* Ek, float* Ep, float* Et, float* Eg, float* Ee, 
  float* Em)
{
  ps.reset();
  for(int n = 0; n < N; n++)
  {
    x1[n] = ps.particles[0].pos.x;
    y1[n] = ps.particles[0].pos.y;
    z1[n] = ps.particles[0].pos.z;

    x2[n] = ps.particles[1].pos.x;
    y2[n] = ps.particles[1].pos.y;
    z2[n] = ps.particles[1].pos.z;

    Ek[n] = ps.getKineticEnergy();
    Ep[n] = ps.getPotentialEnergy();
    Et[n] = ps.getTotalEnergy();

    Eg[n] = ps.getGravitationalPotentialEnergy();
    Ee[n] = ps.getElectricPotentialEnergy();
    Em[n] = ps.getMagneticPotentialEnergy();

    ps.updateState();
  }
}
void particleSystem()
{
  // We simulate a simple system of two particles with unit mass and unit charge to see, if they
  // behave as physically expected (i.e. to see, if the force euqations ae plausible)

  static const int N = 2000; // number of steps in the simulations
  float stepSize = 0.0002f;

  // create and set up the particle system:
  rsParticleSystem<float> ps(2);

  // both particles have unit mass, charge and size:
  ps.particles[0].mass   = 1.f;
  ps.particles[0].charge = 1.f;
  ps.particles[0].size   = 1.f;
  ps.particles[1].mass   = 1.f;
  ps.particles[1].charge = 1.f;
  ps.particles[1].size   = 1.f;


  // place them at (-1,0,0) and (+1,0,0) with zero velocity initially:
  ps.initialPositions[0]  = rsVector3D<float>(-0.5f, +0.0f,  +0.0f);
  ps.initialPositions[1]  = rsVector3D<float>(+0.5f, -0.0f,  +0.0f);
  ps.initialVelocities[0] = rsVector3D<float>( 0.0f, -0.01f, -0.0f);
  ps.initialVelocities[1] = rsVector3D<float>( 0.0f, +0.01f, +0.0f);

  // in a first run, we only let them interact via gravitation - they should attract each other:
  ps.setGravitationalConstant(1.0f);
  ps.setElectricConstant(0.0f);
  ps.setMagneticConstant(0.2f);
  ps.setStepSize(stepSize);
  ps.setForceLawExponent(2.0f); // 2: inverse-square law (asymptotic), 0: force distance-independent
  //ps.setForceLawOffset(1.0);   // 0: non-asymptotic inverse power law

  // record trajectories and energies:
  float x1[N], y1[N], z1[N], x2[N], y2[N], z2[N];  // coordinates
  float Ek[N], Ep[N], Et[N];     // kinetic, potential, total energy
  float Eg[N], Ee[N], Em[N];     // gravitational, electric and magnetic potential energy


  getTwoParticleTrajectories(ps, N, x1, y1, z1, x2, y2, z2, Ek, Ep, Et, Eg, Ee, Em);

  // they initially approach each other, fly through each other and then drift apart forever
  // i suppose, this is due to the singularity, when they are very close (division by zero
  // -> infinite force)


  // test cross-product formula - move to unit-tests:
  rsVector3D<float> test = cross(rsVector3D<float>(2,3,5), rsVector3D<float>(7,11,13)); 
  // bool result = test == rsVector3D<float>(-16,9,1); // ok - cross-product is correct


  GNUPlotter plt;
  float t[N];
  createTimeAxis(N, t, 1/stepSize);
  plt.addDataArrays(N, t, x1, y1, z1, x2, y2, z2);
  //plt.addDataArrays(N, t, x1, x2);
  //plt.addDataArrays(N, t, y1, y2);
  //plt.addDataArrays(N, t, z1, z2);
  //plt.addDataArrays(N, t, Et, Ek, Ep);
  //plt.addDataArrays(N, t, Et, Ek, Eg, Ee, Em);
  //plt.addDataArrays(N, t, Et, Ek, Eg);
  //plt.addDataArrays(N, t, Ek);
  plt.plot();

  // Observations:
  // -when starting at p1=(-1,0.0),p2=(+1,0,0),v1=v2=(0,0.2,0) with only magnetic forces active,
  //  i would expect the x-coordinates converge (initially) and the y-coordinates just linearly 
  //  increase for both partcile (and the z-coordinates stay zero) - but something weird happens
  //  -either the magnetic force computaion is wrong, or the cross-product formula or my 
  //   expectation is wrong
  //  -it seems to depend on the stepsize
  //  -maybe we should apply the stepSize to the velocity update?
}


template<class TArg, class TTol>
inline bool isCloseTo(TArg x, TArg y, TTol tol)
{
  if(rsAbs(x - y) <= tol)
    return true;
  else
    return false;
}
// move to rapt

bool quantumSpin()
{
  // This should be turned into a unit test...

  // create some qubits in pure states:
  typedef rsQuantumSpin<double> QS;
  rsQuantumSpin<double> u, d, l, r, i, o; // maybe use capital letters
  u.prepareUpState();
  d.prepareDownState();
  l.prepareLeftState();
  r.prepareRightState();
  i.prepareInState();
  o.prepareOutState();

  // some work variables:
  std::complex<double> p; // for inner products
  double P;               // for probabilities
  double tol = 1.e-13;    // tolerance for rounding errors
  std::complex<double> one(1,0), zero(0,0), two(2,0), half(.5,0), s(1/sqrt(2.0),0); // i(0,1);
  //double s = 1/sqrt(2.0);
  bool pass = true;

  // check normalization of pure states:
  pass &= (p=u*u) == 1.0;
  pass &= (p=d*d) == 1.0;
  pass &= isCloseTo(p=l*l, one, tol);
  pass &= isCloseTo(p=r*r, one, tol);
  pass &= isCloseTo(p=i*i, one, tol);
  pass &= isCloseTo(p=o*o, one, tol);

  // check orthogonality:
  pass &= (p=u*d) == 0.0;  // (1) Eq 2.3
  pass &= (p=d*u) == 0.0;
  pass &= (p=r*l) == 0.0;
  pass &= (p=l*r) == 0.0;
  pass &= (p=i*o) == 0.0;
  pass &= (p=o*i) == 0.0;

  // check, if Eq 2.8 and 2.9 are satisfied:
  pass &= isCloseTo(p = (o*u)*(u*o), half, tol);
  pass &= isCloseTo(p = (o*d)*(d*o), half, tol);
  pass &= isCloseTo(p = (i*u)*(u*i), half, tol);
  pass &= isCloseTo(p = (i*d)*(d*i), half, tol);
  pass &= isCloseTo(p = (o*r)*(r*o), half, tol);
  pass &= isCloseTo(p = (o*l)*(l*o), half, tol);
  pass &= isCloseTo(p = (i*r)*(r*i), half, tol);
  pass &= isCloseTo(p = (i*l)*(l*i), half, tol);

  // check up-spin probabilities of the various pure spin states:
  pass &= isCloseTo(P = QS::getUpProbability(u), 1.0, tol); // pure up-spin   has P(up) = 1
  pass &= isCloseTo(P = QS::getUpProbability(d), 0.0, tol); // pure down-spin has P(up) = 0
  pass &= isCloseTo(P = QS::getUpProbability(r), 0.5, tol); // all other pure spin states (left, 
  pass &= isCloseTo(P = QS::getUpProbability(l), 0.5, tol); // right, in, out) have up-spin 
  pass &= isCloseTo(P = QS::getUpProbability(i), 0.5, tol); // probability of 1/2
  pass &= isCloseTo(P = QS::getUpProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QS::getRightProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(r), 1.0, tol);
  pass &= isCloseTo(P = QS::getRightProbability(l), 0.0, tol);
  pass &= isCloseTo(P = QS::getRightProbability(i), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QS::getInProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(r), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(l), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(i), 1.0, tol);
  pass &= isCloseTo(P = QS::getInProbability(o), 0.0, tol);

  // set up the random number generator to be used for measurements:
  rsNoiseGenerator<double> prng;
  prng.setRange(0.0, 1.0);


  // test arithmetic operators:
  QS A, B, C;

  A = QS(one, one);               // this is an invalid state because..
  P = QS::getTotalProbability(A); // ..it has total probability 2
  A.normalize();                  // this call normalizes the total probability
  pass &= isCloseTo(P = QS::getTotalProbability(A), 1.0, tol);

  A = u;
  B = d;
  C = s*u + s*d;
  //pass &= (C == r);  // (1) Eq 2.5  ..implement ==
  //pass &= (s*u - s*d == l);


  // maybe check probabilities for some mixed states, i.e. spin-states that are not aligned to
  // any axis

  A.randomizeState(&prng);
  B.randomizeState(&prng);
  pass &= isCloseTo(P = QS::getTotalProbability(A), 1.0, tol);
  pass &= isCloseTo(P = QS::getTotalProbability(B), 1.0, tol);

  //C.randomizeState(); // has unassigned PRNG

  // some test with spin operators:
  typedef rsSpinOperator<double> QSO;
  QSO pauliZ, pauliY, pauliX;
  pauliZ.setToPauliZ();
  pauliY.setToPauliY();
  pauliX.setToPauliX();

  p = pauliZ.getEigenvalue1();  pass &= p == -1.0;
  p = pauliZ.getEigenvalue2();  pass &= p == +1.0;
  A = pauliZ.getEigenvector1(); pass &= A.isCloseTo(d, tol); // "down"
  A = pauliZ.getEigenvector2(); pass &= A.isCloseTo(u, tol); // "up"

  p = pauliX.getEigenvalue1();  pass &= p == -1.0;
  p = pauliX.getEigenvalue2();  pass &= p == +1.0;
  A = pauliX.getEigenvector1(); pass &= A.isCloseTo(l, tol); // "left" - wrong - not normalized
  A = pauliX.getEigenvector2(); pass &= A.isCloseTo(r, tol); // "right"

  p = pauliY.getEigenvalue1();  pass &= p == -1.0;
  p = pauliY.getEigenvalue2();  pass &= p == +1.0;
  A = pauliY.getEigenvector1(); pass &= A.isCloseTo(o, tol); // "out"
  A = pauliY.getEigenvector2(); pass &= A.isCloseTo(i, tol); // "in"

  // test eigenvalue and eigenvector compuation:
  QSO op;
  op.setValues(one, two, two, one);
  std::complex<double> e1 = op.getEigenvalue1(); pass &= e1 == -1.0;
  std::complex<double> e2 = op.getEigenvalue2(); pass &= e2 == +3.0;
  QS E1 = op.getEigenvector1(); // (1, 0)     -> wrong result
  QS E2 = op.getEigenvector2(); // (1,-1) * s





  double r1, r2; // results of 1st and 2nd measurement

  // test:
  A = r; prng.reset(); r1 = A.measureObservable(pauliZ, &prng);
  A = r; prng.reset(); r2 = A.measureSpinZ(&prng);
  // r1 = -1. r2 = +1





  // test spin measurements via Pauli matrices:
  A.prepareDownState();
  p = A.measureObservable(pauliZ, &prng); pass &= p == -1.0;
  p = A.measureObservable(pauliZ, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, d);      pass &= P ==  1.0;
  A.prepareUpState();
  p = A.measureObservable(pauliZ, &prng); pass &= p == +1.0;
  p = A.measureObservable(pauliZ, &prng); pass &= p == +1.0;
  P = QS::getStateProbability(A, u);      pass &= P ==  1.0;

  A.prepareLeftState();
  p = A.measureObservable(pauliX, &prng); pass &= p == -1.0;
  p = A.measureObservable(pauliX, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, l);      pass &= isCloseTo(P, 1.0, tol);
  A.prepareRightState();
  p = A.measureObservable(pauliX, &prng); pass &= p == +1.0;
  p = A.measureObservable(pauliX, &prng); pass &= p == +1.0;
  P = QS::getStateProbability(A, r);      pass &= isCloseTo(P, 1.0, tol);

  A.prepareOutState();
  p = A.measureObservable(pauliY, &prng); pass &= p == -1.0;
  p = A.measureObservable(pauliY, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, o);      pass &= isCloseTo(P, 1.0, tol);
  A.prepareInState();
  p = A.measureObservable(pauliY, &prng); pass &= p == +1.0;
  p = A.measureObservable(pauliY, &prng); pass &= p == +1.0;
  P = QS::getStateProbability(A, i);      pass &= isCloseTo(P, 1.0, tol);


  // test spin measurements via dedicated functions:

  A.prepareDownState();
  p = A.measureSpinZ(&prng);         pass &= p == -1.0;
  p = A.measureSpinZ(&prng);         pass &= p == -1.0;
  P = QS::getStateProbability(A, d); pass &= P ==  1.0;

  A.prepareUpState();
  p = A.measureSpinZ(&prng);         pass &= p == +1.0;
  p = A.measureSpinZ(&prng);         pass &= p == +1.0;
  P = QS::getStateProbability(A, u); pass &= P ==  1.0;


  A.prepareLeftState();
  p = A.measureSpinX(&prng);         pass &= p == -1.0;
  p = A.measureSpinX(&prng);         pass &= p == -1.0;
  P = QS::getStateProbability(A, l); pass &= isCloseTo(P, 1.0, tol);

  A.prepareRightState();
  p = A.measureSpinX(&prng);         pass &= p == +1.0;
  p = A.measureSpinX(&prng);         pass &= p == +1.0;
  P = QS::getStateProbability(A, r); pass &= isCloseTo(P, 1.0, tol);

  A.prepareOutState();
  p = A.measureSpinY(&prng);         pass &= p == -1.0;
  p = A.measureSpinY(&prng);         pass &= p == -1.0;
  P = QS::getStateProbability(A, o); pass &= isCloseTo(P, 1.0, tol);

  A.prepareInState();
  p = A.measureSpinY(&prng);         pass &= p == +1.0;
  p = A.measureSpinY(&prng);         pass &= p == +1.0;
  P = QS::getStateProbability(A, i); pass &= isCloseTo(P, 1.0, tol);


  // now, do some actual measurements:

  int N = 100;   // number of measurements to take
  std::vector<double> spins1(N);
  int n;
  prng.reset();
  for(n = 0; n < N; n++) {
    A = r;                      // initialize state - todo: mayb try different states
    r1 = A.measureSpinZ(&prng); // should have a 50/50 chance to be +1 or -1
    r2 = A.measureSpinZ(&prng); // a 2nd measurement must always give the same result
    pass &= (r1 == r2);
    spins1[n] = r1;
  }
  double mean1 = rsMean(spins1);
  //rsPlotVector(spins1);

  // now do the same thing using the Pauli matrix for the z-component - this is more general but 
  // the computations are also more expensive, the results should be the same
  std::vector<double> spins2(N);
  prng.reset();
  for(n = 0; n < N; n++) {
    A = r;
    r1 = A.measureObservable(pauliZ, &prng);
    r2 = A.measureObservable(pauliZ, &prng);
    pass &= (r1 == r2);
    spins2[n] = r1;
  }
  double mean2 = rsMean(spins2);
  //rsPlotVectors(spins1, spins2);

  // now test, if the matrix-based and dedicated function based measurements give the same results
  // (note that they operate on different states because the first measurement alters the state):
  prng.reset();
  for(n = 0; n < N; n++) {
    A.randomizeState(&prng);
    r1 = A.measureObservable(pauliZ, &prng);
    r2 = A.measureSpinZ(&prng);
    pass &= (r1 == r2);

    A.randomizeState(&prng);
    r2 = A.measureSpinZ(&prng);
    r1 = A.measureObservable(pauliZ, &prng);
    pass &= (r1 == r2);

    A.randomizeState(&prng);
    r1 = A.measureObservable(pauliX, &prng);
    r2 = A.measureSpinX(&prng);
    pass &= (r1 == r2);

    A.randomizeState(&prng);
    r2 = A.measureSpinX(&prng);
    r1 = A.measureObservable(pauliX, &prng);
    pass &= (r1 == r2);

    A.randomizeState(&prng);
    r1 = A.measureObservable(pauliY, &prng);
    r2 = A.measureSpinY(&prng);
    pass &= (r1 == r2);

    A.randomizeState(&prng);
    r2 = A.measureSpinY(&prng);
    r1 = A.measureObservable(pauliY, &prng);
    pass &= (r1 == r2);
  }

  // now do a similar test that ensures, that they operate on the same state:
  rsNoiseGenerator<double> prng1, prng2;
  prng1.setRange(0.0, 1.0);
  prng1.reset();
  prng2.setRange(0.0, 1.0);
  prng2.reset();
  for(n = 0; n < N; n++)
  {
    B.randomizeState(&prng);

    A = B; r1 = A.measureObservable(pauliZ, &prng1);
    A = B; r2 = A.measureSpinZ(&prng2);
    pass &= (r1 == r2);

    A = B; r1 = A.measureObservable(pauliX, &prng1);
    A = B; r2 = A.measureSpinX(&prng2);
    pass &= (r1 == r2);

    A = B; r1 = A.measureObservable(pauliY, &prng1);
    A = B; r2 = A.measureSpinY(&prng2);
    pass &= (r1 == r2);
  }

  // todo: test, if the statistical distribution is as desired - set it into a state
  // au = sqrt(0.8), ad = sqrt(0.2) - we should see roughly 80% "up" measurements and 20% down
 

  // todo: implement quantum gates (and, or, Hadamard, cnot, toffoli)
  // this here:
  // https://www.youtube.com/watch?v=ZN0lhYU1f5Q
  // says: measure, hadamard, phase, T (rotate |1> by pi/4), cnot


  // https://homepages.cwi.nl/~rdewolf/qcnotes.pdf
  // http://mmrc.amss.cas.cn/tlb/201702/W020170224608150244118.pdf

  return pass;
  //GNUPlotter plt;

  // todo: maybe use float instead of double
}
// Notes:
// I think, the relationship to what is called the "wavefunction" in quantum mechanics is as 
// follows: The wavefunction is in general some function from a set S into the complex numbers. 
// For position and momentum of a particle, that set S is the set of real numbers R (for 1D) or 
// R^2 or R^3 for 2D and 3D space. In our case, the set S is just S = { up, down } - which can be
// renamed to S = { 0, 1 } or S = { |0>, |1> } to get the common qubit noation. The "collapse of 
// the wavefunction" occurs when we assign the pure up/down states in the measurement operation. In
// a general state, both of these pure states have a complex number associated with them - the 
// "probability amplitude" - and the square of its magnitude gives the actual probability. When the
// wavefunctions is collapsed due to a measurement one the values becomes 1 and the other 0 - it 
// becomes a delta distribution - although its spikey nature is not really obvious in this simple 
// case where we have only two possible input values into the function.


void tennisRacket()
{
  // Numerically integrates the system of differential equations describing the angular velocities 
  // of a rigid object about its 3 principal axes of rotation. Rotation around the axis with the 
  // intermediate moment of inertia is unstable and flips over periodically whenever there is a 
  // small amount of initial angular velocity along any of the other two axes. It's an unstable
  // equilibrium of the Euler equations.
  // see:
  // https://en.wikipedia.org/wiki/Tennis_racket_theorem
  // https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
  // https://www.youtube.com/watch?v=1VPfZ_XzisU
  // https://arxiv.org/pdf/1606.08237.pdf

  // User parameters:
  int N = 5000;       // number of samples
  double h = 0.01;    // step-size ("delta-t")
  double I1, I2, I3;  // the 3 moments inertia
  double c = 2;       // ratio of I1/I2 and I2/I3
  I1 =  c;
  I2 =  1;
  I3 =  1/c;
  double w1, w2, w3;  // the 3 (initial) angular velocities (with respect to the principal axes)
  w1 = 0.01;
  w2 = 1;
  w3 = 0.0;
  bool renormalize = true; // renormalize rotational energy in each step (counteract numeric drift)

  // coefficients for the creative extra terms:
  double p = 0.0;   // a small negative number (-0.01) leads to some asymmetry (and decay, which 
                    // can probably be compensated by enforcing constant rotational energy and/or
                    // angular momentum)
  // move this into a 2nd implementation with more creative add-ons

  // create time axis and vectors to hold the results (w1,w2,w3 as functions of time):
  std::vector<double> t(N), W1(N), W2(N), W3(N);
  double a1, a2, a3;  // the 3 angular accelerations
  double E, E0;       // rotational energy and its initial value
  double k1, k2, k3;  // multipliers resulting from the moments of inertia
  k1 = (I2-I3)/I1;
  k2 = (I3-I1)/I2;
  k3 = (I1-I2)/I3;
  E0 = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3)/2; // https://en.wikipedia.org/wiki/Moment_of_inertia#Kinetic_energy_2

  // Use forward Euler method to integrate the ODE system of the Euler rotation equations 
  // (...Euler is every-fucking-where):
  for(int i = 0; i < N; i++)
  {
    // compute absolute time and record angular velocities:
    t[i]  = i*h;
    W1[i] = w1;
    W2[i] = w2;
    W3[i] = w3;

    // compute angular accelerations:
    a1 = k1*w2*w3;
    a2 = k2*w3*w1 + p*w2;
    a3 = k3*w1*w2;

    // update angular velocities:
    w1 += h*a1;
    w2 += h*a2;
    w3 += h*a3;

    // optionally renormalize rotational energy:
    E = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3)/2;
    if(renormalize) {
      double r = sqrt(E0/E);
      w1 *= r;
      w2 *= r;
      w3 *= r;
    }
  }

  // Plot w1,w2,w3 as functions of time:
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &W1[0], &W2[0], &W3[0]);
  plt.plot();

  //rosic::writeToMonoWaveFile("TennisRacket.wav", &W2[0], N, 44100);

  // Observations:
  // -if there's a small initial w1 or w3 component, the instability lets the w2 component (blue) 
  //  flip back and forth periodically
  // -the excursions of w3 (green) are greater than those of w1 (black)
  // -the difference in the strength of these excursions depends of the ratio of I1 and I3
  // -in case of an initial w1 component, the w1 excursions are always positive and the w3 
  //  excursions alternate
  // -in case of an initial w3 component, the w1 excursions alternate and the w3 excursions are 
  //  always positive
  // -choosing I1=4,I2=5,I3=1 (i.e. the middle I is greatest), we get a stable rotation, as the 
  //  theory predicts (w2 stays constant) 
  //  -there is some (sinusoidal?) wobble between w1 and w3 - they seem to exchange energy
  // -choosing I1=4,I2=0.5,I3=1 (i.e. the middle I is smallest), w2 wobbles a little bit, w1 and w3
  //  also wobble but with unequal amounts (w3 wobbles more)
  // -i guess, the instability around the intermediate axis is due to k2 being negative while k1,k3
  //  are positive?
  // -using I=(8,4,1) instead of I=(4,2,1) results in an increase of frequency of the flips
  //  -guess: is the frequency proportional to the magnitude to the vector valued moment of inertia 
  //   sqrt(I1^2 + I2^2 + I3^2)
  // -the frequency of the flips seems to go up with the length of the initial angular velocity 
  //  vector but also with the ratio of the initial disturbance to w2 (figure out details)



  // todo:
  // -try to figure out, how the frequency of the flips depends on the various parameters - can we 
  //  find a formula - or even better, a sort of inverse formula that let's us dial in the 
  //  frequency as user parameter?
  //  -nope: I=(4,2,1) has the same frequency as I=(8,4,2), so the absolute numbers seem irrelevant
  //  -I=(16,4,2) has much higher frequency than I=(8,4,2)
  //  -try I = (c,1,1/c) and figure out frequency f as function of c (make plots)
  //   -c=1 -> f=0

  //  ...maybe it depends on the total mass or the norm of the inertia vector?
  // -figure out effects of having initial nonzero values for both, w1 and w3 
  // -figure out effects of the sign(s) of the initial angular velocities
  // -compute angular momentum and rotational energy (as functions of time) - they should remain
  //  constant
  // -figure out the formula for the angular momentum in the fixed "lab" reference frame - maybe 
  //  enforce this to stay constant

  // -according to this https://arxiv.org/pdf/1606.08237.pdf (section 3), the system has as analytic 
  //  solutions the jacobi elliptic functions cn,sn,dn - plot them, too
}

// see also:
// https://en.wikipedia.org/wiki/Moment_of_inertia
// https://en.wikipedia.org/wiki/Poinsot%27s_ellipsoid
// https://en.wikipedia.org/wiki/Polhode

template<class T>
T rotationalEnergy(const rsVector3D<T>& I, const rsVector3D<T>& w)
{
  return T(0.5) * (I.x*w.x*w.x + I.y*w.y*w.y + I.z*w.z*w.z);
  // https://en.wikipedia.org/wiki/Moment_of_inertia#Kinetic_energy_2
}

template<class T>
void plotVectorComponents(std::vector<T>& t, std::vector<rsVector3D<T>>& v)
{
  int N = (int) t.size();
  rsAssert((int)v.size() == N);
  std::vector<T> x(N), y(N), z(N);
  for(int n = 0; n < N; n++) {
    x[n] = v[n].x;
    y[n] = v[n].y;
    z[n] = v[n].z;
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &x[0], &y[0], &z[0]);
  plt.setPixelSize(1000, 300);
  plt.plot();
}
// move to rapt, make input vectors const (requires changes to GNUPlotter)


void tennisRacket2()
{
  // In contrast to tennisRacket, we use vectors here and introduce some additional tweaks in an 
  // attempt to provide a musically viable set of user parameters so this thing can be used as a 
  // waveform generator or resonator.

  // User parameters:
  typedef rsVector3D<double> Vec;
  int    N = 8000;          // number of samples
  double h = 0.01;          // step-size ("delta-t")
  Vec I(   4, 2, 1);        // moments of inertia along principa axes
  Vec w(0.01, 1, 0);        // initial angular velocities
  bool renormalize = false; // keep energy constant

  // experimental parameters: 
  Vec d(0.0, 0.02, 0.0); // damping/dissipation - works only when renormalization is off


  // integrate ODE system via forward Euler method :
  Vec a;                     // angular acceleration
  Vec M(0, 0, 0);            // external torque vector (zero at the moment)
  std::vector<double> t(N);  // time axis
  std::vector<Vec>    W(N);  // vectorial angular velocity as function of time
  double E;                  // rotational energy...
  double E0 = rotationalEnergy(I, w); //...and its initial value
  for(int i = 0; i < N; i++)
  {
    // compute absolute time and record angular velocity vector:
    t[i] = i*h;
    W[i] = w;

    // compute angular acceleration vector:
    a.x = (M.x - (I.z - I.y)*w.y*w.z) / I.x;
    a.y = (M.y - (I.x - I.z)*w.z*w.x) / I.y;
    a.z = (M.z - (I.y - I.x)*w.x*w.y) / I.z;
    // formula from https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)

    // todo: add experimental extra terms to the angular acceleration...
    a.x -= d.x*w.x / I.x;
    a.y -= d.y*w.y / I.y;
    a.z -= d.z*w.z / I.z;


    // update angular velocity vector:
    w += h*a;

    // optionally renormalize rotational energy:
    E = rotationalEnergy(I, w);
    if(renormalize)
      w *= sqrt(E0/E);
  }

  // -when only the 1st axis is damped, it tends towards rotating only around the 3rd axis
  // -when only the 3rd axis is damped, it tends towards rotating only around the 1st axis
  // -when only the 2nd axis is damped, it dies out

  // todo: 
  // -try asymmetrical damping that affects the both half-waves of the "square-wave" differently
  //  ->goal: some sort of pulse-width control

  plotVectorComponents(t, W);
}

void tennisRacket3()
{
  // This produces the same plot as the functions above but uses the rsTennisRacket class

  int N = 8000; 
  double h = 0.01;
  rsTennisRacket<double> tr;
  tr.setInertiaRatio(2);
  tr.setState(0.01, 1, 0);  // initial state
  tr.setStepSize(h);

  std::vector<double> t(N), w1(N), w2(N), w3(N);
  for(int n = 0; n < N; n++) {
    t[n]  = n*h;
    w1[n] = tr.getW1();
    w2[n] = tr.getW2();
    w3[n] = tr.getW3();
    tr.updateState(0, 0, 0);
  }

  GNUPlotter plt;
  plt.plotFunctionTables(N, &t[0], &w1[0], &w2[0], &w3[0]);
}

// todo: make a function tennisRacketFreqs that plots the frequency or period of the flipping
// against the inertia ratio - goal: figure out how the inertia ratio controls the freq

double tennisRacketPeriod(rsTennisRacket<double>& tr)
{
  // todo: estimate period of the tennis racket by measuring the distance between two upward
  // zero crossings...

  int first = -1, second = -1;  // sample instants of zero crossings
  double x1 = tr.getSample(0);
  int n = 0;
  while(second == -1) {
    double x  = tr.getSample(0);
    if(x1 <= 0 && x > 0) // zero crossing found
    {
      if(first == -1)
        first = n;
      else
        second = n;
    }
    x1 = x;
    n++;
  }
  return second - first;

  // we could refine this to estimate zeros with subsample precision - but that may be overkill 
  // here
}

void tennisRacketFreq()
{
  // We plot the frequency of the flipping as function of the inertia ratio


  rsTennisRacket<double> tr;
  tr.setStepSize(0.01);

  std::vector<double> r = { 1.01, 1.125,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 };
  //std::vector<double> r = { 2,4,8,16,32,64 };
  int N = (int) r.size();
  std::vector<double> p(N), f(N); // period and frequency

  for(int i = 0; i < N; i++)
  {
    tr.setInertiaRatio(r[i]);
    tr.setState(0.01, 1, 0);  // initial state
    p[i] = tennisRacketPeriod(tr);
    f[i] = 1/p[i];
  }

  GNUPlotter plt;
  //plt.plotFunctionTables(N, &r[0], &p[0]);
  plt.plotFunctionTables(N, &r[0], &f[0]);

  // It looks similar to a straight line, but not quite...
  // maybe some variant of this?
  // https://en.wikipedia.org/wiki/Logarithmic_integral_function
  // https://de.wikipedia.org/wiki/Integrallogarithmus#/media/Datei:Li10.png

  // -maybe plot the numerical derivative
  // -maybe use the logarithm of the ratio for the x-axis - let it extend also to negative values,
  //  i.e. r < 1 - maybe it's symmetrical?
}