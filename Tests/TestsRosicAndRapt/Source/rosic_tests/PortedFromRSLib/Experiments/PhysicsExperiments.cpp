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



//template<class T>
//inline double rsAbs(rsVector2D<T> v)
//{
//  return (double) sqrt(v.x*v.x + v.y*v.y); // abs(vector) is defined as length(vector)
//}
// that's a bit weird naming ("absolute value of a vector") but is needed in isCloseTo for vector
// arguments - find a better way

template<class TArg, class TTol>
inline bool isCloseTo(TArg x, TArg y, TTol tol)
{
  if(rsAbs(x - y) <= tol)
    return true;
  else
    return false;
}
// move to rapt or test utilities

inline bool isCloseTo(
  rsVector2D<std::complex<double>> v, rsVector2D<std::complex<double>> w, double tol)
{
  bool r = true;
  rsVector2D<std::complex<double>> d = w-v;
  r &= d.x.real() <= tol;
  r &= d.x.imag() <= tol;
  r &= d.y.real() <= tol;
  r &= d.y.imag() <= tol;
  return r;
}
// move to



bool quantumSpinMeasurement2()
{
  // This should be turned into a unit test...
  bool pass = true;   // move to unit tests

  typedef rsQuantumSpinFunctions<double> QF;
  typedef rsVector2D<std::complex<double>>  Vec;
  typedef rsMatrix2x2<std::complex<double>> Mat;
  typedef rsNoiseGenerator<double> PRNG;

  // some work variables:
  Vec A, B, C;            // for generic states
  std::complex<double> p; // for inner products
  double P;               // for probabilities
  double tol = 1.e-13;    // tolerance for rounding errors
  double r1, r2;          // results of 1st and 2nd measurement
  std::complex<double> one(1,0), zero(0,0), two(2,0), half(.5,0), s(1/sqrt(2.0),0);
  PRNG prng;  // randum number genertor for measurements
  prng.setRange(0, 1);  // this is super important and a source for trouble if forgotten
                        // todo: allow ony a specific kind of PRNG that only has range 0..1

  // create some qubits in pure states:
  Vec u, d, l, r, i, o; // maybe use capital letters
  QF::prepareUpState(u);
  QF::prepareDownState(d);
  QF::prepareLeftState(l);
  QF::prepareRightState(r);
  QF::prepareInState(i);
  QF::prepareOutState(o);


  // check normalization of pure states:
  pass &= (p=QF::bracket(u,u)) == 1.0;
  pass &= (p=QF::bracket(d,d)) == 1.0;
  pass &= isCloseTo(p=QF::bracket(l,l), one, tol);
  pass &= isCloseTo(p=QF::bracket(r,r), one, tol);
  pass &= isCloseTo(p=QF::bracket(i,i), one, tol);  // fails
  pass &= isCloseTo(p=QF::bracket(o,o), one, tol);  // fails

  // check orthogonality:
  pass &= (p=QF::bracket(u,d)) == 0.0;  // (1) Eq 2.3
  pass &= (p=QF::bracket(d,u)) == 0.0;
  pass &= (p=QF::bracket(r,l)) == 0.0;
  pass &= (p=QF::bracket(l,r)) == 0.0;
  pass &= (p=QF::bracket(i,o)) == 0.0;  // fails
  pass &= (p=QF::bracket(o,i)) == 0.0;  // fails

  // check, if Eq 2.8 and 2.9 are satisfied:
  pass &= isCloseTo(p = QF::bracket(o,u) * QF::bracket(u,o), half, tol);
  pass &= isCloseTo(p = QF::bracket(o,d) * QF::bracket(d,o), half, tol);
  pass &= isCloseTo(p = QF::bracket(i,u) * QF::bracket(u,i), half, tol);
  pass &= isCloseTo(p = QF::bracket(i,d) * QF::bracket(d,i), half, tol);
  pass &= isCloseTo(p = QF::bracket(o,r) * QF::bracket(r,o), half, tol);
  pass &= isCloseTo(p = QF::bracket(o,l) * QF::bracket(l,o), half, tol);
  pass &= isCloseTo(p = QF::bracket(i,r) * QF::bracket(r,i), half, tol);
  pass &= isCloseTo(p = QF::bracket(i,l) * QF::bracket(l,i), half, tol);

  // check up-spin probabilities of the various pure spin states:
  pass &= isCloseTo(P = QF::getUpProbability(u), 1.0, tol); // pure up-spin   has P(up) = 1
  pass &= isCloseTo(P = QF::getUpProbability(d), 0.0, tol); // pure down-spin has P(up) = 0
  pass &= isCloseTo(P = QF::getUpProbability(r), 0.5, tol); // all other pure spin states (left, 
  pass &= isCloseTo(P = QF::getUpProbability(l), 0.5, tol); // right, in, out) have up-spin 
  pass &= isCloseTo(P = QF::getUpProbability(i), 0.5, tol); // probability of 1/2
  pass &= isCloseTo(P = QF::getUpProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QF::getRightProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QF::getRightProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QF::getRightProbability(r), 1.0, tol);
  pass &= isCloseTo(P = QF::getRightProbability(l), 0.0, tol);
  pass &= isCloseTo(P = QF::getRightProbability(i), 0.5, tol);
  pass &= isCloseTo(P = QF::getRightProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QF::getInProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QF::getInProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QF::getInProbability(r), 0.5, tol);
  pass &= isCloseTo(P = QF::getInProbability(l), 0.5, tol);
  pass &= isCloseTo(P = QF::getInProbability(i), 1.0, tol);
  pass &= isCloseTo(P = QF::getInProbability(o), 0.0, tol);

  // test some stuff:
  A = Vec(one, one);              // this is an invalid state because..
  P = QF::getTotalProbability(A); // ..it has total probability 2
  QF::normalizeState(A);          // this call normalizes the total probability
  pass &= isCloseTo(P = QF::getTotalProbability(A), 1.0, tol);
  A = u;
  B = d;
  C = s*u + s*d;
  pass &= (C == r);               // (1) Eq 2.5 
  pass &= (s*u - s*d == l);
  QF::randomizeState(A, &prng);
  QF::randomizeState(B, &prng);
  pass &= isCloseTo(P = QF::getTotalProbability(A), 1.0, tol);
  pass &= isCloseTo(P = QF::getTotalProbability(B), 1.0, tol);
  pass &= QF::bracket(A,B) == conj(QF::bracket(B,A));


  std::complex<double> e1, e2;
  Vec E1, E2;

  // some test with spin operators:
  Mat pauliZ, pauliY, pauliX;
  QF::setToPauliZ(pauliZ);
  QF::setToPauliY(pauliY);
  QF::setToPauliX(pauliX);

  e1 = pauliZ.eigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliZ.eigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliZ.eigenvector1(); pass &= isCloseTo(E1, d, tol); // "down"
  E2 = pauliZ.eigenvector2(); pass &= isCloseTo(E2, u, tol); // "up"
  // E1 ane E2 are swapped - why?

  e1 = pauliX.eigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliX.eigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliX.eigenvector1(); pass &= isCloseTo(E1, l, tol); // "left" - wrong - not normalized
  E2 = pauliX.eigenvector2(); pass &= isCloseTo(E2, r, tol); // "right"

  e1 = pauliY.eigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliY.eigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliY.eigenvector1(); pass &= isCloseTo(E1, o, tol); // "out"
  E2 = pauliY.eigenvector2(); pass &= isCloseTo(E2, i, tol); // "in"

  // test eigenvalue and eigenvector compuation:
  Mat op;
  op.setValues(one, two, two, one);
  e1 = op.eigenvalue1(); pass &= e1 == -1.0;
  e2 = op.eigenvalue2(); pass &= e2 == +3.0;

  //E1 = op.eigenvector1(); // (1, 0)     -> wrong result
  //E2 = op.eigenvector2(); // (1,-1) * s
  //E1 = pauliZ.eigenvector1(); 
  //E2 = pauliZ.eigenvector2();


  // test spin measurements via Pauli matrices:
  QF::prepareDownState(A);
  p = QF::measureObservable(A, pauliZ, &prng); pass &= p == -1.0;
  p = QF::measureObservable(A, pauliZ, &prng); pass &= p == -1.0;
  P = QF::getStateProbability(A, d);           pass &= P ==  1.0;
  QF::prepareUpState(A);
  p = QF::measureObservable(A, pauliZ, &prng); pass &= p == +1.0;
  p = QF::measureObservable(A, pauliZ, &prng); pass &= p == +1.0;
  P = QF::getStateProbability(A, u);           pass &= P ==  1.0;

  QF::prepareLeftState(A);
  p = QF::measureObservable(A, pauliX, &prng); pass &= p == -1.0;
  p = QF::measureObservable(A, pauliX, &prng); pass &= p == -1.0;
  P = QF::getStateProbability(A, l);           pass &= isCloseTo(P, 1.0, tol);
  QF::prepareRightState(A);
  p = QF::measureObservable(A, pauliX, &prng); pass &= p == +1.0;
  p = QF::measureObservable(A, pauliX, &prng); pass &= p == +1.0;
  P = QF::getStateProbability(A, r);           pass &= isCloseTo(P, 1.0, tol);

  QF::prepareOutState(A);
  p = QF::measureObservable(A, pauliY, &prng); pass &= p == -1.0;
  p = QF::measureObservable(A, pauliY, &prng); pass &= p == -1.0;
  P = QF::getStateProbability(A, o);           pass &= isCloseTo(P, 1.0, tol);
  QF::prepareInState(A);
  p = QF::measureObservable(A, pauliY, &prng); pass &= p == +1.0;
  p = QF::measureObservable(A, pauliY, &prng); pass &= p == +1.0;  // fails
  P = QF::getStateProbability(A, i);           pass &= isCloseTo(P, 1.0, tol);

  // test spin measurements via dedicated functions:
  QF::prepareDownState(A);
  p = QF::measureSpinZ(A, &prng);    pass &= p == -1.0;
  p = QF::measureSpinZ(A, &prng);    pass &= p == -1.0;
  P = QF::getStateProbability(A, d); pass &= P ==  1.0;
  QF::prepareUpState(A);
  p = QF::measureSpinZ(A, &prng);    pass &= p == +1.0;
  p = QF::measureSpinZ(A, &prng);    pass &= p == +1.0;
  P = QF::getStateProbability(A, u); pass &= P ==  1.0;

  QF::prepareLeftState(A);
  p = QF::measureSpinX(A, &prng);    pass &= p == -1.0;
  p = QF::measureSpinX(A, &prng);    pass &= p == -1.0;
  P = QF::getStateProbability(A, l); pass &= isCloseTo(P, 1.0, tol);
  QF::prepareRightState(A);
  p = QF::measureSpinX(A, &prng);    pass &= p == +1.0;
  p = QF::measureSpinX(A, &prng);    pass &= p == +1.0;
  P = QF::getStateProbability(A, r); pass &= isCloseTo(P, 1.0, tol);

  QF::prepareOutState(A);
  p = QF::measureSpinY(A, &prng);    pass &= p == -1.0;
  p = QF::measureSpinY(A, &prng);    pass &= p == -1.0;
  P = QF::getStateProbability(A, o); pass &= isCloseTo(P, 1.0, tol);
  QF::prepareInState(A);
  p = QF::measureSpinY(A, &prng);    pass &= p == +1.0;
  p = QF::measureSpinY(A, &prng);    pass &= p == +1.0;
  P = QF::getStateProbability(A, i); pass &= isCloseTo(P, 1.0, tol);


  // test with random states, if matrix based and dedicated function based measurements do the
  // same thing:
  int n;
  int N = 1000;
  rsNoiseGenerator<double> prng1, prng2;
  prng1.setRange(0.0, 1.0);
  prng1.reset();
  prng2.setRange(0.0, 1.0);
  prng2.reset();
  for(n = 0; n < N; n++)
  {
    QF::randomizeState(B, &prng);

    A = B; r1 = QF::measureObservable(A, pauliZ, &prng1);
    A = B; r2 = QF::measureSpinZ(A, &prng2);
    pass &= (r1 == r2);

    A = B; r1 = QF::measureObservable(A, pauliX, &prng1);
    A = B; r2 = QF::measureSpinX(A, &prng2);
    pass &= (r1 == r2);

    A = B; r1 = QF::measureObservable(A, pauliY, &prng1);
    A = B; r2 = QF::measureSpinY(A, &prng2);
    pass &= (r1 == r2);
  }
  // todo: make additional meausrements after the collapse - they should always give the same 
  // results

  // test, if the statistical distribution is as desired - set it into a state
  // au = sqrt(0.8), ad = sqrt(0.2) - we should see roughly 80% "up" measurements and 20% "down"
  B = Vec(std::complex<double>(sqrt(0.8), 0), std::complex<double>(sqrt(0.2), 0));
  std::vector<double> spins(N);
  for(n = 0; n < N; n++) {
    A = B;
    spins[n] = QF::measureSpinZ(A, &prng);
  }
  double mean = rsMean(spins);
  P = (mean+1)/2;  // -1..+1 -> 0..1  P = 0.81 (with N=1000) - close enough to 0.8


  // compute expectation values for the spin component measurements in state B:
  double Ez, Ex, Ey;  
  Ez = QF::getExpectedMeasurement(pauliZ, B); // should be 0.6 = (0.8*2)-1 and close to mean1
  Ex = QF::getExpectedMeasurement(pauliX, B); // 0.8
  Ey = QF::getExpectedMeasurement(pauliY, B); // 0.0
  P  = Ex*Ex + Ey*Ey + Ez*Ez;  // should be one according to (1) Eq 3.27
  pass &= isCloseTo(P,  1.0, tol);
  pass &= isCloseTo(Ez, 0.6, tol);
  // Ex = 0.8 and Ey = 0 - maybe try to compute these manually

  // OK - states and measurements (chapters 2, 3) are done and seem to work
  // next: evolution of states, i.e. applying unitary operators to modify a state
  // todo: implement quantum gates (and, or, Hadamard, cnot, toffoli)

  // structure the tests:
  // 1: states, probabilities
  // 2: measurements
  // 3: state evolution
  // 4: entanglement


  // this here:
  // https://www.youtube.com/watch?v=ZN0lhYU1f5Q
  // says: measure, hadamard, phase, T (rotate |1> by pi/4), cnot

  // https://homepages.cwi.nl/~rdewolf/qcnotes.pdf
  // http://mmrc.amss.cas.cn/tlb/201702/W020170224608150244118.pdf


  rsAssert(pass);
  return pass;
}

template<class T>
void plotQuantumSpinStateTrajectory(
  const std::vector<rsVector2D<std::complex<T>>>& Psi, T dt)
{
  typedef rsQuantumSpinFunctions<double> QF;
  int N = (int) Psi.size();
  std::vector<T> t(N), dr(N), di(N), ur(N), ui(N);
  for(int n = 0; n < N; n++) {
    t[n] = n*dt;   // time axis
    dr[n] = QF::getDownAmplitude(Psi[n]).real(); // rename to DownAmplitude
    di[n] = QF::getDownAmplitude(Psi[n]).imag();
    ur[n] = QF::getUpAmplitude(Psi[n]).real();
    ui[n] = QF::getUpAmplitude(Psi[n]).imag();
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &dr[0], &di[0], &ur[0], &ui[0]);
  plt.plot();
}

bool quantumSpinEvolution()
{
  bool pass = true;   // move to unit tests - or maybe not

  // Simulates a spin (i.e. a spinning electron) in a magnetic field as explained in (1) pg 116 ff. 
  // The magnetic field is aligned with the z-axis. For a classical spinning charge, the energy of 
  // the system is proportional to the dot product of the spin and the field. In quantum mechanics,
  // the Hamiltonian H is proportional to the z-spin operator (the observable associated with the
  // Pauli matrix sigma_z. ....
  //
  // References:
  //  (1) The Theoretical Minimum - Quantum Mechanics (Leonard Susskind, Art Friedman)


  typedef std::complex<double> Complex;
  typedef rsQuantumSpinFunctions<double> QF;
  typedef rsVector2D<Complex>  Vec;
  typedef rsMatrix2x2<Complex> Mat;


  double w    = 1;
  double hBar = 1; // we use https://en.wikipedia.org/wiki/Planck_units

  // set up the random number generator to be used for measurements:
  rsNoiseGenerator<double> prng;
  prng.setRange(0.0, 1.0);

  // Create an initial spin state:
  Vec Psi;
  QF::randomizeState(Psi, &prng);

  // Define the Hamiltonian (Eq 4.23):
  Mat H;
  QF::setToPauliZ(H);
  H = Complex(hBar*w/2) * H;

  // numerically integrate the time dependent Schrödinger equation using the forward Euler method:
  int n, N = 3000;
  Complex i(0, 1);  // imaginary unit
  std::vector<Vec> stateTrajectory(N); // records the trajectory of our spin state Psi
  double step = 0.01;                  // integration step size
  Complex cStep = Complex(step);       // ...needs to be complexified 
  for(n = 0; n < N; n++) {
    stateTrajectory[n] = Psi;
    Vec dPsi = -i * H * Psi;   // (1) Eq 4.9 or 4.10 (time dependent Schrödinger equation)
    Psi = Psi + cStep * dPsi;  // forward Euler step
    QF::normalizeState(Psi);   // avoid divergence due to error build up
  }
  plotQuantumSpinStateTrajectory(stateTrajectory, step);


  // Notes:
  // -The time dependent Schrödinger equation dPsi = -i * H * Psi actually encapsulates a system of
  //  4 differential equations - the state vector Psi has two elements and each is a complex number 
  //  (todo: write the 4 equations out in full). If the state vector would be N-dimensional, it 
  //  would be a system of 2*N (first order) differential equations. If the state "vector" would be 
  //  a continuous function, it would be a partial differential equation and H would be a continuous
  //  differntial operator (right?).
  // -It seems, we are getting some sort of rotation in 4D space

  // todo: take as observable not the spin along the z-axis but along an arbitrary axis (i think, along 
  // the z-axis, the spin is constant anyway?

  // todo: factor out the computation of a state trajectory, given an initial state, a Hamiltonian
  // a number of steps and a step-size
  // move the plotting into a different function

  // plot not only the state but also the expectation value(s) of some observable(s) as function(s)
  // of time, see if this matches hat Eq 4.17 predicts

  // todo: can we accumulate all the steps into a single big step represented by a matrix U(t) 
  // (the time-development operator, see Eq 4.1) that has accumulated all the little steps up to t?

  return pass;
}


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