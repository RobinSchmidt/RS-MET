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

/*
I think, in order to ensure reproducible results when using chaotic systems, it is necessarry to
ensure that the rounding behavior is always the same (for example on different platforms, 
compilers, etc.). There's a compiler switch for the floating point model which can be set to 
"precise", but i'm not sure, if that's enough. The standard library functions may have different
implementations. In any case, I think there shouldn't be any algebraically equivalent different
formulas between versions of a software that uses such chaotic systems.

Idea: 
-start with a simple linear system with well understood and musically useful behavior, like the 
 (driven) damped harmonic oscillator: 
 m*a = -k*x - r*v + F 
 where m: mass, a: acceleration, k: spring hardness, x: excursion, r: damping/resistence, 
       v: velocity, F: driving force
-add nonlinear terms one by one and study their influence on the system, for example a progressive
 spring hardening term: -k3*x^3 - this should distort the waveform and increase the frequency at
 higher amplitudes
-generally, instead of using just the linear contributions (m*a, k*x, r*v), we could use 
 polynomials in x, v, a
-maybe we could also consider more derivatives like the increase/decrease of acceleration ("jerk"?)

*/
