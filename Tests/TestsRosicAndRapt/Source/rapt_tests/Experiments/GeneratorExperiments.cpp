#include "GeneratorExperiments.h"

void getTwoParticleTrajectories(rsParticleSystemF& ps, int N, float* x1, float* y1, float* z1,
  float* x2, float* y2, float* z2)
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

    ps.updateState();
  }
}
void particleSystem()
{
  // We simulate a simple system of two particles with unit mass and unit charge to see, if they
  // behave as physically expected (i.e. to see, if the force euqations ae plausible)

  static const int N = 500; // number of steps in the simulations

  // create and set up the particle system:
  rsParticleSystemF ps(2);

  // both particles have unit mass and unit charge:
  ps.particles[0].mass   = 1.f;
  ps.particles[0].charge = 1.f;
  ps.particles[1].mass   = 1.f;
  ps.particles[1].charge = 1.f;

  // place them at (-1,0,0) and (+1,0,0) with zero velocity initially:
  ps.initialPositions[0]  = rsVector3DF(-1, 0, 0);
  ps.initialPositions[1]  = rsVector3DF(+1, 0, 0);
  ps.initialVelocities[0] = rsVector3DF( 0, 0, 0);
  ps.initialVelocities[1] = rsVector3DF( 0, 0, 0);

  // in a first run, we only let them interact via gravitation - they should attract each other:
  ps.setGravitationalConstant(1);
  ps.setElectricConstant(0);
  ps.setMagneticConstant(0);
  ps.setStepSize(0.01);

  // record trajectories:
  float x1[N], y1[N], z1[N], x2[N], y2[N], z2[N]; // maybe record kinetic and potential energy

  getTwoParticleTrajectories(ps, N, x1, y1, z1, x2, y2, z2);


  // they initially approach each other, fly through each other and then drift apart forever
  // i suppose, this is due to the singularity, when they are very close (division by zero
  // -> infinite force)

  int dummy = 0;
}

void rayBubble()
{
  rsRayBouncerF rb;
  rb.setEllipseParameters(1, 2, float(PI/4), 0.1f, 0.2f);
  rb.setInitialPosition(0.2f, -0.4f);
  rb.setSpeed(0.02f);
  rb.setLaunchAngle(float(PI) * 30.f / 180.f);

  // create output sequence:
  static const int N = 8000;   // number of output samples
  float x[N], y[N];
  rb.reset();
  for(int n = 0; n < N; n++)
    rb.getSampleFrame(x[n], y[n]);

  GNUPlotter plt;
  plt.setRange(-1.5, +1.5, -1.5, +1.5);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");           // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(N, x, y);
  //plt.addDataArrays(Ne, xe, ye);
  plt.plot();
}
