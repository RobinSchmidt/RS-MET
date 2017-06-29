#include "Experiments.h"

void particleBouncerExperiment()
{
  static const int N = 6000;   // number of output samples

  // create and set up particle bouncer:
  ParticleBouncer bouncer;
  bouncer.setEnclosureEllipseAspectRatio(1.15);
  bouncer.setInitialPosition(0.2, -0.4);
  //bouncer.setInitialPosition(-0.41, 0.8);
  bouncer.setSpeed(0.02);
  bouncer.setLaunchAngle(30.0);

  // create output sequence:
  double x[N], y[N];
  bouncer.reset();
  for(int n = 0; n < N; n++)
    bouncer.getSampleFrame(x[n], y[n]);

  // create enclosing ellipse for reference inside the plot:
  static const int Ne = 100;
  double a = bouncer.getEllipseA();
  double b = bouncer.getEllipseB();
  double xe[Ne], ye[Ne];
  for(int n = 0; n < Ne; n++)
  {
    double t = 2 * PI * double(n) / (Ne-1);
    xe[n] = a * cos(t);
    ye[n] = b * sin(t);
  }
  double areaOverPi = a*b; // actual area is PI*a*b

  // plot sequence:
  GNUPlotter plt;
  plt.setRange(-1.2, +1.2, -1.2, +1.2);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");           // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(Ne, xe, ye);
  plt.plot();
  //plt.plotFunctionTables(N, x, y);
}