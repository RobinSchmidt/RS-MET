#include "Experiments.h"

void particleBouncerExperiment()
{
  static const int N = 100;   // number of output samples

  // create and set up particle bouncer:
  ParticleBouncer bouncer;
  bouncer.setEnclosureEllipseAspectRatio(1.0);
  //bouncer.setInitialPosition(0.2, -0.4);
  bouncer.setInitialPosition(-0.41, 0.8);
  bouncer.setSpeed(0.2);
  bouncer.setAngle(0.0);
  //bouncer.setInitialIncrements(0.2, 0.0);


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
  plt.addCommand("set size square");           // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(Ne, xe, ye);
  plt.plot();
  //plt.plotFunctionTables(N, x, y);
}