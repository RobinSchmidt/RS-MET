#include "Experiments.h"

void particleBouncerExperiment()
{
  static const int N = 100;   // number of output samples

  // create and set up particle bouncer:
  ParticleBouncer bouncer;

  // create output sequence:
  double x[N], y[N];
  for(int n = 0; n < N; n++)
    bouncer.getSampleFrame(x[n], y[n]);



  // plot sequence:
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y);
}