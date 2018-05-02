#include "DelayExperiments.h"

void basicIntegerDelayLine()
{
  static const int N = 20;
  double t[N], h[N];
  rsFillWithIndex(t, N);
  rsBasicDelayLine dl;
  dl.setDelayInSamples(5);
  getImpulseResponse(dl, h, N);
  plotData(N, t, h);
}
