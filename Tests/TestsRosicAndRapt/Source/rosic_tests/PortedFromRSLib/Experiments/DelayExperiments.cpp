#include "DelayExperiments.h"

void basicIntegerDelayLine()
{
  static const int N = 20;
  double t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsBasicDelayLineD dl;
  dl.setDelayInSamples(5);
  RAPT::getImpulseResponse(dl, h, N);
  plotData(N, t, h);
}
