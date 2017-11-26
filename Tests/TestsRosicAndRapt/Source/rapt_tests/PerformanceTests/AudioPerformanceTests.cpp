#include "AudioPerformanceTests.h"
using namespace RAPT;

void ladderPerformance()
{
  int numSamples = 50000;
  typedef double Real;
  vector<Real> x = createNoise(numSamples, 0, Real(-1), Real(+1));

  RAPT::rsLadderFilter<Real, Real> lf;
  lf.setCutoff(1000);
  lf.setResonance(0.5);

  ProcessorCycleCounter counter;
  counter.init();
  for(int n = 0; n < numSamples; n++)
  {
    x[n] = lf.getSample(x[n]);
  }
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("RAPT::rsLadderFilter", cycles / (numSamples));

  // Results:
  // double: with 1 / (1 + x^2) nonlinearity: 88 cycles, linear: 58 cycles, softClipHexic: 102


}