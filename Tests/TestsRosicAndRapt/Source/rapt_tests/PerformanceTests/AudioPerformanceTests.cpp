#include "AudioPerformanceTests.h"
using namespace RAPT;

void ladderPerformance()
{
  int numSamples = 50000;
  //typedef double Real;
  //typedef float Real;
  ProcessorCycleCounter counter;

  // scalar-signal, scalar-coeffs:
  vector<double> xs = createNoise(numSamples, 0, double(-1), double(+1));
  vector<double> ys(numSamples);
  RAPT::rsLadderFilter<double, double> filterSS;
  filterSS.setCutoff(1000);
  filterSS.setResonance(0.5);
  counter.init(); 
  for(int n = 0; n < numSamples; n++) ys[n] = filterSS.getSample(xs[n]);
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsLadderFilter<double double>", cycles / numSamples);

  // vector-signal, scalar-coeffs:



  // Results:
  // double: with 1 / (1 + x^2) nonlinearity: 88 cycles, linear: 58 cycles, softClipHexic: 102
  // for float, it seems to be about the same


}