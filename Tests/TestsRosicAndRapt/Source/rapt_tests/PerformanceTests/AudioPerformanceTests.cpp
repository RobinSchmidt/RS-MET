#include "AudioPerformanceTests.h"
using namespace RAPT;

void ladderPerformance()
{
  int numSamples = 20000;
  ProcessorCycleCounter counter;
  int n;

  // create input and allocate output signals:
  vector<double> xs = createNoise(numSamples, 0, double(-1), double(+1));
  vector<double> ys(numSamples);
  vector<rsFloat64x2> xv(numSamples), yv(numSamples);
  for(n = 0; n < numSamples; n++) 
    xv[n] = xs[n]; 

  // scalar-signal, scalar-coeffs:
  RAPT::rsLadderFilter<double, double> filterSS;
  filterSS.setCutoff(1000);
  filterSS.setResonance(0.5);
  counter.init(); 
  for(n = 0; n < numSamples; n++) ys[n] = filterSS.getSample(xs[n]);
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsLadderFilter<double, double>", cycles/numSamples);

  // vector-signal, vector-coeffs:
  RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2> filterVV;
  filterVV.setCutoff(1000);
  filterVV.setResonance(0.5);
  counter.init(); 
  for(n = 0; n < numSamples; n++) yv[n] = filterVV.getSample(xv[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsLadderFilter<rsFloat64x2, rsFloat64x2>", cycles/numSamples);

  int dummy = 0;

  // Results:
  // double: with 1 / (1 + x^2) nonlinearity: 88 cycles, linear: 58 cycles, softClipHexic: 102
  // for float, it seems to be about the same
  // linear: 60 (scalar), 70 (vector)
  // clip:   70-75 cycles (scalar and vector)
  // div:    90 (scalar), 105 (vector)
}

void engineersFilterPerformance()
{
  int numSamples = 20000;
  int order      = 10;      
  ProcessorCycleCounter counter;
  double cycles;
  int n;

  // create input and allocate output signals:
  vector<double> xs = createNoise(numSamples, 0, double(-1), double(+1));

  rosic::rsEngineersFilterOld filterScalar;
  filterScalar.setPrototypeOrder(order);
  counter.init(); 
  for(n = 0; n < numSamples; n++) filterScalar.getSampleFrameDirect2(&xs[n], &xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterOld", cycles/numSamples);

  xs = createNoise(numSamples, 0, double(-1), double(+1));

  rosic::rsEngineersFilterStereo filterVector;
  filterVector.setPrototypeOrder(order);
  counter.init(); 
  for(n = 0; n < numSamples; n++) filterVector.getSampleFrameStereo(&xs[n], &xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterStereo", cycles/numSamples);
}