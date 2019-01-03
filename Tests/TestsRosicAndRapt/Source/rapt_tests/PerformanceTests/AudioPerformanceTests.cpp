#include "AudioPerformanceTests.h"
using namespace RAPT;

void filterSignConventionPerformance()
{
  // test, how the + or - sign convention for filters affects the performance

  int numSamples = 10000;

  //typedef float Real;  
  typedef double Real;
  //typedef rsFloat64x2 Real;

  // filter coeffs:
  Real b0 = 0.5, b1 = 0.25, b2 =  0.125, a1 = 0.5, a2 = -0.25;
  Real x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  Real g;

  // signals and bookeeping:
  PerformanceCounterTSC counter;
  double cycles;
  vector<Real> x = createNoise(numSamples, double(-1), double(+1));
  vector<Real> y(numSamples);
  int n;

  // biquad, direct form 1, positive sign for a-coeffs:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    y[n] = b0*x[n] + b1*x1 + b2*x2 + a1*y1 + a2*y2;
    x2 = x1;
    x1 = x[n];
    y2 = y1;
    y1 = y[n];
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("biquad, DF1, +:", cycles/numSamples);

  // biquad, direct form 1, negative sign for a-coeffs:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    y[n] = b0*x[n] + b1*x1 + b2*x2 - a1*y1 - a2*y2;
    x2 = x1;
    x1 = x[n];
    y2 = y1;
    y1 = y[n];
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("biquad, DF1, -:", cycles/numSamples);


  // biquad, direct form 2, positive sign for a-coeffs:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    g    = x[n] + a1*y1 + a2*y2;
    y[n] = b0*g + b1*y1 + b2*y2;
    y2 = y1;
    y1 = g;
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("biquad, DF2, +:", cycles/numSamples);

  // biquad, direct form 2, negative sign for a-coeffs:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    g    = x[n] - a1*y1 - a2*y2;
    y[n] = b0*g + b1*y1 + b2*y2;
    y2 = y1;
    y1 = g;
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("biquad, DF2, -:", cycles/numSamples);



  // 1st order, DF1, +:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    y[n] = y1 = b0*x[n] + b1*x1 + a1*y1;
    x1 = x[n];
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("1st order, DF1, +:", cycles/numSamples);

  // 1st order, DF1, -:
  counter.init(); 
  for(n = 0; n < numSamples; n++)
  {
    y[n] = y1 = b0*x[n] + b1*x1 - a1*y1;
    x1 = x[n];
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("1st order, DF1, -:", cycles/numSamples);







  // todo: test 1st order filters, test use of parentheses, test one-pole (without zero)


  // Results:
  // biquad, DF1, +: 18
  // biquad, DF1, -: 15  ..so it's 20% slower (18/15 = 1.2)
}

void ladderPerformance()
{
  int numSamples = 20000;
  PerformanceCounterTSC counter;
  int n;

  // create input and allocate output signals:
  vector<double> xs = createNoise(numSamples, double(-1), double(+1));
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

  // vector signal, scalar coeffs...

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
  int numSamples = 10000;
  int order      = 20;      
  PerformanceCounterTSC counter;
  double cycles;
  int n;

  // create input and allocate output signals:
  vector<double> xs = createNoise(numSamples, double(-1), double(+1));

  /*
  // this is obsolete:
  rosic::rsEngineersFilterOld filterScalar;
  filterScalar.setPrototypeOrder(order);
  filterScalar.setMode(rsInfiniteImpulseResponseDesignerF::BANDPASS);
  counter.init(); 
  for(n = 0; n < numSamples; n++) 
    filterScalar.getSampleFrameDirect1(&xs[n], &xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterOld", cycles/numSamples);
  */

  xs = createNoise(numSamples, double(-1), double(+1));

  rosic::rsEngineersFilterStereo filterVector;
  filterVector.setPrototypeOrder(order);
  filterVector.setMode(rsInfiniteImpulseResponseDesignerF::BANDPASS);
  counter.init(); 
  for(n = 0; n < numSamples; n++) 
    filterVector.getSampleFrameStereo(&xs[n], &xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterStereo", cycles/numSamples);

  // the performance gain is just a factor of 1.4, not 2...maybe try to have the filter coeffs of
  // type rsFloar64x2 too?
}

void turtleGraphicsPerformance()
{
  //rosic::TurtleGraphics tg;
  rosic::LindenmayerRenderer lr;
  //int order = 5;
  std::vector<double> x, y;

  PerformanceCounterTSC counter;
  double cycles;

  counter.init(); 
  lr.getMooreCurve(4, x, y);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Moore Curve 4, 1st time", cycles); // around 500 000

  counter.init(); 
  lr.getMooreCurve(4, x, y);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Moore Curve 4, 2nd time", cycles); // around 350 000
}