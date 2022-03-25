using namespace RAPT;
using namespace std;

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
  std::vector<double> xs = createNoise(numSamples, double(-1), double(+1));
  std::vector<double> ys(numSamples);
  std::vector<rsFloat64x2> xv(numSamples), yv(numSamples);
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
  RAPT::rsLadderFilter<rsFloat64x2, double> filterVS;
  filterVS.setCutoff(1000);
  filterVS.setResonance(0.5);
  counter.init(); 
  for(n = 0; n < numSamples; n++) yv[n] = filterVS.getSample(xv[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsLadderFilter<rsFloat64x2, double>", cycles/numSamples);

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

// maybe generalize this function for use with any kind of filter - we just need to factor out
// the class-specific setup code:
template<class TSig, class TPar>
double getStateVectorFilterCyclesPerSample(TSig sig, TPar par, int numSamples, int numTests)
{
  rsStateVectorFilter<TSig, TPar> flt;
  flt.setupFromBiquad(TPar(1), TPar(2), TPar(0.5), TPar(-0.5), TPar(0.25));
  TSig y;
  double minCycles = DBL_MAX;

  PerformanceCounterTSC counter;
  for(int i = 1; i <= numTests; i++) {
    counter.init();
    for(int n = 0; n < numSamples; n++)
      y = flt.getSample(0);
    double cycles = (double)counter.getNumCyclesSinceInit();
    if(cycles > 0 && cycles < minCycles)
      minCycles = cycles;
  }
  dontOptimize(&y);

  return minCycles/numSamples;
}

void stateVectorFilterPerformance()
{
  int numSamples = 2000;
  int numTests = 10;
  double cycles;

  // double precision:
  cycles = getStateVectorFilterCyclesPerSample(1.0, 1.0, numSamples, numTests);
  printPerformanceTestResult("rsStateVectorFilter<double, double>", cycles);

  cycles = getStateVectorFilterCyclesPerSample(rsFloat64x2(1.0), 1.0, numSamples, numTests);
  printPerformanceTestResult("rsStateVectorFilter<rsFloat64x2, double>", cycles);

  // single precision:
  cycles = getStateVectorFilterCyclesPerSample(1.f, 1.f, numSamples, numTests);
  printPerformanceTestResult("rsStateVectorFilter<float, float>", cycles);

  cycles = getStateVectorFilterCyclesPerSample(rsFloat32x4(1.f), 1.f, numSamples, numTests);
  printPerformanceTestResult("rsStateVectorFilter<rsFloat32x4, float>", cycles);

  cycles = getStateVectorFilterCyclesPerSample(rsFloat32x4(1.f), rsFloat32x4(1.f), numSamples, numTests);
  printPerformanceTestResult("rsStateVectorFilter<rsFloat32x4, rsFloat32x4>", cycles);
}

void engineersFilterPerformance()
{
  int numSamples = 2000;
  int order      = 20;      
  PerformanceCounterTSC counter;
  double cycles;
  int n;

  // create input and allocate output signals:
  vector<double> xs = createNoise(numSamples, double(-1), double(+1));

  // scalar version (mono):
  rosic::rsEngineersFilterMono filterScalar;
  filterScalar.setPrototypeOrder(order);
  filterScalar.setMode(rsInfiniteImpulseResponseDesigner<float>::BANDPASS);
  counter.init(); 
  for(n = 0; n < numSamples; n++) 
    xs[n] = filterScalar.getSample(xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterMono", cycles/numSamples);
  
  xs = createNoise(numSamples, double(-1), double(+1));

  // vector version (stereo):
  rosic::rsEngineersFilterStereo filterVector;
  filterVector.setPrototypeOrder(order);
  filterVector.setMode(rsInfiniteImpulseResponseDesigner<float>::BANDPASS);
  counter.init(); 
  for(n = 0; n < numSamples; n++) 
    filterVector.getSampleFrameStereo(&xs[n], &xs[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsEngineersFilterStereo", cycles/numSamples);

  // Here, the performance gain is indeed a factor of 2 - it works exactly as it should
  // Both, scalar and vector version take about 550 cycles per sample with order=20 and 
  // bandpass setting (which doubles the actual order to 40 biquads). That makes 
  // 550/40 = 13.75 cycles per sample (pair) per biquad. 
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

void samplerEnginePerformance()
{
  using Vec   = std::vector<float>;
  using SE    = rosic::Sampler::rsSamplerEngineTest;
  using OC    = rosic::Sampler::Opcode;
  using OT    = rosic::Sampler::OpcodeType;
  using Shape = rosic::Sampler::WaveshaperCore::Shape;
  using Mode  = rosic::Sampler::ModMode;


  int N = 2000;      // number of samples to produce for the test

  // Create and set up sampler engine. We use a single region with a single-cycle sinewave and 
  // modulate the DC parameter of a waveshape with an LFO:
  SE se;
  se.setSampleRate(44100.f);
  setupForSineWave(&se, 2048);
  Vec outL(N), outR(N);

  PerformanceCounterTSC counter;
  double cycles;
  counter.init(); 
  getSamplerNote(&se, 60, 100, outL, outR);
  cycles = (double) counter.getNumCyclesSinceInit();
  double cyclesPerSample = cycles / N;
  rsPlotVectors(outL, outR);

  int dummy = 0;
}