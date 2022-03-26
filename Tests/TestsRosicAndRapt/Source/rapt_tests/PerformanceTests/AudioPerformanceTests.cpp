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
  using Ev    = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp  = Ev::Type;



  int N = 5000;             // number of samples to produce for the test in each run
  int numRuns = 100;        // number of test runs

  Vec outL(N), outR(N);
  PerformanceCounterTSC counter;


  // Create and set up sampler engine. We use a single region with a single-cycle sinewave:
  SE se;
  se.setSampleRate(44100.f);
  se.preAllocateDspMemory();   // try to avoid the need to call this!


  // Helper functions:

  // Measures the number of cycles that it takes the sampler engine to produce a single sample 
  // frame by producing N sample frames (while writing them into the output buffers).
  auto measureCyclesPerFrame = [&]()
  {
    counter.init(); 
    for(int n = 0; n < N; n++)
      se.processFrame(&outL[n], &outR[n]);
    return (double) counter.getNumCyclesSinceInit() / (double) N;
    // ToDo: 
    // -Maybe don't write the produced samples to the outL, outR buffers but into dummy 
    //  variables instead? This would allow us to take the number N as function parameter. It may 
    //  be a less realistic scenario, though. 
  };

  // Using the measureCyclesPerFrame defined above, it collects data of running this function a 
  // given number of times. This data may be visualized or analyzed statistically in order to draw 
  // any conclusions about the performance. A single run of measureCyclesPerFrame is not conclusive
  // because the result of a single run may vary wildly. Multiple runs and subsequent analysis of 
  // the measurement data is necessary.
  auto collectCyclesPerFrameData = [&](int numDataPoints)
  {
    std::vector<double> data(numDataPoints);
    for(int i = 0; i < numDataPoints; i++)
      data[i] = measureCyclesPerFrame();
    return data;
  };


  auto testSingleNote = [&](const std::string& testName, int key = 60, int vel = 100)
  {
    se.reset();
    se.handleMusicalEvent(Ev(EvTp::noteOn, key, vel));

    // new:
    std::vector<double> data = collectCyclesPerFrameData(numRuns);
    rsPlotVector(data);

    /*
    // old:
    double cycles = measureCyclesPerFrame();
    std::string str = "Key=" + to_string(key); // maybe print also vel
    printPerformanceTestResult(str, cycles);
    */
  };
  // -Instead of printing a result of a single run, use collectCyclesPerFrameData and show a plot.
  //  Maybe write a reusable visualizePerformanceData function taking a vector of cyclesPerOp data
  //  creates a plot that is suitable and meaningful. Maybe plot a histogram or a kernel density 
  //  estimation.
  // -maybe rename to testSingleKey or measureSingleKey

  // Triggers "numNotes" notes starting at the given lowest notes where each note is one semitone 
  // higher than the previous.
  auto testMultiNotes = [&](const std::string& testName, int numNotes = 10, int lowest = 50)
  {
    se.reset();
    for(int i = 0; i < numNotes; i++)
      se.handleMusicalEvent(Ev(EvTp::noteOn, lowest + i, 100));  // trigger the notes

    // new:
    std::vector<double> data = collectCyclesPerFrameData(numRuns) / double(numNotes);
    rsPlotVector(data);

    /*
    // old:
    double cycles = measureCyclesPerFrame() / numNotes;          // per sample and note
    std::string str = "Num=" + to_string(numNotes);              // number of keys
    printPerformanceTestResult(str, cycles);
    */
  };
  // maybe rename to testManyKeys or measureManyKeys, maybe take the event handling out of the 
  // measurement, i.e. drag it to before counter.init. But then we should probably do the same 
  // thing in the single key test - getSamplerNote also contains the event handling

  // Calls both testSingleNote and testMultiNotes:
  auto playTests = [&](const std::string& testName, int loKey = 60, int numKeys = 10)
  {
    std::cout << testName << ":\n";
    testSingleNote(testName, loKey);
    testMultiNotes(testName, numKeys, loKey);
    std::cout << "\n\n";
  };





  // Play the empty patch to figure out CPU load in idle state:
  playTests("Empty");     // 25 / 2

  // Play just one layer of the looped single cycle sample:
  setupForSineWave(&se, 2048);
  playTests("1 region");  // 150 / 125
  //rsPlotVectors(outL, outR);  // just to sanity check the output

  // Modulate the DC parameter of a waveshape with an LFO:
  se.setRegionSetting(   0, 0, OC::distortN_dc, 0.f, 1);
  se.setRegionSetting(   0, 0, OC::lfoN_freq, 200.f, 1);
  se.setRegionModulation(0, 0, OT::FreeLfo, 1, OC::distortN_dc, 1, 0.2f, Mode::absolute);
  playTests("1 region, 1 LFO to DC");    // 330 / 235
  //rsPlotVectors(outL, outR);

  // Modulate the DC parameter by a second LFO:
  se.setRegionSetting(   0, 0, OC::lfoN_freq, 300.f, 2);
  se.setRegionModulation(0, 0, OT::FreeLfo, 2, OC::distortN_dc, 1, 0.1f, Mode::absolute);
  playTests("1 region, 2 LFOs to DC");  // 405 / 350
  //rsPlotVectors(outL, outR); // does the waveshape look right? use high key to see shape better

  // Modulate the DC parameter by a third LFO:
  se.setRegionSetting(   0, 0, OC::lfoN_freq, 400.f, 3);
  se.setRegionModulation(0, 0, OT::FreeLfo, 3, OC::distortN_dc, 1, 0.05f, Mode::absolute);
  playTests("1 region, 3 LFOs to DC");  // 585 / 510
  //rsPlotVectors(outL, outR); 

  // Observations:
  // -It's actually interesting to plot the collected data as time-series. There are patterns in it
  //  which are more easily visible with the simpler patches.
  // -The single note tests have clearly visible flat portions in the data when plotted as time 
  //  series. The multi-note tests are more erratic
  // -Adding another sine LFO to modulate DC seems to increase the per sample cost by roughly
  //  130 cycles. That's for the additional modulation infrastructure and the LFO's signal
  //  processing. The first LFO is more expensive because it triggers the one-time cost of invoking
  //  the modulation infrastructure
  // -When multiple keys are playing, the cost per sample per key seems to go down. We seem to get
  //  a sort of quantity rebate. This is good news!

  // ToDo:
  // -Try to not use the outL/outR arrays and investigate how (or if) that changes the patterns in
  //  the collected data.
  // -For better consistency, maybe run each test M times and take the minimum but with some 
  //  constraints to sort out mismeasurements (sometime we even get negative numbers)
  // -Implement some data collection and visualize and analyze the data.
  // -Add measurements for the cost of starting a new RegionPlayer on noteOn
  // -Measure costs for handling midi-events
  // -Implement block-based processing and measure its cost. It should hopefully be a lot cheaper
  //  than sample-by-sample processing.
  // -Implement EGs and make sure that they don't incur the coeff-recomputation costs when they
  //  output a constant level. That may even be a unit test (but it may be a flaky one).
  // -Maybe collect batches of data, say 1000 datapoints each, sort results, throw away bottom and
  //  top quartiles (outliers) and show visualizations of the rest. then appaly statistical 
  //  analysis on this rest.
  // -Maybe plot results in histograms or using kernel density estimators.
  //  https://en.wikipedia.org/wiki/Kernel_density_estimation
  //  https://de.wikipedia.org/wiki/Kerndichtesch%C3%A4tzer
  // -Maybe just plot the sorted results (maybe cleaned up from outliers before). Maybe we can 
  //  identify a region which has very low slope - if so, then the height there could be a good 
  //  estimate for the value we are trying to measure. Rationale: where the plot has little slope,
  //  we had a lot of very similar results.
  // -For kernel density plots, we can use non-unform filters: 
  //  -sort the array of results
  //  -use this sorted array as t-values
  //  -use as x-array all ones
  //  -apply non-uniform filtering to these t,x arrays (although t does not represent time in this
  //   case but rather cpu-load and x represents the number of runs that had this kind of load)




  int dummy = 0;
}