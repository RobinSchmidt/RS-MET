
bool samplerEngineUnitTest()
{
  bool ok = true;

  using TSig = rsFloat64x2;           // sampler output signal values
  using TPar = double;                // parameter values
  using TSmp = float;                 // sample values, as stored in RAM
  using VecS = std::vector<TSmp>;     // vector of sample values in RAM
  using SE   = rsSamplerEngine<TSig, TPar, TSmp>;
  using RC   = SE::ReturnCode;
  SE se;                              // create the sampler engine object

  // Create a sinewave as example sample:
  double fs = 44100;  // sample rate
  double f  = 440.0;  // frequency of sinewave sample
  double a  = 0.5;    // amplitude of sinewave sample
  int    N  = 500;    // length of sinewave sample
  VecS sample(N);
  double w = 2*PI*f/fs;
  for(int n = 0; n < N; n++)
    sample[n] = (TSmp) sin(w*n);
  //rsPlotVector(sample);

  // Create an array of pointers to the channels and add the sample to the sample pool in the 
  // sampler engine:
  TSmp* pSmp[2];
  pSmp[0] = &sample[0];
  pSmp[1] = nullptr;
  se.addSampleToPool(pSmp, N, 1, fs, "Sine440Hz");
  // should return the sample-index in the sample-pool si

  // Add a region for the sinewave sample to the sampler engine:
  int gi = se.addGroup();     // add new group to instrument definition, gi: group index
  ok &= gi == 0;
  int ri = se.addRegion(gi);  // add new region to group gi, ri: region index
  ok &= ri == 0;
  // addRegion should optionally take a loKey and hiKey parameter because it may actually
  // add entries to the regionsForKey array and if we set up loKey, hiKey later, they may be added
  // just to be removed again shortly thereafter which is wasteful. Of course, it should be 
  // possible to set it up later, but if it can start with the right setting from the beginning, we
  // save a lot of memory operations.

  //int rc = se.setRegionSample(gi, ri, si); 
  //ok &= rc == RC::success;

  // ToDo: 
  // -create a couple of simple samples (maybe sine-waves or something) and assign them to
  //  regions, trigger notes, record output and compare the produced output to what is expected
  // -set up performance tests, too


  rsAssert(ok);
  return ok;
}