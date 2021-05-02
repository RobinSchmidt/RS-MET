
bool samplerEngineUnitTest()
{
  bool ok = true;

  using TSig = rsFloat64x2;           // sampler output signal values
  using TPar = double;                // parameter values
  using TSmp = float;                 // sample values, as stored in RAM
  using VecS = std::vector<TSmp>;     // vector of sample values in RAM
  using SE   = rsSamplerEngine<TSig, TPar, TSmp>;
  SE se;                              // create the sampler engine object

  // Create a sinewave as example sample:
  double fs = 44100;  // sample rate
  double f  = 440.0;  // frequency of sinewave sample
  double a  = 0.5;    // amplitude of sinewave sample
  int    N  = 5000;   // length of sinewave sample
  VecS sample(N);
  double w = 2*PI*f/fs;
  for(int n = 0; n < N; n++)
    sample[n] = (float) sin(w*n);

  // Create an array of pointers to the channels and add the sample to the sample pool in the 
  // sampler engine:
  TSmp* pSmp[2];
  pSmp[0] = &sample[0];
  pSmp[1] = nullptr;
  se.addSampleToPool(pSmp, N, 1, fs, "Sine440Hz");

  // Add a region for the sinewave sample to the sampler engine:
  //se.addRegion


  // ToDo: 
  // -create a couple of simple samples (maybe sine-waves or something) and assign them to
  //  regions, trigger notes, record output and compare the produced output to what is expected
  // -set up performance tests, too
  // -fix the memory leak!


  return ok;
}