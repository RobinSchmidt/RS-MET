
bool samplerEngineUnitTest()
{
  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rsSamplerEngineTest;
  using RC   = SE::ReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  SE se;                              // create the sampler engine object

  // Create a sinewave as example sample:
  double fs = 44100;  // sample rate
  double f  = 440.0;  // frequency of sinewave sample
  double a  = 0.5;    // amplitude of sinewave sample
  int    N  = 500;    // length of sinewave sample
  VecF sample(N);
  double w = 2*PI*f/fs;
  for(int n = 0; n < N; n++)
    sample[n] = (float) sin(w*n);
  //rsPlotVector(sample);

  // Create an array of pointers to the channels and add the sample to the sample pool in the 
  // sampler engine:
  float* pSmp[2];
  pSmp[0] = &sample[0];
  pSmp[1] = nullptr;
  int si = se.addSampleToPool(pSmp, N, 1, fs, "Sine440Hz");
  ok &= si == 0; // should return the sample-index in the sample-pool

  // Add a region for the sinewave sample to the sampler engine:
  int gi = se.addGroup();     // add new group to instrument definition, gi: group index
  ok &= gi == 0;
  int ri = se.addRegion(gi);  // add new region to group gi, ri: region index
  ok &= ri == 0;

  const SE::Region* r = se.getRegion(gi, ri);
  int rc = se.setRegionSample(gi, ri, si); 
  ok &= rc == RC::success;
  se.setRegionSetting(r, PST::PitchKeyCenter, 69.f);


  // Maybe it should be a non-const pointer, so we can actually edit it - or let the edits be done
  // by the sampler-engine, like se.setRegionSetting(PS::RootKey, 69). This way, the SamplerEngine
  // object remains in control over anything that is done to the region, which is a good thing



  // todo:
  //const SE::Group* group1 = se.getGroup(0);
  //const SE::Group* group2 = region.getGroup();
  //ok &= group1 == group2;






  // ToDo: 
  // -create a couple of simple samples (maybe sine-waves or something) and assign them to
  //  regions, trigger notes, record output and compare the produced output to what is expected
  // -set up performance tests, too



  int regionPlayerSize = SE::getRegionPlayerSize();
  // 64 without filters and eq, 512 with - move this into some performance test function

  rsAssert(ok);
  return ok;
}