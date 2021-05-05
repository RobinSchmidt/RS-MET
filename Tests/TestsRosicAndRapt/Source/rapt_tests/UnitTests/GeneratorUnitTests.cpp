
bool samplerEngineUnitTest()
{
  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rsSamplerEngineTest;
  using RC   = SE::ReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  using Ev   = rsMusicalEvent<float>;
  using EvTp = Ev::Type;


  int maxLayers = 8;  
  // Maximum number of layers, actually, according to the sfz spec, this is supposed to be 
  // theoretically infinite, but i guess, in practice, there has to be some (high) limit which is
  // never supposed to be reached

  SE se(maxLayers);                     // create the sampler engine object

  // Create a sinewave as example sample:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of sinewave sample
  float a  = 0.5;    // amplitude of sinewave sample
  int   N  = 500;    // length of sinewave sample
  VecF sample(N);
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
    sample[n] = sinf(w*n);
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

  SE::Region* r = se.getRegion(gi, ri);
  int rc = se.setRegionSample(gi, ri, si); 
  ok &= rc == RC::success;
  se.setRegionSetting(r, PST::PitchKeyCenter, 69.f); // key = 69 is A4 at 440 Hz

  // Now the engine is set up with the sinewave sample in a single region that spans the whole
  // keyboard. We trigger a note at the original PitchCenterKey with maximum velocity and let the
  // engine produce N samples. We expect that the output exactly matches the original sample:
  VecF outL(N), outR(N);
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumIdleLayers()   == maxLayers-1;
  ok &= se.getNumActiveLayers() == 1;
  for(int n = 0; n < N; n++)
  {
    float tmp[2];
    se.processFrame(tmp);  // this api sucks! change it!
    outL[n] = tmp[0];
    outR[n] = tmp[1];
  }




  // todo:
  //const SE::Group* group1 = se.getGroup(0);
  //const SE::Group* group2 = region.getGroup();
  //ok &= group1 == group2;

  // todo: take the settings member of rsSamplerEngine into an instrument call inside
  // rsInstrumentDataSFZ..or maybe don't create an "Instrument" nested class





  // ToDo: 
  // -create a couple of simple samples (maybe sine-waves or something) and assign them to
  //  regions, trigger notes, record output and compare the produced output to what is expected
  // -set up performance tests, too



  int regionPlayerSize = SE::getRegionPlayerSize();
  // 64 without filters and eq, 512 with - move this into some performance test function

  rsAssert(ok);
  return ok;
}