
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
  VecF sin440(N);    // sine wave
  VecF cos440(N);    // cosine wave 
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
  {
    sin440[n] = sinf(w*n);
    cos440[n] = cosf(w*n);
  }
  //rsPlotVector(sample);

  // Create an array of pointers to the channels and add the sample to the sample pool in the 
  // sampler engine:
  float* pSmp[2];
  pSmp[0] = &sin440[0];
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
  // engine produce N samples. We expect that the output in both channels exactly matches the 
  // original sample:
  VecF outL(N), outR(N);
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumIdleLayers()   == maxLayers-1;
  ok &= se.getNumActiveLayers() == 1;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(sin440, outL, outR);
  ok &= outL == sin440 && outR == sin440;
  //se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 0.f));  // is interpreted as note-off
  int i = se.stopAllPlayers();
  ok &= i == 1;
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  // hmm - a note-off should actually just trigger entering the release state. But maybe by 
  // default, that should indeed just cut off the note immediately. Only when an amp-env with 
  // nonzero release is used, the player will remain active for a while after note-off.
  // ok - we are now actually calling stopAllPlayers which is a hard reset and should indeed stop
  // all playback immediately. Let's see, if the output is zero:

  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= rsIsAllZeros(outL);
  ok &= rsIsAllZeros(outR);



  // todo:
  // -set pan of the first region to hard left
  // -add a second sample (maybe a cosine wave of the same frequency)
  // -add a region for the cosine wave
  // -pan the cosine to hard right
  // -check stereo output - left should be the sine, right the cosine
  pSmp[0] = &cos440[0];
  si = se.addSampleToPool(pSmp, N, 1, fs, "Cosine440Hz"); // add cosine sample to pool
  ok &= si == 1; 
  ri = se.addRegion(gi);                                  // add new region for cosine
  ok &= ri == 1;
  r = se.getRegion(gi, ri);                               // retrieve pointer to region
  rc = se.setRegionSample(gi, ri, si);                    // set region sample to cosine
  ok &= rc == RC::success;
  se.setRegionSetting(r, PST::PitchKeyCenter, 69.f);      // cosine has same root key as sine
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumIdleLayers()   == maxLayers-2;
  ok &= se.getNumActiveLayers() == 2;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR);
  ok &= outL == sin440 + cos440;
  ok &= outR == sin440 + cos440;

  // After having read both samples until the end, trying to produce more output thereafter should
  // produce all zeros, even if we don't stop all players. The engine should detect that the end of
  // the sample was reached and stop the players automatically...




  i = se.stopAllPlayers();
  ok &= i == 2;
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  r = se.getRegion(gi, 0);
  se.setRegionSetting(r, PST::Pan, -100.f);               // pan sine to hard left
  r = se.getRegion(gi, 1);
  se.setRegionSetting(r, PST::Pan, +100.f);               // pan cosine to hard right
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR); 
  ok &= outL == (2.f * sin440);
  ok &= outR == (2.f * cos440);
  // ToDo: Check, how the sfz player handles the amp/pan parameters with respect to toal gain. 
  // Should there be a factor of 2 for hard left/right settings or a factor of 0.5 for a center 
  // setting? Make sure to match the behavior of the reference player. I actually would tend to 
  // prefer the former because it implies that with the neutral default settings, samples are 
  // played back as is as opposed to having acquired a gain of 0.5. Maybe if sfz behaves the other
  // way, we could provide both options by introducing additional pan rules.
   





  // this looks like the sum of sine an cosine because pan is not yet implemented - but that is
  // a good test to keep too - just do it before setting the pan values





  // -clear the regions
  // -add sine and cosine samples again, but this time as a single stereo sample
  // -create a region for that and test the output



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