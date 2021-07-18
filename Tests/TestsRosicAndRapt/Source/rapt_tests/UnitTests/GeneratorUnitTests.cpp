bool samplerDataUnitTest()
{
  bool ok = true;

  //using SD = rsSamplerData;
  using SD = rosic::rsSamplerData;

  SD d1;

  ok &= d1.getNumGroups() == 0;

  // todo: add some groups/regions and test copy-assignment we need deep copies and that doesn't 
  // seem to work yet

  int gi, ri;             // gi: group index, ri: region index

  // Add 2 groups:
  gi = d1.addGroup(); ok &= gi == 0;
  gi = d1.addGroup(); ok &= gi == 1;

  // Add2 regions to 1st group and 3 regions to 2nd group:
  ri = d1.addRegion(0); ok &= ri == 0;
  ri = d1.addRegion(0); ok &= ri == 1;
  ri = d1.addRegion(1); ok &= ri == 0;
  ri = d1.addRegion(1); ok &= ri == 1;
  ri = d1.addRegion(1); ok &= ri == 2;

  // Test assignment operator:
  SD d2 = d1; ok &= d2 == d1; // copy assign
  SD d3(d1);  ok &= d3 == d1; // copy construct
  ok &= d2 == d3;             // equality should be transitive



  // this triggers an assert - the assignment operator needs to be implemented differently
  // ..or better: avoid using pointers and new operator for adding regions
  // -in rsSamplerEngine, refer to regions not by a pointer to the region (in setRegionSetting 
  //  etc.), but instead by group index and region index
  // -oh - no - we need to use pointers because of rsSamplerEngine::regionsForKey
  // -actually, we need a pointer array for the groups also because the regions refer back to them
  // 

  rsAssert(ok);
  return ok;
}

bool testSamplerNote(rosic::rsSamplerEngine* se, float key, float vel, 
  const std::vector<float>& targetL, const std::vector<float>& targetR, 
  float tol = 0.f)
{
  using AT   = RAPT::rsArrayTools;
  using Ev   = rosic::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  int N = (int) targetL.size();
  rsAssert((int)targetR.size() == N);
  std::vector<float> outL(N), outR(N);

  se->handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
  for(int n = 0; n < N; n++)
    se->processFrame(&outL[n], &outR[n]);
  float errL = AT::maxDeviation(&outL[0], &targetL[0], N);
  float errR = AT::maxDeviation(&outR[0], &targetR[0], N);
  //rsPlotVectors(outL, outR); // uncomment for debugging
  return errL <= tol && errR <= tol;
};
// maybe have a bool resetBefore that optionally resets the engine before playing...but maybe it's
// better when the caller does the reset directly, if desired - it's not much longer but clearer

bool samplerEngineUnitTest1()
{
  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using AT   = RAPT::rsArrayTools;
  using SD   = rosic::rsSamplerData;
  using SE   = rosic::rsSamplerEngineTest;
  using RC   = rosic::rsReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  using Ev   = rosic::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  int maxLayers = 8;  
  // Maximum number of layers, actually, according to the sfz spec, this is supposed to be 
  // theoretically infinite, but i guess, in practice, there has to be some (high) limit which is
  // never supposed to be reached. However, for the unit tests, we use an unrealistically low value
  // to facilitate also testing the behavior when the limits are indeed reached.

  SE se(maxLayers);                     // create the sampler engine object

  // Create a sine- and cosine wave as example samples:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of (co)sinewave sample
  //float a  = 0.5;    // amplitude of (co)sinewave sample
  int   N  = 500;    // length of (co)sinewave sample
  VecF sin440(N);    // sine wave
  VecF cos440(N);    // cosine wave 
  VecF zeros(N);     // zero samples, silence
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
  {
    sin440[n] = sinf(w*n);
    cos440[n] = cosf(w*n);
    zeros[n]  = 0.f;
  }
  //rsPlotVectors(sin440, cos440);

  // Create an array of pointers to the channels and add the sample to the sample pool in the 
  // sampler engine:
  float* pSmp[2];
  pSmp[0] = &sin440[0];
  pSmp[1] = nullptr;
  int si = se.addSampleToPool(pSmp, N, 1, fs, "Sine440Hz");
  ok &= si == 0; // should return the sample-index in the sample-pool
  // Maybe instead of directly adding the sample-data to the sample-pool, write it to a wave file 
  // and instruct the engine to read the file, which should in turn add it. Or make and additional 
  // function loadSampleToPool(sampleDir, fileName)

  // Add a region for the sinewave sample to the sampler engine:
  int gi = se.addGroup();     // add new group to instrument definition, gi: group index
  ok &= gi == 0;
  int ri = se.addRegion(gi);  // add new region to group gi, ri: region index
  ok &= ri == 0;

  int rc = se.setRegionSample(gi, ri, si);           // rc: return code
  ok &= rc == RC::success;
  rc = se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f); // key = 69 is A4 at 440 Hz
  ok &= rc == RC::success;

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

  // After having played both samples until the end, trying to produce more output thereafter 
  // should produce all zeros, even if we don't stop all players. The engine should detect that 
  // the end of the sample was reached and stop the players automatically:
  ok &= se.getNumActiveLayers() == 0;            // player should have stopped
  ok &= se.getNumIdleLayers()   == maxLayers;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);         // output should be zero
  ok &= rsIsAllZeros(outL);
  ok &= rsIsAllZeros(outR);

  // Test with 2 regions:
  // -add a second sample to the pool (a cosine wave of the same frequency)
  // -add and set up a region for the cosine wave
  // -expected output: both, left and right should be the sum of sine and cosine
  pSmp[0] = &cos440[0];
  si = se.addSampleToPool(pSmp, N, 1, fs, "Cosine440Hz"); // add cosine sample to pool
  ok &= si == 1; 
  ri = se.addRegion(gi);                                  // add new region for cosine
  ok &= ri == 1;
  rc = se.setRegionSample(gi, ri, si);                    // set region sample to cosine
  ok &= rc == RC::success;
  se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f);   // cosine is also at A4 = 440 Hz
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
  ok &= se.getNumActiveLayers() == 0;                     // players should have stopped
  ok &= se.getNumIdleLayers()   == maxLayers;

  // Test panning:
  // -set pan of the first region to hard left
  // -pan the cosine to hard right
  // -expected output: left should be the twice sine, right twice the cosine (the factor two arises
  //  from the pan law)
  se.setRegionSetting(gi, 0, PST::Pan, -100.f);           // pan sine to hard left
  se.setRegionSetting(gi, 1, PST::Pan, +100.f);           // pan cosine to hard right
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR); 
  ok &= outL == (2.f * sin440);
  ok &= outR == (2.f * cos440);
  ok &= se.getNumActiveLayers() == 0;

  // Test handling of noteOff. At the moment, we have no amp envelope yet, so a noteOff should 
  // immediately stop the playback of all layers which it applies to.
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));  // the noteOn, again
  for(int n = 0; n < N/2; n++)                           // play through half the sample
    se.processFrame(&outL[n], &outR[n]);
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 0.f));    // vel=0 is interpreted as note-off
  for(int n = N/2; n < N; n++)                           // play through 2nd half, this should...
    se.processFrame(&outL[n], &outR[n]);                 // ...produce silence
  for(int n = 0; n < N/2; n++) {
    ok &= outL[n] == 2.f * sin440[n];
    ok &= outR[n] == 2.f * cos440[n]; }
  for(int n = N/2; n < N; n++) {
    ok &= outL[n] == 0.f;
    ok &= outR[n] == 0.f; }
  //rsPlotVectors(outL, outR);

  // Test realtime downsampling. Play a note an octave above the root key:
  se.handleMusicalEvent(Ev(EvTp::noteOn, 81.f, 127.f));  // noteOn, 1 octave above root key
  for(int n = 0; n < N/2; n++)                           // play the sample at double speed which
    se.processFrame(&outL[n], &outR[n]);                 // makes it half as long
  //rsPlotVectors(outL, outR);
  ok &= se.getNumActiveLayers() == 0;
  for(int n = 0; n < N/2; n++) {
    ok &= outL[n] == 2.f * sin440[2*n];
    ok &= outR[n] == 2.f * cos440[2*n]; }

  // Test realtime upsampling. Play a note an octave below the root key:
  se.handleMusicalEvent(Ev(EvTp::noteOn, 57.f, 127.f));
  for(int n = 0; n < N; n++)                  // play the sample at half speed
    se.processFrame(&outL[n], &outR[n]);      // this loop goes only through half of it
  ok &= se.getNumActiveLayers() == 2;         // ..so the 2 layers shall remain active
  //rsPlotVectors(sin440, cos440, outL, outR);
  for(int n = 0; n < N/2; n++) {
    ok &= outL[2*n] == 2.f * sin440[n];
    ok &= outR[2*n] == 2.f * cos440[n]; }
  for(int n = 0; n < N; n++)                  // ...now go through the second half
    se.processFrame(&outL[n], &outR[n]);
  ok &= se.getNumActiveLayers() == 0;         // now, the 2 layers should be off
  for(int n = 0; n < N/2; n++) {
    ok &= outL[2*n] == 2.f * sin440[n+N/2];
    ok &= outR[2*n] == 2.f * cos440[n+N/2]; }
  // Other tests to do: set the root-key differently, set the sample-rates for playback and 
  // audiofile differently, test detuning opcodes

  const SD& sfzData = se.getInstrumentData();
  std::string sfzString = sfzData.getAsSFZ();
  SD sfzData2;
  sfzData2.setFromSFZ(sfzString);
  ok &= sfzData2 == sfzData;
  //SE se2;
  //se2.setupFromSFZ(sfzData2);
  //ok &= se2.matchesInstrumentDefinition(se);

  // Remove the 1st region (sin440) -> cos440 should become 1st region, then re-assign the sample
  // for 1st region to sine. Do the removal in the middle of playback - desired bevavior: it still 
  // plays to the end, but the next time the key is triggered, it doesn't start playing anymore. 
  // This works because in RegionPlayer::getFrame/processBlock, the region pointer is not 
  // referenced anymore. It's referenced only when the playback starts...but what about noteOff?
  // At the moment, it also doesn't reference the region pointer, but later we will want to trigger 
  // the release...maybe we can do that without referencing the region pointer, by setting up the 
  // envelopes already at the start...but what if a release sample should get triggered?
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N/2; n++)
    se.processFrame(&outL[n], &outR[n]);
  se.removeRegion(0, 0);
  for(int n = N/2; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR); 
  ok &= outL == (2.f * sin440);
  ok &= outR == (2.f * cos440);
  ok &= se.getNumActiveLayers() == 0;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR); 
  ok &= outL == zeros;
  ok &= outR == (2.f * cos440);
  se.setRegionSample(0, 0, 0);
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= outL == zeros;
  ok &= outR == (2.f * sin440);
  //rsPlotVectors(outL, outR); 


  // Computes a linearly interpolated value from the gievn vector at the given position:
  auto getSampleAt = [](const std::vector<float>& v, float pos)
  {
    if(pos < 0.f) return 0.f;

    int   i = (int) pos;
    float f = pos - (float)i;

    float x0 = 0.f, x1 = 0.f;
    int N = (int)v.size();
    if(i   < N) x0 = v[i];
    if(i+1 < N) x1 = v[i+1];

    return (1.f - f) * x0 + f * x1;
  };
  // move up or maybe move to library as rsInterpolateAt or rsLerpAt
  // hmmm...we assume here that between v[-1] and v[0], there is silence and between v[N-1] and 
  // v[N] we ramp linearly down to zero...maybe we should make that consistent - either ramp up and
  // ramp down or don't ramp at all - but the latter requires special handling of pos == N-1: if
  // it's exactly N-1, we should return the last sample, if it's slightly (e.g. 0.0001) above, 
  // return 0...that seems complicated. ramping up at the start seems more consistent with how the
  // linear interpolator filter kernel looks like...should the sampler engine also behave this way,
  // i.e. when we play a sample that contains a single impulse and play it back at 1/10 of the 
  // normal speed, it would ramp up and down over 10 samples? ...but that's not feasible because it
  // would require the region to start producing output before it was triggered...or wait..maybe 
  // not. It could be feasible, if we replace  if(sampleTime < 0.0)  by  if(sampleTime < 1.0) in
  // rsSamplerEngine::RegionPlayer::getFrame and implement the  stream->getFrameStereo  call in the
  // same way as above. ...but it will work only if delay is used...but such minute details are 
  // perhaps not very important anyway

  // Test delay:
  float delaySamples = 10.75f;
  float delaySeconds = delaySamples / fs;
  se.setRegionSetting(0, 0, PST::Pan, 0.f);              // back to center, makes testing easier
  se.setRegionSetting(0, 0, PST::Delay, delaySeconds);
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));  // the noteOn, again
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  for(int n = 0; n < delaySamples; n++) {
    ok &= outL[n] == 0.f;
    ok &= outR[n] == 0.f;  }
  float tol = 1.e-7f;  // ~= 140 dB SNR
  for(int n = (int) delaySamples; n < N; n++) 
  {
    float tgt = getSampleAt(sin440, n-delaySamples);
    ok &= rsIsCloseTo(outL[n], tgt, tol);
    ok &= rsIsCloseTo(outR[n], tgt, tol);
  }
  //rsPlotVectors(sin440, outL);


  // move this into samplerEngine2UnitTest
  // Test setGroupSettingsOnTop: We set the volume of the group and the region and check the 
  // behavior in both modes:
  float regionAmp = 0.5f;
  float groupAmp  = 0.25f;
  //float instrumentAmplitude  = 0.75f;
  se.setRegionSetting(0, 0, PST::Delay, 0.f);  // Turn delay off again
  se.setRegionSetting(0, 0, PST::Volume, rsAmpToDb(regionAmp));
  se.setGroupSetting( 0,    PST::Volume, rsAmpToDb(groupAmp));


  auto testNote = [&](
    float key, float vel, const VecF& targetL, const VecF& targetR, float tol = 0.f)
  {
    se.handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
    for(int n = 0; n < N; n++)
      se.processFrame(&outL[n], &outR[n]);
    float errL = AT::maxDeviation(&outL[0], &targetL[0], N);
    float errR = AT::maxDeviation(&outR[0], &targetR[0], N);
    return errL <= tol && errR <= tol;
  };
  // ToDo: move up and use it to reduce boilerplate for many other tests as well - maybe make it
  // a free function, taking the engine as reference argument...or a pointer
  // done: testSamplerNote ...use that function in the tests above to reduce the boilerplate


  /*
  // todo - goes into another test, maybe samplerEngineRoutingUnitTest:
  se.setGroupSettingsOnTop(false);
  se.reset();
  ok &= testNote(69.f, 127.f, regionAmp*sin440, regionAmp*sin440);

  se.setGroupSettingsOnTop(true);
  se.reset();
  ok &= testNote(69.f, 127.f, groupAmp*regionAmp*sin440, groupAmp*regionAmp*sin440);

  se.unsetRegionSetting(0, 0, PST::Volume);
  se.reset();
  ok &= testNote(69.f, 127.f, groupAmp*sin440, groupAmp*sin440);

  se.setGroupSettingsOnTop(false);
  se.reset();
  ok &= testNote(69.f, 127.f, groupAmp*sin440, groupAmp*sin440);
  */




  //se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  //for(int n = 0; n < N; n++)
  //  se.processFrame(&outL[n], &outR[n]);
  //ok &= outL == 0.5f * sin440;
  //ok &= outR == 0.5f * sin440;
  // maybe factor this out into a function: 
  //   testNoteOutput(69.f, 127.f, 0.5f * sin440, 0.5f* sin440)
  // which returns a bool. use it to reduce the boilerplate in the tests. it should take the key 
  // and vel and expected left and right signals (and maybe a tolerance)


  // ToDo: We also need a unsetRegionSetting, unsetGroupSetting, etc. Maybe before implementing, 
  // them it would indeed make sense to refactor such that the error conditions are detected in 
  // rsSamplerData...that requires to define the rsReturnCodes somewhere else....


  int dummy = 0;

  // ToDo:
  // -implement the signal flow: regions -> groups -> instrument, introduce 2 switches in the 
  //  engine: groupSettingsAccumulate, instrumentSettingsAccumulate - test, if fallback vs 
  //  accumulate works as intended ...check, if fallback is actually indeed the sfz behavior, i 
  //  just assumed so
  // -implement and test opcodes for key- and vel-tracking for:
  //  pitch, volume, pan, delay
  // -write a performance test for the sampler
  // -switch to an int+float representation of the current sample position and increment and check, 
  //  if this improves performance...even if not, it's still better because it doesn't lose 
  //  precision for later samples
  // -implement opcodes: pos, width, start, loop_start/end, loop_mode


  // -Implement and test loadFromFile/saveToFile
  //  -may have an option whether or not to save the samples, too (and if so, where)
  // -Test with more complex instruments, featuring:
  //  -multiple groups (with their own settings)
  //  -instrument-wide settings
  // Maybe provide comparison functions with different degrees of strictness, with respect to
  // in which order the opcodes in the settings occur, the handling of the "custom" pointer in the 
  // region, etc. Maybe also have a function to "clean up" a data object by keeping only the very 
  // last setting of a particular kind (the last one would overwrite all others anyway). Maybe it 
  // should return, how many settings had to be removed. Ideally, we would like to always keep only
  // at most one setting of each kind.

  /*
  // Test exporting the instrument-definition related state ins an .sfz-file compliant string:
  using Loader = rsSamplerEngineLoaderSFZ;
  //rsSamplerData sfzData = se.getInstrumentData();
  //rsSamplerData& sfzData = se.getInstrumentData();
  const rsSamplerData& sfzData = se.getInstrumentData();
  // This creates actually a copy - this is good to test (deep) copying of these objects, too.
  // Oh - it's empty. OK, yes - it's because the rsSamplerEngine still has this 
  // std::vector<Group> groups; which should actually go into the instrument in rsSamplerData. This 
  // needs to be refactored first. ok - done, but we need to protect some variables again...
  // ...also, i think, we have a shallow copy here now which would actually allow use to change
  // instrument settings behind the back of the engine. Find some way to make that impossible.
  // Maybe implement deep copying for the rsSamplerData class such that all region-pointers in the 
  // groups actually point to new region objects. Or make rsSamplerData non-copyable to force using a
  // const reference like:
  // const rsSamplerData& sfzData = se.getInstrumentData();
  // ok, done - but using a const reference is not enforced yet - either enforce it or indeed 
  // implement deep copying - or both
  std::string sfzFile = Loader::getAsSFZ(sfzData);
  */



  // ToDo: 
  // -implement and test sfz export/import
  //  -maybe retrive the current instrument settings, set up a 2nd engine object with the same
  //   settings and compare, if both engines are in the same state with respect to these instrument
  //   settings, maybe rsSamplerEngineTest should provide such a comparison function
  // -implement and test better realtime resampling (linear interpolation at first, later cubic and
  //  sinc, maybe some sort of "Elephant" interpolation, too - although, they are supposed to work 
  //  with 2x oversampling) 
  //  -we need a double for the sample-time and an increment...but later, that increment shall be
  //   modulated by pitch-env and -lfo, or maybe use and int and a a float to represent sampleTime
  //   and increment -> do benchmarks, which is faster
  //  -should the played note affect the delay?...nah - i don't think so. maybe sfz had an opcode 
  //   for controlling this? in some situations, that may makes sense, in others not so much

  // -add a SFZPlayer-like simple GUI, so we can import and test actual sfz files
  // -Check, how the sfz player handles the amp/pan parameters with respect to total gain. 
  //  Should there be a factor of 2 for hard left/right settings or a factor of 0.5 for a center 
  //  setting? Make sure to match the behavior of the reference player. I actually would tend to 
  //  prefer the former because it implies that with the neutral default settings, samples are 
  //  played back as is as opposed to having acquired a gain of 0.5. Maybe if sfz behaves the other
  //  way, we could provide both options by introducing additional pan rules.

  // todo:
  //const SE::Group* group1 = se.getGroup(0);
  //const SE::Group* group2 = region.getGroup();
  //ok &= group1 == group2;

  // todo: take the settings member of rsSamplerEngine into an instrument call inside
  // rsInstrumentDataSFZ..or maybe don't create an "Instrument" nested class

  // ToDo: 
  // -set up performance tests, too



  int regionPlayerSize = SE::getRegionPlayerSize();
  // 64 without filters and eq, 512 with - move this into some performance test function
  // currently 160, with virtual functions, it had 16 bytes more. Apparently, that's what the 
  // vtable takes

  rsAssert(ok);
  return ok;
}

bool samplerEngine2UnitTest()
{
  // We test the extended functionality of rsSamplerEngine2, in particular, the signal routing 
  // through per group DSP processes.

  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using AT   = RAPT::rsArrayTools;
  using SE   = rosic::rsSamplerEngine2;
  using RC   = rosic::rsReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  using Ev   = rosic::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  // Create a sine wave as example sample:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of wave
  int   N  = 500;    // length of (co)sinewave sample
  VecF sin440(N);    // sine wave
  VecF tgt;          // target output in tests
  for(int n = 0; n < N; n++)
    sin440[n] = sinf((float)(2*PI*f/fs) * n);

  // Create the sampler engine object and set up a region with the sine sample:
  int maxLayers = 8;  
  SE se(maxLayers);
  float* pSmp[2];
  pSmp[0] = &sin440[0];
  pSmp[1] = nullptr;
  int si = se.addSampleToPool(pSmp, N, 1, fs, "Sine440Hz");    ok &= si == 0;
  int gi = se.addGroup();                                      ok &= gi == 0;
  int ri = se.addRegion(gi);                                   ok &= ri == 0;
  int rc = se.setRegionSample(gi, ri, si);                     ok &= rc == RC::success;
  rc = se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;

  // Set up volume opcode for region and group:
  float regionAmp = 0.5f;
  float groupAmp  = 0.25f;
  float instrAmp  = 1.5f;
  se.setRegionSetting(0, 0, PST::Volume, rsAmpToDb(regionAmp));
  se.setGroupSetting( 0,    PST::Volume, rsAmpToDb(groupAmp));
  se.setInstrumentSetting(  PST::Volume, rsAmpToDb(instrAmp));

  // In the default setting, the regionAmp should override the instrument and the group setting,
  // so the produced output should only have the region volume applied:
  se.reset();
  tgt = regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt);
  ok &= se.getNumActiveLayers() == 0;  // rename to getNumActiveRegions or getNumPlayingRegions
  //ok &= se.getNumActiveGroups() == 0;

  // Set up the engine such that the group settings are applied on top of the region settings:
  se.setGroupSettingsOnTop(true);
  tgt = groupAmp*regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt);
  ok &= se.getNumActiveLayers() == 0;




  // ToDo: 
  // -When setting up the region players, we need to take into account, whether the settings should
  //  have override or accumulate behavior
  //  getRegionPlayerFor, rsSamplerEngine::RegionPlayer::setRegionToPlay, prepareToPlay
  //  need to take 2 boolean parameters for groupSettingsOnTop, instrumentSettingsOnTop
  // -do the same test for other parameters like delay, pan, pitch, etc.
  // -we need to also try it with processes that do not just modify, what algo parameters the 
  //  RegionPlayers use (this is the case for volume, delay, pitch) but that actually apply two
  //  independent DSP processes. That would be the case for a filter. For this, we need to 
  //  implement the DSP chain functionality. ...the maybe implement the filter first. it should 
  //  support 1st order lowpass/highpass, 2nd order bandpass etc...maybe implement it using 
  //  rsBiquad first - later switch to rsStateVariableFilter...which may need some modifications.
  //  maybe let the user also use the ladder filter. or maybe use the rsStateVectorFilter
  // -maybe the baseclass should already support these "on top" settings but implement them by
  //  duplicating/triplicating the DSP processes in each RegionPlayer, if the accumulate option is
  //  chosen. The subclass uses the GroupPlayer only for optimization purposes and to change the 
  //  signal flow (mix before fx - that matters only for nonlinear effects
  // -maybe rename the subclass rsSamplerEngineRoutable



  rsAssert(ok, "samplerEngine2UnitTest failed");
  return ok;
}

bool samplerEngineUnitTestFileIO()
{
  // This test also tests the file I/O

  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rosic::rsSamplerEngineTest;
  using RC   = rosic::rsReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  using Ev   = rosic::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  // Create a sine- and cosine wave-file as example samples:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of (co)sinewave sample
  int   N  = 500;    // length of (co)sinewave sample
  VecF sin440(N);    // sine wave
  VecF cos440(N);    // cosine wave 
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
  {
    sin440[n] = sinf(w*n);
    cos440[n] = cosf(w*n);
  }
  rosic::writeToMonoWaveFile("Sin440Hz.wav", &sin440[0], N, (int)fs, 16);
  rosic::writeToMonoWaveFile("Cos440Hz.wav", &cos440[0], N, (int)fs, 16);
  rosic::writeToStereoWaveFile("SinCos440Hz.wav", &sin440[0], &cos440[0], N, (int)fs, 16);
  // Using "Samples/Sin440Hz.wav" works only, iff the "Samples" folder already exists. ToDo: maybe
  // writeToMonoWaveFile should create it, if it doesn't exist already.

  // Create the engine and instruct it to load the just created sample files into its sample pool:
  int maxLayers = 8;
  SE se(maxLayers);
  int si;                                            // sample index
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  si = se.loadSampleToPool("Cos440Hz.wav"); ok &= si == 1;
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == RC::nothingToDo; // already there

  // Add a group and to that group, add regions for the sine and cosine:
  int gi, ri, rc;                                    // group index, region index, return code
  gi = se.addGroup(); ok &= gi == 0;

  // Set up region for sine:
  ri = se.addRegion(0);                                   ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0);                       ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::Pan, -100.f);          ok &= rc == RC::success;

  // Set up region for cosine:
  ri = se.addRegion(0);                                   ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1);                       ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::Pan, +100.f);          ok &= rc == RC::success;

  // Let the engine produce the sine and cosine:
  VecF outL(N), outR(N);
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumIdleLayers()   == maxLayers-2;
  ok &= se.getNumActiveLayers() == 2;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;

  // The produced signals should equal 2 times the original sin440,cos440 up to 16 bit 
  // quantization noise. The factor 2 comes from the pan law.
  VecF errL = 0.5f * outL - sin440;
  VecF errR = 0.5f * outR - cos440;
  float errMax = rsMax(rsMaxAbs(errL), rsMaxAbs(errR));
  float errMaxDb = rsAmpToDb(errMax);
  ok &= errMaxDb < 90.f;
  //rsPlotVectors(errL, errR);
  // The error is around 3 * 10^-5, or around -90.3 dB. Is that within the expectations? 
  // ToDo: work out exact formula for max-error and compare to that. Maybe also compare signal and
  // error powers and from that the signal-to-noise ratio - this should come out aroun -98dB.

  // Save the state of the engine object se into an sfz file, create a new engine object that loads
  // the sfz file and then test, if both engines are indeed in the same state with respect to the
  // instrument definition:
  se.saveToSFZ("SineCosine.sfz");
  SE se2(maxLayers);
  rc = se2.loadFromSFZ("SineCosine.sfz");
  ok &= rc == RC::success;
  si = se2.loadSampleToPool("Sin440Hz.wav"); ok &= si == RC::nothingToDo;
  si = se2.loadSampleToPool("Cos440Hz.wav"); ok &= si == RC::nothingToDo;
  int nl = se2.getNumSamplesLoaded();  ok &= nl == 2;
  int nr = se2.getNumSamplesRemoved(); ok &= nr == 0;
  int nf = se2.getNumSamplesFailed();  ok &= nf == 0;

  // To test, if se2 has really the same instrument definition as se, produce output and compare. 
  // That's not fool-proof though, but anyway:
  VecF outL2(N), outR2(N);
  se2.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se2.getNumIdleLayers()   == maxLayers-2;
  ok &= se2.getNumActiveLayers() == 2;
  for(int n = 0; n < N; n++)
    se2.processFrame(&outL2[n], &outR2[n]);
  ok &= se2.getNumIdleLayers()   == maxLayers;
  ok &= se2.getNumActiveLayers() == 0;
  ok &= outL2 == outL;
  ok &= outR2 == outR;

  // Test, if it also works when the engine already has some of the samples in its pool already. We
  // create a 3rd, fresh engine object and load the sine sample and the load the sfz. The desired
  // behavior is that loading the sfz only triggers the cosine sample to be loaded to the pool. The
  // sine sample is already there and should not be unloaded and then reloaded. We don't really 
  // have a way to check that automatically, though
  SE se3(maxLayers);
  si = se3.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  rc = se3.loadFromSFZ("SineCosine.sfz");
  ok &= rc == RC::success;
  nl = se3.getNumSamplesLoaded();  ok &= nl == 1;
  nr = se3.getNumSamplesRemoved(); ok &= nr == 0;
  nf = se3.getNumSamplesFailed();  ok &= nf == 0;
  se3.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se3.getNumIdleLayers()   == maxLayers-2;
  ok &= se3.getNumActiveLayers() == 2;
  for(int n = 0; n < N; n++)
    se3.processFrame(&outL2[n], &outR2[n]);
  ok &= se3.getNumIdleLayers()   == maxLayers;
  ok &= se3.getNumActiveLayers() == 0;
  ok &= outL2 == outL;
  ok &= outR2 == outR;

  // Create an sfz patch that uses only the cosine sample and load it into an engine that has both 
  // loaded. The desired behavior is that the engine unloads the sine an keeps the cosine. The 
  // output of the right channel should still be the cosine and the left channel should be silent:
  si = se.findSampleIndexInPool("Sin440Hz.wav"); ok &= si == 0;
  int n = se.getNumRegionsUsing("Sin440Hz.wav"); ok &= n == 1;
  n = se.unUseSample("Sin440Hz.wav"); ok &= n == 1;
  n = se.unUseSample("Sin440Hz.wav"); ok &= n == 0;
  //rc = se.removeSample(si);  // should not really matter for the rest of the test
  se.saveToSFZ("Cosine.sfz");
  // The 1st region in the sfz file now has no sample assigned at all. Maybe unUseSample should 
  // optionally(!) remove the regions that have an empty sample now. Maybe we should also have a 
  // "cleanUp" function that removes all empty regions

  rc = se2.loadFromSFZ("Cosine.sfz");              // Load the cosine-only sfz int se2
  ok &= rc == RC::success;
  nl = se2.getNumSamplesLoaded();  ok &= nl == 0;  // requires no samples to be loaed
  nr = se2.getNumSamplesRemoved(); ok &= nr == 1;  // but one should be removed (the sine)
  nf = se2.getNumSamplesFailed();  ok &= nf == 0;
  se2.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se2.getNumIdleLayers()   == maxLayers-1;
  ok &= se2.getNumActiveLayers() == 1;
  for(int n = 0; n < N; n++)
    se2.processFrame(&outL2[n], &outR2[n]);
  ok &= rsIsAllZeros(outL2);
  ok &= outR2 == outR;

  // Test using a stereo sample:
  se.clearInstrument();
  si = se.loadSampleToPool("SinCos440Hz.wav"); ok &= si == 0;
  ok &= se.getNumSamples() == 1;
  gi = se.addGroup();   ok &= gi == 0;
  ri = se.addRegion(0); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumIdleLayers()   == maxLayers-1;
  ok &= se.getNumActiveLayers() == 1;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL2[n], &outR2[n]);
  ok &= se.getNumIdleLayers()   == maxLayers;
  ok &= se.getNumActiveLayers() == 0;
  errL = 0.5f * outL - outL2; ok &= rsIsAllZeros(errL);
  errR = 0.5f * outR - outR2; ok &= rsIsAllZeros(errR);

  // Save instrument to sfz and load it into engine 3:
  se.saveToSFZ("SineCosine2.sfz");
  rc = se3.loadFromSFZ("SineCosine2.sfz");
  ok &= rc == RC::success;
  nl = se3.getNumSamplesLoaded();  ok &= nl == 1;
  nr = se3.getNumSamplesRemoved(); ok &= nr == 2;
  nf = se3.getNumSamplesFailed();  ok &= nf == 0;

  // Test hikey/lokey opcodes by defining 2 regions:
  se.clearInstrument();
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  si = se.loadSampleToPool("Cos440Hz.wav"); ok &= si == 1;
  gi = se.addGroup(); ok &= gi == 0;
  ri = se.addRegion(0, 59, 61); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f); ok &= rc == RC::success;
  ri = se.addRegion(0, 69, 71); ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::PitchKeyCenter, 70.f); ok &= rc == RC::success;
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);

  // ToDo:
  // -maybe make a local function testSaveLoadRoundTrip(se, ...) that saves the state of se and 
  //  loads it into a new instance and produces and compares some output of both engines...it also
  //  needs a list of events that should be used to produce the output...or maybe it should loop
  //  all keys and velocities within some range and produce a couple of samples for each
  //  ..or maybe just compare the sfz data object
  // -Test using a custom sfz and/or sample directory. Maybe the engine needs members
  //  sfzDir, wavDir. If they are empty, the project folder is used by default, but that's not how
  //  it should work in ToolChain. Perhaps, in ToolChain, we need to re-implement the loading 
  //  anyway in order to support more file formats (in particular, flac). For the moment, we can 
  //  only use 16 bit wav. ...or maybe use the TinyWav library? ...but we also want .flac
  // -make SamplePool, AudioStream, SamplerEngine etc. objects or non-copyable - maybe it's
  //  enough to make AudioStream non-copyable - the feature will then propagate to all aggregates

  rsAssert(ok);
  return ok;
}



bool samplerEngineUnitTest()
{
  bool ok = true;

  ok &= samplerDataUnitTest();
  ok &= samplerEngineUnitTest1();
  ok &= samplerEngineUnitTestFileIO();

  ok &= samplerEngine2UnitTest();


  rsAssert(ok, "samplerEngineUnitTest failed");
  return ok;
}
