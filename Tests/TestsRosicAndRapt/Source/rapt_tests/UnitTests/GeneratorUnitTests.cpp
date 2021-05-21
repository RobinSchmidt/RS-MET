bool samplerDataUnitTest()
{
  bool ok = true;

  using SD = rsDataSFZ;

  SD d1;

  ok &= d1.getNumGroups() == 0;

  rsAssert(ok);
  return ok;
}

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
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
  {
    sin440[n] = sinf(w*n);
    cos440[n] = cosf(w*n);
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

  SE::Region* r = se.getRegion(gi, ri);
  int rc = se.setRegionSample(gi, ri, si);           // rc: return code
  ok &= rc == RC::success;
  rc = se.setRegionSetting(r, PST::PitchKeyCenter, 69.f); // key = 69 is A4 at 440 Hz
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
  r = se.getRegion(gi, ri);                               // retrieve pointer to region
  rc = se.setRegionSample(gi, ri, si);                    // set region sample to cosine
  ok &= rc == RC::success;
  se.setRegionSetting(r, PST::PitchKeyCenter, 69.f);      // cosine is also at A4 = 440 Hz
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

  const rsDataSFZ& sfzData = se.getInstrumentData();
  std::string sfzString = sfzData.getAsSFZ();
  rsDataSFZ sfzData2;
  sfzData2.setFromSFZ(sfzString);
  ok &= sfzData2 == sfzData;
  // ToDo:
  // -Implement and test loadFromFile/safeToFile
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


  //SE se2;
  //se2.setupFromSFZ(sfzData2);
  //ok &= se2.matchesInstrumentDefinition(se);


  int dummy = 0;


  /*
  // Test exporting the instrument-definition related state ins an .sfz-file compliant string:
  using Loader = rsSamplerEngineLoaderSFZ;
  //rsDataSFZ sfzData = se.getInstrumentData();
  //rsDataSFZ& sfzData = se.getInstrumentData();
  const rsDataSFZ& sfzData = se.getInstrumentData();
  // This creates actually a copy - this is good to test (deep) copying of these objects, too.
  // Oh - it's empty. OK, yes - it's because the rsSamplerEngine still has this 
  // std::vector<Group> groups; which should actually go into the instrument in rsDataSFZ. This 
  // needs to be refactored first. ok - done, but we need to protect some variables again...
  // ...also, i think, we have a shallow copy here now which would actually allow use to change
  // instrument settings behind the back of the engine. Find some way to make that impossible.
  // Maybe implement deep copying for the rsDataSFZ class such that all region-pointers in the 
  // groups actually point to new region objects. Or make rsDataSFZ non-copyable to force using a
  // const reference like:
  // const rsDataSFZ& sfzData = se.getInstrumentData();
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
  // currently 176

  rsAssert(ok);
  return ok;
}

bool samplerEngineUnitTestFileIO()
{
  // This test also tests the file I/O

  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rsSamplerEngineTest;
  using RC   = SE::ReturnCode;
  using PST  = SE::PlaybackSetting::Type;
  using Ev   = rsMusicalEvent<float>;
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
  // Using "Samples/Sin440Hz.wav" works only, iff the "Samples" folder already exists. ToDo: maybe
  // writeToMonoWaveFile should create it, if it doesn't exist already.

  // Create the engine and instruct it to load the just created sample files into its sample pool:
  int maxLayers = 8;
  SE se(maxLayers);
  int si;                                            // sample index
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  si = se.loadSampleToPool("Cos440Hz.wav"); ok &= si == 1;

  // Add a group and to that group, add regions for the sine and cosine:
  int gi, ri, rc;                                    // group index, region index, return code
  SE::Region* r;
  gi = se.addGroup(); ok &= gi == 0;

  // Set up region for sine:
  ri = se.addRegion(0);                                   ok &= ri == 0;
  r  = se.getRegion(0, 0);                                ok &= r != nullptr;
  rc = se.setRegionSample(0, 0, 0);                       ok &= rc == RC::success;
  rc = se.setRegionSetting(r, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  rc = se.setRegionSetting(r, PST::Pan, -100.f);          ok &= rc == RC::success;

  // Set up region for cosine:
  ri = se.addRegion(0);                                   ok &= ri == 1;
  r  = se.getRegion(0, 1);                                ok &= r != nullptr;
  rc = se.setRegionSample(0, 1, 1);                       ok &= rc == RC::success;
  rc = se.setRegionSetting(r, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  rc = se.setRegionSetting(r, PST::Pan, +100.f);          ok &= rc == RC::success;

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
  // loKey/hiKey settings for the regions are wrong - the sfz recall does not yet include the 
  // loKey/hiKey settings

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
  // This test still fails - i think, when loading an sfz, we also need to update the regionsForKey
  // array


  // ToDo: test, if it also works when the engine already has some of the samples in its pool 
  // already


  // ToDo: 
  // -export the settings to an sfz file.
  // -import the sfz file with a new instance and compare the states
  // -instead of passing the region pointer r to setRegionSetting, pass group/region indices
  //  -that should be more efficient, because it bypasses the findRegion(region, &gi, &ri); call
  //   there


  rsAssert(ok);
  return ok;
}