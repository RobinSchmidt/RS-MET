
// Helper functions:

// Computes a linearly interpolated value from the given vector at the given position:
template<class T>
T getSampleAt(const std::vector<T>& v, T pos)
{
  if(pos < T(0)) return T(0);

  int i = (int) pos;
  T f = pos - (T)i;

  T x0(0), x1(0);
  int N = (int)v.size();
  if(i   < N) x0 = v[i];
  if(i+1 < N) x1 = v[i+1];

  return (T(1) - f) * x0 + f * x1;
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

// New implementation, may make the one above obsolete
// -xL,xR are values that indicate, how the signal is supposed to continue beyond the left and 
//  right boundaries
template<class T>
T getSampleAtNew(const std::vector<T>& v, T pos, T xL = T(0), T xR = T(0))
{
  int N = (int)v.size();
  if(pos <= T(-1)) return xL;
  if(pos >= T(N) ) return xR;

  int i = (int) floor(pos);       // floor of position
  T   f = pos - (T)i;             // fractional part of pos

  if(i   >= 0) xL = v[i];
  if(i+1 <  N) xR = v[i+1];

  return (T(1) - f) * xL + f * xR;
}
// needs verification

bool testGetSampleAt()  // unit test for function above
{
  bool ok = true;
  using Vec = std::vector<float>;

  Vec v;  // empty vector

  float y;

  float xL =  9.f;
  float xR = 11.f;

  // For an empty vector, it ramps from xL to xR for pos in -1..0:
  y = getSampleAtNew(v, -1.25f, xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -1.f,   xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -0.75f, xL, xR); ok &= y == 0.75*xL + 0.25*xR;
  y = getSampleAtNew(v, -0.5f,  xL, xR); ok &= y == (xL+xR)/2.f;
  y = getSampleAtNew(v, -0.25f, xL, xR); ok &= y == 0.25*xL + 0.75*xR;
  y = getSampleAtNew(v,  0.f,   xL, xR); ok &= y == xR;
  y = getSampleAtNew(v, +0.25f, xL, xR); ok &= y == xR;
  y = getSampleAtNew(v, +1.f,   xL, xR); ok &= y == xR;
  y = getSampleAtNew(v, +1.25f, xL, xR); ok &= y == xR;

  // For a 1-element vector with the sole element v[0] = y0, we expect to see a ramp from xL to x0
  // for pos in -1...0 and a ramp from x0 to xR for pos in 0...+1:
  float x0 = 20.0f;  // rename to x0
  v.push_back(x0);
  y = getSampleAtNew(v, -1.25f, xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -1.f,   xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -0.75f, xL, xR); ok &= y == 0.75*xL + 0.25*x0;
  y = getSampleAtNew(v, -0.25f, xL, xR); ok &= y == 0.25*xL + 0.75*x0;
  y = getSampleAtNew(v,  0.f,   xL, xR); ok &= y == x0;
  y = getSampleAtNew(v, +0.25f, xL, xR); ok &= y == 0.75*x0 + 0.25*xR;
  y = getSampleAtNew(v, +0.75f, xL, xR); ok &= y == 0.25*x0 + 0.75*xR;
  y = getSampleAtNew(v, +1.f,   xL, xR); ok &= y == xR;
  y = getSampleAtNew(v, +1.25f, xL, xR); ok &= y == xR;

  // For a 2-element vector, with elements x0,x1, we expect to see a ramp from xL to x0
  // for pos in -1...0, a linear interpolation between x0 and x1 for pos in 0...1 and a ramp from 
  // x1 to xR for pos in 1...2:
  float x1 = 25.0f;
  v.push_back(x1);
  y = getSampleAtNew(v, -1.25f, xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -1.f,   xL, xR); ok &= y == xL;
  y = getSampleAtNew(v, -0.75f, xL, xR); ok &= y == 0.75*xL + 0.25*x0;
  y = getSampleAtNew(v, -0.25f, xL, xR); ok &= y == 0.25*xL + 0.75*x0;
  y = getSampleAtNew(v,  0.f,   xL, xR); ok &= y == x0;
  y = getSampleAtNew(v, +0.25f, xL, xR); ok &= y == 0.75*x0 + 0.25*x1;
  y = getSampleAtNew(v, +0.75f, xL, xR); ok &= y == 0.25*x0 + 0.75*x1;
  y = getSampleAtNew(v, +1.f,   xL, xR); ok &= y == x1;
  y = getSampleAtNew(v, +1.25f, xL, xR); ok &= y == 0.75*x1 + 0.25*xR;
  y = getSampleAtNew(v, +1.75f, xL, xR); ok &= y == 0.25*x1 + 0.75*xR;
  y = getSampleAtNew(v, +2.f,   xL, xR); ok &= y == xR;
  y = getSampleAtNew(v, +2.25f, xL, xR); ok &= y == xR;

  // ToDo: maybe make a little function that can be called like: 
  //   ok &= test(-0.75, 0.75*xL + 0.25*x0); etc. to reduce boilerplate


  return ok;
}

// pan is suppoesed to be in -1...+1
template<class T>
void rsApplyPan(std::vector<T>& L, std::vector<T>& R, T pan, bool constPow = false)
{
  int N = L.size(); rsAssert(R.size() == N);
  using AT = RAPT::rsArrayTools;

  T t = (pan + 1) * 0.5;     // -1...+1  ->   0...1
  AT::scale(&L[0], N, 2 * (1 - t));
  AT::scale(&R[0], N, 2 *      t);
}

// delay is supposed to be in samples - todo allow non-integer delays
template<class T>
void rsApplyDelay(std::vector<T>& x, int delay)
{
  rsArrayTools::shift(&x[0], (int) x.size(), delay);
}

template<class T>
void rsApplyDelay(std::vector<T>& x, T delay)
{
  std::vector<T> tmp = x;  // do we need this?
  for(int i = 0; i < (int) x.size(); i++)
    x[i] = getSampleAt(tmp, T(i-delay));
  // ToDo: use getSampleAtNew and update the unit tests accordingly
  // ...maybe before doing that, move getSampleAtNew into rsArrayTools, mabye as 
  // lerp(TSig* x, int N, TPos n)
}

template<class T>
T rsEstimateMidiPitch(const std::vector<T>& x, T sampleRate)
{
  std::vector<double> xd = rsConvert(x, double());
  rsCycleMarkFinder<double> cmf(sampleRate);

  cmf.setSubSampleApproximationPrecision(0);
  std::vector<double> marks = cmf.findCycleMarks(&xd[0], (int)x.size());

  //rsPlotSignalWithMarkers(&xd[0], (int)xd.size(), &marks[0], (int)marks.size());
  // not very precise, the setSubSampleApproximationPrecision call seems to have no effect
  // ..oh - in the plot we only have a linear interpolation, so maybe our etsimated values are 
  // indeed more precise than eyballing - but with precision setting of zero, we should actually
  // find the zeros of an linear interpolant...hmmm...investigate this at some point

  double T = (marks[1] - marks[0]) / sampleRate;  // period of 1st cycle in seconds
  double f = 1/T;
  double p = rsFreqToPitch(f);
  return p; 
}

template<class T>
std::vector<T> createColoredNoise(int N, T spectralSlope, int seed = 0)
{
  // Create white noise:
  std::vector<T> y(N);
  RAPT::rsNoiseGenerator<T> ng;
  ng.setSeed(seed);
  for(int n = 0; n < N; n++)
    y[n] = ng.getSample();

  // Apply coloring:
  rosic::SlopeFilter flt;
  flt.setSlope(spectralSlope);
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(y[n]);

  return y;
}
// move to test tools

//=================================================================================================

bool samplerDataUnitTest()
{
  bool ok = true;

  //using SD = rsSamplerData;
  using SFZT = rosic::Sampler::SfzCodeBook;
  using SD   = rosic::Sampler::rsSamplerData;
  using PST  = rosic::Sampler::Opcode;

  SFZT::createInstance();
  // Normally, this is supposed to be done in the constructor of rsSamplerEngine and objects of 
  // type rsSamplerData are supposed to live only inside the engine. But here in the test, we 
  // create these data objects with having an engine around, so we must take over the 
  // responsibility for the lifetime of the SfzCodeBook which is used in the data object.

  SD d1;
  ok &= d1.getNumGroups() == 0;

  // todo: add some groups/regions and test copy-assignment we need deep copies and that doesn't 
  // seem to work yet

  int gi, ri;             // gi: group index, ri: region index
  std::string sfz;

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


  /*
  // These still fail:
  // Test parsing of some degenerate sfz strings:
  sfz = "";
  d2.setFromSFZ(sfz);  // should result in an empty instrument
  //ok &= d2.isEmpty();  // implement this

  sfz = "<region>";
  d2.setFromSFZ(sfz); // crashes! the parser gets tripped up by empty regions
  // crash: in rsSamplerData::setFromSFZ, j0 = 18446744073709551615

  // Set volume settings for instrument, group(0), and region(0,0):
  d1.setInstrumentSetting(  PST::Volume, -2.f);
  d1.setGroupSetting(0,     PST::Volume, -3.f);
  d1.setRegionSetting(0, 0, PST::Volume, -5.f);

  // To test sfz generaion and parsing, retrieve the sfz string, set up d2 from the string and 
  // check, if it equal to d1:
  sfz = d1.getAsSFZ();
  d2.setFromSFZ(sfz);
  ok &= d2 == d1;
  */


  // Create an isntument with 1 group containing 1 region:
  d3.clearInstrument();
  gi = d3.addGroup();   ok &= gi == 0;
  ri = d3.addRegion(0); ok &= ri == 0;

  // Set volume settings for instrument, group(0), and region(0,0):
  //d3.setInstrumentSetting(  PST::Volume, -2.f);
  d3.setGroupSetting(0,     PST::Volume, -3.f);
  d3.setRegionSetting(0, 0, PST::Volume, -5.f);

  // Test sfz generation and parsing:
  sfz = d3.getAsSFZ();
  d2.setFromSFZ(sfz);  // crashes, if we uncomment the setInstrumentSetting call
  ok &= d2 == d3;

  // ToDo: make a more complex sfz patch, using more opcodes, samples, etc.



  // Parser bugs:
  // -empty regions lead to crashes
  // -instument settings lead to crashes










  // this triggers an assert - the assignment operator needs to be implemented differently
  // ..or better: avoid using pointers and new operator for adding regions
  // -in rsSamplerEngine, refer to regions not by a pointer to the region (in setRegionSetting 
  //  etc.), but instead by group index and region index
  // -oh - no - we need to use pointers because of rsSamplerEngine::regionsForKey
  // -actually, we need a pointer array for the groups also because the regions refer back to them
  // 

  SFZT::deleteInstance();
  // We must clean up what we have created - otherwise, we'll trigger a memleak.

  rsAssert(ok);
  return ok;
}

/** Fills the outL, outR arrays with the output of the given sampler engine for the given note. */
void getSamplerNote(rosic::Sampler::rsSamplerEngine* se, float key, float vel,
  std::vector<float>& outL, std::vector<float>& outR)
{
  rsAssert(outL.size() == outR.size());
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;
  se->handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
  for(int n = 0; n < (int) outL.size(); n++)
    se->processFrame(&outL[n], &outR[n]);
}
// should we clear the outL/R arrays first? maybe not, if we want instruments to accumuluate 
// their outputs in ToolChain

bool testSamplerNote(rosic::Sampler::rsSamplerEngine* se, float key, float vel, 
  const std::vector<float>& targetL, const std::vector<float>& targetR, 
  float tol = 0.f, bool plot = false)
{
  using AT   = RAPT::rsArrayTools;
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  int N = (int) targetL.size();
  rsAssert((int)targetR.size() == N);
  std::vector<float> outL(N), outR(N);

  se->handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
  for(int n = 0; n < N; n++)
    se->processFrame(&outL[n], &outR[n]);
  float errL = AT::maxDeviation(&outL[0], &targetL[0], N);
  float errR = AT::maxDeviation(&outR[0], &targetR[0], N);
  if(plot)
    rsPlotVectors(targetL, targetR, outL, outR);
  return errL <= tol && errR <= tol;
};
// maybe have a bool resetBefore that optionally resets the engine before playing...but maybe it's
// better when the caller does the reset directly, if desired - it's not much longer but clearer
// factor out a getSamplerOutput(rosic::rsSamplerEngine* se, float key, float vel, 
// const std::vector<float>& targetL, const std::vector<float>& targetR, bool plot) function

bool samplerEngineUnitTest1()
{
  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using AT   = RAPT::rsArrayTools;
  using SD   = rosic::Sampler::rsSamplerData;
  using SE   = rosic::Sampler::rsSamplerEngineTest;
  using RC   = rosic::Sampler::rsReturnCode;
  using PST  = rosic::Sampler::Opcode;
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
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

  // Remove the 1st region (sin440, on left channel) in the middle of playback. Desired bevavior: 
  // the playback of the respective region immediately stops, the rest in unaffected.
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumActiveLayers() == 2;
  for(int n = 0; n < N/2; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= se.getNumActiveLayers() == 2;
  se.removeRegion(0, 0);               // remove a region in the middle of playback 
  ok &= se.getNumActiveLayers() == 1;
  for(int n = N/2; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= se.getNumActiveLayers() == 0;
  //rsPlotVectors(outL, outR); 
  ok &= outR == (2.f * cos440);        // 2nd region (R channel) is unaffected by removeRegion
  for(int n = 0; n < N/2; n++)         // L channel's 1st half is as always
    ok &= outL[n] == 2.f * sin440[n];
  for(int n = N/2; n < N; n++)         // L channel's 2nd half is zero
    ok &= outL[n] == 0.f;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  ok &= se.getNumActiveLayers() == 1;
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= se.getNumActiveLayers() == 0;
  //rsPlotVectors(outL, outR); 
  ok &= outL == zeros;
  ok &= outR == (2.f * cos440);

  // Re-assign the sample for 1st region to sine:
  se.setRegionSample(0, 0, 0);
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  ok &= outL == zeros;
  ok &= outR == (2.f * sin440);
  //rsPlotVectors(outL, outR); 


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
  VecF tgt = sin440;
  rsApplyDelay(tgt, delaySamples);
  //rsPlotVectors(tgt, outL); 


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
  // -Currently at 176
  // -With virtual functions, it had 16 bytes more. Apparently, that's what the vftable take.
  // -Move this into some performance test function

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
  using SE   = rosic::Sampler::rsSamplerEngine2Test;
  using RC   = rosic::Sampler::rsReturnCode;
  using PST  = rosic::Sampler::Opcode;
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;

  // Create a sine wave as example sample:
  float fs = 44100;      // sample rate
  float f  = 440.0;      // frequency of wave
  int   N  = 500;        // length of (co)sinewave sample
  VecF sin440(N);        // sine wave
  VecF tgt;              // target output in tests
  VecF tgtL, tgtR;       // ...for when we need different left and right target signal
  VecF outL(N), outR(N); // for the output signals
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

  //---------------------------------------------------------------------------
  // Test accumulation of amp setting:

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
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 0.f, false);
  ok &= se.getNumActiveLayers() == 0;  // rename to getNumActiveRegions or getNumPlayingRegions
  ok &= se.getNumActiveGroupPlayers() == 0;

  // Set up the engine such that the region settings do not override the group settings but instead
  // region settings and group settings are combined:
  se.setRegionSettingsOverride(false);
  tgt = groupAmp*regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 0.0f, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0;

  // Set up the engine such that the group settings do not override the instrument settings but 
  // instead group settings and instrument settings are combined. The region settings are combined
  // into that, too (from the setting before):
  se.setGroupSettingsOverride(false);
  tgt = instrAmp*groupAmp*regionAmp*sin440;
  float tol = 1.e-7f;  // why can't we use 0 tolerance?
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, tol, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0;

  // Now set the regionSettingsOverride back to true. That renders the question whether or not 
  // group settings override instrument settings irrelevant in cases when a region setting is 
  // available (which is the case here). The region setting will override the combined 
  // instrument + group setting:
  //se.reset();
  se.setRegionSettingsOverride(true);
  tgt = regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 0.0, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0;


  // ToDo: test behavior when there is no region setting available...maybe test also when there is
  // only a region and instrument setting


  //---------------------------------------------------------------------------
  // Test accumulation of pan setting:

  se.clearAllSfzSettings();                                   // remove all the amp settings
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f);  // restore the rootkey setting
  ok &= rc == RC::success;
  float regionPan = 10.f;   // slightly right, pan range is -100...+100
  float groupPan  = 20.f;
  float instrPan  = 30.f;
  se.setRegionSetting(0, 0, PST::Pan, regionPan);
  se.setGroupSetting( 0,    PST::Pan, groupPan);
  se.setInstrumentSetting(  PST::Pan, instrPan);

  // We want to see only the region pan:
  se.setGroupSettingsOverride(true);     // group settings override again
  se.setRegionSettingsOverride(true);    // not required but anyway
  tgtL = tgtR = sin440;
  rsApplyPan(tgtL, tgtR, regionPan/100);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgtL, tgtR, 1.e-6f, false); 

  // Now we want to see region and group pan combined:
  se.setRegionSettingsOverride(false);
  tgtL = tgtR = sin440;
  rsApplyPan(tgtL, tgtR, regionPan/100);
  rsApplyPan(tgtL, tgtR, groupPan /100);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgtL, tgtR, 1.e-6f, false);

  // Now we want to see region, group and instrument pan combined:
  se.setGroupSettingsOverride(false);
  tgtL = tgtR = sin440;
  rsApplyPan(tgtL, tgtR, regionPan/100);
  rsApplyPan(tgtL, tgtR, groupPan /100);
  rsApplyPan(tgtL, tgtR, instrPan /100);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgtL, tgtR, 1.e-6f, false);

  // ToDo: test it also with constant power pan rule

  //---------------------------------------------------------------------------
  // Test delay accumulation:

  se.clearAllSfzSettings();                                   // remove all the amp settings
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f);  // restore the rootkey setting
  int regionDelay = 10;   // in samples - todo: use float
  int groupDelay  = 20;
  int instrDelay  = 40;
  se.setRegionSetting(0, 0, PST::Delay, regionDelay / fs);
  se.setGroupSetting( 0,    PST::Delay, groupDelay  / fs);
  se.setInstrumentSetting(  PST::Delay, instrDelay  / fs);

  // We want to see only the region delay:
  se.setGroupSettingsOverride(true);
  se.setRegionSettingsOverride(true); 
  tgt = sin440;
  rsApplyDelay(tgt, regionDelay);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);
  ok &= se.getNumActiveLayers() == 1;        // it's still playing due to the delay
  ok &= se.getNumActiveGroupPlayers() == 0;  // no group player is/was used due to settings

  // Now we want to see region and group delay combined:
  se.setRegionSettingsOverride(false);
  tgt = sin440;
  rsApplyDelay(tgt, regionDelay);
  rsApplyDelay(tgt, groupDelay);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);

  // Now we want to see region, group and instrument delay combined:
  se.setGroupSettingsOverride(false);
  tgt = sin440;
  rsApplyDelay(tgt, regionDelay);
  rsApplyDelay(tgt, groupDelay);
  rsApplyDelay(tgt, instrDelay);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);


  //---------------------------------------------------------------------------
  // Test offset accumulation:

  float regionOffset = 10;   // "offset" opcode in samples
  float groupOffset  = 20;
  float instrOffset  = 40;
  se.clearAllSfzSettings();
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f);
  se.setRegionSetting(0, 0, PST::Offset, regionOffset);
  se.setGroupSetting( 0,    PST::Offset, groupOffset);
  se.setInstrumentSetting(  PST::Offset, instrOffset);

  // We want to see only the region offset:
  se.setGroupSettingsOverride(true);
  se.setRegionSettingsOverride(true); 
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset);  // offset is like a negative delay
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0; 

  // We want to see region and group offset:
  se.setRegionSettingsOverride(false);
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset); 
  rsApplyDelay(tgt, -groupOffset); 
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);

  // We want to see region, group and instrument offset:
  se.setGroupSettingsOverride(false);
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset); 
  rsApplyDelay(tgt, -groupOffset); 
  rsApplyDelay(tgt, -instrOffset); 
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, false);


  //---------------------------------------------------------------------------
  // Test offset and delay (but only for the region setting):

  se.clearAllSfzSettings();
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f);
  se.setRegionSetting(0, 0, PST::Offset, regionOffset);
  se.setRegionSetting(0, 0, PST::Delay,  regionDelay / fs);
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset); 
  rsApplyDelay(tgt,  regionDelay);
  //ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7, true);
  // Still fails. They are actually in sync but the sampler has skipped the first sample

  // when uncommenting stream->getFrame... stuff RegionPlayer::getFrame, it produces nonzero 
  // values - but they are slightly wrong and this messes up also another test


  int dummy = 0;


  // ToDo: 
  // -test offset together with delay. offset should jump forward in the sample, delay should
  //  delay it
  // -maybe allow floating point offsets, maybe even negative ones


  //---------------------------------------------------------------------------
  // Test pitch accumulation 

  // Set up the transposition and tune opcodes (for coarse and fine detune):
  float regionTrans = 1;   // "transpose" opcode (in semitones, i think)
  float groupTrans  = 2;
  float instrTrans  = 3;
  float regionTune  = 10;  // "tune" opcode (in cents)
  float groupTune   = 20;
  float instrTune   = 30;
  se.clearAllSfzSettings();                               // remove all the amp settings
  se.setGroupSetting( 0,    PST::PitchKeyCenter, 50.f);   // should always be overriden
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f);   // restore the rootkey setting
  se.setRegionSetting(0, 0, PST::Transpose, regionTrans);
  se.setGroupSetting( 0,    PST::Transpose, groupTrans);
  se.setInstrumentSetting(  PST::Transpose, instrTrans);
  se.setRegionSetting(0, 0, PST::Tune,      regionTune);
  se.setGroupSetting( 0,    PST::Tune,      groupTune);
  se.setInstrumentSetting(  PST::Tune,      instrTune);

  tol = 0.02f; 
  // The error of the pitch estimation is around 2 cents...that's quite a large error actually.
  // Try to improve this at some point!I think, the algo should give better results!

  // Test override mode. We expect to see only the region transpose and tune. That's 
  // 69 + 1 + 10/100 = 70.1
  se.setGroupSettingsOverride(true);
  se.setRegionSettingsOverride(true);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVectors(outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 70.1f, tol);

  // Now we want to see region and group transpose combined. That's
  // 69 + 1 + 2 + 10/100 + 20/100 = 72.3
  se.setRegionSettingsOverride(false);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVectors(outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 72.3f, tol);
  // that works! ...i'm actually surprised that it does -> figure out why

  // Now we want to see region, group and instrument transpose combined. That's
  // 69 + 1 + 2 + 3 + 10/100 + 20/100 + 30/100 = 75.6
  se.setGroupSettingsOverride(false);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 75.6f, tol);
  // ok - works, too

  // ToDo: check what happens when region setting and/or group setting and/or instrument setting
  // are removed and if it behaves as expected. we may need functionality to delete particular 
  // opcodes like se.removeRegionSetting(0, 0, PST::Tune) etc.

  // Remove the transpose setting for the region. We expect to see a combination of transpose
  // settings of instrument and group and a combination of the tune settings of all 3, so the pitch
  // should be: 69 + 2 + 3 + 0.1 + 0.2 + 0.3 = 74.6
  ok &= se.removeRegionSetting(0, 0, PST::Transpose) == RC::success;
  ok &= se.removeRegionSetting(0, 0, PST::Transpose) == RC::nothingToDo;
  //se.reset();  // what happens if we don't reset? seems to make no difference
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.6f, tol);

  // Remove tune setting from group. p = 69 + 2 + 3 + 0.1 + 0.3 = 74.4
  ok &= se.removeGroupSetting(0, PST::Tune) == RC::success;
  ok &= se.removeGroupSetting(0, PST::Tune) == RC::nothingToDo;
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVector(outL);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.4f, tol);

  // Remove tune setting from instrument. p = 69 + 2 + 3 + 0.1 = 74.1
  ok &= se.removeInstrumentSetting(PST::Tune) == RC::success;
  ok &= se.removeInstrumentSetting(PST::Tune) == RC::nothingToDo;
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVector(outL);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.1f, tol);

  // Restore the instrument's tune setting and let the group settings override the instrument 
  // settings again. Now we should see for tune the instrument setting combined with the region 
  // setting because the region accumulates and the group has no tune setting. For the transpose,
  // we should just see the group setting because it overrides the instrument's tune and the 
  // region has no tune anymore. So we expect: p = 69 + 2 + 0.1 + 0.3 = 71.4
  se.setInstrumentSetting(PST::Tune, instrTune);
  se.setGroupSettingsOverride(true);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVector(outL);
  float tmp = rsEstimateMidiPitch(outL, fs);
  //ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 71.4f, tol);
  // fails! it's 71.1. seems like when the group settings overrides the instrument setting, the 
  // absence of a tune setting in the group causes the tune to be overriden with 0? i.e. when the 
  // group overrides the instrument, it will also override it with the default value, in case no 
  // value is defined?. why would that happen?






  // with se.setGroupSetting( 0,    PST::PitchKeyCenter, 50.f) uncommented, the test fails



  // we are in accumulate mode for both, so after that removal, we should see a combinations of
  // instrument and group transpose...but still see the region's tune

  //int dummy = 0;
  
  // -For PitchKeyCenter, accumulation makes actually no sense. So maybe, this parameter should 
  //  always work in override mode. Check what happens, when the instrument and or group also
  //  defines a pitch keycenter. We currently wouldn't notice any problems with that because the
  //  keycenter is only defined for the region




  // Maybe we should also have a detuneHz opcode -> check sfz spec, if such a thing exists 
  // (-> nope, not in sfz 1.0 at least)


  // ...wait...the code for the accumulation of parameters will not work as desired because some
  // user parameters are themselves accumulative. Take PitchKeyCenter, DetuneCoarse, DetuneFine: 
  // all accumulate into the increment. Check what happens in the different when a group defines 
  // coarse and fine and the region defines only one of them? or vice versa. I think, we cannot 
  // just accumulate all values into the increment as we currently do. If we are in override mode 
  // and the group has already baked both of its detunes into the increment, we won't get them out
  // again with the current code. But that's what we would need. Maybe we need a temporary 
  // data-structure to hold the accumulated or overriden parameters seperately. Test it with the 
  // following settings: 
  //   Instr:  coarse: 1 semitone,  fine: 10 cents
  //   Group:  coarse: 2 semitones, fine: 20 cents
  //   Region: coarse: 4 semitones, fine: 40 cents
  // 


  // The same problem may arise when we want to support the stereo-width opcode because it 
  // interferes with Pan





  // ToDo: 
  // -do the same test for other parameters like delay, pitch, etc.
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
  // -figure out what should happen when only one of the "onTop" settings is true and check if the 
  //  behavior is as desired



  //rsAssert(ok, "samplerEngine2UnitTest failed");
  return ok;
}

bool samplerEngineUnitTestFileIO()
{
  // This test also tests the file I/O

  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rosic::Sampler::rsSamplerEngineTest;
  using RC   = rosic::Sampler::rsReturnCode;
  using PST  = rosic::Sampler::Opcode;
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
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
  ri = se.addRegion(0);             ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::Pan, -100.f);          ok &= rc == RC::success;

  // Set up region for cosine:
  ri = se.addRegion(0);             ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1); ok &= rc == RC::success;
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

  // Test filter opcodes:
  using FltType = rosic::Sampler::FilterType;
  se.clearInstrument();
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  gi = se.addGroup(); ok &= gi == 0;
  ri = se.addRegion(0); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f); ok &= rc == RC::success;
  se.setRegionSetting(0, 0, PST::filN_type, (float) FltType::lp_6, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,  1000.f, 1);
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);
  se2.saveToSFZ("tmp2.sfz");                        // For manual inspection 
  //ok &= rsAreFilesEqual("tmp.sfz", "tmp2.sfz");   // ToDo: write this function
  // FAILS! 
  // -the signalProcessors array of the region in se2 is not the same as in se
  // i think, in loadFromSFZ, we must make sure to add the appropriate dsp-types to the array

  // Test equalizer opcodes:

  // Set up eq1,2,3 without specifying frequencies or bandwidths (i.e. using the defaults):
  se.setRegionSetting(0, 0, PST::eqN_gain, 1.f, 1);
  se.setRegionSetting(0, 0, PST::eqN_gain, 2.f, 2);
  se.setRegionSetting(0, 0, PST::eqN_gain, 3.f, 3);
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);

  // Set up some more eq bands, still without specifying freqs or widths:
  se.setRegionSetting(0, 0, PST::eqN_gain,  4.f,  4);
  se.setRegionSetting(0, 0, PST::eqN_gain,  5.f,  5);
  se.setRegionSetting(0, 0, PST::eqN_gain,  8.f,  8);
  se.setRegionSetting(0, 0, PST::eqN_gain, 13.f, 13); 
  se.setRegionSetting(0, 0, PST::eqN_gain, 10.f, 10); 
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);

  // Set up frequencies and bandwidths of some of the eq bands:
  se.setRegionSetting(0, 0, PST::eqN_freq, 800.f,  8);
  se.setRegionSetting(0, 0, PST::eqN_bw,     1.8f, 8);
  se.setRegionSetting(0, 0, PST::eqN_freq, 500.f,  5);
  se.setRegionSetting(0, 0, PST::eqN_bw,     1.3f, 3);
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);

  // Clear the region's settings, then set up 3 filters:
  using Region = rosic::Sampler::rsSamplerData::Region;
  Region* r = se.getRegion(0, 0);
  r->clearSettings();
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::hp_12, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,    200.f,                  1);
  se.setRegionSetting(0, 0, PST::resonanceN, 10.f,                   1);
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::lp_12, 2);
  se.setRegionSetting(0, 0, PST::cutoffN,    800.f,                  2);
  se.setRegionSetting(0, 0, PST::resonanceN, 15.f,                   2);
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::lp_6,  3);
  se.setRegionSetting(0, 0, PST::cutoffN,    5000.f,                 3);
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

bool samplerFilterTest()
{
  bool ok = true;

  ok &= testGetSampleAt(); // preliminary - should go into a rapt unit test

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rosic::Sampler::rsSamplerEngineTest;
  using PST  = rosic::Sampler::Opcode;
  using Type = rosic::Sampler::FilterType;

  // Create a pinkish noise as example sample:
  float fs     = 44100.f;    // sample rate
  float cutoff = 1000.f;     // filter cutoff
  float reso   = 0.f;        // filter resonance gain in dB
  float slope  = -3.01;      // spectral slope of the noise
  int   N      = 1000;       // length of sample
  VecF  noise;               // noise sample
  VecF  tgt(N);              // target output in tests
  VecF  outL(N), outR(N);    // for the output signals
  noise = createColoredNoise(N, slope);   // maybe use a sawtooth wave instead
  //rsPlotVector(noise);

  // Create and set up sampler engine:
  SE se;
  float *pSmp = &noise[0];
  se.addSampleToPool(&pSmp, N, 1, fs, "Noise");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter,  60.f);
  se.setRegionSetting(0, 0, PST::cutoffN,         cutoff, 1);
  se.setRegionSetting(0, 0, PST::resonanceN,      reso,   1); // affects only 2nd order modes

  // Test the sampler's 1st order filter modes against the 1-pole-1-zero implementation from RAPT:
  using OPF = RAPT::rsOnePoleFilter<float, float>;
  OPF flt;
  flt.setSampleRate(fs);
  flt.setCutoff(cutoff);
  auto testAgainstOpf = [&](OPF::modes opfMode, Type sfzType, bool plot)
  {
    flt.setMode(opfMode);
    flt.reset();
    for(int n = 0; n < N; n++)
      tgt[n] = flt.getSample(noise[n]);
    se.setRegionSetting(0, 0, PST::filN_type, (float) sfzType, 1);
    return testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-7, plot);

  };
  ok &= testAgainstOpf(flt.LOWPASS_IIT,  Type::lp_6, false);
  ok &= testAgainstOpf(flt.HIGHPASS_MZT, Type::hp_6, false);
  // ToDo: 
  // -Maybe use magnitude matching design formulas, unless they are wildly more expensive (we'll
  //  later want per-sample updates of the freq via an envelope). Or: maybe use cheap polynomial
  //  or rational approximations for the coeff computations (and benchmark(!) them against the 
  //  exact formulas). Or try a fast exp implementation. This will be done in the filter-design
  //  stage.
  // -Add 1st order allpass, low-shelf, high-shelf. These modes are not defined in sfz, maybe 
  //  define an extension - but first figure out, if there are already suitable extensions in use
  //  and if not, consult other sampler engine writers on KVR to find common standard. Maybe ls_6, 
  //  hs_6 and ap_6 could be used - although the 6 would be kinda wrong here because in alpasses,
  //  there's no slope at all and in shelvers there's no asymptotic slope. Maybe lsf_1p, hsf_1p, 
  //  apf_1p would make more sense and be still consistent with the conventions used for 
  //  sfz-defined modes. 
  // -Trying to render a note before setting up the filter type will trigger an error because
  //  the filter's default setting for the type is "Unknown" -> fix this! Use maybe lp_6 as 
  //  default (check, what sfz prescribes as default and use that)

  // Under construction:
  // Test the sampler's 2nd order filter modes against the biquad implementation from RAPT:
  //using BQD = RAPT::rsBiquad<float>;  
  //BQD bqd;
  // ...tbc...
  // ToDo: 
  // -Rename rsBiquad to rsBiquadBase. It has only coeffs and no state. We should have 2 subclasses
  //  that implement DF1/TDF2 and DF2/TDF1 respectively. The DF2/TDF1 variant needs less memory
  //  for the states. Actually, the DF1/TDF2 variant could provide all 4 modes, so maybe that 
  //  should get the name rsBiquad and have 4 getSampleXYZ functions where XYZ is one of DF1,DF2
  //  TDF1,TDF2. The one with lesser memory consumption and lesser flexibility could be called 
  //  rsBiquadCanonical (the implementations with least possible memory are called like that in 
  //  the DSP literature)
  //  https://www.dsprelated.com/freebooks/filters/Four_Direct_Forms.html
  // -Implement a new function in RAPT::rsBiquadDesigner that computes the biquad coeffs
  //  from cutoff-freq (as omega) and resonance gain (either in dB or as raw factor - whatever is
  //  cheaper)

  // Test the sampler's 2nd order filter modes against the SVF implementation from RAPT:
  using SVF = RAPT::rsStateVariableFilter<float, float>;
  using BWC = RAPT::rsBandwidthConverter;
  SVF svf;
  svf.setSampleRate(fs);
  svf.setFrequency(cutoff);
  auto testAgainstSvf = [&](SVF::modes svfMode, Type sfzType, float cutoff, float resoGain, 
    float tol, bool plot)
  {
    float Q = 0.f;
    float resoAmp = RAPT::rsDbToAmp(resoGain);
    switch(svfMode)
    {
    case SVF::LOWPASS:        Q = BWC::lowpassResoGainToQ( resoAmp); break;
    case SVF::HIGHPASS:       Q = BWC::lowpassResoGainToQ( resoAmp); break;
    case SVF::BANDPASS_SKIRT: Q = BWC::bandpassResoGainToQ(resoAmp); break;
    case SVF::BANDREJECT:     Q = BWC::bandpassResoGainToQ(resoAmp); break;
    }
    svf.setMode(svfMode);
    svf.setFrequency(cutoff);
    svf.setGain(Q);
    svf.reset();
    for(int n = 0; n < N; n++)
      tgt[n] = svf.getSample(noise[n]);
    se.setRegionSetting(0, 0, PST::filN_type,  (float) sfzType, 1);
    se.setRegionSetting(0, 0, PST::cutoffN,    cutoff,          1);
    se.setRegionSetting(0, 0, PST::resonanceN, resoGain,        1); 
    return testSamplerNote(&se, 60.f, 127.f, tgt, tgt, tol, plot);
    // Tolerance needs to be a bit higher for 2nd order filters. 1/10^5 corresponds to a relative
    // SNR of 100 dB. It's "relative" in the sense that it is measured against the actual signal 
    // level and not against the maximum possible signal level (i think).
  };
  ok &= testAgainstSvf(svf.LOWPASS,        Type::lp_12,  cutoff, reso, 1.e-5, false);
  ok &= testAgainstSvf(svf.HIGHPASS,       Type::hp_12,  cutoff, reso, 1.e-5, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, cutoff, reso, 1.e-5, false);
  //ok &= testAgainstSvf(svf.BANDREJECT,     Type::br_6_6, true);
  // BRF fails! they look very similar though. Maybe there are different definitions in place for
  // how to intepret the resoGain parameter. It's questionable anyway, if we have implemented
  // to correct behavior as sfz wants it. This needs to be verified! Maybe compare to directly
  // using an RBJ biquad. And/or maybe try using the filters with double precision. Maybe its a 
  // numerical issue - although the error is visible, so that's perhaps a bit too much for roundoff
  // errors.

  // Test numerical stability using a bandpass filter with a high resonance. Bandpasses have 
  // slightly higher Q than low- or highpasses with the same resonance gain. We test it for 
  // frequencies from all the way up to all the way down. For the lower cutoffs, we need higher 
  // tolerances:
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 22050.f,   40.f, 1.e-4, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 20000.f,   40.f, 1.e-4, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 10000.f,   40.f, 1.e-4, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,  1000.f,   40.f, 1.e-4, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,   100.f,   40.f, 1.e-3, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,    10.f,   40.f, 1.e-3, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.1f,  40.f, 1.e-3, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.01f, 40.f, 1.e-3, false);
  //ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.0f,  40.f, 1.e-3, true);
  // todo: 
  // -maybe write a helper function taht goes through these checks als for lowpass and highpass
  // -maybe we should allow cutoff up to 30000 as in the EQ

  // ToDo
  // -Cutoff=0 does not yet work - the svf produces silence and the sampler goes into bypass. Maybe
  //  we should do something different in this case. For a bandpass, it seems to make sense that 
  //  the limiting case is a lowpass. For a highpass, the limiting case should indeed be a bypass.
  //  for a lowpass, we may see the resonance peak at DC?
  // -Try a series connection of waveshaper and filter in both possible orders - check if the order
  //  is indeed determined by the first opcode that applies to the given dsp as it should be
  // -Maybe change the default filter type to make it work even the fil_type opcode is missing.
  //  Currently, we trigger an error when e.g.
  //    se.setRegionSetting(0, 0, PST::FilterType, (float) Type::lp_6);
  //  is not being called before playing a note because then, the filter is in "Unknown" mode.
  // -Reduce the boilerplate by defining a lambda taking e.g. svf.LOWPASS and Type::lp_12
  // -Implement mapping between filter-types from our enum and their sfz-strings.
  // -Use that to allow defining filtering in the sfz.
  // -SFZ 1 supports: lpf_1p, hpf_1p, lpf_2p, hpf_2p, bpf_2p, brf_2p
  // -Support the filter opcodes in the SFZ parser (type and cutoff - later reso). -> add a unit
  //  test to the "...FileIO" test which tests that.
  // -Add more filter modes. Make sure, they behave the same as in sfz+ and/or sfizz with respect
  //  to parameter settings.
  // -I think, we should generally use formulas that match the magnitude of the analog filter. I 
  //  think, the algorithms by Martin Vicanek are quite nice: use impulse-invariance for the poles
  //  and magnitude matching for the zeros. But to get the ball rolling, just use impulse 
  //  invariance.
  // -Try to replicate the behavior of filters in sfz+ or sfizz
  // -Because the resonance is given in dB in sfz files, it could be re-used as gain for 
  //  bell-filters. Then, an additional variable could be the bandwidth (in octaves). but figure 
  //  out, how sfz does it. ...if it has specs for bell-filters, that is. If not, maybe start a 
  //  thread on KVR about extending the sfz specs suitably.
  // -Maybe include a slope filter. Cutoff may be re-interpreted as unit-gain freq, resonance,
  //  which is in db, may be re-interpreted as dB/oct
  // -For the mapping between the resonance parameter in the sfz spec and filter Q in the design
  //  algorithms, see the biquadResoGain() experiment.

  rsAssert(ok);
  return ok;
}

bool samplerWaveShaperTest()
{
  bool ok = true;

  using VecF  = std::vector<float>;
  using SE    = rosic::Sampler::rsSamplerEngineTest;
  using PST   = rosic::Sampler::Opcode;
  //using WS    = rosic::Sampler::rsSamplerWaveShaper;
  using Shape = rosic::Sampler::rsSamplerWaveShaper::Shape;

  // Create a sinewave as example sample:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of (co)sinewave sample
  int   N  = 500;    // length of (co)sinewave sample
  VecF sin440(N);    // sine wave
  VecF tgt(N);       // target output in tests
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
    sin440[n] = sinf(w*n);
  //rsPlotVector(sin440);

  // Waveshaper settings:
  float drive  = 4.0f;
  //float postGain = 0.5f;
  float dcOffset = 0.0;
  Shape shape    = Shape::Tanh;

  // Create target signal:
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(drive * sin440[n]);
  //rsPlotVector(tgt);
  // todo: use a cheap approximation to tanh and/or rsTanh (based on exp), maybe use exp based on
  // 2Dat's code (somewher on the kvr forum)...also implement mystran's "random cheap sigmoid"
  // maybe use the quake alo for the inverse sqrt...but first things first...and the first thing is
  // to put the infrastructure in place

  // Create and set up sampler engine:
  SE se;
  float *pSmp = &sin440[0];
  se.addSampleToPool(&pSmp, N, 1, fs, "Sine440");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f);
  se.setRegionSetting(0, 0, PST::DistShape, float(shape));
  se.setRegionSetting(0, 0, PST::DistDrive, drive);
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-7, false);
  rsAssert(ok);
  // ToDo:
  // -Try different shapes, use different sets of parameters, use DC, postGain, etc.
  // Cosmetics:
  // -Maybe drag out RegionPlayer from rsSamplerEngine
  // -maybe let testSamplerNote take a plotMode parameter which can be: 0: never plot, 1: always 
  //  plot, 2: plot when failed
  // -make a class SignalProcessorPool, let the engine maintain such a pool as member 
  //  -this class should have a function grabProcessor(SignalProcessorType type) that returns
  //   a pointer to a processor of the desired type or a nullptr if no such processor is available
  //   anymore...or maybe a pointer to a dummy-processor?

  // -would be nice, if we could wrap the namespace Sampler around the includes in rosic.h/cpp
  //  but for that, we first need to move all the namespace rosic stuff there, too
  // -maybe try to do this task with a python script, see
  //  https://stackoverflow.com/questions/1120707/using-python-to-execute-a-command-on-every-file-in-a-folder
  //  the answer using pathlib seems to be the simplest way to do it
  // -before that, merge the current update_juce branch to master
  // -create a new branch develop where we can do all this stuff, see also:
  //  https://newbedev.com/how-to-iterate-over-files-in-a-given-directory
  //  https://www.codegrepper.com/code-examples/python/loop+through+all+files+in+a+directory+python


  // ToDo: 
  // -Figure out, if there is already an existing extension to sfz that defines opcodes for 
  //  waveshaping. If so, use these. Otherwise define opcodes: dist_drive, dist_shape, dist_gain
  //  ..or maybe distort_ or distortion_

  rsAssert(ok);
  return ok;
}

bool samplerDspChainTest()
{
  bool ok = true;

  using VecF  = std::vector<float>;     // vector of sample values in RAM
  using SE    = rosic::Sampler::rsSamplerEngineTest;
  using PST   = rosic::Sampler::Opcode;
  using Type  = rosic::Sampler::FilterType;
  using Shape = rosic::Sampler::rsSamplerWaveShaper::Shape;

  // Setup:
  float fs      = 44100.f;  // sample rate
  float slope   = -3.01;    // spectral slope of the noise
  float cutoff1 = 2000.f;   // lowpass cutoff
  float cutoff2 =  500.f;   // highpass cutoff
  int   N       =  500;     // length of sample

  // Create a pinkish noise as example sample:
  VecF  noise;              // noise sample
  VecF  tgt(N);             // target output in tests
  VecF  outL(N), outR(N);   // for the output signals
  noise = createColoredNoise(N, slope);

  // Create target signal using a filter from RAPT:
  using OPF = RAPT::rsOnePoleFilter<float, float>;
  OPF flt;
  flt.setSampleRate(fs);
  flt.setCutoff(cutoff1);
  flt.setMode(flt.LOWPASS_IIT);
  for(int n = 0; n < N; n++)
    tgt[n] = flt.getSample(noise[n]);
  flt.setCutoff(cutoff2);
  flt.setMode(flt.HIGHPASS_MZT);
  flt.reset();
  for(int n = 0; n < N; n++)
    tgt[n] = flt.getSample(tgt[n]);
  //rsPlotVector(tgt);

  // Build a sampler patch using DSP chain containing two filters, a lowpass and a highpass, let it
  // produce a note and check it against the target signal:
  SE se;
  float *pSmp = &noise[0];
  se.addSampleToPool(&pSmp, N, 1, fs, "Noise");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter,  60.f);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::lp_6, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff1, 1);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::hp_6, 2);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff2, 2);
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);

  // Add a waveshaper and after that a 3rd (lowpass) filter into the chain, such that 
  // the chain is now: LPF -> HPF -> WS -> LPF:
  float drive    = 4.0f;
  Shape shape    = Shape::Tanh;
  float cutoff3  = 1000.f; 
  se.setRegionSetting(0, 0, PST::DistShape, float(shape));
  se.setRegionSetting(0, 0, PST::DistDrive, drive);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::lp_6, 3);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff3, 3);

  // Create new target signal and run test:
  flt.setMode(flt.LOWPASS_IIT);
  flt.setCutoff(cutoff3);
  flt.reset();
  for(int n = 0; n < N; n++)
    tgt[n] = flt.getSample(tanh(drive * tgt[n]));
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);

  // Updates the target signal according to new values of cutoff1,2,3 and drive
  auto updateTgt = [&]()
  {
    flt.setMode(flt.LOWPASS_IIT);
    flt.setCutoff(cutoff1);
    flt.reset();
    for(int n = 0; n < N; n++)
      tgt[n] = flt.getSample(noise[n]);
    flt.setMode(flt.HIGHPASS_MZT);
    flt.setCutoff(cutoff2);
    flt.reset();
    for(int n = 0; n < N; n++)
      tgt[n] = flt.getSample(tgt[n]);
    flt.setMode(flt.LOWPASS_IIT);
    flt.setCutoff(cutoff3);
    flt.reset();
    for(int n = 0; n < N; n++)
      tgt[n] = flt.getSample(tanh(drive * tgt[n]));
  };

  // Set the cutoffs of the filters to a different values:
  cutoff2 = 100.f;
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff2, 2);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);
  cutoff1 = 4000.f;
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff1, 1);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);
  cutoff3 = 3000.f;
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff3, 3);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);
  drive = 8.0;
  se.setRegionSetting(0, 0, PST::DistDrive, drive);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);

  // OK - we currently have 3 filters with a waveshaper between filter 2 and filter 3. Now let's
  // check the behavior when we try to add a 5th filter without adding a 4th before. The desired 
  // behavior should be, that the DSP chain actually has 5 filters in it but the 4th one is in its
  // default setting which means that it is effectively bypassed. For the added filter at index 5,
  // we don't specify the filter type in which case it should default to lpf_2p.
  float cutoff5 = 5000.f; 
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff5, 5);
  ok &= se.getRegion(0, 0)->getNumProcessors() == 6;  // 5 filters + 1 waveshaper
  updateTgt();
  using SVF = RAPT::rsStateVariableFilter<float, float>;
  SVF svf;
  svf.setSampleRate(fs);
  svf.setFrequency(cutoff5);
  float G = 1.f / sqrt(2.f);
  svf.setGain(G);
  for(int n = 0; n < N; n++)
    tgt[n] = svf.getSample(tgt[n]);
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);

  // ToDo: 
  // -maybe write a test that creates a random dsp chain programmatically using lots of filters and 
  //  waveshapers, set the filter parameters in random order, like type3, cutoff1, type2, reso4,
  // -test to include filters into an actual sfz instrument for playing, i.e. if the parsing works
  //  right...well - it certainly won't for more than one filter because this is not yet 
  //  implemented
  // -test what happens, if 
  //  -we pass index = 0, -1, 4, -2
  //  -we pass the parameters in interleaved order
  //  -we have not enough DSPs available (the region should not play at all in this case)

  rsAssert(ok);
  return ok;
}

bool samplerEqualizerTest()
{
  // Tests the support of the eqN_freq, eqN_bw, eqN_gain opcodes where N=1,2,3 in the original sfz
  // spec, but we want to allow arbitrary N.

  bool ok = true;

  using Vec    = std::vector<float>;
  using SE     = rosic::Sampler::rsSamplerEngineTest;
  using OC     = rosic::Sampler::Opcode;
  using Type   = rosic::Sampler::FilterType;
  using Region = rosic::Sampler::rsSamplerData::Region;


  using namespace RAPT;

  // Setup:
  float fs    = 44100.f;     // sample rate
  float gain1 =    24.0f;    // gain of eq1
  float freq1 =  1500.f;
  float bw1   =     0.3f;
  float gain2 =    -2.0f;    // gain of eq2
  float gain3 =     5.0f;    // gain of eq3
  int   N     =   500;       // length of sample
  // Don't use bandwidth < 0.25 or > 6 because the reference filter doesn't support them and clips
  // the values

  // Create a pinkish noise as example sample:
  Vec noise;                 // noise sample
  Vec tgt(N);                // target output in tests
  Vec outL(N), outR(N);      // for the output signals
  noise = createColoredNoise(N, -3.01f);

  // Helper function that takes an input x, applies an arbitrary number of eq stages to it and 
  // writes the result into the output y:
  auto applyEqs = [&](const Vec& x, Vec& y, const Vec& gains, const Vec& freqs, const Vec& bws)
  {
    // Input sanity checks:
    rsAssert(rsAreSameSize(gains, freqs, bws));
    rsAssert(rsAreSameSize(x, y));

    // Create filter object (we use only the 1st stage):
    using  TPF = rosic::TwoPoleFilter;
    using DTPF = rosic::DualTwoPoleFilter;
    DTPF flt;
    flt.setSampleRate(fs);
    flt.setMode1(TPF::PEAK);

    // Apply the filter stages one after another:
    int numStages = (int)gains.size();
    rsCopy(x, y);
    for(int i = 0; i < numStages; i++)
    {
      flt.setFrequency1(freqs[i]);
      flt.setGain1(gains[i]);
      flt.setBandwidth1(bws[i]);
      flt.reset();
      for(int n = 0; n < N; n++)
        y[n] = flt.getSample(y[n]);
    }
  };

  // Create the basic sampler patch with no eq yet:
  SE se;
  Region* r;
  float *pSmp = &noise[0];
  se.addSampleToPool(&pSmp, N, 1, fs, "Noise");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, OC::PitchKeyCenter,  60.f);

  // Test the simplest case: use only eq1 and specify all 3 parameters:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq1, 1);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw1, 1);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 1;
  applyEqs(noise, tgt, { gain1 }, { freq1 }, { bw1 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-5, false);
  // Tolerance needs to be even higher than in the filter tests

  // Now we set only the gain of eq1. The freq and bandwidth should default to 50 Hz and 1 oct:
  r->clearSettings();
  se.setRegionSetting(0, 0, OC::PitchKeyCenter, 60.f);
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  ok &= r->getNumProcessors() == 1;
  applyEqs(noise, tgt, { gain1 }, { 50 }, { 1 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-3, false);
  // Whoa - here we need a really high tolerance! OK, we have a moderately high Q and high gain and
  // a low frequency. Settings which are prone to numeric errors. But still, this is actually quite 
  // bad. Maybe we need to use double precision for filters and equalizers indeed. Or maybe use 
  // double for the coefficient calculation at least. Wait - no - the Q is actually rather low. And 
  // sfz allows very high Q settings

  // Now we set only the gain of eq3. The desired behavior is that we actually get 3 filter stages 
  // in the dsp chain but the first two are in neutral setting. We don't specify the center 
  // frequency or bandwidth. Therefore, the default values should be used which are 5 kHz and 1 
  // octave:
  r->clearSettings();
  se.setRegionSetting(0, 0, OC::PitchKeyCenter, 60.f);
  se.setRegionSetting(0, 0, OC::eqN_gain, gain3, 3);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 3;
  applyEqs(noise, tgt, { gain3 }, { 5000.f }, { 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6, false);

  // Add band 2. This has a default freq of 500Hz:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain2, 2); 
  applyEqs(noise, tgt, { gain2, gain3 }, { 500.f, 5000.f }, { 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4, false);

  // Add band 1. This has a default freq of 50Hz:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { 50.f, 500.f, 5000.f }, { 1.f, 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-3, false);
  // Here, we also need a high tolerance

  // Change settings of band 1 to something more benign from a numeric point of view:
  freq1 = 1500;
  gain1 = 5.f;
  bw1   = 0.5f;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq1, 1);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw1, 1);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, 500.f, 5000.f }, { bw1, 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4, false);

  // Change settings of band 2:
  gain2 = -6;
  float freq2 = 2500;
  float bw2   = 0.3;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain2, 2);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq2, 2);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw2, 2);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, freq2, 5000.f }, { bw1, bw2, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4, false);

  // Change settings of band 3:
  gain3 = 7;
  float freq3 = 8000;
  float bw3   = 2.2;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain3, 3);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq3, 3);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw3, 3);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, freq2, freq3 }, { bw1, bw2, bw3 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4, false);

  // Add a 5th band with no 4th in between, so we get a dspChain with 5 eqs but the 4th is neutral.
  // The new eq should be at its default freq of 1000:
  float gain5 = -3;
  float bw5   = 1.2;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain5, 5);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw5, 5);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 5;
  applyEqs(noise, tgt, { gain1, gain2, gain3, gain5 }, { freq1, freq2, freq3, 1000 }, 
    { bw1, bw2, bw3, bw5 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4, false);

  // ToDo:
  // -Check that it all works also when we have one or more filters (via the actual filter opcode)
  //  in the chain, i.e. make sure that there isn't any sort of bad interference between equalizers
  //  nad filters. Check also, what happens, if the eq opcodes are defined before the filter 
  //  opcodes.
  // -We currently use the Orfanidis design, but that shouldn't be the final word. First of all,
  //  we may have to support center frequencies above the Nyquist limit because sfz specifies the
  //  frequency range to be 0..30000. Maybe the Vicanek design could be a good choice as well 
  //  especially, when we want to make the eq modulatable later (because the formulas are 
  //  simpler). ...but these questions can be postponed to a filter design stage.
  // -Implment and test saving and parsing of eq parameters

  rsAssert(ok);
  return ok;
}


bool samplerProcessorsTest()
{
  bool ok = true;

  // For inspection in the debugger. We want to keep the sizes of these DSP objects small because 
  // we'll potentially have to pre-allocate a lot of them when a patch is loaded:
  int size;
  size = sizeof(rosic::Sampler::rsSamplerFilter);       // 64
  size = sizeof(rosic::Sampler::rsSamplerWaveShaper);   // 24
  //size = sizeof(SP::Filter);
  //size = sizeof(SP::WaveShaper);
  // later move this into a (yet to be written) benchmark testbed

  ok &= samplerFilterTest();      // tests the different filter modes
  ok &= samplerWaveShaperTest();  // tests the wvashaping DSP module
  ok &= samplerDspChainTest();    // uses multiple filters and a waveshaper in between
  ok &= samplerEqualizerTest();

  // ToDo:
  // -Implement more DSP modules: echo, vibrato, flanger, phaser, chorus, etc., 
  //  ...delay based algorithms could become a memory-hog when we need to pre-allocate many of 
  //  them. Maybe let's stay away from them for the moment. How about freq-shifting?

  rsAssert(ok);
  return ok;
}

bool samplerEngineUnitTest()
{
  bool ok = true;

  // The new test that is currently under construction:
  //ok &= samplerEqualizerTest();

  // The tests, that already pass and are supposed to continue to do so:
  ok &= samplerDataUnitTest();
  ok &= samplerEngineUnitTest1();
  ok &= samplerEngineUnitTestFileIO();
  ok &= samplerProcessorsTest();
  ok &= samplerEngine2UnitTest(); 

  // ToDo:
  // -Figure out and implement the correct filter design formulas for lpf_2p, hpf_2p, bpf_2p and 
  //  the equalizers. It's rather interesting that for the filters, the range is 0..fs/2 and for
  //  equalizers it's 0..30kHz
  // -Test filters with the most extreme Q settings, also at low frequnecies. That's where we 
  //  should expect numerical precision issues, especially when using single precision. Maybe using
  //  an SVF implementation can help against this.



  //rsAssert(ok, "samplerEngineUnitTest failed");
  return ok;
}

/*

ToDo:
-maybe wrap all the sampler-related code into a a sub-namespace rosic::Sampler
 -the, the names of the classes may get rid of the "Sampler" part - that would then be redundant, 
  i.e. classes are named Data, Engine, SignalProcessor, Modulator, etc.
 -the SignalProcessor and Modulator class could be dragged out of the Engine
 -maybe the RegionPlayer class, too

 resources:
 http://sfzformat.com/legacy/           basic opcode reference
 http://drealm.info/sfz/plj-sfz.xhtml   explanantions and some unofficial opcodes

maybe open a thread at kvr..something like SFZ: interpreting and extending the format specification

*/