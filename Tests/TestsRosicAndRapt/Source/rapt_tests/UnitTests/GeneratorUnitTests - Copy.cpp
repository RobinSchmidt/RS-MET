
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
  int N = L.size(); rsAssert((int)R.size() == N);
  using AT = RAPT::rsArrayTools;

  T t = (pan + 1) * T(0.5);    // -1...+1  ->   0...1
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
std::vector<T> rsApplyResampling(std::vector<T>& x, T readSpeed)
{
  int Nx= (int)x.size();
  int Ny= (int) ceil(Nx / readSpeed);
  std::vector<T> y(Ny);
  for(int n = 0; n < Ny; n++)
    y[n] = getSampleAt(x, T(readSpeed*n));
  return y;
}

std::vector<float> rsApplyFilter(std::vector<float>& x, rosic::Sampler::FilterType type, 
  float cutoff, float sampleRate, float resonance)
{
  int N = (int)x.size();
  std::vector<float> y(N);
  using FilterDsp  = rosic::Sampler::Filter;
  using FilterCore = rosic::Sampler::FilterCore;
  FilterCore flt;
  float omega = float(2*PI) * cutoff / sampleRate;
  flt.setupCutRes(FilterDsp::convertTypeEnum(type), omega, resonance);
  for(int n = 0; n < N; n++) {
    y[n] = x[n]; float dummy;
    flt.processFrame(&y[n], &dummy); }
  return y;
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
    y[n] = (float)flt.getSample(y[n]);

  return y;
}
// move to test tools

// Helper function to add a single region for the given sample to the engine. The region is added
// to the first group, which is added if not already there:
void addSingleSampleRegion(rosic::Sampler::rsSamplerEngine* se,
  const std::vector<float>& sample, float keyCenter = 60.f, double sampleRate = 44100)
{
  using PST  = rosic::Sampler::Opcode;
  const float *pSmp = &sample[0];
  int si = se->addSampleToPool((float**) &pSmp, (int)sample.size(), 1, sampleRate, "Sample");
  if(se->getNumGroups() == 0)
    se->addGroup();
  int ri = se->addRegion(0);
  se->setRegionSample( 0, ri, si);
  se->setRegionSetting(0, ri, PST::PitchKeyCenter, keyCenter, -1);
  // ToDo: try to get rid of casting ways the const in addSampleToPool((float**) &pSmp,...). 
  // addSampleToPool does not modify anything - make it const correct...but that my need ugly
  // and confusing syntax in the function declaration
}

// Helper funtion to set up the sampler engine with single a sinewave region:
void setupForSineWave(rosic::Sampler::rsSamplerEngine* se)
{
  int N = 2048;
  std::vector<float> sineWave(N);
  for(int n = 0; n < N; n++)
    sineWave[n] = (float) sin(2.0*PI*n/N);
  se->clearInstrument();
  addSingleSampleRegion(se, sineWave, 21, 56320);
}

//=================================================================================================

bool samplerDataTest()
{
  bool ok = true;

  //using SD = SfzInstrument;
  using SFZT = rosic::Sampler::SfzCodeBook;
  using SD   = rosic::Sampler::SfzInstrument;
  using PST  = rosic::Sampler::Opcode;

  SFZT::createInstance();
  // Normally, this is supposed to be done in the constructor of rsSamplerEngine and objects of 
  // type SfzInstrument are supposed to live only inside the engine. But here in the test, we 
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
  // crash: in SfzInstrument::setFromSFZ, j0 = 18446744073709551615

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
  d3.setGroupSetting(0,     PST::volumeN, -3.f, 1);
  d3.setRegionSetting(0, 0, PST::volumeN, -5.f, 1);

  // Test sfz generation and parsing:
  sfz = d3.getAsSFZ();
  d2.setFromSFZ(sfz);  // crashes, if we uncomment the setInstrumentSetting call
  ok &= d2 == d3;

  // Set up a feq settings that involve DSPs:
  d3.setRegionSetting(0, 0, PST::cutoffN,  1000.f, 1);
  d3.setRegionSetting(0, 0, PST::eqN_gain,    6.f, 3);  // should prodcue 3 eq bands
  sfz = d3.getAsSFZ(); d2.setFromSFZ(sfz); ok &= d2 == d3;

  d1 = d3;
  ok &= d1 == d3;





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

bool testSamplerOutput(rosic::Sampler::rsSamplerEngine* se,
  const std::vector<float>& targetL, const std::vector<float>& targetR,
  float tol = 0.f, bool plot = false)
{
  int N = (int) targetL.size();
  rsAssert((int)targetR.size() == N);
  std::vector<float> outL(N), outR(N);
  for(int n = 0; n < N; n++)
    se->processFrame(&outL[n], &outR[n]);
  using AT   = RAPT::rsArrayTools;
  float errL = AT::maxDeviation(&outL[0], &targetL[0], N);
  float errR = AT::maxDeviation(&outR[0], &targetR[0], N);
  if(plot)
    rsPlotVectors(targetL, targetR, outL, outR, targetL-outL, targetR-outR);
  return errL <= tol && errR <= tol;
}

bool testSamplerNote(rosic::Sampler::rsSamplerEngine* se, float key, float vel, 
  const std::vector<float>& targetL, const std::vector<float>& targetR, 
  float tol = 0.f, bool plot = false)
{
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;
  se->handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
  return testSamplerOutput(se, targetL, targetR, tol, plot);
};
// maybe have a bool resetBefore that optionally resets the engine before playing...but maybe it's
// better when the caller does the reset directly, if desired - it's not much longer but clearer
// factor out a getSamplerOutput(rosic::rsSamplerEngine* se, float key, float vel, 
// const std::vector<float>& targetL, const std::vector<float>& targetR, bool plot) function

bool samplerRegionPlayerTest()
{
  // Tests the basic functionality of RegionPlayer

  bool ok = true;

  using VecF = std::vector<float>;     // vector of sample values in RAM
  using AT   = RAPT::rsArrayTools;
  using SD   = rosic::Sampler::SfzInstrument;
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
  rc = se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f, -1); // key = 69 is A4 at 440 Hz
  ok &= rc == RC::success;

  // Try retrieving Group and Region pointers. The pointers need to be const because we are not 
  // allowed to modify these things behind the back of the SamplerEngine. All mutations have to use
  // member functions of the engine:
  const SE::Group*  grp = se.getGroup(0);
  const SE::Region* rgn = se.getRegion(0, 0);
  ok &= rgn->getGroup() == grp;
  //rgn->setSetting(SD::PlaybackSetting(PST::LoKey, 0.f, -1); // should not compile

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
  se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f, -1);   // cosine is also at A4 = 440 Hz
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
  se.setRegionSetting(gi, 0, PST::panN, -100.f, 1);           // pan sine to hard left
  se.setRegionSetting(gi, 1, PST::panN, +100.f, 1);           // pan cosine to hard right
  //se.setRegionSetting(gi, 0, PST::positionN, -100.f, 1);        // pan sine to hard left
  //se.setRegionSetting(gi, 1, PST::positionN, +100.f, 1);        // pan cosine to hard right
  se.handleMusicalEvent(Ev(EvTp::noteOn, 69.f, 127.f));
  for(int n = 0; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);
  //rsPlotVectors(outL, outR); 
  ok &= outL == (2.f * sin440);
  ok &= outR == (2.f * cos440);
  ok &= se.getNumActiveLayers() == 0;
  // ToDo:
  // -Do the same test with positionN instead of panN
  // -Figure out, if 100% is really the right default value for the width. The sfz spec says so, 
  //  but that would be a bad default. It should be 100 because that leaves the sample unchanged.

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
  se.setRegionSetting(0, 0, PST::panN,   0.f, 1);        // back to center, makes testing easier
  se.setRegionSetting(0, 0, PST::Delay,  delaySeconds, -1);
  VecF tgt = sin440;
  float tol = 1.e-7f;  // ~= 140 dB SNR
  rsApplyDelay(tgt, delaySamples);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, tol, false); 

  // Test some corner cases:

  // Test playing a sample of length 1:
  float one = 1.f; float* pOne = &one; float** ppOne = &pOne;
  se.clearInstrument();
  se.addSampleToPool(ppOne, 1, 1, 44100.f, "RAM");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample(0, 0, 0);
  float L, R;

  // The unit impulse played at its original pitch should just produce a single sample of 1:
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= se.getNumActiveLayers() == 1;
  se.processFrame(&L, &R);
  ok &= L == 1.f && R == 1.f;
  ok &= se.getNumActiveLayers() == 0;
  se.processFrame(&L, &R);
  ok &= L == 0.f && R == 0.f;

  // Playing it two octaves below it's keycenter should produce 4 samples: 1.0, 0.75, 0.5, 0.25 due 
  // to the linear interpolation:
  se.handleMusicalEvent(Ev(EvTp::noteOn, 36.f, 127.f));
  ok &= se.getNumActiveLayers() == 1;
  se.processFrame(&L, &R);  ok &= L == 1.0f  && R == 1.0f;
  se.processFrame(&L, &R);  ok &= L == 0.75f && R == 0.75f;
  se.processFrame(&L, &R);  ok &= L == 0.5f  && R == 0.5f;
  se.processFrame(&L, &R);  ok &= L == 0.25f && R == 0.25f;
  ok &= se.getNumActiveLayers() == 0;
  se.processFrame(&L, &R);  ok &= L == 0.f && R == 0.f;

  // Test playing a sample of length 0:
  // ...or should we? is this something we should be able to handle? What should actually be the
  // desired state of the samplePool after adding a sample of size 0? should we indeed add an
  // object with zero data or should we just don't add anything? ..I'm not yet sure...but perhaps
  // 



  // ToDo:
  // -write a performance test for the sampler
  // -switch to an int+float representation of the current sample position and increment and check, 
  //  if this improves performance...even if not, it's still better because it doesn't lose 
  //  precision for later samples
  // -implement and test better realtime resampling (linear interpolation at first, later cubic and
  //  sinc, maybe some sort of "Elephant" interpolation, too - although, they are supposed to work 
  //  with 2x oversampling) 
  //  -we need a double for the sample-time and an increment...but later, that increment shall be
  //   modulated by pitch-env and -lfo, or maybe use and int and a a float to represent sampleTime
  //   and increment -> do benchmarks, which is faster
  //  -should the played note affect the delay?...nah - i don't think so. maybe sfz had an opcode 
  //   for controlling this? in some situations, that may makes sense (for creating comb-filter 
  //   effects), in others not so much


  rsAssert(ok);
  return ok;
}

bool samplerBusModeTest()
{
  // We test the extended functionality of rsSamplerEngine2, in particular, the signal routing 
  // through per group DSP processes.
  // rename to samplerBusModeTest

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
  rc = se.setRegionSetting(gi, ri, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;

  //---------------------------------------------------------------------------
  // Test accumulation of amp setting:

  // Set up volume opcode for region and group:
  float regionAmp = 0.5f;
  float groupAmp  = 0.25f;
  float instrAmp  = 1.5f;
  float tol       = 1.e-7f;  // why can't we use 0 tolerance?
  se.setRegionSetting(0, 0, PST::volumeN, rsAmpToDb(regionAmp), 1);
  se.setGroupSetting( 0,    PST::volumeN, rsAmpToDb(groupAmp),  1);
  se.setInstrumentSetting(  PST::volumeN, rsAmpToDb(instrAmp),  1);

  // In the default setting, the regionAmp should override the instrument and the group setting,
  // so the produced output should only have the region volume applied:

  // In bus-mode, we want to see all 3 settings applied:
  se.setBusMode(true);
  tgt = instrAmp*groupAmp*regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, tol, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0;

  // In default mode (i.e. non-busMode), we want to see only the region setting applied:
  se.setBusMode(false);
  tgt = regionAmp*sin440;
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 0.0, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0;

  //---------------------------------------------------------------------------
  // Test accumulation of pan setting:

  se.clearAllSfzSettings();                                   // remove all the amp settings
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1);  // restore the rootkey setting
  ok &= rc == RC::success;
  float regionPan = 10.f;   // slightly right, pan range is -100...+100
  float groupPan  = 20.f;
  float instrPan  = 30.f;
  se.setRegionSetting(0, 0, PST::panN, regionPan, 1);
  se.setGroupSetting( 0,    PST::panN, groupPan,  1);
  se.setInstrumentSetting(  PST::panN, instrPan,  1);

  // We want to see only the region pan:
  se.setBusMode(false);
  tgtL = tgtR = sin440;
  rsApplyPan(tgtL, tgtR, regionPan/100);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgtL, tgtR, 1.e-6f, false); 

  // Now we want to see region, group and instrument pan combined:
  se.setBusMode(true);
  tgtL = tgtR = sin440;
  rsApplyPan(tgtL, tgtR, regionPan/100);
  rsApplyPan(tgtL, tgtR, groupPan /100);
  rsApplyPan(tgtL, tgtR, instrPan /100);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgtL, tgtR, 1.e-6f, false);

  // ToDo: test it also with constant power pan rule

  //---------------------------------------------------------------------------
  // Test delay accumulation:

  se.clearAllSfzSettings();                                   // remove all the amp settings
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1);  // restore the rootkey setting
  int regionDelay = 10;   // in samples - todo: use float
  int groupDelay  = 20;
  int instrDelay  = 40;
  se.setRegionSetting(0, 0, PST::Delay, regionDelay / fs, -1);
  se.setGroupSetting( 0,    PST::Delay, groupDelay  / fs, -1);
  se.setInstrumentSetting(  PST::Delay, instrDelay  / fs, -1);

  // We want to see only the region delay:
  se.setBusMode(false);
  tgt = sin440;
  rsApplyDelay(tgt, regionDelay);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7f, false);
  ok &= se.getNumActiveLayers() == 1;        // it's still playing due to the delay
  ok &= se.getNumActiveGroupPlayers() == 0;  // no group player is/was used due to settings

  // Now we want to see region, group and instrument delay combined:
  se.setBusMode(true);
  tgt = sin440;
  rsApplyDelay(tgt, regionDelay);
  rsApplyDelay(tgt, groupDelay);
  rsApplyDelay(tgt, instrDelay);
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7f, false);

  //---------------------------------------------------------------------------
  // Test offset accumulation:

  float regionOffset = 10;   // "offset" opcode in samples
  float groupOffset  = 20;
  float instrOffset  = 40;
  se.clearAllSfzSettings();
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1);
  se.setRegionSetting(0, 0, PST::Offset, regionOffset, -1);
  se.setGroupSetting( 0,    PST::Offset, groupOffset, -1);
  se.setInstrumentSetting(  PST::Offset, instrOffset, -1);

  // We want to see only the region offset:
  se.setBusMode(false);
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset);  // offset is like a negative delay
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7f, false);
  ok &= se.getNumActiveLayers() == 0;
  ok &= se.getNumActiveGroupPlayers() == 0; 

  // We want to see region, group and instrument offset:
  se.setBusMode(true);
  tgt = sin440;
  rsApplyDelay(tgt, -regionOffset); 
  rsApplyDelay(tgt, -groupOffset); 
  rsApplyDelay(tgt, -instrOffset); 
  ok &= testSamplerNote(&se, 69.f, 127.f, tgt, tgt, 1.e-7f, false);

  //---------------------------------------------------------------------------
  // Test offset and delay (but only for the region setting):

  se.clearAllSfzSettings();
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1);
  se.setRegionSetting(0, 0, PST::Offset, regionOffset, -1);
  se.setRegionSetting(0, 0, PST::Delay,  regionDelay / fs, -1);
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
  //regionTune = groupTune = instrTune = 0;  // for debug
  se.clearAllSfzSettings();                               // remove all the amp settings
  se.setGroupSetting( 0,    PST::PitchKeyCenter, 50.f, -1);   // should always be overriden
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1);   // restore the rootkey setting
  se.setRegionSetting(0, 0, PST::Transpose, regionTrans, -1);
  se.setGroupSetting( 0,    PST::Transpose, groupTrans, -1);
  se.setInstrumentSetting(  PST::Transpose, instrTrans, -1);
  se.setRegionSetting(0, 0, PST::Tune,      regionTune, -1);
  se.setGroupSetting( 0,    PST::Tune,      groupTune, -1);
  se.setInstrumentSetting(  PST::Tune,      instrTune, -1);

  float pitch;  
  tol = 0.02f; 
  // The error of the pitch estimation is around 2 cents...that's quite a large error actually.
  // Try to improve this at some point!I think, the algo should give better results!

  // Test override mode. We expect to see only the region transpose and tune. That's 
  // 69 + 1 + 10/100 = 70.1
  se.setBusMode(false);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVectors(outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 70.1f, tol);
  // its 69.08 ~= 69.1 so it's about 1 semitone too low


  // Now we want to see region, group and instrument transpose combined. That's
  // 69 + 1 + 2 + 3 + 10/100 + 20/100 + 30/100 = 75.6
  se.setBusMode(true);
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  pitch = rsEstimateMidiPitch(outL, fs);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 75.6f, tol);
  // pitch = 74.58 - again 1 semitone too low

  // ToDo: check what happens when region setting and/or group setting and/or instrument setting
  // are removed and if it behaves as expected. we may need functionality to delete particular 
  // opcodes like se.removeRegionSetting(0, 0, PST::Tune) etc.

  // Remove the transpose setting for the region. We expect to see a combination of transpose
  // settings of instrument and group and a combination of the tune settings of all 3, so the pitch
  // should be: 69 + 2 + 3 + 0.1 + 0.2 + 0.3 = 74.6
  ok &= se.removeRegionSetting(0, 0, PST::Transpose, -1) == RC::success;
  ok &= se.removeRegionSetting(0, 0, PST::Transpose, -1) == RC::nothingToDo;
  //se.reset();  // what happens if we don't reset? seems to make no difference
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.6f, tol);

  // Remove tune setting from group. p = 69 + 2 + 3 + 0.1 + 0.3 = 74.4
  ok &= se.removeGroupSetting(0, PST::Tune, -1) == RC::success;
  ok &= se.removeGroupSetting(0, PST::Tune, -1) == RC::nothingToDo;
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVector(outL);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.4f, tol);

  // Remove tune setting from instrument. p = 69 + 2 + 3 + 0.1 = 74.1
  ok &= se.removeInstrumentSetting(PST::Tune, -1) == RC::success;
  ok &= se.removeInstrumentSetting(PST::Tune, -1) == RC::nothingToDo;
  getSamplerNote(&se, 69.f, 127.f, outL, outR);
  //rsPlotVector(outL);
  ok &= rsIsCloseTo(rsEstimateMidiPitch(outL, fs), 74.1f, tol);

  /*
  // Test behavior of pitch_keycenter when it's only defined on the instrument level:
  se.clearAllSfzSettings();
  //se.setInstrumentSetting(PST::PitchKeyCenter, 48.f, -1); 
  getSamplerNote(&se, 60.f, 127.f, outL, outR);
  pitch = rsEstimateMidiPitch(outL, fs);    // 69
  rsPlotVector(outL);
  */





  // Test behavior of pitch_keycenter when it's only defined on the group level:

  // Test behavior of pitch_keycenter when it's only defined on the instrument and group and level. 
  // In this case, the group setting should override the instrument setting:

  // Test behavior of pitch_keycenter when it's defined on all 3 levels. The region setting should
  // override both:




  /*
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
  */

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
  // -make tests where triggering a new region doe not start a new GroupPlayer
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

bool samplerSaveLoadTest()
{
  // This test also tests the file I/O

  bool ok = true;

  using namespace rosic::Sampler;
  using VecF = std::vector<float>;     // vector of sample values in RAM
  using SE   = rsSamplerEngineTest;
  using RC   = rsReturnCode;
  using PST  = Opcode;
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
  system("mkdir tmpwav");
  rosic::writeToStereoWaveFile("tmpwav/SinCos440Hz.wav", &sin440[0], &cos440[0], N, (int)fs, 16);
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
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::panN, -100.f, 1);      ok &= rc == RC::success;

  // Set up region for cosine:
  ri = se.addRegion(0);             ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::panN, +100.f, 1);      ok &= rc == RC::success;

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
  si = se.loadSampleToPool("tmpwav/SinCos440Hz.wav"); ok &= si == 0;
  ok &= se.getNumSamples() == 1;
  gi = se.addGroup();   ok &= gi == 0;
  ri = se.addRegion(0); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;
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

  // Helper function to save the state of se into a file and load it into se2 and compare their 
  // states. Returns true when both states are equal which indicates that the load/save roundtrip
  // has worked as expected:
  auto testSaveLoadRoundtrip = [&]()
  {
    se.saveToSFZ("tmp.sfz");
    se2.loadFromSFZ("tmp.sfz");
    return se2.isInSameStateAs(se);
  };

  // Test hikey/lokey opcodes by defining 2 regions:
  se.clearInstrument();
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  si = se.loadSampleToPool("Cos440Hz.wav"); ok &= si == 1;
  gi = se.addGroup(); ok &= gi == 0;
  ri = se.addRegion(0, 59, 61); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1); ok &= rc == RC::success;
  ri = se.addRegion(0, 69, 71); ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::PitchKeyCenter, 70.f, -1); ok &= rc == RC::success;
  ok &= testSaveLoadRoundtrip();

  // Test filter opcodes:
  using FltType = rosic::Sampler::FilterType;
  se.clearInstrument();
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  gi = se.addGroup(); ok &= gi == 0;
  ri = se.addRegion(0); ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1); ok &= rc == RC::success;
  se.setRegionSetting(0, 0, PST::filN_type, (float) FltType::lp_6, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,  1000.f, 1);
  se.saveToSFZ("tmp.sfz");
  se2.loadFromSFZ("tmp.sfz");
  ok &= se2.isInSameStateAs(se);
  se2.saveToSFZ("tmp2.sfz");                        // For manual inspection 
  //ok &= rsAreFilesEqual("tmp.sfz", "tmp2.sfz");   // ToDo: write this function

  // Test equalizer opcodes:

  // Set up eq1,2,3 without specifying frequencies or bandwidths (i.e. using the defaults):
  se.setRegionSetting(0, 0, PST::eqN_gain, 1.f, 1);
  se.setRegionSetting(0, 0, PST::eqN_gain, 2.f, 2);
  se.setRegionSetting(0, 0, PST::eqN_gain, 3.f, 3);
  ok &= testSaveLoadRoundtrip();

  // Set up some more eq bands, still without specifying freqs or widths:
  se.setRegionSetting(0, 0, PST::eqN_gain,  4.f,  4);
  se.setRegionSetting(0, 0, PST::eqN_gain,  5.f,  5);
  se.setRegionSetting(0, 0, PST::eqN_gain,  8.f,  8);
  se.setRegionSetting(0, 0, PST::eqN_gain, 13.f, 13); 
  se.setRegionSetting(0, 0, PST::eqN_gain, 10.f, 10); 
  ok &= testSaveLoadRoundtrip();

  // Set up frequencies and bandwidths of some of the eq bands:
  se.setRegionSetting(0, 0, PST::eqN_freq, 800.f,  8);
  se.setRegionSetting(0, 0, PST::eqN_bw,     1.8f, 8);
  se.setRegionSetting(0, 0, PST::eqN_freq, 500.f,  5);
  se.setRegionSetting(0, 0, PST::eqN_bw,     1.3f, 3);
  ok &= testSaveLoadRoundtrip();

  // Clear the region's settings, then set up 3 filters:
  se.clearRegionSettings(0, 0);
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::hp_12, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,    200.f,                  1);
  se.setRegionSetting(0, 0, PST::resonanceN, 10.f,                   1);
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::lp_12, 2);
  se.setRegionSetting(0, 0, PST::cutoffN,    800.f,                  2);
  se.setRegionSetting(0, 0, PST::resonanceN, 15.f,                   2);
  se.setRegionSetting(0, 0, PST::filN_type,  (float) FltType::lp_6,  3);
  se.setRegionSetting(0, 0, PST::cutoffN,    5000.f,                 3);
  ok &= testSaveLoadRoundtrip();

  // Test save and recall of loop settings:
  se.clearRegionSettings(0, 0);
  se.setRegionSetting(0, 0, PST::LoopMode,  (float) LoopMode::loop_continuous, -1);
  se.setRegionSetting(0, 0, PST::LoopStart, 0.f, -1);
  se.setRegionSetting(0, 0, PST::LoopEnd,   9.f, -1);
  ok &= testSaveLoadRoundtrip();


  // ToDo:
  // -Test save/load of loop settings
  // -Maybe make a local function testSaveLoadRoundTrip(se, ...) that saves the state of se and 
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


bool samplerPreProcessorTest()
{
  // Tests the pre-processing of the sfz parser, i.e. stripping off comments, etc.

  bool ok = true;

  // Create a string starting with a comment and then has a group and region and strip off the 
  // comment.

  using namespace rosic::Sampler;



  //auto stripComments = [](std::string& str) { rsRemoveLineComments(str, '/'); };

  // Helper function that takes a sfz-string with comments and the target string without comments
  // that should result from stripping the comments from the former string:
  auto check = [](std::string withComments, const std::string& withoutComments) 
  { 
    rsRemoveLineComments(withComments, '/');
    return withComments == withoutComments;
  };

  ok &= check("/Some comment\n<group>\n<region>",   "\n<group>\n<region>");
  ok &= check("/Some comment \n<group>\n<region>",  "\n<group>\n<region>");
  ok &= check("/ Some comment\n<group>\n<region>",  "\n<group>\n<region>");
  ok &= check("/ Some comment \n<group>\n<region>", "\n<group>\n<region>");

  ok &= check("<group>/Some comment\n<region>",     "<group>\n<region>");
  ok &= check("<group>  /Some comment\n<region>",   "<group>  \n<region>");
  ok &= check("<group>  / Some comment \n<region>", "<group>  \n<region>");

  ok &= check("/ Some = comment \n<group>\n<region>", "\n<group>\n<region>");


  rsAssert(ok);
  return ok;
}


bool samplerParserTest()
{
  // Tests the sfz parser by manually creating some sfz-strings and throwing them at the sampler
  // engine. We deliberately introduce some potential stumbling blocks such as repititions of 
  // seperator characters, etc.

  using SE  = rosic::Sampler::rsSamplerEngineTest;
  using RC  = rosic::Sampler::rsReturnCode;
  using PST = rosic::Sampler::Opcode;

  bool ok = true;

  ok &= samplerPreProcessorTest();

  // Create the engine and instruct it to load the just created sample files into its sample 
  // pool. These files are supposed to exist because samplerEngineUnitTestFileIO has created them
  // in some other unit test that ran before this one:
  int maxLayers = 8;
  SE se(maxLayers);
  int si;            // sample index
  si = se.loadSampleToPool("Sin440Hz.wav"); ok &= si == 0;
  si = se.loadSampleToPool("Cos440Hz.wav"); ok &= si == 1;

  // Add a group and to that group, add regions for the sine and cosine:
  int gi, ri, rc;                                    // group index, region index, return code
  gi = se.addGroup(); ok &= gi == 0;

  // Set up region for sine:
  ri = se.addRegion(0);             ok &= ri == 0;
  rc = se.setRegionSample(0, 0, 0); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 0, PST::panN, -100.f, 1);      ok &= rc == RC::success;

  // Set up region for cosine:
  ri = se.addRegion(0);             ok &= ri == 1;
  rc = se.setRegionSample(0, 1, 1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::PitchKeyCenter, 69.f, -1); ok &= rc == RC::success;
  rc = se.setRegionSetting(0, 1, PST::panN, +100.f, 1);      ok &= rc == RC::success;

  std::string sfzStr;
  SE se2(maxLayers);


  // Some tests that previously have triggered asserts:
  auto test = [&](const std::string& str)
  {
    return se2.setFromSFZ(sfzStr) == RC::success;
  };
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35   pan=79");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35  pan=79");
  ok &= test("<group>\n<region>sample=Sin440Hz.wav\n<region>sample=Cos440Hz.wav");
  ok &= test("<group>\n<region>sample=Sin440Hz.wav");
  ok &= test("<group>\n<region>\nsample=Sin440Hz.wav");
  ok &= test("<group>\n<region>sample=Sin440Hz.wav pan=79");
  ok &= test("<group>\n<region>sample=Sin440Hz.wav volume=35 pan=79");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35 pan=79");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35 cutoff=1234 pan=79 \n");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35 cutoff=1234 pan=79  ");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35 cutoff=1234 pan=79\n");
  ok &= test("<group>\n<region> sample=Sin440Hz.wav volume=35 cutoff=1234 pan=79 ");
  ok &= test("<group> <region>sample=Sin440Hz.wav pan=79");
  ok &= test(" <group> <region>sample=Sin440Hz.wav pan=79"); // triggers assert


  // ToDo: test a patch that doesn't define a group - the preprocessor should perhaps add one

  // Test reading an sfz-string where each opcode in on one line. This is the string that would be
  // generated and written into a file by a call to se.saveToSFZ("SinCos.sfz"); Then, a 2nd engine
  // tries to set itself up according to sfzStr. If all works as it should, se2 should aftwards be
  // in the same state as se:
  sfzStr = "\
<group>\n\
<region>\n\
sample=Sin440Hz.wav\n\
pitch_keycenter=69.000000\n\
pan=-100.000000\n\
<region>\n\
sample=Cos440Hz.wav\n\
pitch_keycenter=69.000000\n\
pan=100.000000";
  rc = se2.setFromSFZ(sfzStr); ok &= rc == RC::success;
  ok &= se2.isInSameStateAs(se);

  // Now we try it with spaces instead of newlines between the region opcodes:
  sfzStr = "\
<group>\n\
<region> \
sample=Sin440Hz.wav \
pitch_keycenter=69.000000 \
pan=-100.000000\n\
<region> \
sample=Cos440Hz.wav \
pitch_keycenter=69.000000 \
pan=100.000000";
  rc = se2.setFromSFZ(sfzStr); ok &= rc == RC::success;
  ok &= se2.isInSameStateAs(se);

  // Now with multiple spaces:
  sfzStr = "\
<group>\n\
<region> \
sample=Sin440Hz.wav  \
pitch_keycenter=69.000000  \
pan=-100.000000\n\
<region> \
sample=Cos440Hz.wav \
pitch_keycenter=69.000000    \
pan=100.000000";
  rc = se2.setFromSFZ(sfzStr); ok &= rc == RC::success;
  ok &= se2.isInSameStateAs(se);

  // Now with comments:
  sfzStr = "\
<group>\n\
<region>\n\
sample=Sin440Hz.wav / some comment\n\
pitch_keycenter=69.000000\n\
pan=-100.000000\n\
/ another comment   \n\
<region>\n\
sample=Cos440Hz.wav\n\
pitch_keycenter=69.000000\n\
pan=100.000000";
  rc = se2.setFromSFZ(sfzStr); ok &= rc == RC::success;
  ok &= se2.isInSameStateAs(se);

  // ToDo: test with a sample path that contains a '/', i.e. goes into a subdirectory. In this 
  // case, the '/' within that path should not be mistaken for a comment








  rsAssert(ok);
  return ok;
}

bool samplerAmplifierCoreTest()
{
  bool ok = true;

  rosic::Sampler::AmplifierCore amp;

  // Helper function to check, if the gain matrix coeffs of the amp have the desired values:
  auto checkCoeffs = [&](float a, float b, float c, float d, float tol = 0.f)
  {
    float A, B, C, D;
    amp.getChannelMixCoeffs(&A, &B, &C, &D);
    bool ok = true;
    ok &= rsIsCloseTo(a, A, tol);
    ok &= rsIsCloseTo(b, B, tol);
    ok &= rsIsCloseTo(c, C, tol);
    ok &= rsIsCloseTo(d, D, tol);
    return ok;
  };

  // Helper function to check, if the givne set of parameters leads to the given set of desired
  // mix coeffs:
  auto testAmp = [&](float vol, float pan, float width, float pos, 
    float a, float b, float c, float d, float tol = 0.f)
  {
    amp.setup(vol, pan, width, pos);
    return checkCoeffs(a,b,c,d);
  };

  // Test panorama and position:
  ok &= testAmp(0.f,    0.f,  100.f,    0.f,  1,0,0,1);  // neutral
  ok &= testAmp(0.f, -100.f,  100.f,    0.f,  2,0,0,0);  // pan hard left
  ok &= testAmp(0.f, +100.f,  100.f,    0.f,  0,0,0,2);  // pan hard right
  ok &= testAmp(0.f,    0.f,  100.f, -100.f,  2,0,0,0);  // pos hard left
  ok &= testAmp(0.f,    0.f,  100.f, +100.f,  0,0,0,2);  // pos hard right

  // Test width:
  float s = sqrt(0.5f);
  ok &= testAmp(0.f,    0.f, -100.f,    0.f,  0,1,1,0);  // swap L/R
  ok &= testAmp(0.f,    0.f,    0.f,    0.f,  s,s,s,s);  // zero width
  // It works but i think nevertheless the formula is still wrong - maybe we need something based
  // on sin/cos. Maybe sideGain should go like a sine segment from -sqrt(2) to +sqrt(2) in x=-2..+2
  // and midGain should go like a cosine segment from 0 to sqrt(2) to 0 in -2...+2. both should
  // meet at (1,1)

  // ToDo:
  // -Figure out if our parameterization makes indeed the fully general set of 2x2 matrices 
  //  reachable - and if so, how.
  // -Try to realize the following settings: invert L, invert R, invert both, invert M, invert S. 
  //  Combine all of them with a channel swap. Document the rules and settings, how to achieve 
  //  those.
  // -Figure out, if stereo signals should just ignore width and pos - i think so, so maybe we
  //  need a flag to switch behavior.
  // -I think, for stereo signals, the "pan" parameter "selects" between the two input channels and
  //  the "position" distributes the signal to the two output channels. A bit like gather and 
  //  scatter. So the mental model is: scale -> gather -> mix/shuffle -> scatter. But the sfz doc
  //  suggests more that "pan" is the scatter operation because pan applies to mono and stereo 
  //  signal alike whereas position applies only to stereo signals

  return ok;
}

bool samplerAmplifierTest()
{
  bool ok = true;

  ok &= samplerAmplifierCoreTest();

  int N = 200;

  using VecF = std::vector<float>;
  using SE   = rosic::Sampler::rsSamplerEngine2Test;
  using PST  = rosic::Sampler::Opcode;

  // Create and set up sampler engine:
  VecF noise = createColoredNoise(N, -6.02f);  
  SE se;
  addSingleSampleRegion(&se, noise, 60.f);

  // Test panning:
  ok &= testSamplerNote(&se, 60.f, 127.f, noise, noise, 0.0, false);
  se.setRegionSetting(0, 0, PST::panN, +100.f, 1);  // hard right
  ok &= testSamplerNote(&se, 60.f, 127.f, 0.f*noise, 2.f*noise, 0.f, false);
  se.setRegionSetting(0, 0, PST::panN, -100.f, 1);  // hard left
  ok &= testSamplerNote(&se, 60.f, 127.f, 2.f*noise, 0.f*noise, 0.f, false);
  se.setRegionSetting(0, 0, PST::panN, 0.f, 1);     // back to center

  // Set amplifier settings for region, group, instrument to 1,2,3 respectively and test it in both
  // modes (busMode and normal). In busMode, the gain should be 1+2+3 = 6dB and in normal mode, it 
  // should be 1dB:
  float tol = 1.e-6f;

  // gains for region, group and instrument
  float rVol = 3.f;
  float gVol = 6.f;
  float iVol = 0.f;   // for the first test, we'll leave that at zero
  se.setRegionSetting(0, 0, PST::volumeN, rVol, 1);
  se.setGroupSetting( 0,    PST::volumeN, gVol, 1);
  se.setInstrumentSetting(  PST::volumeN, iVol, 1);
  se.setBusMode(false);                // it should actually already be in that mode but anyway
  float g1 = RAPT::rsDbToAmp(rVol);    // only region gain counts
  ok &= testSamplerNote(&se, 60.f, 127.f, g1*noise, g1*noise, 0.f, false);

  float g2 = RAPT::rsDbToAmp(rVol + gVol + iVol);  // gains accumulate
  se.setBusMode(true); 
  ok &= testSamplerNote(&se, 60.f, 127.f, g2*noise, g2*noise, tol, false);

  // Now also with instrument-wide gain:
  iVol = 12.f;
  g2 = RAPT::rsDbToAmp(rVol + gVol + iVol);
  //se.reset();
  se.setInstrumentSetting(  PST::volumeN, iVol, 1);
  ok &= testSamplerNote(&se, 60.f, 127.f, g2*noise, g2*noise, tol, false);


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
  float slope  = -3.01f;     // spectral slope of the noise
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
  se.setRegionSetting(0, 0, PST::PitchKeyCenter,  60.f, -1);
  //addSingleSampleRegion(&se, noise, 60.f); // replace code above by that call

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
    return testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-7f, plot);

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
  ok &= testAgainstSvf(svf.LOWPASS,        Type::lp_12,  cutoff, reso, 1.e-5f, false);
  ok &= testAgainstSvf(svf.HIGHPASS,       Type::hp_12,  cutoff, reso, 1.e-5f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, cutoff, reso, 1.e-5f, false);
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
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 22050.f,   40.f, 1.e-4f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 20000.f,   40.f, 1.e-4f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6, 10000.f,   40.f, 1.e-4f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,  1000.f,   40.f, 1.e-4f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,   100.f,   40.f, 1.e-3f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,    10.f,   40.f, 1.e-3f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.1f,  40.f, 1.e-3f, false);
  ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.01f, 40.f, 1.e-3f, false);
  //ok &= testAgainstSvf(svf.BANDPASS_SKIRT, Type::bp_6_6,     0.0f,  40.f, 1.e-3, true);
  // todo: 
  // -maybe write a helper function taht goes through these checks als for lowpass and highpass
  // -maybe we should allow cutoff up to 30000 as in the EQ

  /*
  // Test impulse response ...or maybe not:
  float one = 1.f; float* pOne = &one; float** ppOne = &pOne;
  se.clearInstrument();
  se.addSampleToPool(ppOne, 1, 1, fs, "RAM");
  se.addGroup();   
  se.addRegion(0);
  se.setRegionSample(0, 0, 0);
  RAPT::rsFill(outL, 0.f);
  RAPT::rsFill(outR, 0.f);
  */


  // ToDo
  // -Cutoff=0 does not yet work - the svf produces silence and the sampler goes into bypass. Maybe
  //  we should do something different in this case. For a bandpass, it seems to make sense that 
  //  the limiting case is a lowpass. For a highpass, the limiting case should indeed be a bypass.
  //  for a lowpass, we may see the resonance peak at DC? ..which gets narrower when the cutoff 
  //  goes down and in the limit, nothing passes anymore?
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
  using Shape = rosic::Sampler::WaveshaperCore::Shape;

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
  float drive1 = 2.0f;
  Shape shape1 = Shape::tanh;
  //float postGain = 0.5f;
  //float dcOffset = 0.0;

  // Create target signal:
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(drive1 * sin440[n]);
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
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1);
  se.setRegionSetting(0, 0, PST::distortN_shape, float(shape1), 1);
  se.setRegionSetting(0, 0, PST::distortN_drive, drive1, 1);
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-7f, false);

  // Set up one region within one group and add a waveshaper to the group. When two notes are being 
  // played, the waveshaper should be applied to the sum of both notes:
  se.clearInstrument();  // clears the sample pool as well
  se.addSampleToPool(&pSmp, N, 1, fs, "Sine440");
  se.addGroup();   ok &= se.getNumGroups()   == 1;
  se.addRegion(0); ok &= se.getNumRegions(0) == 1;
  se.setRegionSample(0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1);
  se.setGroupSetting(0, PST::distortN_shape, float(shape1), 1);
  se.setGroupSetting(0, PST::distortN_drive, drive1, 1);

  // The class rsSamplerEngine should treat the group settings as fallback for when there is no
  // region setting and the DSP should be applied to each region separately:
  for(int n = 0; n < N; n++) 
    tgt[n] += tanh(drive1 * getSampleAt(sin440, 0.5f*n));
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= testSamplerOutput(&se, tgt, tgt, 1.e-13f, false);

  // Let the region override the group setting for drive. The shape setting should still come from
  // the group. We again play two notes at 60 and 48:
  float drive2 = 2*drive1;
  se.setRegionSetting(0, 0, PST::distortN_drive, drive2, 1);
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(drive2 * sin440[n]) + tanh(drive2 * getSampleAt(sin440, 0.5f*n));
  se.reset();
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= testSamplerOutput(&se, tgt, tgt, 1.e-13f, false);

  // Remove the group setting for the shape - now the waveshaper should use the default shape which
  // is linear, i.e. no shaping at all. 
  se.removeGroupSetting(0, PST::distortN_shape, -1);
  for(int n = 0; n < N; n++)
    tgt[n] = (drive2 * sin440[n]) + (drive2 * getSampleAt(sin440, 0.5f*n));
  se.reset();
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= testSamplerOutput(&se, tgt, tgt, 1.e-13f, false);

  // Add another region to the group without giving it a distortion setting - it should keep 
  // falling back to the group setting (linear, drive2) whereas the old region should keep using
  // the overriden setting (linear, drive1):
  se.addRegion(0); ok &= se.getNumRegions(0) == 2;
  se.setRegionSample(0, 1, 0);
  se.setRegionSetting(0, 1, PST::PitchKeyCenter, 60.f, -1);
  for(int n = 0; n < N; n++)
    tgt[n] += (drive1 * sin440[n]) + (drive1 * getSampleAt(sin440, 0.5f*n));
  se.reset();
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= testSamplerOutput(&se, tgt, tgt, 1.e-5f, false);
  // Why do we need such a high tolerance here? This seems wrong! figure out!

  // Start fresh, define the shape on the insrument level, the drive at the group level and the
  // DC-offset at the region level. Make sure that only 1 waveshaper is present in the chain.
  float dc1 = 0.125f;
  se.clearInstrument();
  se.addSampleToPool(&pSmp, N, 1, fs, "Sine440");
  se.addGroup();   ok &= se.getNumGroups()   == 1;
  se.addRegion(0); ok &= se.getNumRegions(0) == 1;
  se.setRegionSample(0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1);
  se.setInstrumentSetting(  PST::distortN_shape,  float(shape1), 1);
  se.setGroupSetting( 0,    PST::distortN_drive,  drive1, 1);
  se.setRegionSetting(0, 0, PST::distortN_dc, dc1, 1);
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(drive1 * sin440[n] + dc1);  // maybe include a gain1, too
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-7f, false);

  // ToDo:
  // -Try different shapes, use different sets of parameters, use DC, postGain, etc.
  // -Try the accumulate mode with both SamplerEngine and SamplerEngine2 with a region and group
  //  waveshaper. They should behave differently: 1 should apply both waveshapers and then mix 
  //  whereas 2 should first mix the shaped region output and then apply the group shaper to the 
  //  mix of both. Maybe engine 2 should provide the mix-before-apply mode as 3rd mode. Yes maybe
  //  that makes sense...hmmm...or maybe 1 should only provide fallback mode and both other types 
  //  of modes should be implemented by 2? That seems to make more sense because 1 is supposed to 
  //  implement straight sfz - all other "routing modes" shall be delegated to 2. Let there be 3 
  //  modes: fallback (default, sfz-like), accumulate-before-mix, accumulate-after-mix...maybe the
  //  3 modes should be selectable separately for group and instrument setting...but maybe that 
  //  would be too complex to implement
  // Cosmetics:
  // -Maybe drag out RegionPlayer from rsSamplerEngine
  // -maybe let testSamplerNote take a plotMode parameter which can be: 0: never plot, 1: always 
  //  plot, 2: plot when failed
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
  //  ..or maybe distort_ or distortion_ ...there is something like eg2drive_shape in sfz2...not 
  //  sure what that is


  rsAssert(ok);
  return ok;
}

bool samplerWaveShaperTest2()
{
  // Tests the additional routing options of rsSamplerEngine2 using waveshapers

  bool ok = true;

  using VecF  = std::vector<float>;
  using SE2   = rosic::Sampler::rsSamplerEngine2; 
  using PST   = rosic::Sampler::Opcode;
  using Shape = rosic::Sampler::DistortShape;
  using Ev    = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp  = Ev::Type;

  // Create a sinewave as example sample:
  float fs = 44100;  // sample rate
  float f  = 440.0;  // frequency of (co)sinewave sample
  int   N  = 500;    // length of (co)sinewave sample
  VecF sin440(N);    // sine wave
  VecF tgt(N);       // target output in tests
  float w = (float)(2*PI*f/fs);
  for(int n = 0; n < N; n++)
    sin440[n] = sinf(w*n);
  float *pSmp = &sin440[0];

  // Set up a sampler engine with 1 region within one group where the group has some waveshaper
  // settings. Play two notes in the 3 possible modes and compare with expected output signals:
  SE2 se;
  float driveG = 2.0f;
  Shape shapeG = Shape::tanh;
  se.addSampleToPool(&pSmp, N, 1, fs, "Sine440");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, PST::PitchKeyCenter, 60.f, -1);
  se.setGroupSetting(0, PST::distortN_shape, float(shapeG), -1);
  se.setGroupSetting(0, PST::distortN_drive, driveG, -1);

  // Default setting (fallback mode):
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(driveG * sin440[n]) + tanh(driveG * getSampleAt(sin440, 0.5f*n));
  se.reset();
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  ok &= testSamplerOutput(&se, tgt, tgt, 1.e-13f, false);
  // maybe factor this out into testSamplerNotes2

  // Mix-and-accumulate mode:
  se.setBusMode(true);
  for(int n = 0; n < N; n++)
    tgt[n] = tanh(driveG*sin440[n] + driveG*getSampleAt(sin440, 0.5f*n));
  se.reset();
  se.handleMusicalEvent(Ev(EvTp::noteOn, 48.f, 127.f));
  se.handleMusicalEvent(Ev(EvTp::noteOn, 60.f, 127.f));
  //ok &= testSamplerOutput(&se, tgt, tgt, 1.e-13, true);
  // Fails! The signal produced by the sampler has no saturation applied at all. That's halfway 
  // right: that the region player do not apply distortion anymore is correct - but now we need the
  // group player to apply it instead. this requires updating rsSamplerEngine2::startGroupPlayerFor

  // ToDo:
  // -Drag out the class RegionPlayer from the SamplerEngine and the class GroupPlayer from
  //  SamplerEngine2 - move them into a file SamplerPlayers.h/cpp
  // -Factor out a baseclass with the common stuff (dspChain, etc.)

  // Set up a 2nd waveshaper for the region and do the same test:
  // Set up a 3rd waveshaper for the instrument and do the same test:

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
  using Shape = rosic::Sampler::DistortShape;

  // Setup:
  float fs      = 44100.f;  // sample rate
  float slope   = -3.01f;   // spectral slope of the noise
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
  se.setRegionSetting(0, 0, PST::PitchKeyCenter,  60.f, -1);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::lp_6, 1);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff1, 1);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::hp_6, 2);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff2, 2);
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);

  // Add a waveshaper and after that a 3rd (lowpass) filter into the chain, such that 
  // the chain is now: LPF -> HPF -> WS -> LPF:
  float drive    = 4.0f;
  Shape shape    = Shape::tanh;
  float cutoff3  = 1000.f; 
  se.setRegionSetting(0, 0, PST::distortN_shape, float(shape), -1);
  se.setRegionSetting(0, 0, PST::distortN_drive, drive, -1);
  se.setRegionSetting(0, 0, PST::filN_type, (float)Type::lp_6, 3);
  se.setRegionSetting(0, 0, PST::cutoffN,   cutoff3, 3);

  // Create new target signal and run test:
  flt.setMode(flt.LOWPASS_IIT);
  flt.setCutoff(cutoff3);
  flt.reset();
  for(int n = 0; n < N; n++)
    tgt[n] = flt.getSample(tanh(drive * tgt[n]));
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);

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
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);
  cutoff1 = 4000.f;
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff1, 1);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);
  cutoff3 = 3000.f;
  se.setRegionSetting(0, 0, PST::cutoffN, cutoff3, 3);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);
  drive = 8.0;
  se.setRegionSetting(0, 0, PST::distortN_drive, drive, -1);
  updateTgt();
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);

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
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);

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
  using Region = rosic::Sampler::SfzInstrument::Region;


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
        y[n] = (float)flt.getSample(y[n]);
    }
  };

  // Create the basic sampler patch with no eq yet:
  SE se;
  const Region* r;
  float *pSmp = &noise[0];
  se.addSampleToPool(&pSmp, N, 1, fs, "Noise");
  se.addGroup();
  se.addRegion(0);
  se.setRegionSample( 0, 0, 0);
  se.setRegionSetting(0, 0, OC::PitchKeyCenter,  60.f, -1);

  // Test the simplest case: use only eq1 and specify all 3 parameters:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq1, 1);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw1, 1);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 1;
  applyEqs(noise, tgt, { gain1 }, { freq1 }, { bw1 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-5f, false);
  // Tolerance needs to be even higher than in the filter tests

  // Now we set only the gain of eq1. The freq and bandwidth should default to 50 Hz and 1 oct:
  //r->clearSettings();
  se.clearRegionSettings(0, 0);
  se.setRegionSetting(0, 0, OC::PitchKeyCenter, 60.f, -1);
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  ok &= r->getNumProcessors() == 1;
  applyEqs(noise, tgt, { gain1 }, { 50 }, { 1 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-3f, false);
  // Whoa - here we need a really high tolerance! OK, we have a moderately high Q and high gain and
  // a low frequency. Settings which are prone to numeric errors. But still, this is actually quite 
  // bad. Maybe we need to use double precision for filters and equalizers indeed. Or maybe use 
  // double for the coefficient calculation at least. Wait - no - the Q is actually rather low. And 
  // sfz allows very high Q settings

  // Now we set only the gain of eq3. The desired behavior is that we actually get 3 filter stages 
  // in the dsp chain but the first two are in neutral setting. We don't specify the center 
  // frequency or bandwidth. Therefore, the default values should be used which are 5 kHz and 1 
  // octave:
  //r->clearSettings();
  se.clearRegionSettings(0, 0);
  se.setRegionSetting(0, 0, OC::PitchKeyCenter, 60.f, -1);
  se.setRegionSetting(0, 0, OC::eqN_gain, gain3, 3);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 3;
  applyEqs(noise, tgt, { gain3 }, { 5000.f }, { 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-6f, false);

  // Add band 2. This has a default freq of 500Hz:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain2, 2); 
  applyEqs(noise, tgt, { gain2, gain3 }, { 500.f, 5000.f }, { 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4f, false);

  // Add band 1. This has a default freq of 50Hz:
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { 50.f, 500.f, 5000.f }, { 1.f, 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-3f, false);
  // Here, we also need a high tolerance

  // Change settings of band 1 to something more benign from a numeric point of view:
  freq1 = 1500;
  gain1 = 5.f;
  bw1   = 0.5f;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain1, 1);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq1, 1);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw1, 1);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, 500.f, 5000.f }, { bw1, 1.f, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4f, false);

  // Change settings of band 2:
  gain2 = -6;
  float freq2 = 2500;
  float bw2   = 0.3f;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain2, 2);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq2, 2);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw2, 2);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, freq2, 5000.f }, { bw1, bw2, 1.f });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4f, false);

  // Change settings of band 3:
  gain3 = 7;
  float freq3 = 8000;
  float bw3   = 2.2f;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain3, 3);
  se.setRegionSetting(0, 0, OC::eqN_freq, freq3, 3);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw3, 3);
  applyEqs(noise, tgt, { gain1, gain2, gain3 }, { freq1, freq2, freq3 }, { bw1, bw2, bw3 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4f, false);

  // Add a 5th band with no 4th in between, so we get a dspChain with 5 eqs but the 4th is neutral.
  // The new eq should be at its default freq of 1000:
  float gain5 = -3;
  float bw5   = 1.2f;
  se.setRegionSetting(0, 0, OC::eqN_gain, gain5, 5);
  se.setRegionSetting(0, 0, OC::eqN_bw,     bw5, 5);
  r = se.getRegion(0, 0);
  ok &= r->getNumProcessors() == 5;
  applyEqs(noise, tgt, { gain1, gain2, gain3, gain5 }, { freq1, freq2, freq3, 1000 }, 
    { bw1, bw2, bw3, bw5 });
  ok &= testSamplerNote(&se, 60.f, 127.f, tgt, tgt, 1.e-4f, false);

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

bool samplerOverloadTest()
{
  bool ok = true;

  using Vec    = std::vector<float>;
  using SE2    = rosic::Sampler::rsSamplerEngine2Test;
  using OC     = rosic::Sampler::Opcode;
  using Region = rosic::Sampler::SfzInstrument::Region;
  using Ev     = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp   = Ev::Type;

  int N = 200;
  Vec noise = createColoredNoise(N, -3.01f);

  // Create an engine which has only 8 filters available and add a region using 3 filters:
  SE2 se;
  se.setMaxNumFilters(8);
  addSingleSampleRegion(&se, noise, 60.f);
  se.setRegionSetting(0, 0, OC::cutoffN, 1000.f, 1);
  se.setRegionSetting(0, 0, OC::cutoffN, 2000.f, 2);
  se.setRegionSetting(0, 0, OC::cutoffN, 3000.f, 3);

  // Helper function that triggers a note and returns true, iff after triggering, the number of 
  // active RegionPlayers in the engine is equal to expectedLayers and the number of used filters
  // is equal to expectedFilters:
  auto testNote = [&](int note, int expectedLayers, int expectedFilters)
  {
    se.handleMusicalEvent(Ev(EvTp::noteOn, (float) note, 127));
    bool ok = true;
    ok &= se.getNumActiveLayers() == expectedLayers;
    ok &= se.getNumUsedFilters()  == expectedFilters;
    return ok;
  };
  // Maybe rename to testLayersAndFilters and write also helpers like testGroupsAndLayers, 
  // testVoicesAndLayers etc...or maybe just add 2 further parameters expectedGroups, 
  // expectedVoices to the function...yeah - always testing everything gives greater test coverage.
  // But expectedGroups would always be zero unless in busMode

  // Try to trigger notes. The first and second should play but the third should not due to not 
  // enough filters available:
  ok &= testNote(60, 1, 3);
  ok &= testNote(61, 2, 6);
  ok &= testNote(62, 2, 6);

  // Now put also cutoff settings for 3 filters into the group. In normal mode, this should make
  // no difference because the group merely defines a fallback value:
  se.reset();
  se.setGroupSetting(0, OC::cutoffN, 1100.f, 1);
  se.setGroupSetting(0, OC::cutoffN, 1200.f, 2);
  se.setGroupSetting(0, OC::cutoffN, 1300.f, 3);
  ok &= testNote(60, 1, 3);
  ok &= testNote(61, 2, 6);
  ok &= testNote(62, 2, 6);

  // Now switch to busMode in which case the group sub-bus needs an additional filter so the
  // first not already consumes 6 filters. A second note should not play:
  se.setBusMode(true);  // calls reset
  ok &= testNote(60, 1, 6);
  ok &= testNote(61, 1, 6);

  // Restrict the key-range for the region to 60..69 and add a second group and add a region to 
  // this group. The second region responds to key-range 70..79. Both group and region have two
  // filters. When we play one note in this 2nd region, 4 filters are consumed. If we then try to
  // trigger a 3rd note in the 1st region, it fails because it would need another 6 filters:
  se.setRegionSetting(0, 0, OC::LoKey, 60, -1);
  se.setRegionSetting(0, 0, OC::HiKey, 69, -1);
  se.addGroup();
  se.addRegion(1, 70, 79);
  se.setRegionSample(1, 0, 0);
  se.setGroupSetting( 1,    OC::cutoffN, 110.f, 1);
  se.setGroupSetting( 1,    OC::cutoffN, 220.f, 2);
  se.setRegionSetting(1, 0, OC::cutoffN, 100.f, 1);
  se.setRegionSetting(1, 0, OC::cutoffN, 200.f, 2);
  se.reset();
  ok &= testNote(70, 1, 4);  // 2nd group
  ok &= testNote(60, 1, 4);  // 1st group - can't play
  ok &= testNote(71, 2, 6);  // 2nd group - can play bcs uses less filters
  se.reset();

  // Remove 2 filters from the first group and do the same test. Now, it should work because the 
  // 2nd note needs only the 3 filters for the RegionPlayer and 1 for the GroupPlayer:
  se.removeGroupSetting(0, OC::cutoffN, 2);
  se.removeGroupSetting(0, OC::cutoffN, 3);
  ok &= testNote(70, 1, 4);  // uses 2/2 filters
  ok &= testNote(60, 2, 8);  // uses 3/1 filters
  se.reset();
  ok &= testNote(60, 1, 4);  // uses 3/1 filters
  ok &= testNote(61, 2, 7);  // uses 6/1 filters

  // OK - let's start fresh with regard to the settings. Then give each region just one filter
  // and each group again 3 filters:
  se.clearAllSfzSettings(); 
  se.setRegionSetting(0, 0, OC::cutoffN, 100.f,  1);
  se.setRegionSetting(1, 0, OC::cutoffN, 200.f,  1);
  se.setGroupSetting( 0,    OC::cutoffN, 1100.f, 1);
  se.setGroupSetting( 0,    OC::cutoffN, 1200.f, 2);
  se.setGroupSetting( 0,    OC::cutoffN, 1300.f, 3);
  se.setGroupSetting( 1,    OC::cutoffN, 110.f,  1);
  se.setGroupSetting( 1,    OC::cutoffN, 120.f,  2);
  se.setGroupSetting( 1,    OC::cutoffN, 130.f,  3);
  se.reset();
  ok &= testNote(60, 1, 4);  // uses 3/1 filters
  ok &= testNote(61, 2, 5);  // also 3/1 but only the 1 is new
  ok &= testNote(70, 2, 5);  // also 3/1 but all would be new

  // Test running out of GroupPlayers:
  se.clearAllSfzSettings();        // no filters anymore!
  se.setMaxNumLayers(5);
  se.setMaxNumGroups( 2);
  se.addGroup();
  se.addRegion(2, 80, 89);
  se.setRegionSample(2, 0, 0);
  se.reset();
  ok &= testNote(60, 1, 0);
  ok &= testNote(70, 2, 0);
  ok &= testNote(80, 2, 0);   // we don't have enough GroupPlayers

  // Test running out of RegionPlayers. Allow a maximum of 8 layers, create a Region with 3 layers,
  // trigger 3 notes. The 3rd note should not trigger because 6 RegionPlayers are already used up 
  // by the other two notes and the 3rd would need 3 more but only 2 are left:
  se.clearInstrument();                    // start fresh with clean slate
  se.setMaxNumLayers(8);
  se.setBusMode(false);                    // group limitations shouldn't play a role here
  addSingleSampleRegion(&se, noise, 60.f); // add sample and 1 region
  se.addRegion(0);                         // add 2nd region
  se.setRegionSample(0, 1, 0);             // ...and set up its sample
  se.addRegion(0);                         // add 3rd region
  se.setRegionSample(0, 2, 0);             // ...and set up its sample
  ok &= testNote(60, 3, 0);
  ok &= testNote(61, 6, 0);
  ok &= testNote(62, 6, 0);

  // ToDo: 
  // 

  rsAssert(ok);
  return ok;
}

bool samplerLoopTest()
{
  // Tests the looping feature of the sampler engine

  bool ok = true;

  using namespace rosic::Sampler;
  using Vec = std::vector<float>;
  using SE  = rosic::Sampler::rsSamplerEngineTest;
  using OC  = rosic::Sampler::Opcode;

  // Create a sinewave with 3 cycles as example sample:
  int cycleLength = 101;    // length of the sienwave cycle
  int numCycles   = 3;      // number of cycles in the sample
  Vec sinTable(numCycles*cycleLength); // sine wave sample 
  double w = 2*PI/cycleLength; 
  for(size_t n = 0; n < sinTable.size(); n++)
    sinTable[n] = (float)sin(w*n);
  //rsPlotVector(sinTable);

  // Playback settings:
  float fs = 44100;  // playback sample rate
  float f  = 440.0;  // frequency of sinewave to generate
  int   N  = 2000;   // number of samples to generate

  // Set up the engine:
  SE se;
  double f0 = fs/cycleLength;  // fundamental freq of the sample
  double rootKey = RAPT::rsFreqToPitch(f0);
  addSingleSampleRegion(&se, sinTable, (float)rootKey);
  se.setRegionSetting(0,0, OC::LoopMode,  (float)LoopMode::loop_continuous, -1);
  se.setRegionSetting(0,0, OC::LoopStart, 0,               -1);
  se.setRegionSetting(0,0, OC::LoopEnd,   0 + (float)cycleLength, -1);

  // Helper function to produce the error signal between the ideal sinewave of given frequency and
  // the signal produced by the sampler engine for given key:
  auto getError = [&](int key, double freq)
  {
  // Produce target signal:
    Vec tgt(N);
    w = 2*PI*freq/fs;
    for(int n = 0; n < N; n++)
      tgt[n] = (float)sin(w*n);

    // Produce output and return error:
    Vec outL(N), outR(N);
    se.reset();
    se.handleNoteOn(key, 127);
    for(int n = 0; n < N; n++)
      se.processFrame(&outL[n], &outR[n]);
    rsAssert(outR == outL);
    return tgt - outL;
  };

  // We expect some error due to the linear interpolation:  
  float tol = 1.e-3;
  Vec err1 = getError(69, 440);
  ok &= rsMaxAbs(err1) <= tol;


  // Set the loop length to 2 cycles - this should make no difference up to roundoff:  
  se.setRegionSetting(0,0, OC::LoopEnd, float(2*cycleLength), -1);
  Vec err2 = getError(69, 440);
  ok &= rsMaxAbs(err2) <= tol;

  // Set the loop length to all 3 cycles and therefore equal to the total length of the sample:
  se.setRegionSetting(0,0, OC::LoopEnd, (float)sinTable.size(), -1);
  Vec err3 = getError(69, 440);
  ok &= rsMaxAbs(err3) <= tol;

  // Plot error signals:
  //rsPlotVectors(err1, err2, err3);
  // Error is between +-0.0005. That's an SNR of about 66 dB. Actually not that bad for not even
  // trying to imlement decent interpolation. But it will probably get a lot worse for higher 
  // frequencies. The error has strange discontinuities in the derivative - why? I tried without 
  // loop and the corners are still there. This seems to be a feature of the linear interpolator
  // and is not related to the loop implementation.

  // Let's have a look ad the difference of the error. We expect it to be nonzero due to different
  // roundoff errors when one is in the loop while the other one isn't:
  //rsPlotVectors(err1 - err2);
  //rsPlotVectors(err1 - err3);
  //rsPlotVectors(err2 - err3);
  // The error difference is zero for 100 samples, the nonzero at a level of about 10^-6 for the 
  // rest of the longer loop, then zero again when they are again "in phase", etc. That's quite 
  // interesting actually. Shouldn't the 2nd and 3rd cycle contain exactly the same data as the 
  // 1st?

  // To test the one-shot mode, we create a sample of an exponential decay and trigger the same 
  // note 3 times with note-offs in between the note-ons. The note-offs should be ignored and the 
  // sample should alway play until the end is reached:
  int L      = 500;    // length of decay-sample
  int start2 = 200;    // start time of 2nd note
  int start3 = 300;    // start time of 3rd note

  // Render and set up "shot" sample:
  Vec decay(L);
  for(int n = 0; n < L; n++)
    decay[n] = exp(-0.01*n);
  //rsPlotVector(decay);
  se.clearInstrument();
  addSingleSampleRegion(&se, decay, 60);
  se.setRegionSetting(0,0, OC::LoopMode,  (float)LoopMode::one_shot, -1);

  // Produce the output:
  Vec outL(N), outR(N);
  se.handleNoteOn(60, 100);
  for(int n = 0; n < start2; n++)
    se.processFrame(&outL[n], &outR[n]);
  se.handleNoteOn(60, 0);   // noteOff
  se.handleNoteOn(60, 100);
  for(int n = start2; n < start3; n++)
    se.processFrame(&outL[n], &outR[n]);
  se.handleNoteOn(60, 0); 
  se.handleNoteOn(60, 100);
  for(int n = start3; n < N; n++)
    se.processFrame(&outL[n], &outR[n]);

  // Produce target output:
  Vec tgt(N);
  using AT = RAPT::rsArrayTools;
  AT::addInto(&tgt[0], N-0,      &decay[0], L, 0);
  AT::addInto(&tgt[0], N-start2, &decay[0], L, start2);
  AT::addInto(&tgt[0], N-start3, &decay[0], L, start3);
  ok &= outL == tgt && outR == tgt;
  //rsPlotVectors(tgt, outL, outR);


  // ToDo:
  // -add reverse playback mode...maybe this should be one of the loop_modes? or do we need an extra 
  //  opcode for that?
  // -Test what happens when loopEnd is <= loopStart. I guess, it jumps forward by loopLength 
  //  every sample after reaching loopEnd. We should probably ensure that this doesn't happen
  //  maybe use loopStart = min(loopStart, loopEnd). Then, the loop would just be ignored. That may
  //  actually be the most reasonable behavior in such a case.
  // -Check, if the behavior is correct with respect to sustain, note-off, etc. I think, in 
  //  loop_sustain modewe should do the warp around conditionally if the note is being held
  // -Later when better interpolation is implemented, do similar tests with less tolerance.
  // -Implement the default_dir (or whatever it is called) opcode and write some patches using
  //  the CyclePack. Maybe wite a patch using mip-mapped cycles? ...but we don't have any in the 
  //  CyclePack - maybe add Saw_1024, Saw_512, Saw_256, Saw_128, ... Maybe make a SuperSaw patch
  //  from 7 of them using one filter for the whole group...should also have some highpass.
  //  ...use key-crossfading...oh - and the files are in flac format which is not (yet?) supported.
  //  a freq of 22.5 Hz at fs = 44100 gives a clean cycleLength of 1960. Maybe use that for
  //  single-cycle waves for ths sampler. 1960 is divisble by 8 but not 16 it's then also 
  //  divisible by 5 and 7
  //  goals: we want to achieve an frequency that is exactly representable as integer midi-key 
  //  without detune - so we must choose an A
  //  tableSize = L = 2048, rootKey is an A so we don't need any detune, fs = 44000, freq=fs/L
  //  maybe we want 2048 represent 27.5Hz - what sampleRate do wo need: fs = 56320 which 
  //  corresponds exaclty to midi key 21

  rsAssert(ok);
  return ok;
}

bool samplerKeyVelTrackTest()
{
  // Tests the key- and velocity tracking opcodes.

  bool ok = true;

  using Vec = std::vector<float>;
  using SE  = rosic::Sampler::rsSamplerEngineTest;
  using OC  = rosic::Sampler::Opcode;
  using FT  = rosic::Sampler::FilterType;

  int N = 1000;
  SE se;
  Vec noise = createColoredNoise(N, -10.f);
  addSingleSampleRegion(&se, noise);

  // Set up velocity tracking of the volume:
  float vel = 100;            // velocity for note to play
  float vel_track = -100.f;   // in %
  se.setRegionSetting(0, 0, OC::ampN_veltrack, vel_track, 1);

  // Compute target amplitude for given settings of velocity and vel_track:
  float dB  = 40 * log10(127.f/vel);    // change of volume at given velocity
  float amp = RAPT::rsDbToAmp(0.01f * vel_track * dB);
  ok &= testSamplerNote(&se, 60, vel, amp*noise, amp*noise, 0.f, false);
  // unit is percent, formula is dB = 20 log (127^2 / Velocity^2) ..i think, that's the change in
  // dB at 100%? range is -100...+100 ..so it's 40 log (127/Velocity)

  // Reset veltrack and set up keytrack:
  char key = 65;
  float key_track = -2.f;
  dB  = key_track * (key-60);
  amp = RAPT::rsDbToAmp(dB);
  se.setRegionSetting(0, 0, OC::ampN_veltrack,  0.f,        1);  // no vletrack anymore
  se.setRegionSetting(0, 0, OC::ampN_keytrack,  key_track,  1);  // -1 dB/key (range: -96..+12)
  se.setRegionSetting(0, 0, OC::ampN_keycenter, 60,         1);  // neutral at A4
  se.setRegionSetting(0, 0, OC::PitchKeyCenter, key,       -1);  // to avoid transposition
  ok &= testSamplerNote(&se, key, vel, amp*noise, amp*noise, 0.f, false);

  // Test pitch-keytracking by setting it to 50% (or 50 cents per key), triggering a note 2 octaves
  // above the pitch_keycenter and verify that it gets transposed only by one octave:
  se.clearRegionSettings(0, 0);  // resets keycenter to the default of 60, removes veltrack stuff
  key_track = 50.f;              // 50 cents per key
  se.setRegionSetting(0, 0, OC::PitchKeyTrack, key_track, -1);
  Vec tgt = rsApplyResampling(noise, 2.f);
  ok &= testSamplerNote(&se, 84, vel, tgt, tgt, 0.f, false);  // 84 = 60 + 2*12

  // Test filter keytracking by creating a sort filter-whistle patch:
  se.clearRegionSettings(0, 0); 
  se.setRegionSetting(0, 0, OC::filN_type, (float)FT::bp_6_6, 1);
  se.setRegionSetting(0, 0, OC::PitchKeyTrack,    0.f, -1);  // pitch shall not track key
  se.setRegionSetting(0, 0, OC::filN_keytrack,  100.f,  1);  // cutoff shall track key 100%
  se.setRegionSetting(0, 0, OC::filN_keycenter,  69.f,  1);  // A4 is the neutral key
  se.setRegionSetting(0, 0, OC::cutoffN,        440.f,  1);  // at A4, cutoff is 440
  se.setRegionSetting(0, 0, OC::resonanceN,      40.f,  1);  // we use high resonance
  float fs = (float)se.getOutputSampleRate();
  tgt = rsApplyFilter(noise, FT::bp_6_6, 880.f, fs, 40.f); 
  ok &= testSamplerNote(&se, 81, vel, tgt, tgt, 2.e-5f, false); // 81 = 69 + 12

  // Test vel-tracking:
  vel_track = -1200.f;  // -12 semitones reduction at min-vel, i.e. vel=1
  se.setRegionSetting(0, 0, OC::filN_veltrack, vel_track,  1); 
  tgt = rsApplyFilter(noise, FT::bp_6_6, 440.f, fs, 40.f); 
  ok &= testSamplerNote(&se, 69, 127, tgt, tgt, 4.e-5f, false);  // at key=69, keytrack should be neutral
  tgt = rsApplyFilter(noise, FT::bp_6_6, 220.f, fs, 40.f); 
  ok &= testSamplerNote(&se, 69, 1, tgt, tgt, 0.00015f, false);
  // We need quite high tolerances here. I'm not sure about the veltrack formula - it's just a 
  // guess based on what i think, the behavior should be. I think, at vel = 127, the cutoff should 
  // be unmodified and at vel=1, the cutoff should be reduced by 1200 cents, i.e. 12 semitones, 
  // i.e. 1 octave

  // ToDo:
  // -Verify the formula for vel_tracking of amplitude against some reference implemenetation (sfz+)

  rsAssert(ok);
  return ok;
}

bool samplerEffectsTest()
{
  bool ok = true;

  // For inspection in the debugger. We want to keep the sizes of these DSP objects small because 
  // we'll potentially have to pre-allocate a lot of them when a patch is loaded:
  using namespace rosic::Sampler;
  int size;

  // Sizes of DSP cores:
  size = sizeof(AmplifierCore);          // 16
  size = sizeof(FilterCore);             // 64
  size = sizeof(WaveshaperCore);         // 24
  size = sizeof(rsSamplerEnvGen);        // 36
  size = sizeof(rsSamplerLowFreqOsc);    // 16
  // FilterCore is quite large - try to reduce it - maybe by defining different kinds of filters 
  // because many filter types do not use all variables. See below

  // Sizes of infrastructural classes:
  size = sizeof(RegionPlayer);           // 184
  size = sizeof(Effect);                 //  48
  size = sizeof(Parameter);              //   8

  // Sizes of some basic underlying data structures:
  size = sizeof(std::vector<Parameter>); // 32
  //std::cout << sizeof(std::vector<Parameter>) << '\n'; // 24 in release mode
  // Whoa! 32? that's twice as much as i expected! Two pointers of 64 bit (8 byte) size should take 
  // only 16 byte! Maybe when we optimize the memory usage later, switch to a hand-rolled vector 
  // replacement that has smaller overhead use a using Vector = std::vector directive somewhere 
  // where we can switch between the two (std::vector is nice for debugging). Hmm..in release mode,
  // it's only 24 - but that's still 8 too much

  // Sizes of some DSP classes from RAPT:
  size = sizeof(RAPT::rsFirstOrderFilterBase<float, float>);  // 20
  // It takes less than one third of the memory (0.3125x as much) of what FilterCore takes and 
  // first order filters are used a lot in my patches, so this may have a big impact to use a 
  // different class for simpler filters. First order filters are good for processing the layers
  // statically. The full blown filter is good for processing the mix of layers

  //size = sizeof(SP::Filter);
  //size = sizeof(SP::WaveShaper);
  // later move this into a (yet to be written) benchmark testbed

  // -Move this into some performance test function


  ok &= samplerFilterTest();      // tests the different filter modes
  ok &= samplerWaveShaperTest();  // tests the waveshaping DSP module
  ok &= samplerDspChainTest();    // uses multiple filters and a waveshaper in between
  ok &= samplerEqualizerTest();
  ok &= samplerAmplifierTest();
  ok &= samplerWaveShaperTest2(); // tests the different rotung options using waveshaping


  // ToDo:
  // -Implement more DSP modules: echo, vibrato, flanger, phaser, chorus, etc., 
  //  ...delay based algorithms could become a memory-hog when we need to pre-allocate many of 
  //  them. Maybe let's stay away from them for the moment. How about freq-shifting?

  rsAssert(ok);
  return ok;
}

bool samplerModulationsTest()
{
  // ToDo:
  // -Set up a sampler engine with a sample that is just DC and apply an amplitude envelope
  // -Envelopes in sfz have the parameters: delay, start, attack, hold, decay, sustain, release
  // -Maybe implement a linear shape first.
  // -Figure out what shapes sfz+ produces.
  // -Maybe we should also define shape parameters for the segments, i.e. attack_shape, 
  //  decay_shape, release_shape

  bool ok = true;

  rsAssert(ok);
  return ok;
}


bool samplerEngineUnitTest()
{
  bool ok = true;

  // The new test that is currently under construction:
  //ok &= samplerKeyVelTrackTest();
  //ok &= samplerModulationsTest();
  ok &= samplerParserTest();

  // The tests, that already pass and are supposed to continue to do so:
  ok &= samplerDataTest();           // datastructure for representing an sfz
  ok &= samplerRegionPlayerTest();   // basic playback of regions
  ok &= samplerBusModeTest();        // basic playback in busMode
  ok &= samplerSaveLoadTest();       // saving and loading of sfz files
  ok &= samplerParserTest();         // uses some files created by "..FileIO" -> order matters!
  ok &= samplerEffectsTest();        // effect chain
  ok &= samplerModulationsTest();    // modulation system
  ok &= samplerOverloadTest();       // behavior in overload conditions
  ok &= samplerKeyVelTrackTest();    // key- and velocity tracking
  ok &= samplerLoopTest();           // loop modes


  // ToDo:
  // -implement key/vel crossfade
  // -implement cutoff_ccN ..or actually cutoffN_ccX
  //  -extend the PlayStatus class to a MidiStatus class that provides the info about 
  //   controllers, pitch-wheel, etc.
  // -implement loop_mode=loop_sustain, 
  // -implement and test loop_end beyond the last sample - it should append an appropriate amount 
  //  of silence before looping back. rationale: together with delay, the feature can be used to 
  //  create repititve hits (of a drum sample, say) to create looped rhythms
  //  Here (under "Playing part of the sample repeatedly"):
  //    http://www.drealm.info/sfz/plj-sfz.xhtml
  //  it is said that when no loop start point is given, it doesn't default to 0 but to offset.
  //  There's an unofficial opcode delay_beats described there (at the bottom). That could be 
  //  useful for programming drumloops directly in sfz (but where does the bpm info come from?
  //  from an opcode or from the host?)
  // -maybe keep tests for the basic SamplerEngine and samplerEngine2 (supporting busMode) in 
  //  seperate functions. We may want to be able to easily seperate out the basic implementation
  //  along with its test, if needed.
  // -clean up codebase and fix warnings
  // -Test to define the sample opcode on instrument and group level. We can have various regions
  //  in a group that use the same sample, test what happens when a region has no sample defined
  //  (it should be silent, of course)
  // -Add an overload test that simulates conditions when the engine is running out of resources 
  //  such as DSPs, players, memory, etc. Make sure that things get cleaned up correctly in cases 
  //  where a partially assembled RegionPlayer must be rolled back due to lack of resources.
  // -Figure out and implement the correct filter design formulas for lpf_2p, hpf_2p, bpf_2p and 
  //  the equalizers. It's rather interesting that for the filters, the range is 0..fs/2 and for
  //  equalizers it's 0..30kHz (done?)
  // -Test filters with the most extreme Q settings, also at low frequnecies. That's where we 
  //  should expect numerical precision issues, especially when using single precision. Maybe using
  //  an SVF implementation can help against this. OK - up to resonance=40 (the limit set by the 
  //  sfz spec), we don't seem to get any problems with float, even for subsonic cutoffs.



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


 https://github.com/sfz/opcode-suggestions/issues

maybe open a thread at kvr..something like SFZ: interpreting and extending the format specification


a nice list of all the various kinds of envelopes (ADSR, DAHDSR, ...):
https://synth.fandom.com/de/wiki/H%C3%BCllkurve

*/