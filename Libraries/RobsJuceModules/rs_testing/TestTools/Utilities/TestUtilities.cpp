#include "TestUtilities.h"

void reportUnitTestSuccess(const std::string& name = std::string())
{ 
  std::cout << name << "OK\n"; 
}

void reportUnitTestFailure(const std::string& name = std::string())
{ 
  std::cout << name << "!!!!----> F A I L E D <----!!!!\n"; 
}

bool runUnitTest(bool (*test)(), const std::string& name)
{
  std::cout << name + ": ";
  bool ok = test();
  //rsAssert(ok); // break, if test fails
  if(ok) reportUnitTestSuccess();
  else   reportUnitTestFailure();
  return ok;
}

/*
bool detectMemoryLeaks()
{
  #ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
  #else
  return false;
  #endif
}
*/

std::vector<double> rsLinearRangeVector(int N, double min, double max)
{
  std::vector<double> v(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&v[0], N, min, max);
  return v;
}

std::vector<double> rsExponentialRangeVector(int N, double min, double max)
{
  std::vector<double> v(N);
  RAPT::rsArrayTools::fillWithRangeExponential(&v[0], N, min, max);
  return v;
}

std::vector<double> rsRandomVector(int N, double min, double max, int seed)
{
  std::vector<double> v(N);
  RAPT::rsArrayTools::fillWithRandomValues(&v[0], N, min, max, seed);
  return v;
}

std::vector<double> rsRandomIntVector(int N, int min, int max, int seed)
{
  std::vector<double> v(N);
  rsNoiseGenerator<double> ng;
  ng.setSeed(seed);
  for(int i = 0; i < N; i++)
  {
    unsigned long raw = ng.getSampleRaw();
    int iVal = raw % (max-min) + min;
    v[i] = (double) iVal;
  }
  return v;
}

std::vector<double> rsApplyFunction(const std::vector<double>& v, double p, 
  double (*f) (double, double))
{
  std::vector<double> r(v.size());
  for(size_t i = 0; i < r.size(); i++)
    r[i] = f(v[i], p);
  return r;
}


#if defined(_MSC_VER)
// works only on MSVC and we need this function only in the performance tests, which we do only
// with the MSVC compiler anyway:
std::string toString(int n)
{
  //return std::to_string((_Longlong)n);
  return std::to_string((long long)n);
}
#endif

/*
double rsSquare(double x)
{
  return x*x;
}
*/

using namespace std;

#undef min

bool areNumbersEqual(double x, double y, double relativeTolerance)
{
  // nan == nan:
  double tmp = RS_NAN(double);
  if( memcmp(&x, &tmp, sizeof(double)) == 0 && memcmp(&y, &tmp, sizeof(double)) == 0 )
    return true;

  // inf == inf:
  tmp = RS_INF(double);
  if( x == tmp && y == tmp )
    return true;

  // -inf == -inf:
  tmp = -tmp;
  if( x == tmp && y == tmp )
    return true;

  // catch case where one or bothe of x, y are denormals - in this case, the 
  // relativeTolerance*rsMax(fabs(x), fabs(y)) yields a zero absolute tolerance whereas the 
  // difference on the left hand side is a very small (denormal) nonzero number
  //double denormThresh = RS_MIN(double);
  double denormThresh = std::numeric_limits<double>::min();
  if( fabs(x) <= denormThresh && fabs(y) <= denormThresh )
    return true;

  // x == y, if the absolute difference is below a relative tolerance:
  return fabs(x-y) <= relativeTolerance * max(fabs(x), fabs(y));
}

RAPT::rsWindowFunction::WindowType stringToWindowType(const std::string& wt)
{
  typedef RAPT::rsWindowFunction::WindowType WT;
  if(wt == "rc") return WT::rectangular;
  if(wt == "hn") return WT::hanningZZ;
  if(wt == "hm") return WT::hamming;
  if(wt == "bm") return WT::blackman;
  if(wt == "bh") return WT::blackmanHarris;
  if(wt == "dc") return WT::dolphChebychev;
  RAPT::rsError("Unknown window type");
  return WT::rectangular;
}


void addSingleSampleRegion(rosic::Sampler::rsSamplerEngine* se,
  const std::vector<float>& sample, float keyCenter, double sampleRate)
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

void setupForSineWave(rosic::Sampler::rsSamplerEngine* se, int N)
{
  std::vector<float> sineWave(N);
  for(int n = 0; n < N; n++)
    sineWave[n] = (float) sin(2.0*PI*n/N);
  se->clearInstrument();
  addSingleSampleRegion(se, sineWave, 21.f, 56320.f);

  // Set up loop settings:
  using namespace rosic::Sampler;
  se->setRegionSetting(0, 0, Opcode::LoopMode, (float) LoopMode::loop_continuous, -1);
  se->setRegionSetting(0, 0, Opcode::LoopStart, 0, -1);
  se->setRegionSetting(0, 0, Opcode::LoopEnd,   N, -1);
}

void setupForLoopedDC(rosic::Sampler::rsSamplerEngine* se, int N)
{ 
  using OC  = rosic::Sampler::Opcode;
  std::vector<float> dc(N);
  rsFill(dc, 1.f);
  se->clearInstrument();

  addSingleSampleRegion(se, dc);
  // What about the root-key and sample-rate? We are depending on the defaults here (60, 44100)
  // but that may not be the best thing to do. We see looping artifacts in samplerEnvTest

  se->setRegionSetting(0, 0, OC::LoopMode, (float)rosic::Sampler::LoopMode::loop_continuous, 1);
  se->setRegionSetting(0, 0, OC::LoopStart, 0.f,      1);
  se->setRegionSetting(0, 0, OC::LoopEnd,  (float) N, 1);
  // code almost the same as for sine wave -> get rid of duplication
}

void getSamplerNote(rosic::Sampler::rsSamplerEngine* se, float key, float vel,
  std::vector<float>& outL, std::vector<float>& outR)
{
  rsAssert(outL.size() == outR.size());
  using Ev   = rosic::Sampler::rsMusicalEvent<float>;
  using EvTp = Ev::Type;
  se->handleMusicalEvent(Ev(EvTp::noteOn, key, vel));
  for(int n = 0; n < (int) outL.size(); n++)
    se->processFrame(&outL[n], &outR[n]);
  // Should we clear the outL/R arrays first? Maybe not, if we want instruments to accumuluate 
  // their outputs in ToolChain
}

