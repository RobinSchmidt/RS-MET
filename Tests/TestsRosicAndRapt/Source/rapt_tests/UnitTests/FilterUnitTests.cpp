
//#include "../../../../Libraries/ThirdParty/Nayuki/SlidingWindowMinMax.hpp"
#include "../../../../Libraries/ThirdParty/Nayuki/SlidingWindowMinMax.hpp"

using namespace std;

// maybe move somewhere else for sharing
bool isCloseTo(complex<float> x, complex<float> y, float tol)
{
  return abs(x-y) <= tol;
}

bool prototypeDesignUnitTest()
{
  // shorthands:
  typedef rsPrototypeDesigner<float> PD;
  typedef complex<float> CF;
  float inf = std::numeric_limits<float>::infinity();

	bool r = true;      // test result
  //float tol = 0.f;  // zero tolerance - float comparisons are currently exact
  float tol = 1.e-7f;
  CF p[3], z[3];      // arrays for retrieved poles
  PD pd;              // prototype designer object

  // Papoulis filters:
  pd.setApproximationMethod(PD::PAPOULIS);
  pd.setPrototypeMode(PD::LOWPASS_PROTOTYPE);
  pd.setOrder(5);                                               // 5th order lowpass
  pd.getPolesAndZeros(&p[0], &z[0]);              
  r &= isCloseTo(p[0], CF(-0.153586745f, 0.968145967f), tol);
  r &= isCloseTo(p[1], CF(-0.388139844f, 0.588632464f), tol);
  r &= isCloseTo(p[2], CF(-0.468089849f, 0.f),          tol);
  r &= z[0] == CF(inf, 0.f);
  r &= z[1] == CF(inf, 0.f);
  r &= z[2] == CF(inf, 0.f);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 0;
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);                  // 5th order low-boost
  pd.setGain(+6);
  pd.getPolesAndZeros(&p[0], &z[0]);              
  r &= isCloseTo(p[0], CF(-0.147757247f, 0.931399286f), tol);
  r &= isCloseTo(p[1], CF(-0.373407722f, 0.566290498f), tol);
  r &= isCloseTo(p[2], CF(-0.450323164f, 0.f),          tol);
  r &= isCloseTo(z[0], CF(-0.185926169f, 0.996051848f), tol);
  r &= isCloseTo(z[1], CF(-0.481426001f, 0.611809969f), tol);
  r &= isCloseTo(z[2], CF(-0.590867400f, 0.f),          tol);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 5;
  pd.setGain(-6);                                               // 5th order low-cut
  pd.getPolesAndZeros(&p[0], &z[0]); 
  r &= isCloseTo(p[0], CF(-0.185926169f, 0.996051848f), tol);
  r &= isCloseTo(p[1], CF(-0.481426001f, 0.611809969f), tol);
  r &= isCloseTo(p[2], CF(-0.590867400f, 0.f),          tol);
  r &= isCloseTo(z[0], CF(-0.147757247f, 0.931399286f), tol);
  r &= isCloseTo(z[1], CF(-0.373407722f, 0.566290498f), tol);
  r &= isCloseTo(z[2], CF(-0.450323164f, 0.f),          tol);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 5;
  pd.setPrototypeMode(PD::LOWPASS_PROTOTYPE);
  pd.setOrder(6);                                               // 6th order lowpass
  pd.getPolesAndZeros(&p[0], &z[0]);
  r &= isCloseTo(p[0], CF(-0.115192428f, 0.977922320f), tol);
  r &= isCloseTo(p[1], CF(-0.308961034f, 0.698167443f), tol);
  r &= isCloseTo(p[2], CF(-0.438901573f, 0.239981338f), tol);
  r &= z[0] == CF(inf, 0.f);
  r &= z[1] == CF(inf, 0.f);
  r &= z[2] == CF(inf, 0.f);
  r &= pd.getNumFinitePoles() == 6;
  r &= pd.getNumFiniteZeros() == 0;
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);                  // 6th order low-boost
  pd.setGain(+6);
  pd.getPolesAndZeros(&p[0], &z[0]);  
  r &= isCloseTo(p[0], CF(-0.111895919f, 0.949936688f), tol);
  r &= isCloseTo(p[1], CF(-0.300119370f, 0.678187668f), tol);
  r &= isCloseTo(p[2], CF(-0.426341325f, 0.233113691f), tol);
  r &= isCloseTo(z[0], CF(-0.139246285f, 0.999132454f), tol);
  r &= isCloseTo(z[1], CF(-0.380611271f, 0.720496416f), tol);
  r &= isCloseTo(z[2], CF(-0.534630060f, 0.254959404f), tol);

  // 6th order Bessel low-cut:
  pd.setApproximationMethod(PD::BESSEL);
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);
  pd.setGain(-6);
  pd.getPolesAndZeros(&p[0], &z[0]);
  r &= isCloseTo(p[0], CF(-0.659143269f, 1.40492857f),  tol);
  r &= isCloseTo(p[1], CF(-1.15471625f,  1.03821945f),  tol);
  r &= isCloseTo(p[2], CF(-1.53899157f,  0.387341291f), tol);
  r &= isCloseTo(z[0], CF(-0.750605762f, 1.34034932f),  tol);
  r &= isCloseTo(z[1], CF(-1.11451340f,  0.783524334f), tol);
  r &= isCloseTo(z[2], CF(-1.26746202f,  0.258809835f), tol);
	
	return r;
}

bool filterSpecUnitTest()
{
  bool r = true;      // test result

  typedef RAPT::rsFilterSpecificationBA<double>  BA;
  typedef RAPT::rsFilterSpecificationZPK<double> ZPK;
  typedef std::complex<double> Complex;

  // create example ZPK-filter with 3 zeros, 2 poles and gain = 3
  Complex q1(-1.0, +2.0), q2(-1.0, -2.0), q3(-4.0, 0.0);
  Complex p1(-0.5, +0.5), p2(-0.5, -0.5);
  Complex k = 4.0;
  double inf = RS_INF(double);
  ZPK zpk32({ q1, q2, q3 }, { p1, p2}, k, inf); // sampleRate = inf -> analog filter

  // Analog case:
  //             (s-q1)*(s-q2)*(s-q3)     B0 + B1*s + B2*s^2 + B3*s^3
  // H(s) = k * ---------------------- = -----------------------------
  //             (s-p1)*(s-p2)            A0 + A1*s + A2*s^2
  //
  // multiplying out the zpk representation gives:
  Complex B0 = -k*q1*q2*q3, B1 = k*(q1*q2+q1*q3+q2*q3), B2 = -k*(q1+q2+q3), B3 = k;
  Complex A0 = p1*p2, A1 = -(p1+p2), A2 = 1;

  // we check now, if the built-in conversion-function gives the desired result:
  BA ba32 = zpk32.toBA();
  r &= ba32.b[0] == B0;
  r &= ba32.b[1] == B1;
  r &= ba32.b[2] == B2;
  r &= ba32.b[3] == B3;
  r &= ba32.a[0] == A0;
  r &= ba32.a[1] == A1;
  r &= ba32.a[2] == A2;

  // now, we convert back from ba to zpk and check, if we get our original zpk specifiction
  // properly reconstructed:
  ZPK zpkTmp = ba32.toZPK();
  r &= zpkTmp.equals(zpk32); 

  // Digital case:
  //             (1-q1/z)*(1-q2/z)*(1-q3/z)     b0 + b1/z + b2/z^2 + b3/z^3
  // H(z) = k * ---------------------------- = -----------------------------
  //             (1-p1/z)*(1-p2/z)              a0 + a1/z + a2/z^2
  // multiplying out the zpk representation gives:
  Complex b0 = k, b1 = -k*(q1+q2+q3), b2 = k*(q1*q2+q1*q3+q2*q3), b3 = -k*q1*q2*q3;
  Complex a0 = 1, a1 = -(p1+p2), a2 = p1*p2;
  // The coeffs are the same as in the analog case but in reverse order because in the digital 
  // domain, we multiply inverse powers of z (instead of regular powers of s in the analog domain)

  // We re-interpret the zpk32 specification as a digital one by setting the sample-rate to 
  // some finite value (1 in this case):
  zpk32.sampleRate = 1;
  // ...and now do the same tests as we did in the analog case:

  // ZPK -> BA:
  ba32 = zpk32.toBA();
  r &= ba32.b[0] == b0;
  r &= ba32.b[1] == b1;
  r &= ba32.b[2] == b2;
  r &= ba32.b[3] == b3;
  r &= ba32.a[0] == a0;
  r &= ba32.a[1] == a1;
  r &= ba32.a[2] == a2;

  // BA -> ZPK:
  zpkTmp = ba32.toZPK();
  r &= zpkTmp.equals(zpk32); 

  // test of conversions is done - now we evaluate the transfer-function at a couple of randomly
  // selected values for s or z an see, if both representations (ZPK and BA) give the same results:
  int numValues = 50;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(-2.0, +2.0);
  Complex z, s, H_zpk, H_ba, d;
  double tol = 1.e-10; // seems like we need a quite high tolerance - check for numeric issues
  int i;

  // digital transfer function computation:
  for(i = 0; i < numValues; i++) {
    Complex z = Complex(prng.getSample(), prng.getSample());
    H_zpk = zpk32.transferFunctionAt(z);
    H_ba  = ba32.transferFunctionAt( z);
    d     = H_zpk - H_ba;
    r &= abs(d) <= tol;
    //rsAssert(r);
  }

  // analog transfer function computation:
  zpk32.sampleRate = inf;
  //ba32.sampleRate  = inf; // this doesn't work because it doesn't reverse the coeff-arrays
  ba32 = zpk32.toBA();
  for(i = 0; i < numValues; i++) {
    Complex z = Complex(prng.getSample(), prng.getSample());
    H_zpk = zpk32.transferFunctionAt(z);
    H_ba  = ba32.transferFunctionAt( z);
    d     = H_zpk - H_ba;
    r &= abs(d) <= tol;
    //rsAssert(r);
  }

  return r;
}


// (naive) reference implementation of moving maximum "filter":
#undef min
#undef max 
std::vector<int> movingMax(const std::vector<int>& x, int L)
{
  std::vector<int> r(x.size());
  for(int i = 0; i < (int)r.size(); i++) {
    int max = x[i];
    for(int k = 1; k <= L; k++) {
      if(i-k < 0) 
        break;
      if(x[i-k] > max)
        max = x[i-k]; }
    r[i] = max; }
  return r;
}
bool testMovingMaxFilter(rsMovingMaximumFilter<int>& flt, const std::vector<int>& x, int L)
{
  std::vector<int> target = movingMax(x, L);

  std::vector<int> result(x.size());
  flt.setLength(L);
  flt.reset();
  for(size_t n = 0; n < x.size(); n++) 
    result[n] = flt.getSample(x[n]);  // remove the "Naive" later

  return result == target;
}


bool movingMaximumUnitTest()
{
  bool r = true;   

  //std::vector<int> v = { 1,2,6,8,3,2,7,3,2,6,7,2,5,8,8,4,3,1,5,7,8,4,3,2,5,7,8,6 };
  std::vector<int> v = { 1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1 };
  std::vector<int> vMax1 = movingMax(v, 1);
  std::vector<int> vMax2 = movingMax(v, 2);
  std::vector<int> vMax3 = movingMax(v, 3);
  std::vector<int> vMax4 = movingMax(v, 4);
  std::vector<int> vMax5 = movingMax(v, 5);

  rsMovingMaximumFilter<int> flt(6);
  r &= testMovingMaxFilter(flt, v, 0); 
  r &= testMovingMaxFilter(flt, v, 1);
  r &= testMovingMaxFilter(flt, v, 2);
  r &= testMovingMaxFilter(flt, v, 3);
  r &= testMovingMaxFilter(flt, v, 4);
  r &= testMovingMaxFilter(flt, v, 5);
  r &= testMovingMaxFilter(flt, v, 6);
  //r &= testMovingMaxFilter(flt, v, 7);
  //r &= testMovingMaxFilter(flt, v, 8);
  // 8 doesn't work - maybe it needs to be strictly less than capacity
  // ...maybe write a loop for these tests

  std::vector<int> vMax3_Nayuki = computeSlidingWindowMinOrMax(v, 3, true);
  // seems to compute a maximum where the window is centered over the current datapoint



  size_t N = v.size();
  size_t n;
  std::vector<int> tmp(N);

  flt.setLength(1);
  flt.reset();
  for(n = 0; n < N; n++) tmp[n] = flt.getSample(v[n]);
  r &= tmp == vMax1;


  //// test if masking works also for negative integers:
  ////int tmp, mask = 7;
  //unsigned int tmp, mask = 7;
  //tmp = 10 & mask;  // 2
  //tmp = -6 & mask;  // 2
  //tmp = 11 & mask;  // 3
  //tmp = -5 & mask;  // 3
  //// yes
  //// test, if it also works when neagtive integers are forced to unsigned int -> yes
  //// we may use unsigned int (size_t) and masking in rsDoubleEndedQueue

  return r;
}
// todo: movingMedian, movingMinimum, movingQuantile
// efficient O(log(N)) implementation: insert the new incoming value into a sorted list
// hmmm....but in a linked list, we don't have random access...maybe with memmove, it's not
// so inefficient to insert into a sorted array?
/*
                   10                  20                  30
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0  2n digit of index
1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3     v 
1 2 3 4 5 6 7 8 9 9 8 7 6 5 4 3 2 2 3 4 5 6 7 8 9 9 8 7 6 5 4     vMax2
1 2 3 4 5 6 7 8 9 9 9 9 8 8 7 5 4 3

// see here:
https://stackoverflow.com/questions/8905525/computing-a-moving-maximum
https://stackoverflow.com/questions/4802038/implement-a-queue-in-which-push-rear-pop-front-and-get-min-are-all-consta

https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm
https://www.nayuki.io/res/sliding-window-minimum-maximum-algorithm/SlidingWindowMinMax.hpp

*/


bool testMovingQuantileCore(int maxLength, int smallLength, int largeLength, int numSamples, 
  int seed = 0)
{
  rsAssert(maxLength >= smallLength + largeLength);

  bool r = true;

  rsMovingQuantileFilterCore<double> fltH;      // H for heap-based implementation
  fltH.setMaxLength(maxLength);
  fltH.setLength(smallLength + largeLength);
  fltH.setReadPosition(smallLength);

  rsMovingQuantileFilterNaive<double> fltN(smallLength, largeLength); // N for naive implementation

  // Create output signals of the naive and heap based implementation using as input signal random
  // numbers between 0 and 99 and along the way, check, if both outputs match
  using Vec = std::vector<double>;
  Vec x = rsRandomIntVector(numSamples, 0, 99, seed);
  Vec y(numSamples), z(numSamples);
  for(int n = 0; n < numSamples; n++)  {
    double q = y[n] = fltH.getSample(x[n]);
    double p = z[n] = fltN.getSample(x[n]);
    r &= p == q; }

  //rsPlotVectors(y, z);  // uncomment to see the result
  return r;
}

bool testMovingQuantileModulation()
{
  bool r = true;

  int maxLength = 15;
  int N = 200;           // number of samples

  // create vector of settings, each with a timestamp:
  struct Settings
  {
    int time;
    int nS;
    int nL;
  };
  //std::vector<Settings> settings ={ {0, 5, 7}, {40, 7, 5}, /*{80, 5, 7},*/ {120, 7, 7}, {160, 5, 5} };
  //std::vector<Settings> settings ={ {0, 7, 5}, {40, 5, 7}, {80, 7, 5}, {120, 7, 7}, {160, 5, 5} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {50, 7, 5}, {150, 5, 7} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {50, 6, 6}, {150, 5, 7} };
  std::vector<Settings> settings ={ {0, 3, 1}, {20, 2, 2}, {40, 1, 3}, {60, 3, 1} };



  rsMovingQuantileFilterNaive<double> fltN(maxLength, maxLength);
  rsMovingQuantileFilterCore<double>  fltH;      // H for heap-based implementation
  fltH.setMaxLength(maxLength);
  fltH.setModulatable(true);

  using Vec = std::vector<double>;
  Vec x = rsRandomIntVector(N, 0, 99);
  Vec y(N), z(N);
  int i = 0; // index of the settings that we choose next
  for(int n = 0; n < N; n++)
  {
    // switch settings, if desired:
    if(i < settings.size() && n == settings[i].time) {
      int nS = settings[i].nS;
      int nL = settings[i].nL;
      fltN.setLengths(nS, nL);
      fltH.setLengthAndReadPosition(nS+nL, nS);
      i++;  
    }
    y[n] = fltH.getSample(x[n]);
    z[n] = fltN.getSample(x[n]);
    int dummy = 0;
  }

  // heap-based version shows the resets to zero at the switches

  // sample 40 is still correct, at sample 41 it gets wrong and stays wrong up to sample 84 at 
  // which point it becomes correct again - unless we comment out the switch at sample 80 - in this
  // case, it stays wrong even after that - how is this possible? this far after the switch, any
  // influence of it should have died away anyway

  // todo: 
  // -check manually correctness of results, specially at the the switch samples
  // -use a heap-based filter and compare the results
  // -try longer filters with more drastic switches, say: med, min, max, very-short, very long

  // 20: 7,89,7,17,18,76,70,68,60,7,83,52 -> 7,7,7,17,18,52,60,68,70,76,83,89 -> 68

  rsPlotVectors(y, z);  // uncomment to see the result
  return r;
}


bool movingQuantileUnitTest()
{
  bool r = true;

  // Notation: nS: length of small heap, nL: length of large heap, mL: max length, L: length, 
  // q: quantile


  r &= testMovingQuantileModulation();

  int N = 500;  // number of samples
  r &= testMovingQuantileCore(64, 32, 32, N);
  r &= testMovingQuantileCore(64, 50, 14, N);  // nS + nL = 50 + 14 = 64
  r &= testMovingQuantileCore(64, 30, 14, N);
  r &= testMovingQuantileCore(64, 63,  1, N);
  r &= testMovingQuantileCore(64,  1, 63, N);
  r &= testMovingQuantileCore(64, 40,  1, N);
  r &= testMovingQuantileCore(64,  1, 40, N);

  // these do not work - the heap size of either the small or the large heap goes to zero and we
  // get index-out-of-range errors:
  //r &= testMovingQuantile(64, 40,  0, N);
  //r &= testMovingQuantile(64,  0, 40, N);
  //r &= testMovingQuantile(64, 64,  0, N);
  //r &= testMovingQuantile(64,  0, 64, N);





  return r;


  /*
  // maybe move this into an experiment later:

  int nS = 8;
  int nL = 9;
  flt.setLengths(nS, nL);

  // re-activate later:
  rsMovingQuantileFilterNaive<double> fltN(nS, nL);

  using Vec = std::vector<double>;

  int N = 200;  // number of samples

  Vec y(N), z(N), t(N);

  Vec x = rsRandomVector(N, -1.0, +1.0);

  flt.reset();
  for(int n = 0; n < N; n++)
  {
    y[n] = flt.getSample(x[n]);
    z[n] = fltN.getSample(x[n]);
  }

  for(int n = nS; n < N-nL; n++)
    t[n+nS] = rsArrayTools::median(&x[n-nS], nS+nL); 
    // why t[n+nS]?


  //rsPlotVectors(x, t, z);
  //rsPlotVectors(t, z); // match from sample 16 onwards...ok

  // now the big task is to make y match z...
  rsPlotVectors(y, z); // ..of course, this does not yet work
  */


  return r;
}