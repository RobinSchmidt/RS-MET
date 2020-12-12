
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
  double tol = 1.e-14;
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
  r &= zpkTmp.equals(zpk32, tol);

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
  r &= zpkTmp.equals(zpk32, tol);

  // test of conversions is done - now we evaluate the transfer-function at a couple of randomly
  // selected values for s or z an see, if both representations (ZPK and BA) give the same results:
  int numValues = 50;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(-2.0, +2.0);
  Complex z, s, H_zpk, H_ba, d;
  tol = 1.e-10; // seems like we need a quite high tolerance - check for numeric issues
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
    result[n] = flt.getSample(x[n]);

  return result == target;
}

bool movingMaximumUnitTest1()
{
  bool r = true;

  //std::vector<int> v = { 1,2,6,8,3,2,7,3,2,6,7,2,5,8,8,4,3,1,5,7,8,4,3,2,5,7,8,6 };
  std::vector<int> v = { 1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1 };
  std::vector<int> vMax1 = movingMax(v, 1);
  std::vector<int> vMax2 = movingMax(v, 2);
  std::vector<int> vMax3 = movingMax(v, 3);
  std::vector<int> vMax4 = movingMax(v, 4);
  std::vector<int> vMax5 = movingMax(v, 5);

  //rsMovingMaximumFilter<int> flt(6);
  rsMovingMaximumFilter<int> flt(7);
  r &= testMovingMaxFilter(flt, v, 0);
  r &= testMovingMaxFilter(flt, v, 1);
  r &= testMovingMaxFilter(flt, v, 2);
  r &= testMovingMaxFilter(flt, v, 3);
  r &= testMovingMaxFilter(flt, v, 4);
  r &= testMovingMaxFilter(flt, v, 5);
  r &= testMovingMaxFilter(flt, v, 6);

  //r &= testMovingMaxFilter(flt, v, 7);
  // this triggers an assert and fails - pushes on a full deque - why?

  //r &= testMovingMaxFilter(flt, v, 8);
  // 8 doesn't work - maybe it needs to be strictly less than capacity
  // ...maybe write a loop for these tests
  // 7 and 8 are supposed to not work because the capacity is only 6. However, with 6, there's
  // something strange going on: when we implement rsDoubleEndedQueue::isFull() as
  // getLength() >= getMaxLength(), which seems to be correct, we trigger an assert to try to
  // push onto a full deque. when we use getLength() > getMaxLength() (which is supposed to be
  // wrong), the test passes just fine - hmm. it's not surprising that pushing and popping to/from
  // the deque still works - we actually can use one memory slot more than the nominal maximum
  // length - the head/tail arithmetic doesn't care - it's just that when we do this, the length
  // computation will compute wrong results


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

template<class T>
std::vector<T> filter(const std::vector<T>& x, rsMovingMaximumFilter<T>& flt, int length)
{
  flt.setLength(length);
  flt.reset();
  int N = (int) x.size();
  std::vector<T> y(N);
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(x[n]);
  return y;
}
// later make length a type T

bool movingMaximumLengthModulation()
{
  bool r = true;

  // Under construction...
  // todo: test non-integer length and length modulation (i.e. switching the length to a new value 
  // while the filter is running):

  using Real = double;
  using Vec  = std::vector<Real>;

  //Vec x = { 1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1 }; // input
  Vec x = { 1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2 }; // input
  x = rsConcatenate(x, x);
  x = rsConcatenate(x, x);
  int N = (int) x.size();
  Vec y(N);


  rsMovingMaximumFilter<Real> flt(7);  // todo: implement and use a standard constructor without argument

  // produce max-filtered signal with different integer filter lengths:
  Vec y0 = filter(x, flt, 0);
  Vec y1 = filter(x, flt, 1);
  Vec y2 = filter(x, flt, 2);
  Vec y3 = filter(x, flt, 3);
  Vec y4 = filter(x, flt, 4);
  Vec y5 = filter(x, flt, 5);

  // produce max-filtered signal with a length switch from 2 to 5 in the middle:
  flt.setLength(2);
  flt.reset();
  Vec y_2_5(N);
  for(int n = 0; n < N/2; n++)
    y_2_5[n] = flt.getSample(x[n]);
  flt.setLength(5);
  for(int n = N/2; n < N; n++)
    y_2_5[n] = flt.getSample(x[n]);
  // this is still wrong because it resets the buffer on the switch

  // todo: do an opposite switch from 5 to 2

  // produce max-filtered signal with filter length of 4.25:
  //...


  //rsPlotVectors(x, y0, y1, y2, y3, y4, y5);
  //rsPlotVectors(x, y2, y5, y_2_5);
  return r;
}

bool movingMaximumUnitTest()
{
  bool r = true;

  r &= movingMaximumUnitTest1();
  r &= movingMaximumLengthModulation();

  return r;
}


bool testQuantileCore(int maxLength, int smallLength, int largeLength, int numSamples,
  int seed = 0)
{
  rsAssert(maxLength >= smallLength + largeLength);

  bool r = true;

  rsQuantileFilterCore<double> fltH;      // H for heap-based implementation
  fltH.setMaxLength(maxLength);
  fltH.setLength(smallLength + largeLength);
  fltH.setReadPosition(smallLength);

  rsQuantileFilterNaive<double> fltN(smallLength, largeLength); // N for naive implementation

  // Create output signals of the naive and heap based implementation using as input signal random
  // numbers between 0 and 99 and along the way, check, if both outputs match
  using Vec = std::vector<double>;
  //Vec x = rsRandomIntVector(numSamples, 0, 99, seed);
  //Vec x = rsLinearRangeVector(numSamples, 1, numSamples);
  Vec x = rsLinearRangeVector(numSamples, -1, -numSamples);
  Vec y(numSamples), z(numSamples);
  for(int n = 0; n < numSamples; n++)  {
    double q = y[n] = fltH.getSample(x[n]);
    double p = z[n] = fltN.getSample(x[n]);
    r &= p == q; }


  //rsPlotVectors(y, z);  // uncomment to see the result
  return r;
}

bool testQuantileModulation()
{
  bool r = true;

  int maxLength = 12;
  int N = 200;           // number of samples - make parameter

  // create vector of settings, each with a timestamp:
  struct Settings
  {
    int time;
    int nS;
    int nL;
  };
  //std::vector<Settings> settings ={ {0, 5, 7}, {40, 7, 5}, {80, 5, 7}, {120, 6, 4}, {160, 7, 5} };
  std::vector<Settings> settings ={ {0, 3, 2}, {40, 5, 7}, {80, 1, 2}, {120, 7, 5}, {160, 5, 3} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {50, 7, 5}, {150, 5, 7} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {50, 6, 6}, {150, 5, 7} };
  //std::vector<Settings> settings ={ {0, 3, 1}, {20, 2, 2}, {40, 1, 3}, {60, 3, 1} };
  //std::vector<Settings> settings ={ {0, 1, 3}, {6, 2, 2} };
  //std::vector<Settings> settings ={ {0, 2, 4}, {9, 3, 3} };
  //std::vector<Settings> settings ={ {0, 4, 4}, {9, 3, 5} };
  //std::vector<Settings> settings ={ {0, 1, 4}, {40, 2, 3} };
  //std::vector<Settings> settings ={ {0, 1, 2}, {20, 2, 1} };
  //std::vector<Settings> settings ={ {0, 1, 4}, {20, 3, 2} };
  //std::vector<Settings> settings ={ {0, 1, 4}, {20, 2, 3}/*, {40, 3, 1} , {60, 3, 1} */ };
  //std::vector<Settings> settings ={ {0, 2, 4}, {20, 3, 3}, {40, 3, 1} /*, {60, 3, 1} */ };
  //std::vector<Settings> settings ={ {0, 1, 2}, {20, 2, 1} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {40, 6, 6} };
  //std::vector<Settings> settings ={ {0, 5, 7}, {40, 7, 5}, {60, 5, 7} };
  //std::vector<Settings> settings ={ {0, 4, 5}, {9, 3, 2} };
  //std::vector<Settings> settings ={ {0, 3, 2} };

  // todo: use a longer list of settings covering more cases - should include edge cases as well

  rsDelayBuffer<double> rngBuf;
  rngBuf.setCapacity(maxLength);

  rsQuantileFilterNaive<double> fltN(maxLength, maxLength);
  rsQuantileFilterCore<double>  fltH;      // H for heap-based implementation
  fltH.setMaxLength(maxLength);
  fltH.setDelayBuffer(&rngBuf);   // for artifact free modulation of length

  using Vec = std::vector<double>;
  Vec x = rsRandomIntVector(N, 0, 99);
  //Vec x = rsLinearRangeVector(N, 1, N);
  //Vec x = rsLinearRangeVector(N, -1, -N);
  Vec y(N), z(N);
  int i = 0; // index of the settings that we choose next
  for(int n = 0; n < N; n++)
  {
    // switch settings, if desired:
    if(i < (int)settings.size() && n == settings[i].time) {
      int nS = settings[i].nS;
      int nL = settings[i].nL;
      fltN.setLengths(nS, nL);
      fltH.setLengthAndReadPosition(nS+nL, nS);
      rngBuf.setLength(nS+nL);
      i++;
    }
    rngBuf.getSample(x[n]);       // we must drive the modulation buffer
    y[n] = fltH.getSample(x[n]);
    z[n] = fltN.getSample(x[n]);
  }

  // todo:
  // -test edge cases
  // -do randomized tests - maybe make a function that takes a vector of settings - maybe it should
  //  be a lmabda function implemented directly here
  // -try longer filters with more drastic switches, say: med, min, max, very-short, very long

  r &= y == z;

  //rsPlotVectors(y, z);  // uncomment to see the result
  return r;
}

// Tests using a length L filter to produce an output that would have been produced by a length
// L+1 filter (that's a preliminary for supporting non-integer lengths by crossfading between a
// length L and L+1 filter - we don't want to literally run 2 filters):
bool testQuantileElongation(int L, int N)
{
  bool r = true;

  // we test the filter with these quantiles:
  using Vec = std::vector<double>;
  using QF  = RAPT::rsQuantileFilterCore<double>;
  double eps = RS_EPS(double);
  double tol = 1000*eps;
  Vec quantiles({0.0, eps, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0-eps, 1.0});
  //Vec quantiles({0.0}); // for debug - 0.5 goes int p1==p-1 branch, 0.6 into p1==p branch


  QF fltR; // reference filter
  fltR.setMaxLength(L+1);
  fltR.setLength(L+1, true);

  QF fltE; // elongated filter
  rsDelayBuffer<double> delayLine;
  delayLine.setCapacity(L+1);
  fltE.setMaxLength(L);
  fltE.setLength(L, true);
  fltE.setDelayBuffer(&delayLine);

  QF fltS; // shortened filter
  fltS.setMaxLength(L+2);
  fltS.setLength(L+2, true);
  fltS.setDelayBuffer(&delayLine);


  // compute outputs and compare them:
  Vec x = rsRandomIntVector(N, 0, 99);  // input
  Vec yR(N), yE(N), yS(N);               // reference, elongated and shortened output
  for(size_t i = 0; i < quantiles.size(); i++)
  {
    // compute output of reference filter - this filter actually is one sample longer than our
    // nominal L:
    fltR.setLengthAndQuantile(L+1, quantiles[i], true); // true: hard reset
    for(int n = 0; n < N; n++)
      yR[n] = fltR.getSample(x[n]);

    // compute output of elongated filter - this filter is set to length L and computes the output
    // of a length L+1 filter by additional trickery
    fltE.setLengthAndQuantile(L, quantiles[i], true);
    for(int n = 0; n < N; n++) {
      delayLine.getSample(x[n]);           // feed delayline (output irrelevant)
      yE[n] = fltE.getSample(x[n]);        // feed filter (output irrelevant)
      yE[n] = fltE.getElongatedOutput();   // ...this is what we are interested in
    }

    // compute output of shortened filter - this filter is set to length L+2 and computes the
    // output of a length L+1 filter by additional trickery
    fltS.setLengthAndQuantile(L+2, quantiles[i], true);
    for(int n = 0; n < N; n++) {
      yS[n] = fltS.getSample(x[n]);
      yS[n] = fltS.getShortenedOutput();
    }

    r &= rsAreVectorsEqual(yE, yR, tol);
    r &= rsAreVectorsEqual(yS, yR, tol);

    //rsPlotVectors(yR, yS);
    //rsPlotVectors(yR, yE, yS);
    //rsPlotVectors(x, yR, yE, yS);
  }

  return r;

  // Notes:
  // -We need a tolerance probably because we recompute the quantile q from L,p,w 
  //  rsQuantileFilterCore::getQuantile whereas in the reference filter, we directly use the 
  //  original value from our array in the computations for. So, that means, the w1 value in 
  //  rsQuantileFilterCore2 may differ slightly from the w value in rsQuantileFilterCore.
  // ToDo:
  // -make sure to cover all branches in readOutputWithOneMoreInput with all lengths
  // -maybe make a higher level test that continuously sweeps L and/or q and compares the result
  //  to a filter that literally uses two cores with length L and L+1
}

bool testQuantileDelay(double L, double q, int N)
{
  // tests, if the delay computation is correct by using an rsQuantileFilter instance (the 
  // high-level convenience class) with lowGain = highGain = 1 - this should result in a pure 
  // delay, if everything is right.

  bool r = true;

  int maxLength = (int) ceil(L);
  rsDelayBuffer<double> dly(maxLength);
  rsQuantileFilter<double> flt;
  flt.setSampleRateAndMaxLength(1.0, maxLength);
  flt.setFrequency(1.0/L);
  flt.setQuantile(q);
  flt.setLowpassGain(1.0);
  flt.setHighpassGain(1.0);
  flt.updateInternals();

  double d = flt.getDelayInSamples();
  r &= d == 0.5*(L-1);

  using Vec = std::vector<double>;
  Vec x(N); createWaveform(&x[0], N, 0, 1./100, 1.0);
  Vec yD(N), yF(N);              // delayed and filtered output
  for(int n = 0; n < N; n++) {
    yF[n] = flt.getSample(x[n]);
    dly.getSample(x[n]); yD[n] = dly[d]; }

  //rsPlotVectors(yD, yF);
  r &= rsIsCloseTo(yD, yF, 1.e-13);
  return r;
}

bool testQuantileSmallLengths(int N)
{
  // Tests, if rsQuantileFilterCore2 does the right thing for very small filter lengths, i.e. 
  // lengths < 2. This is an atypical case that we treat specially by crossfading between a 
  // length 2 quantile filter and the input. The length 2 quantile filter itself outputs a 
  // weighted sum between the minimum and maximum of the current and previous sample with weights
  // given by 1-q and q, where q is the quantile. So, in the case of the median (q=0.5), it 
  // becomes a 2-point moving average.

  bool r = true;
  using Vec = std::vector<double>;
  Vec x = rsRandomIntVector(N, 0, 99);
  Vec y(N), t(N);                                   // filter output and target values
  Vec quantiles({ 0.0, 0.25, 0.5, 0.75, 1.0 });
  Vec lengths(  { 1.0, 1.25, 1.5, 1.75, 1.9 });     // all must be strictly less than 2
  rsQuantileFilterCore2<double> flt;
  for(size_t i = 0; i < quantiles.size(); i++) {
    double q = quantiles[i];
    for(size_t j = 0; j < lengths.size(); j++) {
      flt.setLengthAndQuantile(lengths[j], q);
      flt.reset();
      double f   = lengths[j] - floor(lengths[j]);  // fractional part of length
      double min = rsMin(0.0, x[0]); 
      double max = rsMax(0.0, x[0]);
      t[0] = (1-q)*min  + q*max;
      t[0] = (1-f)*x[0] + f*t[0];
      y[0] = flt.getSample(x[0]);
      for(int n = 1; n < N; n++)  {
        min  = rsMin(x[n-1], x[n]); 
        max  = rsMax(x[n-1], x[n]);
        t[n] = (1-q)*min  + q*max;    // 2-value quantile filter output with quantile q...
        t[n] = (1-f)*x[n] + f*t[n];   // ...blended with input via fractional part of length
        y[n] = flt.getSample(x[n]); }
      r &= rsIsCloseTo(t, y, 1.e-13);  }}
  return r;
}

bool movingQuantileUnitTest()
{
  bool r = true;

  // Notation: nS: length of small heap, nL: length of large heap, mL: max length, L: length,
  // q: quantile

  int N = 500;  // number of samples for the tests - small value for plots - bump up later


  // test the general operation of the core:
  r &= testQuantileCore(64, 10, 10, N);
  r &= testQuantileCore(64, 20, 20, N);
  r &= testQuantileCore(64, 31, 31, N);
  r &= testQuantileCore(64, 31, 32, N);
  r &= testQuantileCore(64, 32, 31, N);
  r &= testQuantileCore(64, 32, 32, N);
  r &= testQuantileCore(64, 50, 14, N);  // nS + nL = 50 + 14 = 64
  r &= testQuantileCore(64, 30, 14, N);
  r &= testQuantileCore(64, 63,  1, N);
  r &= testQuantileCore(64,  1, 63, N);
  r &= testQuantileCore(64, 40,  1, N);
  r &= testQuantileCore(64,  1, 40, N);

  // test modulatability (setting new parameters during operation):
  r &= testQuantileModulation();

  // test the read out of a filter one sample longer than nominal length:
  //r &= testQuantileElongation(7, N);
  r &= testQuantileElongation(2, N);
  r &= testQuantileElongation(3, N);
  r &= testQuantileElongation(4, N);
  r &= testQuantileElongation(5, N);
  r &= testQuantileElongation(6, N);
  r &= testQuantileElongation(7, N);
  r &= testQuantileElongation(8, N);


  double L = 51.0; // try: 50, 51, 60, 100
  //L = 50.5;
  L = 120;  // lowpass part looks strange
  //L = 20*PI;  // results look very different
  r &= testQuantileDelay(L, 0.0, N);
  r &= testQuantileDelay(L, 0.1, N);
  r &= testQuantileDelay(L, 0.2, N);
  r &= testQuantileDelay(L, 0.4, N);
  r &= testQuantileDelay(L, 0.5, N);
  r &= testQuantileDelay(L, 0.6, N);
  r &= testQuantileDelay(L, 0.7, N);
  r &= testQuantileDelay(L, 0.8, N);
  r &= testQuantileDelay(L, 0.9, N);
  r &= testQuantileDelay(L, 1.0, N);



  r &= testQuantileSmallLengths(N);


  // try to extract the maximum over the last 8 samples:
  using Vec = std::vector<double>;
  Vec x = Vec({ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });  // input
  Vec t = Vec({ 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 });  // target output
  N = (int) x.size();
  Vec y(N);
  rsQuantileFilterCore<double>  flt;
  flt.setMaxLength(8);
  flt.setLengthAndReadPosition(8, 7);  // should read the maximum
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(x[n]);
  r &= y == t;

  // test computing the quantile from the internal algo parameters:
  auto testQuantileComputation = [&](int L, int p, double w, double q)->bool
  { 
    flt.setLengthAndReadPosition(L, p); 
    flt.setRightWeight(w);
    double qc = flt.getQuantile();
    return q == qc;
  };
  r &= testQuantileComputation(8, 1, 0.0,  0.0);  // minimum
  r &= testQuantileComputation(8, 2, 0.75, 0.25); // lower quartile
  r &= testQuantileComputation(8, 4, 0.5,  0.5);  // median
  r &= testQuantileComputation(8, 6, 0.25, 0.75); // upper quartile
  r &= testQuantileComputation(8, 7, 1.0,  1.0);  // minimum
  // the weights for the quartiles seem counterintuitive...hmmm

  return r;
}
