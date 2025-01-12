
#include "../../../../../Libraries/ThirdParty/Nayuki/SlidingWindowMinMax.hpp"

using namespace std;

// maybe move somewhere else for sharing - or get rid and use rsIsCloseTo instead
bool isCloseTo(complex<float> x, complex<float> y, float tol)
{
  return abs(x-y) <= tol;
}


bool onePoleFilterUnitTest()
{
  bool ok = true;

  using Real   = float;
  using Filter = rsOnePoleFilter<Real, Real>;
  using Vec    = std::vector<Real>;

  int N = 100;

  Filter flt;
  flt.setSampleRate(44100);
  flt.setCutoff(1000);
  flt.setMode(Filter::modes::HIGHSHELV_BLT);
  flt.setShelvingGain(0.6);

  Vec h = impulseResponse(flt, N, Real(1));

  // Let's see, if we can invert the filter. This should bring us back to a unit impulse:
  flt.invert();
  Vec h2 = filterResponse(flt, N, h);
  ok &= rsIsUnitImpulse(h2, 1.e-7f);

  return ok;
}

bool directFormFilterUnitTest()
{
  bool ok = true;

  using Real   = double;
  using Filter = rsDirectFormFilter<Real, Real>;
  using Vec    = std::vector<Real>;
  using AT     = rsArrayTools;

  int maxOrder = 6;
  int order    = 4;
  int N        = 50;           // Number of samples to generate

  // Create some random vectors for the filter coefficients:
  int maxLength = maxOrder+1;
  Vec a(maxLength), b(maxLength);
  a = rsRandomVector(maxLength, -0.5, +0.5, 0);
  b = rsRandomVector(maxLength, -3.0, +3.0, 1);
  a[0] = 1.0;                  // Respect the convention!

  // Create a random signal and filter it with the rnadom coefficients using rsArrayTools::filter
  // to produce the target output:
  Vec x = rsRandomVector(N, -1.0, +1.0, 2);
  Vec yt(N);
  AT::filter(&x[0], N, &yt[0], N, &b[0], order, &a[0], order);

  // Do the same filtering with rsDirectFormFilter and compare the results:
  Filter flt(maxOrder);
  flt.setCoefficients(&a[0], &b[0], order);
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(x[n]);
  ok &= rsIsCloseTo(y, yt, 1.e-14);
  //rsPlotVectors(x, yt, y);

  // Check inversion:
  flt.invert();
  flt.reset();
  Vec z(N);
  for(int n = 0; n < N; n++)
    z[n] = flt.getSample(y[n]);
  ok &= rsIsCloseTo(z, x, 1.e-4);   // We need a very high tolerance here!
  //rsPlotVectors(x, z);   // They do look the same visually
  //rsPlotVectors(x - z);  // ..but there's an error that explodes exponentially
  // I guess maybe the inverse filter is unstable and therefore the error explodes? I guess, we 
  // should do this test with a more sensible filter - not one with random coeffs. But even then - 
  // up to a tolerance of 10^-4, it still works! So, even an unstable filter can be used for 
  // inverting a given filter for the first few samples - until the instability takes over.

  return ok;
}

bool basicFiltersUnitTests()
{
  bool ok = true;

  ok &= onePoleFilterUnitTest();
  ok &= directFormFilterUnitTest();

  return ok;
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
  r &= zpkTmp.isCloseTo(zpk32, tol);

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
  r &= zpkTmp.isCloseTo(zpk32, tol);

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

bool ladderUnitTest()
{
  bool ok = true;

  using TSig = double;  // todo: maybe test it with double, float, rsFloat64x2
  using TPar = double;
  //using LDR  = RAPT::rsLadderFilter<TSig, TPar>;
  using LDR  = rsLadderTest<TSig, TPar>;
  using RF   = RAPT::rsRationalFunction<TPar>;
  using Mode = LDR::Mode;
  using TCmp = rsComplex<TPar>;

  // ToDo: set it up and compare the results of getTransferFunctionAt, getTransferFunction and
  // getTransferFunctionOld. They should be equal up to roudoff error. This test is for making sure
  // that the creating of the transfer function works properly..

  TPar sampleRate = 44100;
  TPar nyquist    = sampleRate/2;

  LDR ldr;
  ldr.setSampleRate(sampleRate);
  ldr.setCutoff(1000);
  ldr.setResonance(0.9);
  ldr.setMode(Mode::LP_24);
  ldr.setB1(0.0);                    // allpole (no zero), cheap special case
  //ldr.setB1(0.5);                    // bilinear (zero at z = -1)
  //ldr.setB1(0.23);                   // sweet spot  


  // Compares the results of getTransferFunctionAt, getTransferFunction, getTransferFunctionOld:
  auto testTransferFunc = [&](TPar cutoff, TPar reso, Mode mode, TPar B1, TPar testFreq, 
    TPar tol1, TPar tol2)
  {
    bool ok = true;
    ldr.setup(cutoff, reso, mode, B1);
    TCmp i(0, 1);                               // imaginary unit
    TPar w    = 2*PI*testFreq / sampleRate;     // normalized radian frequency
    TCmp z    = rsExp(i*w);                     // evaluation point in the z-plane
    TCmp Hz   = ldr.getTransferFunctionAt(z);   // H(z) at z, our reference value
    RF   H1   = ldr.getTransferFunctionOld();   // H(z) via RF arithmetic...
    TCmp H1z  = H1(z);                          // ...evaluated at z
    TCmp err1 = Hz - H1z;
    ok &= rsAbs(err1) <= tol1;
    RF   H2   = ldr.getTransferFunction();      // H(z) via formulas...
    TCmp H2z  = H2(z);                          // ...evaluated at z
    TCmp err2 = Hz - H2z;
    ok &= rsAbs(err2) <= tol2;
    return ok;
  };

  // Does the above test with a couple of pre-selected settings for cutoff, reso, testFreq. The 
  // caller can still choose the mode, B1 and the tolerance.
  auto testTransferFuncs1 = [&](Mode mode, TPar B1, TPar tol1, TPar tol2)
  {
    bool ok = true;
    ok &= testTransferFunc(1000, 0.9, mode, B1,       0, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1,     100, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1,     500, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1,    1000, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1,    2000, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1,   10000, tol1, tol2);
    ok &= testTransferFunc(1000, 0.9, mode, B1, nyquist, tol1, tol2);
    return ok;
  };

  //TPar tol = RS_EPS(TPar) * 1.e9;
  // We need a high tolerance for these tests because getTransferFunctionOld is very imprecise 
  // numerically. ToDo: use two separate tolerances for old and new implementation

  // Allpole mode:
  ok &= testTransferFuncs1(Mode::FLAT,  0.0, 1.e-6, 1.e-10); // whoa! that's a large tol1!!!
  ok &= testTransferFuncs1(Mode::LP_6,  0.0, 1.e-6, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_12, 0.0, 1.e-6, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_18, 0.0, 1.e-6, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_24, 0.0, 1.e-6, 1.e-10);

  // General mode:
  ok &= testTransferFuncs1(Mode::FLAT,  0.2, 1.e-1, 1.e-10); // whoa! that's a large tol1!!!
  ok &= testTransferFuncs1(Mode::LP_6,  0.2, 1.e-1, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_12, 0.2, 1.e-1, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_18, 0.2, 1.e-1, 1.e-10);
  ok &= testTransferFuncs1(Mode::LP_24, 0.2, 1.e-1, 1.e-10);

  //ok &= testTransferFuncs1(Mode::LP_24, 0.0, 1.e-6, 1.e-10);
  //ok &= testTransferFuncs1(Mode::LP_24, 0.5, 1.e-1, 1.e-10);

  // todo: test for other filter modes and other values of B1: test at least 0.0, 0.5, 0.23

  return ok;
}


template<class Real>
bool stateVariableFilterUnitTest1(Real tol)
{
  bool ok = true;

  //using Real = double;
  using Vec  = std::vector<Real>;
  using SVF  = RAPT::rsStateVariableFilterOld<Real, Real>;
  using FDF  = rsFilterDesignFormulas;

  // Helper function to test if an SVF can faithfully emulate a biquad with the given set of 
  // coefficients:
  auto testBiquadCoeffs = [&](Real b0, Real b1, Real b2, Real a1, Real a2)
  {
    int N = 200;

    std::vector<Real> x(N), yBqd(N), ySvf(N);
    x[0] = 1;
    rsBiquadResponse(&x[0], &yBqd[0], N, b0, b1, b2, a1, a2); 

    SVF svf;
    svf.setupFromBiquad(b0, b1, b2, a1, a2);
    for(int n = 0; n < N; n++)
      ySvf[n] = svf.getSample(x[n]);

    std::vector<Real> err = yBqd - ySvf;
    Real maxErr = rsMaxAbs(err);
    bool ok = maxErr <= tol;

    if(!ok) {
      rsError("stateVariableFilterUnitTest1 failed");
      rsPlotVectors(yBqd, ySvf, err);
      rsPlotVectors(err); }

    return ok;
  };

  // Designs biquads according to specifications given by the arrays fc and Q (with given sample 
  // rate fs) and checks, if an SVF can emulate such a desiged biquad:
  auto testBiquadDesigns = [&](Real fs, const Vec& fc, const Vec& Q)
  {
    bool ok = true;
    //Real tol = 1.e-13;
    for(size_t i = 0; i < fc.size(); i++)
    {
      for(size_t j = 0; j < Q.size(); j++)
      {
        Real wc = Real(2*PI)*fc[i]/fs;
        Real b0, b1, b2, a1, a2;

        FDF::mvLowpassSimple(wc, Q[j], &b0, &b1, &b2, &a1, &a2);
        ok &= testBiquadCoeffs(b0, b1, b2, a1, a2);

        FDF::mvHighpassSimple(wc, Q[j], &b0, &b1, &b2, &a1, &a2);
        ok &= testBiquadCoeffs(b0, b1, b2, a1, a2);

        FDF::mvBandpassSimple(wc, Q[j], false, &b0, &b1, &b2, &a1, &a2);
        ok &= testBiquadCoeffs(b0, b1, b2, a1, a2);
      }
    }
    return ok;
  };


  // Test biquads with hand-picked coefficients:              // vars in SVF::setupFromBiquad
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, -0.5f, -0.5f);  // u1 =  0, u2 = -1  ->  s = inf
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, +0.5f, +0.5f);  // u1 = -2, u2 = -1  ->  s = -0.707
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, -0.5f, +0.5f);  // u1 = -1, u2 = -2  ->  s = -0.707
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, +0.5f, -0.5f);  // u1 = -1, u2 =  0  ->  s = inf
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, -1.5f,  0.0f);  // u1 = 0.5, u2 =-2.5 -> s =inf

  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, +0.2f, +0.5f);  // u1 = -1.7, u2 = -1.3
  ok &= testBiquadCoeffs(+4.0f, +0.5f, +2.0f, -0.8f, +0.9f);  // u1 = -1.1, u2 = -2.7

  ok &= testBiquadCoeffs( 0.0f, +4.0f,  0.0f, -0.8f, +0.9f);
  ok &= testBiquadCoeffs( 0.0f,  0.0f, +4.0f, -0.8f, +0.9f);
  ok &= testBiquadCoeffs(+4.0f,  0.0f, +4.0f, -0.8f, +0.9f);

  // Test designed biquads:
  ok &= testBiquadDesigns(44100.0f, 
    Vec({ 10.0f, 100.0f, 1000.0f, 10000.0f }), 
    Vec({ 0.1f, 1.0f, 10.0f, 100.0f, Real(sqrt(0.5)) }) );
    // It seems like lower cutoff frequencies require higher tolerances - which is typical for IIR
    // filters.

  // Observations:
  // -With (b0,b1,b2, a1,a2) = (4,0,0, -0.8,+0.9), the SVF output looks as if it is exactly delayed
  //  by 2 samples with respect to the reference output
  // -With (b0,b1,b2, a1,a2) = (0,0,4, -0.8,+0.9), it seems to be the other way around: the SVF is 
  //  two samples earlier
  // -With (b0,b1,b2, a1,a2) = (0,4,0, -0.8,+0.9), the output seem to match exactly. This leads to
  //  cB=0 in the SVF, i.e. the bandpass output is zero
  // -With (b0,b1,b2, a1,a2) = (4,0,4, -0.8,+0.9), we also get a perfect match, cB=0 as well
  // -Whenever (b0 - b2) = 0, cB becomes zero and in these cases, it looks good
  // -Multiplying cB by -1 seems to fix it in case 1

  // ToDo:
  // -Try extreme settings like < 0, 0, fs/2, > fs/2, > fs for the frequency in the designs. But that
  //  should püerhaps go into a different unit test - one for the rsFilterDesignFormulas class
  // -Try more tests with different settings
  // -Cover cases where u1,u2 in setupFromBiquad are ++, +-, -+, --
  //  --: -0.8,+0.9; -+: 2.0, -0.25 (unstable), ...todo...maybe write a helper function that takes
  //  the biquad coeffs and performs the test.
  // -Try single precision floats
  // -Try unstable filters (use a relative tolerance)

  return ok;
}

bool stateVariableFilterUnitTest2()
{
  // We test if the Simper-SVF implementation produces the same results as the RBJ cookbook
  // filters.

  bool ok = true;

  using Real = double;
  using Vec  = std::vector<Real>;
  using Mode = rsStateVariableFilterSimper<Real, Real>::Mode;

  // Setup:
  int  N          =   128;    // Number of samples to produce for each test case
  Real sampleRate = 44100;

  // Helper function to run a single test:
  auto runTest = [&](Mode mode, Real freq, Real Q, Real gainDb, Real tol)
  {
    // Produce the RBJ impulse response as reference signal:
    rosic::CookbookFilter cbf;
    Vec yCbf(N);
    cbf.setSampleRate(sampleRate);
    cbf.setNumStages(1); 
    cbf.setFreq(freq);
    cbf.setQ(Q);
    cbf.setGain(gainDb);
    cbf.setMode(mode);                     // Dirty! See comment below.
    getImpulseResponse(cbf, &yCbf[0], N);
    // The mode parameter is actually from the wrong enum. But the enums are compatible. That's 
    // very dirty, though!

    // Produce the SVF impulse response and compare against reference:
    rsStateVariableFilterSimper<Real, Real> svf;
    Vec ySvf(N);
    Real w = 2*PI*freq/sampleRate;
    Real A = pow(10, gainDb/40);
    svf.setup(mode, w, Q, A);
    getImpulseResponse(svf, &ySvf[0], N);
    ok &= rsIsCloseTo(ySvf, yCbf, tol);
    rsAssert(ok);

    // Plot impulse responses of SVF and RBJ and their difference in case fo failure:
    if(!ok)
      rsPlotVectors(ySvf, yCbf, yCbf - ySvf);
    // Can be uncommented when the test fails to see what's going on


    // Under construction:
    // Produce my older SVF's impulse response and check it against the reference:
    //rsStateVariableFilter<Real, Real> svf2;
    //int modeSvf2 = convertEnumMode_RBJ_to_SVF(mode);
    //svf2.setMode(modeSvf2);
    //svf2.setSampleRate(sampleRate);
    //svf2.setFrequency(freq);
    //svf2.setGain(Q);
    // ...TBC...we can re-use the ySvf array for the output
  };

  // Test different settings:
  runTest(Mode::Bypass,        1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::Lowpass,       1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::Highpass,      1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::BandpassSkirt, 1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::BandpassPeak,  1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::Notch,         1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::Allpass,       1000, 4.0, 0.0, 1.e-15);
  runTest(Mode::Bell,          1000, 4.0, 5.0, 1.e-14);
  runTest(Mode::LowShelf,      1000, 4.0, 5.0, 1.e-14);
  runTest(Mode::HighShelf,     1000, 4.0, 5.0, 1.e-14);
  //runTest(Mode::Peak,          1000, 4.0, 0.0);  // Not yet available in RBJ filter

  return ok;


  // ToDo: 
  //
  // - Test other SVF implementation as well! It's under construction. But maybe we should create
  //   a separate helper function for that. The older SVF implementation is not quite so compatible
  //   in terms of its parametrization. But maybe instead of accomodating for this incompatibility,
  //   we should actiually make it compatible!
}

bool stateVariableFilterUnitTest3()
{
  // We test the instantiation of rsStateVariableFilterSimper with different datatypes for the 
  // template parameters. Of special interest is the case where the signal type is a vector
  // version of the parameter type.

  bool ok = true;

  // Setup:
  using TPar = double;                    // Scalars for parameters
  using TSig = RAPT::rsVector2D<double>;  // Vectors for signals (i.e. multichannel)
  int   N    = 100;                       // Number of samples for test signals

  // Arrays for signals:
  std::vector<double> xL(N), xR(N), tL(N), tR(N), yL(N), yR(N);

  // Create stereo noise input signal:
  RAPT::rsNoiseGenerator<double> prng;
  for(int n = 0; n < N; n++)
  {
    xL[n] = prng.getSample();
    xR[n] = prng.getSample();
  }

  // Create target output by two independendent mono filters:
  using SVF_MN = rsStateVariableFilterSimper<TPar, TPar>;
  using ModeMn = SVF_MN::Mode;
  SVF_MN svfMnL, svfMnR;
  svfMnL.setup(ModeMn::Bell, 2.0*PI*1000.0/44100.0, 4.0, 5.0);
  svfMnR.setup(ModeMn::Bell, 2.0*PI*1000.0/44100.0, 4.0, 5.0);
  for(int n = 0; n < N; n++)
  {
    tL[n] = svfMnL.getSample(xL[n]);
    tR[n] = svfMnR.getSample(xR[n]);
  }

  // Create stereo output by a single vectorized filter:
  using SVF_ST = rsStateVariableFilterSimper<TSig, TPar>;
  using ModeSt = SVF_ST::Mode;
  SVF_ST svfSt;
  svfSt.setup(ModeSt::Bell, 2.0*PI*1000.0/44100.0, 4.0, 5.0);
  for(int n = 0; n < N; n++)
  {
    TSig x(xL[n], xR[n]);
    TSig y = svfSt.getSample(x);
    yL[n] = y.x;
    yR[n] = y.y;
  }

  // Check if the result when using one single vector/stereo filter is the same as when
  // using two independent scalar/mono filters:
  ok &= yL == tL;
  ok &= yR == tR;
  if(!ok)
  {
    rsPlotVectors(tL, yL, yL - tL);
    rsPlotVectors(tR, yR, yR - tR);
  }

  return ok;

  // ToDo:
  //
  // - Instead of using rsVector<double> use a proper SIMD vector type. We just use a normal vector
  //   type to simulate the SIMD operation at the moment. Using rosic::rsFloat64x2 doesn't compile
  //   here. I guess, rosic is not included or something.
  //
  // - Test it in full SIMD mode, i.e. with both TSig and TPar being SIMD-vector types.
}

bool stateVariableFilterUnitTest4()
{
  // We test the class rsStateVariableFilterMystran2 here.

  bool ok = true;

  using Real  = double;
  using Vec   = std::vector<Real>;
  using Mode  = rsStateVariableFilterSimper<Real, Real>::Mode;
  using ModeM = rsStateVariableFilter2<Real, Real>::Mode;

  // Setup:
  int  N          =   512;    // Number of samples to produce
  Real sampleRate = 44100;

  // Helper function to compare the outputs of the two implementations and report if thy match: 
  auto runImpRespTest = [&](Mode mode, Real freq, Real Q, Real gainDb, Real tol)
  {
    // Compute normalized radian frequency omega and linear gain:
    Real w = 2*PI*freq/sampleRate;
    Real A = pow(10, gainDb/40);

    // Create and set up filters:
    rsStateVariableFilterSimper<Real, Real>  svf_s;
    svf_s.setup(mode, w, Q, A);

    rsStateVariableFilter2<Real, Real> svf_m;
    ModeM mode_m = (ModeM)(int)mode;
    svf_m.setup(mode_m, w, Q, A);

    // Compare impusle responses:
    Vec h_s = getImpulseResponse(svf_s, N, Real(1));
    Vec h_m = getImpulseResponse(svf_m, N, Real(1));
    bool ok = rsIsCloseTo(h_s, h_m, tol);
    if(!ok)
    {
      rsPlotVectors(h_s, h_m, h_m - h_s);
      rsPlotVectors(h_s, h_m, h_m / h_s);
    }

    return ok;
  };


  // Test if outputs match the Simper SVF:
  Real tol = 1.e-13;
  ok &= runImpRespTest(Mode::Bypass,        1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::Lowpass,       1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::Highpass,      1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::BandpassSkirt, 1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::BandpassPeak,  1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::Notch,         1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::Allpass,       1000.0, 5.0, 0.0, tol);
  ok &= runImpRespTest(Mode::Bell,          1000.0, 5.0, 8.0, tol);
  ok &= runImpRespTest(Mode::LowShelf,      1000.0, 5.0, 8.0, tol);
  ok &= runImpRespTest(Mode::HighShelf,     1000.0, 5.0, 8.0, tol);


  // Test the inquiry functions
  Real freq = 1000;
  Real Q    = 7.0;
  Real w    = 2*PI*freq/sampleRate;
  Real A    = 3.0;

  Real res;
  tol = 1.e-15;

  rsStateVariableFilter2<Real, Real> svf;

  svf.setupLowpass(w, Q);
  ok &= svf.isLowpass()       == true;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupHighpass(w, Q);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == true;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupBandpassSkirt(w, Q);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == true;
  ok &= svf.isBandpassSkirt() == true;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupBandpassPeak(w, Q);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == true;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == true;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupBandstop(w, Q);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == true;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupAllpass(w, Q);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == true;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);

  svf.setupBell(w, Q, A);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == true;
  ok &= svf.isShelf()         == false;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);
  res = svf.getBellGain();       ok &= rsIsCloseTo(res, A, tol);

  svf.setupLowShelf(w, Q, A);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == true;
  ok &= svf.isLowShelf()      == true;
  ok &= svf.isHighShelf()     == false;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);
  res = svf.getShelfGain();      ok &= rsIsCloseTo(res, A, tol);

  svf.setupHighShelf(w, Q, A);
  ok &= svf.isLowpass()       == false;
  ok &= svf.isHighpass()      == false;
  ok &= svf.isBandpass()      == false;
  ok &= svf.isBandpassSkirt() == false;
  ok &= svf.isBandpassPeak()  == false;
  ok &= svf.isBandstop()      == false;
  ok &= svf.isAllpass()       == false;
  ok &= svf.isBell()          == false;
  ok &= svf.isShelf()         == true;
  ok &= svf.isLowShelf()      == false;
  ok &= svf.isHighShelf()     == true;
  res = svf.getOmega();          ok &= rsIsCloseTo(res, w, tol);
  res = svf.getQualityFactor();  ok &= rsIsCloseTo(res, Q, tol);
  res = svf.getShelfGain();      ok &= rsIsCloseTo(res, A, tol);



  // Helper function to test then conversion from svf coeffs to biquad coeffs: 
  auto runBiquadConversionTest = [&](Mode mode, Real freq, Real Q, Real gainDb, Real tol)
  {
    // Compute normalized radian frequency omega and linear gain:
    Real w = 2*PI*freq/sampleRate;
    Real A = pow(10, gainDb/40);

    // Create and set up filter:
    rsStateVariableFilter2<Real, Real> svf;
    ModeM mode_m = (ModeM)(int)mode;
    svf.setup(mode_m, w, Q, A);

    // Convert svf coeffs to biquad coeffs:
    Real b0, b1, b2, a1, a2;
    svf.convertToBiquad(&b0, &b1, &b2, &a1, &a2);

    // Produce biquad impulse response:
    Vec h_s = getImpulseResponse(svf, N, Real(1));
    Vec h_b(N);
    RAPT::rsBiquadDF1<Real, Real> bqd;
    bqd.setCoefficients(b0, b1, b2, -a1, -a2);  // ToDo: switch sign convention in RAPT::rsBiquad
    h_b[0] = bqd.getSample(1.0);
    for(int n = 1; n < N; n++)
      h_b[n] = bqd.getSample(0.0);

    // Compare both impulse responses:
    bool ok = rsIsCloseTo(h_s, h_b, tol);
    if(!ok)
      rsPlotVectors(h_s, h_b, h_b - h_s);

    // Test conversion from biquad coeffs back to svf-coeffs:
    rsStateVariableFilter2<Real, Real> svf2;
    svf2.setupFromBiquad(b0, b1, b2, a1, a2);
    ok &= svf2.hasSameCoeffsAs(svf, tol);


    return ok;
  };

  tol = 1.e-13;
  ok &= runBiquadConversionTest(Mode::Lowpass, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::Highpass, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::BandpassSkirt, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::BandpassPeak, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::Notch, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::Allpass, 1000, 8.0, 0.0, tol);
  ok &= runBiquadConversionTest(Mode::Bell, 1000, 8.0, 6.0, tol);
  ok &= runBiquadConversionTest(Mode::LowShelf, 1000, 8.0, 6.0, tol);
  ok &= runBiquadConversionTest(Mode::HighShelf, 1000, 8.0, 6.0, tol);





  int dummy = 0;

  // Helper function to design a cookbook biquad. We need this to produce the reference magnitude 
  // responses
  auto designBiquad = [&](Mode mode, Real freq, Real Q, Real gainDb,
    Real& b0, Real& b1, Real& b2, Real& a1, Real& a2)
  {
    Real A = rsDbToAmp(gainDb);
    using BD = rosic::BiquadDesigner;
    Real fsr = 1/sampleRate;
    switch(mode)
    {
    case Mode::Lowpass:
      BD::calculateCookbookLowpassCoeffs(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    case Mode::Highpass:
      BD::calculateCookbookHighpassCoeffs(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    case Mode::BandpassSkirt:
      BD::calculateCookbookBandpassConstSkirtCoeffsViaQ(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    //case Mode::BandpassPeak: 
    //  BD::calculateCookbookBandpassConstPeakCoeffsViaQ(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    //  ...such a function does not exist - maybe we should add it.
    case Mode::Notch:
      BD::calculateCookbookBandrejectCoeffsViaQ(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    case Mode::Allpass:
      BD::calculateCookbookAllpassCoeffs(b0, b1, b2, a1, a2, fsr, freq, Q); break;
    case Mode::Bell:
      BD::calculateCookbookPeakFilterCoeffsViaQ(b0, b1, b2, a1, a2, fsr, freq, Q, A); break;
    };
  };

  // Helper function to compare the calculated transfer functions between the SVF and some 
  // reference filter
  auto runTransferFuncTest = [&](Mode mode, Real freq, Real Q, Real gainDb, Real tol)
  {
    // Compute normalized radian frequency omega and linear gain:
    Real w = 2*PI*freq/sampleRate;
    Real A = pow(10, gainDb/40);

    // Create and set up filter:
    rsStateVariableFilter2<Real, Real> svf;
    ModeM mode_m = (ModeM)(int)mode;
    svf.setup(mode_m, w, Q, A);

    // Create coeffs of reference filter:
    Real b0, b1, b2, a1, a2;
    designBiquad(mode, freq, Q, gainDb, b0, b1, b2, a1, a2);

    // Compute magnitude responses of SVF and reference biquad:
    Vec ws = rsLinearRangeVector(N, 0, PI);
    Vec mag_svf(N), mag_bqd(N);
    for(int k = 0; k < N; k++)
    {
      mag_svf[k] = svf.getMagnitudeAt(ws[k]);
      mag_bqd[k] = rosic::BiquadDesigner::getBiquadMagnitudeAt(
        b0, b1, b2, a1, a2, ws[k]/(2*PI), 1.0);
    }

    // Check if they are the same. If not, we may want to look at a plot to spot the problem:
    bool ok = rsIsCloseTo(mag_svf, mag_bqd, tol);
    if(!ok)
    {
      Vec err = mag_svf - mag_bqd;
      //rsPlotVectorsXY(ws, err);
      rsPlotVectorsXY(ws, mag_svf, mag_bqd, err);
      rsPlotVectorsXY(ws, mag_svf);
      rsPlotVectorsXY(ws, mag_bqd);
    }
    return ok;
  };

  N = 1024;
  tol = 1.e-9;
  ok &= runTransferFuncTest(Mode::Lowpass, 1000, 8.0, 0.0, 1.e-9);
  ok &= runTransferFuncTest(Mode::Highpass, 1000, 8.0, 0.0, 1.e-6);
  ok &= runTransferFuncTest(Mode::BandpassSkirt, 1000, 8.0, 0.0, 1.e-9);
  //ok &= runTransferFuncTest(Mode::BandpassPeak,  1000, 8.0, 0.0, tol); // Not implemented
  ok &= runTransferFuncTest(Mode::Notch, 1000, 8.0, 0.0, 1.e-10);
  ok &= runTransferFuncTest(Mode::Allpass, 1000, 8.0, 0.0, 1.e-12);

  //ok &= runTransferFuncTest(Mode::Bell,          1000, 8.0, 6.0, 1.e-10); // Fails!
  // One of the filters has the wrong gain. I think, it's the old one. Is this parameterized
  // differently? But that seems strange.

  // We need a rather high tolerance for some types. The error is greatest around the resonance 
  // peak for the lowpass. For the highpass, the error is big at DC.




  // Temporary throw-away code (can be thrown away when code like that has been integrated into one
  // of the runTest() functions):

  // Test the conversion to biquad coeffs:
  Real b0, b1, b2, a1, a2;
  Real b0s, b1s, b2s, a1s, a2s;
  tol = 1.e-12;

  designBiquad(Mode::Lowpass, 1000, 5, 0, b0, b1, b2, a1, a2);
  svf.setupLowpass(2*PI*1000/sampleRate, 5);
  svf.convertToBiquad(&b0s, &b1s, &b2s, &a1s, &a2s);
  ok &= rsIsCloseTo(b0, b0s, tol);
  ok &= rsIsCloseTo(b1, b1s, tol);
  ok &= rsIsCloseTo(b2, b2s, tol);
  ok &= rsIsCloseTo(-a1, a1s, tol);
  ok &= rsIsCloseTo(-a2, a2s, tol);
  // The a-coeffs have a different sign but that's ok. The old code uses the other sign convention.

  designBiquad(Mode::BandpassSkirt, 1000, 5, 0, b0, b1, b2, a1, a2);
  svf.setupBandpassSkirt(2*PI*1000/sampleRate, 5);
  svf.convertToBiquad(&b0s, &b1s, &b2s, &a1s, &a2s);
  ok &= rsIsCloseTo(b0, b0s, tol);
  ok &= rsIsCloseTo(b1, b1s, tol);
  ok &= rsIsCloseTo(b2, b2s, tol);
  ok &= rsIsCloseTo(-a1, a1s, tol);
  ok &= rsIsCloseTo(-a2, a2s, tol);

  designBiquad(Mode::Highpass, 1000, 5, 0, b0, b1, b2, a1, a2);
  svf.setupHighpass(2*PI*1000/sampleRate, 5);
  svf.convertToBiquad(&b0s, &b1s, &b2s, &a1s, &a2s);
  ok &= rsIsCloseTo(b0, b0s, tol);
  ok &= rsIsCloseTo(b1, b1s, tol);
  ok &= rsIsCloseTo(b2, b2s, tol);
  ok &= rsIsCloseTo(-a1, a1s, tol);
  ok &= rsIsCloseTo(-a2, a2s, tol);

  designBiquad(Mode::Notch, 1000, 5, 0, b0, b1, b2, a1, a2);
  svf.setupBandstop(2*PI*1000/sampleRate, 5);
  svf.convertToBiquad(&b0s, &b1s, &b2s, &a1s, &a2s);
  ok &= rsIsCloseTo(b0, b0s, tol);
  ok &= rsIsCloseTo(b1, b1s, tol);
  ok &= rsIsCloseTo(b2, b2s, tol);
  ok &= rsIsCloseTo(-a1, a1s, tol);
  ok &= rsIsCloseTo(-a2, a2s, tol);


  //designBiquad(Mode::Bell, 1000, 5.0, 6.0, b0, b1, b2, a1, a2);
  //svf.setupBell(2*PI*1000/sampleRate, 5.0, 6.0);
  //svf.getBiquadCoeffs(&b0s, &b1s, &b2s, &a1s, &a2s);
  //ok &= rsIsCloseTo( b0, b0s, tol);  // !!!FAILS!!!
  //ok &= rsIsCloseTo( b1, b1s, tol);
  //ok &= rsIsCloseTo( b2, b2s, tol);
  //ok &= rsIsCloseTo(-a1, a1s, tol);
  //ok &= rsIsCloseTo(-a2, a2s, tol);
  //// Fails! Even the a-coeffs are wrong! Well, we had problems with the peak/bell mode before. The
  //// old code may have a bug or my use a different parametrization


  // Helper function to test stability of a biquad with denominator coeffs a1, a2:
  auto isBiquadStable = [](Real a1, Real a2)
  {
    using Poly = rsPolynomial<Real>;
    bool stable = Poly::areRootsOnOrInsideUnitCircle(a2, a1, 1.0);
    return stable;

    // ToDo: Document why the coeffs must be passed in reverse order. It's because the polynomial 
    // is in z^-1 rather than z, I think.
  };

  
  /*
  // Test biquad roundtrip with carefully chosen coeffs. Problems occur when 
  // T = (a1*a1 - a2*a2 - 2*a2 - 1) is positive. So, let's choose a1 = 2, a2 = 0. 
  // Then T = 4 - 0 - 0 - 1 = 3. For comparison, we also use the old implementation where I have 
  // implemented the formulas from the Wishnick paper
  rsStateVariableFilterOld<Real, Real> svf2;
  Real T;
  b0 =  1; 
  b1 =  0; 
  b2 =  0;
  //a1 =  0.0;
  //a2 = -1.01;
  a1 =  0.5;
  a2 = -0.75;
  T  = (a1*a1 - a2*a2 - 2*a2 - 1);
  //svf.setupFromBiquad( b0, b1, b2, a1, a2);
  //svf2.setupFromBiquad(b0, b1, b2, a1, a2);  // This implements the Wishnick formulas
  a1 = 2.0; 
  a2 = 0.0;
  //svf.setupFromBiquad( b0, b1, b2, a1, a2);
  //svf2.setupFromBiquad(b0, b1, b2, a1, a2);
  a1 =  0.8; 
  a2 = -0.4;
  bool stable = isBiquadStable(a1, a2); // is false -> unstable!
  svf.setupFromBiquad( b0, b1, b2, a1, a2);
  //svf2.setupFromBiquad(b0, b1, b2, a1, a2);
  dummy = 0;
  // ...This hints that the workability of the formulas may have something to do with stability,
  // after all.
  */
  
  // If a1 = 0, the function T(a2) is the parabola -x^2 - 2x - 1  which touches the x-axis at 
  // x = -1. https://www.desmos.com/calculator/ycadkwni7d  So, with a1 = 0, T can never become 
  // positive. With a2 = -1, we get T = 0. This is the boundary case between what works and what
  // doesn't work. Maybe let's try to approach the boundary. With a2 = -0.99, we are in the good 
  // range. With a2 = 1.01 we are actually also in the good case.

  // The expression for T can be factored as T = (a1 - a2 - 1) * (a1 + a2 + 1). We want that to be
  // less than zero. If we take x = a1 and y = a2, the allowed range is a double cone in the 
  // xy-plane: https://www.desmos.com/calculator/hakfgeujet
  //
  // Some "nice" pairs outside the cone:        (a1,a2) = (0.5,-0.75), (0.8,-0.4)
  // Some directly on the boundary of the cone: (a1,a2) = (0.8,-0.2)
  //
  // If we assume that a2 = 0, that means we require (a1 - 1) * (a1 + 1) < 0. I think, that means 
  // we need to have -1 < a1 < +1. That seems to correspond to the condition of stability. If we 
  // assume a1 = 0, we require (-a2 - 1) * (a2 + 1) < 0 which seems to always hold unless a2 = -1 
  // in which case we get exactly 0. With A := a2 + 1, this expression basically means -A * A < 0.
  // That's always the case except when A = 0. With that definition, the general condition for the
  // formula to work can be written as:  (a1 - A) * (a1 + A) < 0  or  (A - a1) * (A + a1) > 0
  //
  // Let's try  H(z) = (1/(1 - 0.8*d)) *  (1/(1 - 0.5*d))  where d = z^-1. That's a chain of two
  // 1-pole filters. It expands to
  // H(z) = 1 / (1 - 0.8*d - 0.5*d + 0.4*d^2) = 1 / (1 - 1.3*d + 0.4*d^2), so we have 
  // a1 = -1.3, a2 = 0.4. Therefore, T = (-1.3 - 0.4 - 1) * (-1.3 + 0.4 + 1) = -2.7 * 0.1 = -0.27
  // which is < 0, so it should work.
  //
  // Let's try  H(z) = (1/(1 + 0.8*d)) * (1/(1 + 0.5*d)). It expands to 
  // H(z) = 1 / (1 + 0.8*d + 0.5*d + 0.4*d^2) = 1 / (1 + 1.3*d + 0.4*d^2), so we have 
  // a1 = 1.3, a2 = 0.4. Therefore, T = (1.3 - 0.4 - 1) * (1.3 + 0.4 + 1) = -0.1 * 2.7 = -0.27 
  // which the same, so it's also < 0, so that should work, too.
  //
  // Let's try  H(z) = (1/(1 + 0.8*d)) * (1/(1 - 0.5*d)). It expands to 
  // H(z) = 1 / (1 + 0.8*d - 0.5*d - 0.4*d^2) = 1 / (1 + 0.3*d - 0.4*d^2), so we have 
  // a1 = 0.3, a2 = -0.4. Therefore, T = (0.3 + 0.4 - 1) * (0.3 - 0.4 + 1) = -0.3 * 0.9 = -0.27.
  // Again, the same result.
  //
  // Let's write down the stability condition: The poles, i.e. solutions of 1 + a1/z + a2/z^2 = 0
  // must be inside the unit circle. We can write this as z^2 + a1*z + a2 = 0 which has solutions
  // -a1/2 +- sqrt(a1^2/4 - a2)
  //
  // (-x/2 + sqrt(x^2/4 - y)) * (-x/2 - sqrt(x^2/4 - y))  <  1
  //
  // OK - It seems like old implementation with the Wishnick formulas also doesn't work and I think
  // The conditions for when it does and doesn't work are the same. 


  
  // Test biquad roundtrip conversions with random coeffs:
  int numTests = 1000;
  RAPT::rsNoiseGenerator<Real> prng;
  prng.setRange(-2.5, +2.5);
  int numStable  = 0;
  int numAllowed = 0;
  for(int n = 0; n < numTests; n++)
  {
    // Create random biquad coeffs:
    b0 = prng.getSample();
    b1 = prng.getSample();
    b2 = prng.getSample();
    a1 = prng.getSample();
    a2 = prng.getSample();

    // Check that stability implies that the formulas work, i.e. if it's stable, it should be 
    // allowed but not necessarily the other way around, i.e. some unstable biquads may be 
    // allowed, too:
    bool stable  = isBiquadStable(a1, a2);
    bool allowed = (a1*a1 - a2*a2 - 2*a2 - 1) < 0.0;
    if(stable)
    {
      numStable++;
      ok &= allowed;                    // If it's stable, it should always be allowed.
    }
    
    if(allowed)
    {
      // Set up an SVF from the biquad coeffs and then retrieve the biquad coeffs again, i.e. make
      // a DF -> SVF -> DF roundtrip and check that it worked:
      svf.setupFromBiquad( b0,   b1,   b2,   a1,   a2 );
      svf.convertToBiquad(&b0s, &b1s, &b2s, &a1s, &a2s);
      ok &= rsIsCloseTo(b0, b0s, tol);
      ok &= rsIsCloseTo(b1, b1s, tol);
      ok &= rsIsCloseTo(b2, b2s, tol);
      ok &= rsIsCloseTo(a1, a1s, tol);
      ok &= rsIsCloseTo(a2, a2s, tol);
      numAllowed++;
    }
  }
  // It seems the formulas always work for stable biquads. For unstable biquads, they may or may 
  // not work, I think. I think, stability is a sufficient but not necessary condition for the 
  // formulas to work. Try to find a mathematical argument why T >= 0 implies instability. Try
  // to figure out when an unstable filter is realizable by the SVF. Maybe it has to do with the
  // poles being real? Compare the formulas to the Wishnick formulas. I once tried it with 
  // numTests = 1000000000, i.e. one billion random biquads. The unit test does only 1000 because 
  // it needs to be fast, but I once did the test with a billion and it still passed. So, we have 
  // strong empirical evidence that stability implies convertibility. A mathematical proof would be
  // better, though. So, try to find one!
  //
  // In the biquadStability() experiment, it turned out that stable biquads lie in a triangle. We 
  // use a random range of -2.5...+2.5 that contains this triangle and has a bit of margin. With
  // numTests = 1000, we end up with numStable = 168 and numAllowed = 528. So we get 528 allowed
  // filters and 168 stable biquads with these settings.


  // Compare the two ways of evaluating H(z):

  tol = 1.e-13;
  svf.setupLowpass(0.5, 4.0);
  rsComplex<Real> j(0,1);
  rsComplex<Real> z = rsExp(j*1.5);
  rsComplex<Real> H1 = svf.getTransferFunctionAtOld(z);
  rsComplex<Real> H2 = svf.getTransferFunctionAt(z);
  ok &= rsAbs(H2 - H1) <= tol;

  svf.setupHighpass(0.5, 4.0);
  H1 = svf.getTransferFunctionAtOld(z);
  H2 = svf.getTransferFunctionAt(z);
  ok &= rsAbs(H2 - H1) <= tol;

  svf.setupBandpassSkirt(0.5, 4.0);
  H1 = svf.getTransferFunctionAtOld(z);
  H2 = svf.getTransferFunctionAt(z);
  ok &= rsAbs(H2 - H1) <= tol;


  // Try to make a sine oscillator using filters with infinite Q:
  Real inf = RS_INF(Real);
  svf.setupBandpassSkirt(0.1, inf);
  //plotImpulseResponse(svf, 500, 1.0);
  // Looks ok but has a gain of w = 0.1 rather than 1.0, so the output needs to be scaled by 1/w.
  // What about the start phase? The output starts at 0 but then immediately jumps up to a cosine
  // phase. Maybe to achieve a desired start phase, we need to init the integrator states 
  // accordingly. Try switching the frequency in the middle of the signal.

  // Try the roundtrip between SVF and DF with the oscillator:
  svf.convertToBiquad(&b0, &b1, &b2, &a1, &a2);  // a2 == 1
  svf.setupFromBiquad( b0,  b1,  b2,  a1,  a2);


  rsAssert(ok);
  return ok;

  // ToDo:
  //
  // - Figure out why the runTransferFuncTest(Mode::Bell, ..) fails. Test also shelving filters.
  //
  // - Try creating filters with infinite Q and roundtrip the coeffs through a biquad. Done with
  //   bandpass. Try also lowpass, highpass.
}

bool stateVariableFilterUnitTest()
{
  bool ok = true;

  ok &= stateVariableFilterUnitTest1<float>( 1.e-5f);
  ok &= stateVariableFilterUnitTest1<double>(1.e-13);
  ok &= stateVariableFilterUnitTest2();
  ok &= stateVariableFilterUnitTest3();
  ok &= stateVariableFilterUnitTest4();

  return ok;
}

bool engineersFilterUnitTest()
{
  bool ok = true;

  using Real = double;
  using EF   = RAPT::rsEngineersFilter<Real, Real>;
  using AM   = RAPT::rsPrototypeDesigner<Real>::approximationMethods;


  // We test the access violation bug that formerly affected FrequencyShifter. In certain 
  // situations, EngineersFilter wrote into memory beyond where it is allowed to:
  int order = 24;  // 24 is the highest possible value that still works

  // Create two filter objects. That will also allocate some heap memory. The behavior we want to
  // trigger is that some function calls on ef1 will mess up the heap-allocated memory of the ef2
  // object:
  EF ef1, ef2;
  // Hmm...or well, I think, it doesn't necessarily directly corrupt the heap memory of ef2 but
  // rather the a1-member of ef2 which *points* to heap memory but lives itself on the stack.

  // Retrieve the address of a1 in ef2. This seems to be some thing that gets messed up - under 
  // certain conditions:
  Real* addressPre = ef2.getAddressA1();

  // Call the problematic setup method of ef1. The hope is that if it does write beyond its owned 
  // memory, we may see this in the buffer. These calls formerly messeed with the a1 member in ef2.
  ef1.setPrototypeOrder(order);
  ef1.setApproximationMethod(AM::ELLIPTIC);
  ef1.setApproximationMethod(AM::BUTTERWORTH);

  // Check that the a1 member in ef2 is still at the same value:
  Real* addressPost = ef2.getAddressA1();
  ok &= addressPre == addressPost;
  rsAssert(ok);
  // [OLD:] OK - this triggers. Actually, all 3 calls on ef1 mess up the a1 variable in ef2. It's 
  // different after each call.

  return ok;

  // Notes:
  // -OK - it seems to work now. The bug was fixed.
}

bool hilbertFilterUnitTest()
{
  bool ok = true;

  // ToDo: Wrap this into a helper function taking numTaps as parameter and try for various numbers
  // of taps including 1,2,3:

  using WT = RAPT::rsWindowFunction::WindowType; 

  int numTaps    = 100;               // Number of taps for Hilbert filter.
  WT  window     = WT::rectangular;
  int numSamples = 2000;

  // Design the Hilbert filter:
  int M = numTaps;
  using Vec = std::vector<double>;
  Vec h(M);
  using WFD = rsWindowedFilterDesigner;
  WFD::hilbert(&h[0], numTaps, window);  // use M
  //rsPlotVectors(h); 

  // Create input signal:
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 1, 441.0, 44100.0);

  // Obtain Hilbert transform by direct convolution of x and h:
  using AT = rsArrayTools;
  int Ny = N+M-1;                // Length of convolution result y
  Vec y(Ny);                             
  rsArrayTools::convolve(&x[0], N, &h[0], numTaps, &y[0]); 
  //rsPlotVectors(x, y);

  // Obtain Hilbert transform by rsHilbertFilter:
  Vec z(Ny);
  RAPT::rsHilbertFilter<double, double> hlbFlt;
  hlbFlt.setNominalLength(M);
  hlbFlt.setWindow(window);
  for(int n = 0; n < N; n++)
    z[n] = hlbFlt.getSample(x[n]);
  for(int n = N; n < Ny; n++)
    z[n] = hlbFlt.getSample(0);
  Vec err = z-y;
  ok &= rsIsAllZeros(err);
  //rsPlotVectors(x, y, z, err);

  // Obtain complex analytic signal with proper compesation delay for the real part:
  rsComplexifier<double, double> complexifier;
  complexifier.setMaxLength(M);
  complexifier.setLength(M);
  complexifier.setSmoothing(false);
  Vec re(N), im(N);
  re = x;
  for(int n = 0; n < N; n++)
    complexifier.processSampleFrame(&re[n], &im[n]);


  //rsPlotVectors(re, im);
  // Looks correct! ToDo: Include programmatic test. Maybe that can be best done with a sine input

  // ToDo:
  // -Check the behavior with regard to delay compensation


  return ok;
}

bool delayLineUnitTest()
{
  bool ok = true;

  int delay = 5;

  rsBasicDelayLine<double> dl;
  dl.setMaxDelayInSamples(16);
  // Will be rounded up to the next power of two minus 1, i.e. if we pass 16, the max delay will
  // actually be 31. If we pass 15, it will be used as is, if we pass a value between 16 and 31,
  // 21 will be used, etc. 
  dl.setDelayInSamples(delay);

  int N = 100;  // Number of samples


  for(int i = 0; i < N; i++)
  {
    // Test computing the delay time from the taps:
    ok &= dl.getDelayInSamples() == delay;

    // Test getSample(). This triggers an update of the tapIn/tapOut pointers:
    double y = dl.getSample(double(i));
    if(i < delay)
      ok &= y == 0.0;
    else
      ok &= y == double(i-delay);

    // Test readOutput(). This should just read the same output again without triggering any other 
    // action, so we can test the exact same condition afterwards:
    y = dl.readOutput();
    if(i < delay)
      ok &= y == 0.0;
    else
      ok &= y == double(i-delay);

    // Test readOutputAt(). This should read the delayline at an arbitrary delay without triggering 
    // any further action:
    int readDelay = 8;
    y = dl.readOutputAt(readDelay);
    if(i < readDelay)
      ok &= y == 0.0;
    else
      ok &= y == double(i-readDelay);

    int dummy = 0;
  }

  return ok;
}

bool allpassChainUnitTest()
{
  // We compare impulse responses produced by a literal chain of rsAllpassDelayNaive with those of
  // class rsAllpassDelayChain. The latter should produce the same result with less memory usage 
  // and with a more convenient API.

  bool ok = true;

  using Real = double;
  using VecR = std::vector<Real>;
  using VecI = std::vector<int>;

  int numSamples = 256;

  VecI delays = { 1,    2,    3,    5,    7   };
  VecR coeffs = { 0.9, -0.8, +0.7, -0.6, +0.5 };

  // Create and set up a literal chain of allpass filters:
  int numStages = (int) delays.size();
  std::vector<rsAllpassDelayNaive<Real, Real>> literalAllpassChain(numStages);
  for(int i = 0; i < numStages; i++)
  {
    literalAllpassChain[i].setMaxDelayInSamples(delays[i]);
    literalAllpassChain[i].setDelayInSamples(   delays[i]);
    literalAllpassChain[i].setAllpassCoeff(     coeffs[i]);
  }

  // Produce the impulse response of the literal allpass chain:
  int N = numSamples;
  VecR d(N), h(N);
  d[0] = 1;
  for(int n = 0; n < N; n++)
  {
    Real tmp = d[n];
    for(int i = 0; i < numStages; i++)
      tmp = literalAllpassChain[i].getSample(tmp);
    h[n] = tmp;
  }
  ok &= isAllpass(h, 1.e-7);
  //rsPlotVectors(h);

  // Now let's see if we can produce the same result with the class rsAllpassDelayChain:
  rsAllpassDelayChain<Real, Real> allpassChain;
  allpassChain.setMaxNumStages(numStages);
  allpassChain.setNumStages(   numStages);
  for(int i = 0; i < numStages; i++)
  {
    allpassChain.setMaxDelayInSamples(i, delays[i]);
    allpassChain.setDelayInSamples(   i, delays[i]);
    allpassChain.setAllpassCoeff(     i, coeffs[i]);
  }
  VecR h2(N);
  for(int n = 0; n < N; n++)
    h2[n] = allpassChain.getSample(d[n]);
  ok &= rsIsCloseTo(h2, h, 1.e-15);
  //rsPlotVectors(h, h2);
  //rsPlotVectors(h - h2);  // Plot error 

  return ok;

  // Notes:
  //
  // - I think, we cannot expect an exact match because one implementation use direct form 1 and 
  //   the other direct form 2 which introduce different rounding errors.
}

bool nestedAllpassUnitTest()
{
  // We compare impulse responses of rsAllpassDelayNestedL1/2/3 which implement special cases in a
  // naive way with those of rsAllpassDelayNested which implements the general case in a sensible 
  // way.


  bool ok = true;

  using Real = double;
  using Vec  = std::vector<Real>;

  int numSamples = 512;

  // Create Dirac delta function d[n] to be used as input and h[n] to be used for the impulse 
  // response outputs:
  int N = numSamples;
  Vec d(N), h(N);
  d[0] = 1;

  // Set up a nested allpass with one level of nesting and get its impulse response:
  rsAllpassDelayNestedL1<Real, Real> nested1;
  nested1.setMaxDelayInSamples(17);
  nested1.setAllpassCoeff(  0, +0.8);
  nested1.setAllpassCoeff(  1, -0.9);
  nested1.setDelayInSamples(0, 11);
  nested1.setDelayInSamples(1, 17);
  Vec h1 = impulseResponse(nested1, numSamples, 1.0);
  //rsPlotVectors(h1);


  // Set up the general implementation in such a way that it produces the same result. We test 3 
  // different getSample functions. The low level implementations for the general N-stage case and
  // the unrolled implementation for the 2-stage case and finally the high-level dispatcher 
  // function which should dispatch to the unrolled implementation:
  rsAllpassDelayNested<  Real, Real> nested;
  nested.setMaxNumStages(4);
  nested.setMaxDelayInSamples(31);
  nested.setNumStages(2);                         // 2 stages means a nesting level of 1
  nested.setAllpassCoeff(  0, +0.8);
  nested.setAllpassCoeff(  1, -0.9);
  nested.setDelayInSamples(0, 11);
  nested.setDelayInSamples(1, 17);
  h = impulseResponse(nested, numSamples, 1.0);   // Uses high level dispatcher getSample()
  ok &= h1 == h;
  // I think, it may be important to call setMaxNumStages() *before*  calling 
  // setMaxDelayInSamples(). ToDo: Fix this: both calling orders should work the same way!


  nested.reset();
  for(int n = 0; n < N; n++)
    h[n] = nested.getSampleNStages(d[n]);         // General N-stage implementation
  ok &= h1 == h;

  nested.reset();
  for(int n = 0; n < N; n++)
    h[n] = nested.getSample2Stages(d[n]);         // Unrolled 2-stage implementation
  ok &= h1 == h;


  // Now with a nesting level of 2:
  rsAllpassDelayNestedL2<Real, Real> nested2;
  nested2.setMaxDelayInSamples(23);
  nested2.setAllpassCoeff(  0, +0.8);
  nested2.setAllpassCoeff(  1, -0.9);
  nested2.setAllpassCoeff(  2, +0.7);
  nested2.setDelayInSamples(0, 11);
  nested2.setDelayInSamples(1, 17);
  nested2.setDelayInSamples(2, 23);
  Vec h2 = impulseResponse(nested2, numSamples, 1.0);
  //rsPlotVectors(h2);


  // For the general filter, we only need to ramp up the number of stages by one and set the 
  // parameters for the new stage because the first two stages of the 2-level nested filter above
  // use exactly the same lengths and coeffs for the first two stages as in the previous test:
  nested.setNumStages(3);
  nested.setAllpassCoeff(  2, +0.7);
  nested.setDelayInSamples(2, 23);
  h = impulseResponse(nested, numSamples, 1.0);
  ok &= h2 == h;
  //rsPlotVectors(h2, hN);

  nested.reset();
  for(int n = 0; n < N; n++)
    h[n] = nested.getSampleNStages(d[n]);         // General N-stage implementation
  ok &= h2 == h;

  nested.reset();
  for(int n = 0; n < N; n++)
    h[n] = nested.getSample3Stages(d[n]);         // Unrolled 3-stage implementation
  ok &= h2 == h;

  // Now with a nesting level of 3:
  rsAllpassDelayNestedL3<Real, Real> nested3;
  nested3.setMaxDelayInSamples(29);
  nested3.setAllpassCoeff(  0, +0.8);
  nested3.setAllpassCoeff(  1, -0.9);
  nested3.setAllpassCoeff(  2, +0.7);
  nested3.setAllpassCoeff(  3, -0.6);
  nested3.setDelayInSamples(0, 11);
  nested3.setDelayInSamples(1, 17);
  nested3.setDelayInSamples(2, 23);
  nested3.setDelayInSamples(3, 29);
  Vec h3 = impulseResponse(nested3, numSamples, 1.0);
  //rsPlotVectors(h3);

  nested.setNumStages(4);
  nested.setAllpassCoeff(  3, -0.6);
  nested.setDelayInSamples(3, 29);
  h = impulseResponse(nested, numSamples, 1.0);
  ok &= h3 == h;
  //rsPlotVectors(h3, hN);

  // Now we make a test to ensure that the resulting impulse response is actually allpass. For 
  // this we need to use smaller coefficients such that the response decays down more quickly. 
  // Otherwise, the truncation artifacts would severly disturb the test:
  nested.setAllpassCoeff(  0, +0.08);
  nested.setAllpassCoeff(  1, -0.09);
  nested.setAllpassCoeff(  2, +0.07);
  nested.setAllpassCoeff(  3, -0.06);
  h = impulseResponse(nested, numSamples, 1.0);
  ok &= isAllpass(h, 1.e-4);
  //rsPlotVectors(h);

  return ok;


  // ToDo:
  //
  // - Test also the functions getSample2Stages/getSample3Stages of the nestedN object. They are 
  //   unrolled versions for these special cases. Maybe we shouldn't use impulseResponse() to 
  //   produce the outputs of the nestedN object. We want more control over which "getSample" 
  //   function is being called - especially later when we turn getSample into a dispatcher 
  //   function that dispatches to the unrolled (i.e. optimized?) variants form small N and 
  //   defaults to the general implementation for large N.
  //
  // - Maybe add tests that verify the the filters are actually allpasses. Do this by using an FFT
  //   and checking that the frequency response is flat. Maybe we can make a helper function
  //   isAllpass similar to impulseResponse. It should take a number of samples and a tolerance.
  //   Tolerance is needed because we cut off the infinite impulse response. So, when using less 
  //   samples we expect to need higher tolerances. These tests can be applied to all sorts of 
  //   allpass filters.
  //
  // - Add benchmarks to compare the performance of the general and specialized unrolled 
  //   implementations. I assume that the unrolled ones will be more performant but did not yet 
  //   measure it.
}

bool allpassDisperserUnitTest()
{
  // We compare the new implementation rsAllpassDisperser with rosic::rsFlatZapper. The 
  // rsAllpassDisperser is meant to provide a lower level API and gets away with less member
  // variables.

  bool ok = true;

  using Vec = std::vector<double>;

  int    numSamples =   512;
  int    numStages  =    30;
  double sampleRate = 44100;
  double fLo        =  1000;
  double fHi        =  8000;
  double fShape     =     0.7;
  double Q          =     2;

  // Set up a rosic::rsFlatZapper and let it generate the reference signal:
  using Mode = rosic::rsFlatZapper::Mode;
  rosic::rsFlatZapper flatZapper;
  flatZapper.setMode(Mode::biquad);
  flatZapper.setNumStages(numStages);
  flatZapper.setSampleRate(sampleRate);
  flatZapper.setLowFreq( fLo);
  flatZapper.setHighFreq(fHi);
  flatZapper.setFreqShape(fShape);
  flatZapper.setLowQ( Q);
  flatZapper.setHighQ(Q);
  flatZapper.setQShape(0.0);        // Irrelevant when qLo and qHi are the same
  Vec ht = impulseResponse(flatZapper, numSamples, 1.0);
  ok &= isAllpass(ht, 1.e-5);

  // Set up an rsAllpassDisperser and see, if it can produce the same output:
  rsAllpassDisperser<double, double> disperser;
  double wLo = 2*PI*fLo/sampleRate;
  double wHi = 2*PI*fHi/sampleRate;
  disperser.setMaxNumStages(100);
  disperser.setupWithTwoPoles(numStages, wLo, wHi, fShape, Q);
  Vec h = impulseResponse(disperser, numSamples, 1.0);
  //rsPlotVectors(ht, h);
  //rsPlotVectors(ht - h);
  ok &= rsIsCloseTo(ht, h, 1.e-13);

  return ok;
}

bool twoPoleAllpassDelayUnitTest()
{
  // We test the class rsTwoPoleAllpassDelay against the naive prototype implementation 
  // rsTwoPoleAllpassDelayNaive which is horribly wasteful with memory but can be verified
  // more easily by inspection. We also verify that the classes actually produce an allpass output.

  bool ok = true;

  using Vec = std::vector<double>;
  int numSamples = 512;

  rsTwoPoleAllpassDelayNaive<double, double> naive;
  naive.setMaxDelayInSamples(16);
  naive.setDelayInSamples(10);
  naive.setAllpassCoeffs(-0.7, +0.5);
  Vec ht = impulseResponse(naive, numSamples, 1.0);
  Vec mags = rsSpectralMagnitudes(ht);
  //rsPlotVectors(ht);
  //rsPlotVectors(mags);
  ok &= isAllpass(ht, 1.e-7);

  rsTwoPoleAllpassDelay<double, double> optimized;
  optimized.setMaxDelayInSamples(16);
  optimized.setDelayInSamples(10);
  optimized.setAllpassCoeffs(-0.7, +0.5);
  Vec h = impulseResponse(optimized, numSamples, 1.0);
  ok &= rsIsCloseTo(h, ht, 1.e-13);
  //rsPlotVectors(ht, h);

  return ok;

  // ToDo:
  //
  // - Figure out, if it's possible to negate the sign of the output by using negated coeffs. I 
  //   think maybe not because the coeffs are used in both feedforward and feedback path.
}

bool multiPoleAllpassDelayUnitTest()
{
  // We test the class rsMultiPoleAllpassDelay against rsMultiPoleAllpassDelayProto and 
  // rsTwoPoleAllpassDelay.

  bool ok = true;

  using Vec = std::vector<double>;
  int numSamples = 512;

  Vec c({1.0, -0.7, +0.5});
  // The first coeff 1 is a dummy that is always one. It's the a0 coeffient in biquad notation. We 
  // inlcude it to have the array indices match the math notation in the prototype implementation.
  // The production implementation doesn't need this, though.

  rsTwoPoleAllpassDelay<double, double> twoPole;
  twoPole.setMaxDelayInSamples(16);
  twoPole.setDelayInSamples(10);
  twoPole.setAllpassCoeffs(c[1], c[2]);
  Vec ht = impulseResponse(twoPole, numSamples, 1.0);
  //rsPlotVectors(ht);

  rsMultiPoleAllpassDelayProto<double, double> multiPoleProto;
  multiPoleProto.setMaxDelayInSamples(16);
  multiPoleProto.setDelayInSamples(10);
  multiPoleProto.setAllpassCoeffs(c);
  Vec hp = impulseResponse(multiPoleProto, numSamples, 1.0);
  ok &= rsIsCloseTo(hp, ht, 1.e-16);
  //rsPlotVectors(ht, hp); // Looks good!
  //rsPlotVectors(hp-ht);
  // ok &= hp == ht;  
  // I actually expected an exact match because I thought the algos should be exactly equivalent, 
  // but there seems to be a numerical difference. Figure out why!

  rsMultiPoleAllpassDelay<double, double> multiPole;
  multiPole.setMaxDelayInSamples(16);
  multiPole.setDelayInSamples(10);
  multiPole.setAllpassCoeffs(&c[1], 2);
  Vec h = impulseResponse(multiPole, numSamples, 1.0);
  ok &= h == hp;
  //rsPlotVectors(h, hp);

  // Now compare the multiPoleProto with the multiPole for a setup with a higher order prototype:
  c = Vec({1.0, -0.3, +0.2, -0.4, +0.1, -0.1});  // 5th order
  multiPoleProto.setAllpassCoeffs(c);
  multiPole.setAllpassCoeffs(&c[1], 5);
  hp = impulseResponse(multiPoleProto, numSamples, 1.0);
  h  = impulseResponse(multiPole, numSamples, 1.0);
  ok &= h == hp;
  ok &= isAllpass(h, 1.e-5);
  //rsPlotVectors(h, hp);

  return ok;

  // ToDo:
  //
  // - Figure out, why we need a tolerance for the first test. I actually assumed that in this case
  //   the algorithms should also be exactly equivalent, i.e. the two pole case just a special case
  //   of the multi pole case. Maybe some order of operations is different or something?
}


bool dampedCombAllpassUnitTest1()
{
  bool ok = true;

  using Real      = double;
  using Vec       = std::vector<Real>;
  using CombNaive = rsDampedCombAllpassNaive<Real, Real>;
  using Comb      = rsDampedCombAllpass<Real, Real>;
  using Comb_1p   = rsDampedCombAllpass_1p<Real, Real>;

  // Test parameters:
  int  d =   50;       // Main delay roundtrip length in samples
  int  N = 8192;       // Number of samples to generate
  Real w =  0.1;       // Normalized radian frequency of the low shelf for feedback damping
  Real g =  0.7;       // Linear high freq damping gain
  Real k =  0.9;       // Feedback gain factor

  // Create an instance of the naive prototype implemention, generate its impulse response and 
  // check that it is allpass in nature:
  CombNaive naive;
  naive.setMaxDelayInSamples(d);
  rsSetupHighDamp(naive, d, k, w, g, true);
  Vec h = impulseResponse(naive, N, 1.0);
  ok &= isAllpass(h, 1.e-5);

  // Now try to generate the same output with the production version of the general kind, i.e. the
  // one that allows for high order damping filters:
  Comb comb;
  comb.setMaxDelayInSamples(d);
  rsSetupHighDamp(comb, d, k, w, g, true);
  Vec h2 = impulseResponse(comb, N, 1.0);
  ok &= rsIsCloseTo(h, h2, 1.e-15);

  // Now produce the result with the production version of the trimmed down kind, i.e. the one that 
  // allows only 1st order damping filters:
  Comb_1p comb_1p;
  comb_1p.setMaxDelayInSamples(d);
  rsSetupHighDamp(comb_1p, d, k, w, g, true);
  Vec h3 = impulseResponse(comb_1p, N, 1.0);
  ok &= rsIsCloseTo(h, h2, 1.e-15);

  
  // Now do the same test again for the other mode of operation, i.e. the one without predelay:
  rsSetupHighDamp(naive,   d, k, w, g, false);
  rsSetupHighDamp(comb,    d, k, w, g, false);
  rsSetupHighDamp(comb_1p, d, k, w, g, false);
  h  = impulseResponse(naive,   N, 1.0);
  h2 = impulseResponse(comb,    N, 1.0);
  h3 = impulseResponse(comb_1p, N, 1.0);
  ok &= rsIsCloseTo(h, h2, 1.e-15);
  ok &= rsIsCloseTo(h, h3, 1.e-15);
  ok &= isAllpass(h, 1.e-5); 

  // Check the spacing of the spikes of the comb without correction and with unit feedback 
  // settings (i.e. no decay, no damping). This should produce an alternating spike train with the
  // distance between the spikes given by our desired delay. The first spike occurs at sample index  
  // n == delay - 1. 
  rsSetupHighDamp(comb, d, 1.0, 0.5, 1.0, true);
  comb.reset();
  h[0] = comb.getSampleComb(1.0);      // We use getSampleComb() - that's why impulseResponse()
  for(int n = 1; n < N; n++)           // ...can't be used
    h[n] = comb.getSampleComb(0.0);     
  for(int n = 0; n < N; n++) 
  {
    if((n+1) % d == 0)                 // Triggers at 49, 99, 149, 199, ... if delay == 50
    {
      int a = (n+1) / d;
      int b = a % 2;
      if(b == 1)
        ok &= h[n] == +1.0;
      else
        ok &= h[n] == -1.0;
    }
    else
      ok &= h[n] == 0.0;
  }

  // Now the same test without predelay and a negative k. This is much simpler because there's no 
  // alternation and weird first spike location. It's just a unipolar train of spikes at multiples
  // of the delay:
  comb.reset();
  rsSetupHighDamp(comb, d, -1.0, 0.5, 1.0, false);
  h[0] = comb.getSampleComb(1.0);
  for(int n = 1; n < N; n++) 
    h[n] = comb.getSampleComb(0.0);
  for(int n = 0; n < N; n++) 
  {
    if(n % d == 0)                     // Triggers at 0, 50, 100, 150, ... if delay == 50
      ok &= h[n] == +1.0;
    else
      ok &= h[n] ==  0.0;
  }

  return ok;

  // ToDo:
  //
  // - Testa also rsDampedCombAllpass_1p

}

bool dampedCombAllpassUnitTest2()
{
  // We test rsDampedCombAllpass with higher order damping filters. We do this by starting with
  // coefficient arrays of first order shelving filters and to get higher order filters, we 
  // iteratively convolve the coeff array of the current filter with the coeff array of the 1st 
  // order shelver. This way, the higher order filters represent series connections of shelvers 
  // of the same kind, i.e. more and more "aggressive" shelvers because the shelving gains and 
  // slopes accumulate.

  bool ok = true;

  using Real = double;
  using Vec  = std::vector<Real>;
  using Comb = rsDampedCombAllpass<Real, Real>;
  using AT   = rsArrayTools;

  // Test parameters:
  int  d =   50;       // Main delay roundtrip length in samples
  int  N = 8192;       // Number of samples to generate
  Real w =  0.5;       // Normalized radian frequency of the low shelf for feedback damping
  Real g =  0.7;       // Linear high freq damping gain
  Real k =  0.9;       // Feedback gain factor


  // Create the coeff arrays of the prototype 1st order shelver:
  Real a1[2], b1[2];
  rsMake1stOrderHighShelf(w, g, &b1[0], &b1[1], &a1[1]);
  a1[0] = 1;

  // Create arrays for filter coeffs for higher order filters and initialize them to [1 0 0 0 ...].
  // The actual coeff arrays for order i will be obtained by convolving those of order i-1 with the
  // one-pole shelver arrays. The initialization corresponds to i = 0, i.e. zeroth order.
  int maxDampOrder = Comb::getMaxDampingOrder();
  int length       = maxDampOrder+2;
  Vec a(length); a[0] = 1;
  Vec b(length); b[0] = 1;

  // Create impulse responses of damped allpasses with damping orders from 1 up to maxDampOrder
  // and check that we obtain an allpass filter. We also check that the getTransferFunctionAt()
  // member function computes the correct result:
  Comb flt;
  flt.setMaxDelayInSamples(d);
  rsComplex<Real> z(0.6, 0.8);     // On unit circle
  //rsComplex<Real> z(0.5, 0.7);   // Inside unit circle -> z^-n diverges
  //rsComplex<Real> z(0.7, 0.9);   // Outside unit circle  ->  z^-n converges to zero
  for(int i = 0; i <= maxDampOrder; i++)
  {
    // Create impulse respone of allpass comb with feedback damping order i without predelay and
    // check that we get an allpass:
    flt.setup(d, k, i, &b[0], &a[0], false);
    Vec h = impulseResponse(flt, N, 1.0);
    ok &= isAllpass(h, 1.e-7);
    ok &= testTransferFunction(flt, z, N, 1.e-8);
    //rsPlotVectors(h);

    // Now do the same test with predelay:
    flt.setup(d, k, i, &b[0], &a[0], true);
    h = impulseResponse(flt, N, 1.0);
    ok &= isAllpass(h, 1.e-7);
    ok &= testTransferFunction(flt, z, N, 1.e-7);

    // In-place convolve the current a,b arrays with the first order a1,b1 arrays:
    AT::convolve(&a[0], i+1, a1, 2, &a[0]);
    AT::convolve(&b[0], i+1, b1, 2, &b[0]);
  }

  return ok;

  // ToDo:
  //
  // - I'd really like rsDampedCombAllpass::getMaxDampingOrder(); be a static member function and
  //   then use  maxDampOrder = Comb::.getMaxDampingOrder();  rather than 
  //   maxDampOrder = flt.getMaxDampingOrder();  but it seems, I can't combine static with 
  //   constexpr. Figure out, if this is possible! But maybe we should eventually let the max 
  //   damping order also be a dynamic variable and allocate the memory for the filter coeffs and
  //   states on the heap, e.g. use std::vector for them. Then we'll need to use std::vector for 
  //   a,b, here. Not sure, if that's better. Maybe put some benchmakrs in place and measure if it 
  //   makes a difference performance wise. But each vector would then redundantly store the size
  //   (current damping order) and capacity (maximum damping order) - that just feels wrong to me.
  //   I don't know.
  //
  // - Try it with feedback filters with different numbers of poles and zeros, i.e. where either a
  //   or b has a tail of zeros
  //
  // - Integrate tests for the transfer function computations. See dampedCombAllpassTransFunc() in
  //   DelayExperiments.cpp. Also, do all the tests also with the mode with predelay.
}

// Maybe rename it to dampedCombAllpassTransFuncUnitTest
// Give it a boolean parameter that we can pass on to ap.setup() so we can do the test in both
// opertational modes

bool dampedCombAllpassUnitTest3(bool withPreDelay)
{
  // We test the computation of the transfer function in rsDampedCombAllpass, i.e. the 
  // getTransferFunctionAt(complex z) etc. methods. We use the same setup as in dampedCombAllpass4

  // Define types to be used:
  using Real    = double;
  using Vec     = std::vector<Real>;
  using Complex = rsComplex<Real>;
  using Allpass = rsDampedCombAllpass<Real, Real>;

  int  N        = 8192;
  int  delay    = 50;
  Real feedback = 0.9;

  // Design a wideband dip filter to be used as damping filter:
  Real b[3];
  Real a[3];
  a[0] = 1;
  rsStateVariableFilter<Real, Real> svf;  // Only used for designing the feedback filter
  svf.setupBell(0.2, 0.3, 0.5);           // A wideband dip filter
  svf.convertToBiquad(&b[0], &b[1], &b[2], &a[1], &a[2]);

  // Set up the filter and compute the impulse responses of the comb section, the compensation 
  // filter and the overall filter. We use a setup with predelay here:
  Allpass ap;
  ap.setMaxDelayInSamples(delay);
  ap.setup(delay, feedback, 2, b, a, withPreDelay);
  Vec d(N); d[0] = 1;                         // Unit impulse aka Dirac delta function d[n]
  Vec u(N), c(N), h(N);                       // Imp-resps of comb, corrector and allpass
  for(int n = 0; n < N; n++)
  {
    u[n] = ap.getSampleComb( d[n]);
    h[n] = ap.applyCorrector(u[n]);
  }
  ap.reset();
  for(int n = 0; n < N; n++)
    c[n] = ap.applyCorrector(d[n]);
  bool ok = isAllpass(h, 1.e-3);
  //rsPlotVectors(u, c, h);

  // Define our z-value at which we want to evaluate H(z) and compute the sequence z-^n that is 
  // needed in the z-transform:
  Complex z(0.6, 0.8);              // z = 0.6 + 0.8i is on the unit circle.
  std::vector<Complex> zn(N);       // zn[n] = z^-n = pow(z, -n)
  for(int n = 0; n < N; n++)
    zn[n] = rsPow(z, Complex(-n));

  // Evaluate the transfer functions the hard way, i.e. via the definition of the z-transform. See:
  // https://en.wikipedia.org/wiki/Z-transform#Definition  We can start the sum at n = 0 because
  // x[n] = 0 for n < 0 because the impulses responses are causal. The signal x[n] in the formula 
  // is replaced by our impulse responses u[n], c[n], h[n] in this case:
  Complex Ut = 0, Ht = 0, Ct = 0;   // The t stands for "target"
  for(int n = 0; n < N; n++)
  {
    Ut += u[n] * zn[n];             // Comb transfer function
    Ct += c[n] * zn[n];             // Corrector transfer function
    Ht += h[n] * zn[n];             // Overall allpass transfer function
  }

  // Sanity check: The magnitude of H should be 1 because z is on the unit circle and our filter 
  // is an allpass:
  Real Ha = rsAbs(Ht);
  ok &= rsIsCloseTo(Ha, 1.0, 1.e-9);

  // Compute the transfer functions using the respective methods and check if the results match the
  // naively computed target values:
  Complex Uz = ap.getCombTransferFunctionAt(z);      ok &= rsIsCloseTo(Uz, Ut, 1.e-8);
  Complex Cz = ap.getCorrectorTransferFunctionAt(z); ok &= rsIsCloseTo(Cz, Ct, 1.e-8);
  Complex Hz = ap.getTransferFunctionAt(z);          ok &= rsIsCloseTo(Hz, Ht, 1.e-8);


  // Now do the same test for the mode without predelay:
  ap.setup(delay, feedback, 2, b, a, false);
  ap.reset();
  for(int n = 0; n < N; n++)
  {
    u[n] = ap.getSampleComb( d[n]);
    h[n] = ap.applyCorrector(u[n]);
  }
  ap.reset();
  for(int n = 0; n < N; n++)
    c[n] = ap.applyCorrector(d[n]);
  ok &= isAllpass(h, 1.e-3);

  Ut = 0, Ht = 0, Ct = 0;
  for(int n = 0; n < N; n++)
  {
    Ut += u[n] * zn[n];
    Ct += c[n] * zn[n];
    Ht += h[n] * zn[n];
  }
  Ha = rsAbs(Ht);
  ok &= rsIsCloseTo(Ha, 1.0, 1.e-9);

  Uz = ap.getCombTransferFunctionAt(z);      ok &= rsIsCloseTo(Uz, Ut, 1.e-8);
  Cz = ap.getCorrectorTransferFunctionAt(z); ok &= rsIsCloseTo(Cz, Ct, 1.e-8);
  Hz = ap.getTransferFunctionAt(z);          ok &= rsIsCloseTo(Hz, Ht, 1.e-8);


  // Now do a test using our testTransferFunction() helper function. This will only test the 
  // overall getTransferFunctionAt() function not the separate partial functions 
  // getCombTransferFunctionAt(), getCorrectorTransferFunctionAt(). 
  ok &= testTransferFunction(ap, z, N, 1.e-8);

  // Test retrieving and evaluating the full transfer functions:
  rsSparseDigitalTransferFunction<Real> U, C, H;
  U = ap.getCombTransferFunction();
  C = ap.getCorrectorTransferFunction();
  H = ap.getTransferFunction();
  ok &= rsIsCloseTo(Uz, U(z), 1.e-13);
  ok &= rsIsCloseTo(Cz, C(z), 1.e-13);
  ok &= rsIsCloseTo(Hz, H(z), 1.e-13);

  // Try retrieving the transfer functions with the non-allocation methods:
  rsSparseDigitalTransferFunction<Real> Un, Cn, Hn;
  ap.getCombTransferFunction(&Un);
  ap.getCorrectorTransferFunction(&Cn);
  // ap.getTransferFunction(&Hn);  // this is yet to be written
  ok &= Un.isCloseTo(U, 1.e-13);
  ok &= Cn.isCloseTo(C, 1.e-13);
  //ok &= Hn.isCloseTo(H, 1.e-13);




  // Create and set up a rsSparseFilter object from H2 and produce its impulse response:
  rsSparseFilter<Real, Real> sp;
  sp.setMaxDelayInSamples(H.getFilterOrder());
  sp.setup(H);
  Vec h2 = impulseResponse(sp, N, 1.0);
  ok &= rsIsCloseTo(h, h2, 1.e-13);
  //rsPlotVectors(h2-h);


  return ok;
}

bool dampedSchroederAllpassUnitTest()
{
  // We test rsDampedSchroederAllpass for different feedback filter orders. This structure only
  // admits FIR filters in the feedback path. We start with the trivial filter b = { 1 } and 
  // iteratively convolve that array with a prototype filter b1 = { 0.75, 0.25 }. 

  bool ok = true;

  // Define types to be used:
  using Real     = double;
  using Complex  = rsComplex<Real>;
  using Vec      = std::vector<Real>;
  using Allpass  = rsDampedSchroederAllpass<Real, Real>;
  using AllpassN = rsDampedSchroederAllpassNaive<Real, Real>;

  // User parameters:
  int  delay      =    8;     // Delay roundtrip length in samples.
  int  numSamples =  512;     // Number of samples to generate
  Real feedback   =    0.7;   // Feedback gain factor
  int  maxOrder   =    8;     // Maximum order for feedback filter

  // For convenience:
  int  N = numSamples;
  int  M = delay;
  Real k = feedback;

  // Arrays for prototype filter b1 and the current filter b:
  Vec b1 = { 0.75, 0.25 };
  Vec b(maxOrder+1);
  b[0] = 1;

  // Create allpass filter objects:
  Allpass  ap;
  AllpassN apN;
  ap.setMaxDelayInSamples( M);
  apN.setMaxDelayInSamples(M);

  // Obtain impulse responses for various feedback filter orders and check that they are allpass in
  // nature. Check also that both implementations produce the same result. Check also that the
  // getTransferFunctionAt() functions return the correct result:
  Complex z(0.6, 0.8);
  for(int order = 0; order < maxOrder; order++)
  {
    ap.setup( M, k, order, &b[0]);
    apN.setup(M, k, order, &b[0]);
    Vec h  = impulseResponse(ap,  N, 1.0);
    Vec hN = impulseResponse(apN, N, 1.0);

    // Naive and optimized version should give same results and those results should be allpass in
    // nature:
    ok &= rsIsCloseTo(h, hN, 1.e-15);
    ok &= isAllpass(h, 1.e-7);
    //if(!ok)
    //{
    //  Vec mags = rsSpectralMagnitudes(h);
    //  rsPlotVectors(h, hN);
    //  rsPlotVectors(mags);
    //}

    // The getTransferFunctionAt() member functions should return correct results:
    ok &= testTransferFunction(ap,  z, N, 1.e-10);
    ok &= testTransferFunction(apN, z, N, 1.e-10);

    // Prepare b-array of filter coeffs for next iteration:
    rsArrayTools::convolve(&b[0], order+1, &b1[0], 2, &b[0]);
  }

  return ok;

  // ToDo:
  //
  // - Test cases where order > delay. I think, the current implementation will have problems 
  //   with that. Does such a setup even make sense? OK - test done. Indeed, when order > delay, 
  //   the response deviates from allpass. order == delay still seems to work fine, though. Maybe
  //   the class should enforce and assert this limit.
}

bool dampedAllpassBiCombUnitTest()
{
  bool ok = true;

  // Define types to be used:
  using Real         = double;
  using Complex      = rsComplex<Real>;
  using Vec          = std::vector<Real>;
  using Allpass      = rsDampedAllpassBiComb_1p<Real, Real>;
  using SparseFilter = rsSparseFilter<Real, Real>;

  // Test parameters:
  int  N      = 2048;    // Number of samples to render
  int  delay1 = 23;
  int  delay2 = 29;
  Real g1     = 0.6;
  Real g2     = 0.7;
  Real k1     = 0.9;
  Real k2     = 0.8;

  // Compute the feedback filter coeffs:
  Real b10, b11, a11; rsMake1stOrderHighShelf(0.5, 0.9, &b10, &b11, &a11);
  Real b20, b21, a21; rsMake1stOrderHighShelf(0.7, 0.8, &b20, &b21, &a21);

  // Create and set up the allpass:
  Allpass ap;
  ap.setMaxDelayInSamples(delay2);
  ap.setup(delay1, g1, k1, b10, b11, a11,
           delay2, g2, k2, b20, b21, a21);

  // Produce the impulse response of the two parallel comb filters:
  Vec hc(N);
  hc[0] = ap.getSampleCombs(1.0);
  for(int n = 1; n < N; n++)
    hc[n] = ap.getSampleCombs(0.0);

  // Let the ap convert itself into a sparse direct form filter and check that this converted 
  // filter has the same impulse response:
  SparseFilter sf;
  sf.setMaxDelayInSamples(ap.getCombSumOrder());
  ap.convertCombSumToDirectForm(&sf);
  ok &= ap.getCombSumOrder() == sf.getFilterOrder();
  Vec hc2 = impulseResponse(sf, N, 1.0);
  ok &= rsIsCloseTo(hc, hc2, 1.e-14);
  //rsPlotVectors(hc, hc2);

  // Test transfer function computation:
  Complex z(0.7, 0.8);                          // z outside unit circle - z^-n goes to 0
  //Complex H = ap.getCombTransferFunctionAt(z);
  Complex H = sf.getTransferFunctionAt(z);
  Complex Ht = 0;
  for(int n = 0; n < N; n++)
    Ht += hc[n] * rsPow(z, Complex(-n));
  ok &= rsIsCloseTo(H, Ht, 1.e-12);

  // Test the whole filter, i.e. the comb-sum with corrector applied:
  Vec h = impulseResponse(ap, N, 1.0);
  ok &= isAllpass(h, 1.e-4);

  return ok;


  // ToDo:
  //
  // - Test the different modes of operation. Maybe wrap all the tests into a helper function that
  //   takes the mode as parameter and then call that with all the modes.
  //
  // - Implement rsDampedAllpassBiComb_1p::getTransferFunctionAt() and test it.
  //
  // - Apply the inverse comb to a unit impulse just for curiosity
}

bool allpassUnitTest()
{
  bool ok = true;

  ok &= delayLineUnitTest();
  ok &= allpassChainUnitTest();
  ok &= nestedAllpassUnitTest();
  ok &= allpassDisperserUnitTest();
  ok &= twoPoleAllpassDelayUnitTest();
  ok &= multiPoleAllpassDelayUnitTest();
  ok &= dampedCombAllpassUnitTest1();
  ok &= dampedCombAllpassUnitTest2();
  ok &= dampedCombAllpassUnitTest3(false);
  ok &= dampedCombAllpassUnitTest3(true);
  ok &= dampedSchroederAllpassUnitTest();
  ok &= dampedAllpassBiCombUnitTest();

  return ok;
}





bool phonoFilterUnitTest()
{
  bool ok = true;

  using Real   = double;
  using Filter = rsPhonoFilter<Real, Real>;
  using Vec    = std::vector<Real>;

  int N = 50;

  // Generate the impulse respone of the phono filter:
  Filter flt;
  Vec h = impulseResponse(flt, N, 1.0);

  // Now switch the mode from pre-emphasis to de-emphasis. This produces the inverse filter:
  flt.setMode(Filter::DE_EMPHASIS);

  // Now check, if applying the de-emphasis filter to the impulse response of the pre-emphasis
  // filter gives us back our original unit impulse. A pre-emphasis/de-emphasis roundtrip should be
  // an identity opration:
  Vec h2 = filterResponse(flt, N, h);
  ok &= rsIsUnitImpulse(h2, 1.e-14);
  //rsPlotVectors(h, h2);

  return ok;
}

bool sparseFilterUnitTest()
{
  bool ok = true;

  using Real = double;
  using Vec  = std::vector<Real>;
  //using Poly = rsPolynomial<Real>;
  using FltD = rsDirectFormFilter<Real, Real>;   // Dense filter type
  using FltS = rsSparseFilter<Real, Real>;       // Sparse filter type

  int N = 128;

  // Define the filter coefficient arrays to be used. The sparse implementation will only store the
  // nonzero coeffs:
  Vec b({ 0.5, 0.0, -0.7, 0.0, 0.0,  0.0, 0.3});
  Vec a({ 1.0, 0.0,  0.0, 0.6, 0.1, -0.1     });

  // Helper function to set up a dense filter. We need this because the dense filter does not (yet)
  // support different orders for numerator and denominator, so we need to zero-pad the shorter
  // coeff array:
  auto setupDenseFilter = [](FltD& flt, const Vec& b, const Vec& a)
  {
    size_t size = std::max(b.size(), a.size());
    Vec B = b; B.resize(size);
    Vec A = a; A.resize(size);
    flt.setCoefficients(&A[0], &B[0], (int)size - 1);
  };

  // Create, set up and produce impulse response of dense filter:
  FltD df;
  setupDenseFilter(df, b, a);
  Vec hd = impulseResponse(df, N, 1.0);



  // Create, set up and produce impulse response of sparse filter:
  FltS sf;
  sf.setMaxDelayInSamples((int)rsMax(b.size()-1, a.size()-1));  // Verify the -1!
  sf.setupFromDenseCoeffs(b, a, 0.0);
  Vec hs = impulseResponse(sf, N, 1.0);

  // Check, if both impulse responses match:
  ok &= rsIsCloseTo(hd, hs, 1.e-14);
  //rsPlotVectors(hd, hs);

  // Check density calculation:
  double density;
  density = sf.getTransferFunction().getNumeratorDensity();
  ok &= density == double(3) / double (7);

  density = sf.getTransferFunction().getDenominatorDensity();
  ok &= density == double(3) / double (5);

  density = sf.getTransferFunction().getSeparatedDensity();
  ok &= density == double(6) / double (12);

  density = sf.getTransferFunction().getCombinedDensity();
  ok &= density == double(6) / double (13);



  // Invert the dense filter and check if applying the inverse filter to the impulse response of
  // the original filter gives back a unit impulse:
  df.invert();
  Vec y = filterResponse(df, N, hd);
  ok &= rsIsUnitImpulse(y, 1.e-14);
  //rsPlotVectors(y);

  // Try computing the inverse filter of the sparse filter without inverting the actual coeffs of 
  // the filter:
  sf.reset();
  for(int n = 0; n < N; n++)
    y[n] = sf.getSampleInverse(hs[n]);
  ok &= rsIsUnitImpulse(y, 1.e-14);
  //rsPlotVectors(y);

  // Now actually invert the filter and check that is indeed the inverse of the original one:
  sf.invert();
  y = filterResponse(sf, N, hs);
  ok &= rsIsUnitImpulse(y, 1.e-14);
  //rsPlotVectors(y);

  // Invert the filter again and see if we get the original filter back:
  sf.invert();
  hs = impulseResponse(sf, N, 1.0);
  ok &= rsIsCloseTo(hd, hs, 1.e-14);

  // Get the "phased" impulse response of the sparse filter, i.e. the one with the numerator array
  // reversed. This transformation should only affect the phase response and leave the magnitude 
  // response as is:
  Vec hps(N);
  hps[0] = sf.getSamplePhased(1.0);
  for(int n = 1; n < N; n++)
    hps[n] = sf.getSamplePhased(0.0);
  //rsPlotVectors(hs, hps);
  Vec mags  = rsSpectralMagnitudes(hs);
  Vec magsP = rsSpectralMagnitudes(hps);
  ok &= rsIsCloseTo(mags, magsP, 1.e-5); 
  // A big tolerance needed here due to truncation of the impulse response at N samples. 
  //rsPlotVectors(mags, magsP); 
  //rsPlotVectors(mags - magsP); 

  // Now do the phase transform on the actual filter coeffs and check the result:
  sf.reflectZeros();
  Vec hps2 = impulseResponse(sf, N, 1.0);
  ok &= rsIsCloseTo(hps, hps2, 1.e-6);
  //rsPlotVectors(hps, hps2, hps - hps2);
  // Why do we need such a big tolerance here? That is weird! Commenting out the num.reverse() 
  // call in reflectZeros doesn't seem to help.


  // Try inversion when b0 != 0 by introducing a predelay. In this case, we can only invert up to a
  // delay:
  int preDelay = 10;
  sf.setMaxDelayInSamples(sf.getMaxDelayInSamples() + preDelay);
  sf.setupFromDenseCoeffs(b, a, 0.0);      // Start fresh
  sf.addPreDelay(preDelay);
  hs = impulseResponse(sf, N, 1.0);

  // Check, if the newly computed hs matches our earlier hd - but with the predelay. We don't seem
  // to need a tolerance here (at least not with the Microsoft compiler):
  for(int n = 0; n < preDelay; n++)
    ok &= hs[n] == 0.0;
  for(int n = preDelay; n < N; n++)
    ok &= hs[n] == hd[n-preDelay];
  //rsPlotVectors(hd, hs);

  // Now invert the filter and apply it to the impulse response of the sparse filter hs. The result
  // should be a shifted unit impulse:
  sf.invert();
  y = filterResponse(sf, N, hs);
  ok &= rsIsShiftedUnitImpulse(y, preDelay, 1.e-14);
  //rsPlotVectors(y);


  // Check also getSampleInverse(), getSamplePhased, reflectZeros, etc. How does that work with
  // predelay?

  return ok;

  // ToDo:
  //
  // - Implement and test functions for reflection of zeros and poles in the unit circle and also 
  //   getSample functions that perform the desired transformation on the fly in the sparse and 
  //   dense implementation. Compare magnitude responses of both.
  //
  // - Test inversion when b0 == 0. In this case, we should produce a filter that inverts the 
  //   original up to a delay.
}


bool miscFiltersUnitTest()
{
  bool ok = true;

  ok &= phonoFilterUnitTest();
  ok &= sparseFilterUnitTest();

  return ok;
}

