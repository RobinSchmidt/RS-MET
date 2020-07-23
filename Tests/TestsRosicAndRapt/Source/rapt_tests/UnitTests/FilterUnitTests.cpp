
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



bool movingPercentileUnitTest()
{
  bool r = true;

  double q;

  rsMovingQuantileFilter<double> flt;
  flt.setMaxLength(4);
  flt.setLength(4);
  flt.setQuantile(2);

  q = flt.getSample( -1); r &= q ==  0;
  q = flt.getSample( -2); r &= q ==  0;
  q = flt.getSample( -3); r &= q == -1;
  q = flt.getSample( -4); r &= q == -2;
  q = flt.getSample( -5); r &= q == -3;
  q = flt.getSample( -6); r &= q == -4;
  q = flt.getSample( +1); r &= q == -4;
  q = flt.getSample( +2); r &= q == +1;
  q = flt.getSample( +3); r &= q == +2;

  // todo: go through a couple of hand-calculations and compare with results here - test also with
  // different quantiles and try to make the quantile modulatable also with 0 and L-1 which should 
  // give moving min/max filters


  // compare output against naive version for random inputs:

  int nS = 32;
  int nL = 32;
  flt.setMaxLength(nS+nL);
  //flt.setMaxLength(100);  // test - nope - produces garbage
  //flt.setMaxLength(128);  // same garbage
  //flt.setMaxLength(512);
  flt.setLength(nS+nL);
  flt.setQuantile(nS);

  // it works when nS = nL = 2^k and we use setMaxLength(nS+nL) - otherwise, we get garbage results
  // it probably has to do with the wraparounds of the circular buffer - or maybe with
  // rsRingBuffer::getIndexFromOldest - make a unit test for that
  // or could is have something to do with the heaps not being used fully? ...but i don't think so
  // no: even when we use heaps of very different sizes, it still works as long as the sum of their
  // lengths is a power of two - that points to a problem with the delay-buffer wraparound wait:
  // getIndexFromOldest returns an index *directly* into that data-buffer - not an index that can 
  // be fed to the [] operator of the class - but the values probably coincide in case of a 
  // power-of-2 length? i think so, because in this case the read and write pointers arec euqal?
  // maybe provide a buf.swapValues(int i, int j) function - yep: leftIndex == rightIndex
  // when we don't use getIndexFromOldest buf the bufIndex directly, we also get garbage - but 
  // different garbage for different values of maxLength - maybe instead of advancePointers use
  // retractPoitners going a step back

  // maybe try a std::vector instead of rsRingBuffer and use modulo



  rsMovingQuantileFilterNaive<double> fltN(nS, nL);

  double p;

  int N = 500;  // number of samples
  using Vec = std::vector<double>;
  Vec x = rsRandomIntVector(N, 0, 99, 0);
  Vec y(N), z(N), t(N);
  flt.reset();
  fltN.reset();
  for(int n = 0; n < N; n++)
  {
    q = y[n] = flt.getSample(x[n]);
    p = z[n] = fltN.getSample(x[n]);
    r &= p == q;
    int dummy = 0;
  }
  rsPlotVectors(y, z);

  // try it with unqual lengths for nS, nL that add up to a power of two
  nS = 50;  // nS+nL should be a power of two, such as 64
  nL = 14;
  flt.setMaxLength(nS+nL);
  flt.setLength(nS+nL);
  flt.setQuantile(nS);
  flt.reset();
  fltN.setLengths(nS, nL);
  fltN.reset();
  for(int n = 0; n < N; n++)
  {
    q = y[n] = flt.getSample(x[n]);
    p = z[n] = fltN.getSample(x[n]);
    r &= p == q;
    int dummy = 0;
  }
  rsPlotVectors(y, z);
  // ok - that works also - for nS+nL = 64 and nS > 0 and nL > 0, so (1,63) up to (63,1) is ok, but
  // (0,64) or (64,0) is not ok
  // Observation: when nS and nL are very different, the variance of the filter output gets less: 
  // the median has highest variance but quantiles closer to 0 or 1 have lower variance - that 
  // could be a useful feature in signal processing applications


  /*
  // now with a non-power of two:
  nS = 20;
  nL = 20;
  flt.setLength(nS+nL);
  flt.setQuantile(nS);
  fltN.setLengths(nS, nL);
  flt.reset();
  fltN.reset();
  for(int n = 0; n < N; n++)
  {
    q = y[n] = flt.getSample(x[n]);
    p = z[n] = fltN.getSample(x[n]);
    r &= p == q;
    //rsAssert(p == q); // triggers at sample 6
    int dummy = 0;
  }
  rsPlotVectors(y, z);
  */




  return r;

  // stuff below may be obsolete soon:

  /*
  q = flt.getSample( -3); r &= q ==  0;
  q = flt.getSample( -4); r &= q ==  0;
  q = flt.getSample( -5); r &= q == -1;
  q = flt.getSample( -6); r &= q == -2;
  q = flt.getSample( -7); r &= q == -3;
  q = flt.getSample( -8); r &= q == -4;
  q = flt.getSample( -9); r &= q == -5;
  q = flt.getSample( +1); r &= q == -4;

  // here it gets false - but the error may be in the test code:
  q = flt.getSample( +2); r &= q == -3;
  q = flt.getSample( +3); r &= q == -2;
  q = flt.getSample( +4); r &= q == -3;
  q = flt.getSample( +5); r &= q == -2;
  q = flt.getSample( +6); r &= q == -1;
  */



  /*
  rsMovingQuantileFilter<double> flt(3, 4);  // length 7
                                        // heaps (desired)
  q = flt.getSample( -1); r &= q == 0;  //  0  0 -1 | 0 0 0 0
  q = flt.getSample( -2); r &= q == 0;  //  0 -1 -2 | 0 0 0 0 
  q = flt.getSample( -3); r &= q == 0;  // -1 -2 -3 | 0 0 0 0
  q = flt.getSample( -4); r &= q == -1;

  // triggers assert - 
  q = flt.getSample( -5); r &= q == -2;
  q = flt.getSample( -6); r &= q == -3;
  q = flt.getSample( -7); r &= q == -4;
  q = flt.getSample( -8); r &= q == -5;
  */





  /*
  flt.setLengths(8, 9);
  flt.reset();

  q = flt.getSample(-1); r &= q == 0;
  // the -1 should float down to the bottom of the max-heap of small values

  q = flt.getSample(-2); r &= q == 0;
  // gets inserted into the larger heap and floats down - it should actually remain in the pole
  // position for..ah no - it is in the pole position

  q = flt.getSample( -3); r &= q == 0;
  q = flt.getSample( -4); r &= q == 0;
  q = flt.getSample( -5); r &= q == 0;
  q = flt.getSample( -6); r &= q == 0;
  q = flt.getSample( -7); r &= q == 0;
  q = flt.getSample( -8); r &= q == 0;

  q = flt.getSample( -9); r &= q == -1;     // fails here
  q = flt.getSample(-10); r &= q == -2;     // triggers assert here
  q = flt.getSample(-11); r &= q == -3;
  q = flt.getSample(-12); r &= q == -4;
  q = flt.getSample(-13); r &= q == -5;
  q = flt.getSample(-14); r &= q == -6;
  q = flt.getSample(-15); r &= q == -7;
  q = flt.getSample(-16); r &= q == -8;
  q = flt.getSample(-17); r &= q == -9;
  */








  /*

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