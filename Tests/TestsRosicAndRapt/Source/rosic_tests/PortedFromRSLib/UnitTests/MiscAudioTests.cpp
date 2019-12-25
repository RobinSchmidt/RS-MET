//#include "MiscAudioTests.h"


  // write functions for converting between various bandwidth and related variables
  // let fl: lower bandedge frequency, fu: upper bandedge frequency,
  // fc=sqrt(fl*fu): center frequency, bh=fu-fl: bandwidth in Hz, Q=bh/fc: Q-value, 
  // bo = fu/(2*fl): bandwidth in octaves
  // In "The Physics of Musical Instruments", 2nd Ed. p.12 and p.22, there are also relations 
  // given between Q and damping constant alpha and time constant tau (but check, if definition of 
  // Q matches ours)
  // Maybe the sheer mass of possible conversion functions makes it appropriate to write a class
  // rsBandwidthConversion with static functions

  // It might also be interesting to have a function 
  // double rsBiquadRingoutTime(double a1, double a2, double threshold = 1/e)
  // t = threshold, let p be the radius of a pole, then p^n = t -> n = logp(t) is the ringout
  // time for this pole, if N poles are in series (such as in a direct form filter), the total 
  // ringout time i given by the sum of the ringout times of the poles, if they are in parallel,
  // the total ringout time is given by the maximum of the individual ringout times


//void rsBandedgesTonterFreqAndBandwidthHz


bool testBandwidthConversions()
{
  bool testResult = true;

  double fl =  125;    // lower bandedge frequency   fl = fc/k 
  double fu = 8000;    // upper bandedge frequency   fu = fc*k
  double fc;           // center frequency           fc = sqrt(fl*fu) = 1000
  double ba;           // absolute bandwidth in Hz   ba = fu-fl       = 7875
  double br;           // relative bandwidth         br = ba/fc       = 7.875
  double Q;            // Q-value                    Q  = 1/br        = 0.12698412698412698
  double bo;           // bandwidth in octaves       bo = log2(fu/fl) = 6.0
  double k;            // bandedge factor            k  = br/2 + sqrt(br^2/4 + 1) = 8.0
  double tol = 1.e-14; // error tolerance

  fc = rsBandwidthConverter::bandedgesToCenterFrequency(fl, fu);
  testResult &= fc == 1000.0;

  ba = rsBandwidthConverter::bandedgesToAbsoluteBandwidth(fl, fu);
  testResult &= ba == 7875.0;

  br = rsBandwidthConverter::bandwidthAbsoluteToRelative(ba, fc);
  testResult &= br == 7.875;

  Q = rsBandwidthConverter::absoluteBandwidthToQ(ba, fc);
  testResult &= rsIsCloseTo(Q, 1/br, tol);

  k = rsBandwidthConverter::relativeBandwidthToBandedgeFactor(br);
  testResult &= k == 8.0;

  bo = rsBandwidthConverter::absoluteBandwidthToOctaves(ba, fc);
  testResult &= rsIsCloseTo(bo, 6.0, tol);

  return testResult;
}

bool testSincInterpolationAtPickedInstants(double *x, int N, int L, double s, double tol)
{
  bool testResult = true;
  double y, yt;

  yt = signalValueViaSincAt(             x, N,   2.0,  L, s);
  y  = rsResamplerDD::signalValueViaSincAt(x, N,   2.0,  L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;
  yt = signalValueViaSincAt(             x, N,   2.25, L, s);
  y  = rsResamplerDD::signalValueViaSincAt(x, N,   2.25, L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;
  yt = signalValueViaSincAt(             x, N, 100.0,  L, s);
  y  = rsResamplerDD::signalValueViaSincAt(x, N, 100.0,  L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;
  yt = signalValueViaSincAt(             x, N, 100.25, L, s);
  y  = rsResamplerDD::signalValueViaSincAt(x, N, 100.25, L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;
  yt = signalValueViaSincAt(             x, N, 197.0,  L, s); 
  y  = rsResamplerDD::signalValueViaSincAt(x, N, 197.0,  L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;
  yt = signalValueViaSincAt(             x, N, 197.25, L, s); 
  y  = rsResamplerDD::signalValueViaSincAt(x, N, 197.25, L, s);
  testResult &= rsAbs(yt-y) / yt <= tol;

  return testResult;
}

bool testSincInterpolation()
{
  bool testResult = true;

  static const int N = 200;                   // number of samples in our buffer
  double x[N];

  int L = 10;           // length of the sinc function
  double t;             // time instant to read from
  double s;             // stretch factor
  double yt;            // target value
  double y;             // computed value
  double e;             // relative error
  double eMax = 0.0;    // maximum relative error
  double tol  = 1.e-10; // error tolerance


  // create a swatooth wave with period of N samples:
  for(int n = 0; n < N; n++)
    x[n] = rsSawWave((2*PI*n)/N);

  // check error of production-code implementation with respect to reference implementation at 
  // some hand-picked time instants with various stretch values for the sinc:
  testResult &= testSincInterpolationAtPickedInstants(x, N, L, 1.0, tol);
  testResult &= testSincInterpolationAtPickedInstants(x, N, L, 1.7, tol);
  testResult &= testSincInterpolationAtPickedInstants(x, N, L, 2.0, tol);
  testResult &= testSincInterpolationAtPickedInstants(x, N, L, 3.0, tol);

  // read noise signal at random time-instants with random stretch factors
  int numTests = 100; // number of samples to be read
  L = 64;
  RAPT::rsArrayTools::fillWithRandomValues(x, N, -1.0, +1.0, 0); // create noise input signal
  for(int i = 1; i <= numTests; i++)
  {
    t  = rsRandomUniform(0.0, N);
    s  = rsRandomUniform(1.0, 10.0);
    yt = signalValueViaSincAt(x, N, t, L, s);
    y  = rsResamplerDD::signalValueViaSincAt(x, N, t, L, s);
    e  = rsAbs(yt-y) / yt;
    if( e > eMax )
      eMax = e;
  }

  testResult &= eMax <= tol;

  return testResult;
}

bool testSineParameters()
{
  bool testResult = true;

  double w = 0.1;      // normalized radian frequency
  double p = -PI;      // initial phase
  double a = 1.2;      // amplitude
  double y0, y1, y2;   // 3 successive sample values of the sine
  double w2, p2, a2;   // retrieved parameters
  double tol = 1.e-13; // error tolerance

  while( w < 2*PI )
  {
    while( p < PI )
    {
      // compute sine values:
      y0 = a * sin(p);
      y1 = a * sin(p+w);
      y2 = a * sin(p+2*w);

      // retrieve sine parameters:
      w2 = rsSineFrequency(y0, y1, y2);
      rsSineAmplitudeAndPhase(y0, y1, w2, &a2, &p2);

      // check, if retrieved parameters match original ones:
      testResult &= fabs(w-w2) <= tol;
      testResult &= fabs(p-p2) <= tol;
      testResult &= fabs(a-a2) <= tol;

      p += 0.1;
    }
    w += 0.1;
  }

  // test some special cases:
  rsSineAmplitudeAndPhase(0.0, 0.0, 0.0, &a2, &p2);   // samples zero with zero frequency
  testResult &= a2 == 0.0;
  testResult &= p2 == 0.0;
  rsSineAmplitudeAndPhase(0.0, 0.0, 0.1, &a2, &p2);   // samples zero with nonzero frequency
  testResult &= a2 == 0.0;
  testResult &= p2 == 0.0;

  // just to see, what atan2 produces in this case:
  w = atan2(0.0,  0.0);  
  testResult &= w == 0.0;

  // these values might cause an internal division by zero:
  y0 = 0.25;
  w  = 0.25;
  y1 = y0*cos(w);
  rsSineAmplitudeAndPhase(y0, y1, w, &a2, &p2); 
  testResult &= a2 == 0.0;
  testResult &= p2 == 0.0;

  return testResult;
}


bool testZeroCrossingFinder()
{
  bool testResult = true;

  static const int N = 200;  // number of samples
  double x[N];               // test signal (a sinusoid)
  double w = 0.1;            // normalized radian frequency for the sinusoid

  // create test signal:
  for(int n = 0; n < N; n++)
    x[n] = sin(w*n);


  // the zero crossings are at n = 0, 31 (downward), 62 (upward), 94, 

  int n0;
  n0 = rsFindZeroNear(x, N, 15, +1, false, true);
  testResult &= n0 == 31;
  n0 = rsFindZeroNear(x, N, 15, +1, true, true);
  testResult &= n0 == 31;
  n0 = rsFindZeroNear(x, N, 15, +1, true, false);
  testResult &= n0 == 62;

  n0 = rsFindZeroNear(x, N, 80, -1, false, true);
  testResult &= n0 == 31;
  n0 = rsFindZeroNear(x, N, 80, -1, true, true);
  testResult &= n0 == 62;
  n0 = rsFindZeroNear(x, N, 80, -1, true, false);
  testResult &= n0 == 62;

  return testResult;
}