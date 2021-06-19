using namespace rotes;
//#include "rosic/rosic.h"  // superfluous?
using namespace rosic;
using namespace RAPT;

/*
from https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform
def fwht(a) -> None:
    """In-place Fast Walsh–Hadamard Transform of array a."""
    h = 1
    while h < len(a):
        for i in range(0, len(a), h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2
*/
template<class T>
void fght(T* A, int N, T a, T b, T c, T d)
{
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 2*h) {
      for(int j = i; j < i+h; j++) {
        T x = A[j];
        T y = A[j+h];
        A[j]   = a*x + b*y;
        A[j+h] = c*x + d*y;  }}
    h *= 2;  }
}
// maybe make N a template parameter such that the compiler can unroll the loops

bool rotes::testFastGeneralizedHadamardTransform()
{
  bool ok = true;

  // 4D vector:
  double x4[4] = {4, -8, 12, -4};
  double y4[4];
  double work[8];  // workspace

  typedef rosic::FeedbackDelayNetwork FDN;
  typedef RAPT::rsArrayTools AT;

  // 4-point FWHT:
  AT::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work);
  ok &= y4[0] ==   4;
  ok &= y4[1] ==  28;
  ok &= y4[2] == -12;
  ok &= y4[3] == - 4;

  // 4-point FGHT:
  AT::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work, 2, 3, 5, 7);
  ok &= y4[0] ==  4;
  ok &= y4[1] == 24;
  ok &= y4[2] ==  4;
  ok &= y4[3] == 44;

  // New implementation:
  AT::copy(x4, y4, 4);
  fght(y4, 4, 2., 3., 5., 7.);
  ok &= y4[0] ==  4;
  ok &= y4[1] == 24;
  ok &= y4[2] ==  4;
  ok &= y4[3] == 44;


  // 8D vector:
  double x8[8] = {1, 4, -2, 3, 0, 1, 4, -1};
  double y8[8];

  // 8-point FWHT:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work);
  ok &= y8[0] ==  10;
  ok &= y8[1] == - 4;
  ok &= y8[2] ==   2;
  ok &= y8[3] == - 4;
  ok &= y8[4] ==   2;
  ok &= y8[5] == -12;
  ok &= y8[6] ==   6;
  ok &= y8[7] ==   8;

  // 8-point FGHT:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, 7);
  ok &= y8[0] ==   149;
  ok &= y8[1] ==   357;
  ok &= y8[2] ==   360;
  ok &= y8[3] ==   862;
  ok &= y8[4] ==   362;
  ok &= y8[5] ==   866;
  ok &= y8[6] ==   875;
  ok &= y8[7] ==  2092;

  // New implementation:
  AT::copy(x8, y8, 8);
  fght(y8, 8, 2., 3., 5., 7.);
  ok &= y8[0] ==   149;
  ok &= y8[1] ==   357;
  ok &= y8[2] ==   360;
  ok &= y8[3] ==   862;
  ok &= y8[4] ==   362;
  ok &= y8[5] ==   866;
  ok &= y8[6] ==   875;
  ok &= y8[7] ==  2092;


  // 8D Forward/backward trafo - check if input is reconstructed:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(       y8, 8, 3, work, 2, 3, 5, -7);
  FDN::fastInverseGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, -7);
  ok &= fabs(RAPT::rsArrayTools::maxDeviation(x8, y8, 8)) < 1.e-15;

  return ok;

  // ToDo: 
  // -move the FGHT to rapt and the test to the rapt tests
  // -compare results of bigger sizes with explicit matrix multiplication of matrices created by 
  //  the Sylvester construction (use the Kronecker-product in rsMatrix for this)
}

bool rotes::testFeedbackDelayNetwork()
{
  bool result = true;  // get rid! this is not a unit test!

  //FeedbackDelayNetwork16 *fdn16 = new FeedbackDelayNetwork16;

  double amplitude = 0.5;       // amplitude of the input impulse
  double diffusion = 100.0;     // diffusion parametr in percent


  using AT = RAPT::rsArrayTools;

  FeedbackDelayNetwork fdn;
  fdn.setDiffusion(diffusion);

  static const int N = 100000;
  //double hL[N], hR[N];
  double *hL = new double[N];
  double *hR = new double[N];


  AT::fillWithImpulse(hL, N, amplitude);
  AT::fillWithImpulse(hR, N, amplitude);
  for(int n = 0; n < N; n++)
    fdn.processFrame(&hL[n], &hR[n]);

  double t[N];
  AT::fillWithIndex(t, N);
  //Plotter::plotData(N, t, hL, hR);
  //Plotter::plotData(N, t, hL);


  //rosic::writeToStereoWaveFile("d:\\TmpData\\FDNImpulseResponse.wav", hL, hR, N, 44100, 16);
  rosic::writeToStereoWaveFile("FDNImpulseResponse.wav", hL, hR, N, 44100, 16);

  delete[] hL;
  delete[] hR;
  return result;

  // Observations:
  // -With diffusion set to 100, we indeed get some sort of exponentially decaying white noise, as
  //  it should be
  // -The output seems to contain the input impulse, i.e. it's 100% dry + 100% wet...is that 
  //  correct?
  // -The diffusion parameter seems to need a nonlinear mapping that gives more precision towards
  //  higher values - between D=60 and D=100, there is not so much difference

  // ToDo:
  // -Figure out the distribution of the noise (we need an infinite decay-time for this). It is 
  //  written in the literature that exponentially decaying Gaussian white noise sounds best for
  //  a reverb impulse response. Figure out, if it is indeed Gaussian. Maybe it's an Irwin-Hall 
  //  distribution of order equal to the number of the delaylines? That would seem to make some 
  //  sense.
  // -Maybe we can render Gaussian white noise impulse responses with time-variant slope filters
  //  whose slope increases over time. Maybe we can make a convolution reverb that internally 
  //  renders impulse response according to that idea.
}

template<class Effect>
bool testInOutEqual(Effect& eff, int numSamples, double tolerance)
{
  bool result = true;
  RAPT::rsNoiseGenerator<double> ng;
  double xL, xR, yL, yR;
  for(int n = 0; n < numSamples; n++)
  {
    yL = xL = ng.getSample();
    yR = xR = ng.getSample();
    eff.getSampleFrameStereo(&yL, &yR);
    result &= RAPT::rsIsCloseTo(yL, xL, tolerance);
    result &= RAPT::rsIsCloseTo(yR, xR, tolerance);
    rsAssert(result == true);
  }
  return result;
}

bool rotes::testMultiComp()
{  
  // We check the multiband band compressor with neutral compressor settings for each band and 
  // various splitting configurations. In any case, the output signal should equal the input
  // signal.

  bool result = true;

  int N = 200;  // number of samples

  rosic::rsMultiBandCompressor mbc;
  double tol = 1.e-14;

  result &= mbc.getNumberOfBands() == 1;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 8000.0);  // maybe pass only the frequency, not the index, 
                              // maybe rename addSplit, maybe define splitters in terms of the 
                              // highpass freq (allows to have 0 freq, when numBands == 1)

  result &= mbc.getNumberOfBands() == 2;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 4000.0);
  result &= mbc.getNumberOfBands() == 3;
  result &= testInOutEqual(mbc, N, tol);


  mbc.insertBand(0, 2000.0);
  result &= mbc.getNumberOfBands() == 4;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 1000.0);
  result &= mbc.getNumberOfBands() == 5;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 500.0);
  result &= mbc.getNumberOfBands() == 6;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 250.0);
  result &= mbc.getNumberOfBands() == 7;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 125.0);
  result &= mbc.getNumberOfBands() == 8;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 62.5);
  result &= mbc.getNumberOfBands() == 9;
  result &= testInOutEqual(mbc, N, tol);

  // todo: 
  // -use loop for adding bands 
  // -test removing bands (also using a loop)




  return result;
}

