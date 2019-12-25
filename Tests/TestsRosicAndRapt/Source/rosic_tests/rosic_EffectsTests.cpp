//#include "rosic_EffectsTests.h"
using namespace rotes;

#include "rosic/rosic.h"
using namespace rosic;

using namespace RAPT;

bool rotes::testFastGeneralizedHadamardTransform()
{
  bool result = true;

  // 4-point FGWHT:
  double x4[4] = {4, -8, 12, -4};
  double y4[4];
  double work[8];  // workspace

  typedef rosic::FeedbackDelayNetwork FDN;


  RAPT::rsArrayTools::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work);
  result &= y4[0] ==   4;
  result &= y4[1] ==  28;
  result &= y4[2] == -12;
  result &= y4[3] == - 4;

  RAPT::rsArrayTools::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work, 2, 3, 5, 7);
  result &= y4[0] ==  4;
  result &= y4[1] == 24;
  result &= y4[2] ==  4;
  result &= y4[3] == 44;


  // 8-point FWHT:
  double x8[8] = {1, 4, -2, 3, 0, 1, 4, -1};
  double y8[8];

  RAPT::rsArrayTools::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work);
  result &= y8[0] ==  10;
  result &= y8[1] == - 4;
  result &= y8[2] ==   2;
  result &= y8[3] == - 4;
  result &= y8[4] ==   2;
  result &= y8[5] == -12;
  result &= y8[6] ==   6;
  result &= y8[7] ==   8;

  RAPT::rsArrayTools::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, 7);
  result &= y8[0] ==   149;
  result &= y8[1] ==   357;
  result &= y8[2] ==   360;
  result &= y8[3] ==   862;
  result &= y8[4] ==   362;
  result &= y8[5] ==   866;
  result &= y8[6] ==   875;
  result &= y8[7] ==  2092;

  // forward/backward trafo - check if input is reconstructed:
  RAPT::rsArrayTools::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(       y8, 8, 3, work, 2, 3, 5, -7);
  FDN::fastInverseGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, -7);
  result &= fabs(RAPT::rsArrayTools::maxDeviation(x8, y8, 8)) < 1.e-15;

  return result;
}

bool rotes::testFeedbackDelayNetwork()
{
  bool result = true;

  //FeedbackDelayNetwork16 *fdn16 = new FeedbackDelayNetwork16;

  FeedbackDelayNetwork fdn;

  static const int N = 100000;
  //double hL[N], hR[N];
  double *hL = new double[N];
  double *hR = new double[N];


  RAPT::rsArrayTools::fillWithImpulse(hL, N);
  RAPT::rsArrayTools::fillWithImpulse(hR, N);
  for(int n = 0; n < N; n++)
    fdn.processFrame(&hL[n], &hR[n]);

  double t[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  //Plotter::plotData(N, t, hL, hR);
  //Plotter::plotData(N, t, hL);


  //rosic::writeToStereoWaveFile("d:\\TmpData\\FDNImpulseResponse.wav", hL, hR, N, 44100, 16);
  rosic::writeToStereoWaveFile("FDNImpulseResponse.wav", hL, hR, N, 44100, 16);

  delete[] hL;
  delete[] hR;
  return result;
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

