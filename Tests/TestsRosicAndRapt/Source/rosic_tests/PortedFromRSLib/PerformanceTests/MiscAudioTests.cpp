#include "MiscAudioTests.h"

/*
std::string toString(int n)
{
  return std::to_string((_Longlong)n);
}
*/

void testSincInterpolator(std::string &reportString, double ratio, int sincLength)
{
  int    xN = 10000;                // number of samples in input signal
  int    yN = (int) ceil(xN/ratio); // number of output samples
  double *x = new double[xN];       // input signal
  double *y = new double[yN];       // output signal

  ::PerformanceCounterTSC counter;
  counter.init();
  rsResamplerDD::transposeSinc(x, xN, y, yN, ratio, sincLength, true);
  double cycles = (double) counter.getNumCyclesSinceInit();
  double cyclesPerSample = cycles / yN;

  std::string s;
  s += "Sinc";
  s += to_string(sincLength);
  s += " Interpolation";
  printPerformanceTestResult(s, cyclesPerSample);
}

void testSincInterpolator(std::string &reportString)
{
  int    xN = 10000;            // number of samples in input signal
  double r  = 2.1;              // resampling ratio
  //double r  = 1.000;              // resampling ratio
  int    yN = (int) ceil(xN/r); // number of output samples
  double *x = new double[xN];   // input signal
  double *y = new double[yN];   // output signal

  // create random noise as input signal:
  RAPT::rsArrayTools::fillWithRandomValues(x, xN, -1.0, +1.0, 0);

  // measure different interpolation effiencies:
  ::ProcessorCycleCounter counter;
  counter.init();
  rsResamplerDD::transposeLinear(x, xN, y, yN, r);
  double cycles = (double) counter.getNumCyclesSinceInit();
  double cyclesPerSample = cycles / yN;
  printPerformanceTestResult("Linear Interpolation", cyclesPerSample);

  testSincInterpolator(reportString, r, 10);    //   2700 -> 1000, 1000
  testSincInterpolator(reportString, r, 100);   //  24400 -> 2500, 4000
  testSincInterpolator(reportString, r, 1000);  // 233000 -> 28000, 33000
  testSincInterpolator(reportString, r, 16);
  testSincInterpolator(reportString, r, 32);
  testSincInterpolator(reportString, r, 64);
  testSincInterpolator(reportString, r, 128);
  testSincInterpolator(reportString, r, 256);
  testSincInterpolator(reportString, r, 512);
  testSincInterpolator(reportString, r, 1024);
  testSincInterpolator(reportString, r, 2048);

  delete[] x;
  delete[] y;
}

