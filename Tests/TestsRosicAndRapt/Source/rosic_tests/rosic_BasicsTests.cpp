//#include "rosic_BasicsTests.h"
using namespace rotes;

//#include "../Shared/Plotting/rosic_Plotter.h"
//#include "rosic/rosic.h"
using namespace rosic;

bool rotes::testBinomialCoefficients()
{
  bool ok = true;
  unsigned int nMax = 20;
  unsigned int c1, c2;
  for(unsigned int n = 0; n <= nMax; n++)
  {
    for(unsigned int k = 0; k <= n; k++)
    {
      c1 = RAPT::rsBinomialCoefficient(      n, k);
      c2 = RAPT::rsBinomialCoefficientUpTo20(n, k);
      ok &= c1 == c2;
    }
  }
  return ok;
}

bool rotes::testMathFunctions()
{
  bool ok = true;

  double y = rosic::besselI0(30.0);
  // https://www.wolframalpha.com/input/?i=besseli%5B0%2C30%5D
  // 7.81672297823977489717389816705295005444944253977947003347688 × 10^11
  // ...it's not very accurate but in the right ballpark

  //GNUPlotter plt; plt.plotFunctions(501, 0.0, 20.0, &besselI0);
  // looks like exponential growth - is that correct?


  using Cmp = std::complex<double>;
  using Vec = std::vector<Cmp>;

  Vec roots({Cmp( 2, 0), Cmp( 3, 0), Cmp(-1, 0)}); // p(x) = (x-2)*(x-3)*(x+1) = x^3-4*x^2+x+6
  Vec coeffs = RAPT::rsPolynomial<double>::rootsToCoeffs(roots);
  ok &= coeffs == Vec({6.0,1.0,-4.0,1.0});

  return ok;
}

void rotes::testWindowFunctions()
{
  // This is an experiment, not a unit test

  using AT = RAPT::rsArrayTools;

  static const int windowLength = 51;
  double windowIndices[windowLength];
  AT::fillWithIndex(windowIndices, windowLength);
  double window[windowLength];

  // test summing of time-shifted cosine-power windows:
  WindowDesigner::getWindow(window, windowLength, WindowDesigner::COSINE_SQUARED);
  //WindowDesigner::getCosinePowerWindow(window, windowLength,  8.0);

  plotData(windowLength, windowIndices, window);

  static const int hopSize      = windowLength/2;
  static const int numWindows   = 16;
  static const int signalLength = numWindows*windowLength - (numWindows-1)*(windowLength-hopSize);
  double signalIndices[signalLength];
  AT::fillWithIndex(signalIndices, signalLength);

  // scale the window such that the sum of the windows gives unity:
  //double normalizer = (2.0 * (double) hopSize / (double) windowLength);
  //scale(window, windowLength, normalizer);

  double shiftedWindows[numWindows][signalLength];
  for(int w=0; w<numWindows; w++)
  {
    AT::fillWithZeros(shiftedWindows[w], signalLength);
    int start = w*hopSize;
    for(int n=0; n<windowLength; n++)
      shiftedWindows[w][start+n] = window[n];
  }

  double windowSum[signalLength];
  AT::fillWithZeros(windowSum, signalLength);
  for(int w=0; w<numWindows; w++)
    AT::add(windowSum, shiftedWindows[w], windowSum, signalLength);
  plotData(signalLength, signalIndices, shiftedWindows[0], shiftedWindows[1], shiftedWindows[2], shiftedWindows[3], windowSum);

  // Test Kaiser window:
  WindowDesigner::getKaiserWindow(window, windowLength,  1.5);
  plotData(windowLength, windowIndices, window);  // looks parabolic - is that correct?

  int dummy = 0;
}


void testDerivativeComputation()
{
  double yPadded[14] = {0, 0, 0, 0, 1, 0.5, 0.1, -0.2, -0.7, -0.1, 0.2, 0.3, 0.7, 0.4};
  double *y = &yPadded[4];

  double yd;
  for(int i = 0; i < 10; i++)
    yd = RAPT::getDelayedSampleAsymmetricHermiteM(0.5, &y[i], 3);

  int dummy = 0;
}


void testCompensatedLinearInterpolator()
{
  double compensationAmount = 0.995; // 0: no compensation, 1:full compensation (but borderline (un)stable), 0.995 seems reasonable

  static const int N = 1024;
  double w[N];
  RAPT::rsArrayTools::fillWithRangeLinear(w, N, 0.0, PI);
  double mag[N];

  double d   = 0.5;
  double b0  = 1.0-d;
  double b1  = d;
  double a1;
  if( d <= 0.5 )
    a1 = d;
  else
    a1 = 1.0-d;
  a1  = 1.0 / (1.0/a1 - 1);
  a1 *= compensationAmount;

  double g  = 1.0+a1;

  // for a delay of d samples, interpolate via
  // w[n] = b0*x[n] + b1*x[n-1] - a1*w[n-1]; y[n] = g * w[n]  or
  // y[n] = g * (b0*x[n] + b1*x[n-1]) - a1*y[n-1]
  // -> test which version behaves best under modulation (witch d between 0.0 and 0.5 - these are the extremes for a1)
  // maybe try this stuff in DelayLine module for Liberty

  rosic::rsFilterAnalyzerD::getBiquadMagnitudeResponse(b0, b1, 0.0, a1, 0.0, w, mag, N, false);
  RAPT::rsArrayTools::scale(mag, mag, N, g);
  for(int n = 0; n < N; n++)
    mag[n] = RAPT::rsAmpToDb(mag[n]);
  RAPT::rsArrayTools::clip(mag, N, -60.0, 20.0);
  plotData(N, w, mag);
}

void rotes::testInterpolation()
{
  //testHermiteTwoPoint1();
  //testHermiteTwoPoint2();
  //testHermiteTwoPoint3();
  //testHermiteTwoPointM();

  int M[5]    = {0, 1, 2, 3, 4};                 // numbers of controlled derivatives for the 5 curves
  //double d[5] = {0.125, 0.25, 0.5, 0.75, 0.875}; // delay values for polyphase responses

  //double d[5] = {0.0, 0.05, 0.1, 0.15, 0.2};
  double d[5] = {0.25, 0.3, 0.345, 0.4, 0.45};  // problem with the 1-derivative  interpolator at d = 0.345
  //double d[5] = {0.25, 0.3, 0.35, 0.4, 0.45};     // problem with the 2-derivatives interpolator at d = 0.35
  //double d[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
  //double d[5] = {0.6, 0.7, 0.8, 0.9, 1.0};

  //plotOneSidedInterpolatorContinuousResponses(M, 1.0);
  plotOneSidedInterpolatorPolyphaseResponses(1, 1.0, d);

  //testCompensatedLinearInterpolator();
}

void rotes::testHermiteTwoPoint1()
{
  double a [4];
  double y0[2];
  double y1[2];

  y0[0] = +1; y0[1] = +1; // value and derivative at x = 0
  y1[0] =  0; y1[1] =  0; // value and derivative at x = 1

  getHermiteCoeffs1(y0, y1, a);
  getHermiteCoeffsM(y0, y1, a, 1);

  static const int N = 100;
  double x[N], y[2*N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.0);

  // evaluate polynomial and 3 derivatives:
  for(int n = 0; n < N; n++)
    RAPT::rsPolynomial<double>::evaluateWithDerivatives(x[n], a, 3, &y[2*n], 1);
  RAPT::rsArrayTools::deInterleave(y, N, 2);

  plotData(N, x, y);
  //plotData(N, x, y, &y[N]);
  //plotData(N, x, y, &y[N], &y[2*N]);
}


void rotes::testHermiteTwoPoint2()
{
  double a [6];
  double y0[3];
  double y1[3];

  y0[0] = +1; y0[1] = +2; y0[2] = -1;  // value and derivatives at x = 0
  y1[0] = -1; y1[1] = +1; y1[2] =  0;  // value and derivatives at x = 1

  getHermiteCoeffs2(y0, y1, a);
  getHermiteCoeffsM(y0, y1, a, 2);

  static const int N = 100;
  double x[N], y[3*N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.0);

  // evaluate polynomial and 3 derivatives:
  for(int n = 0; n < N; n++)
    RAPT::rsPolynomial<double>::evaluateWithDerivatives(x[n], a, 5, &y[3*n], 2);
  RAPT::rsArrayTools::deInterleave(y, N, 3);

  plotData(N, x, y);
  //plotData(N, x, y, &y[N]);
  //plotData(N, x, y, &y[N], &y[2*N]);
}

void rotes::testHermiteTwoPoint3()
{
  double a [8];
  double y0[4];
  double y1[4];

  y0[0] = +1; y0[1] = +1; y0[2] = +1; y0[3] = +1;  // value and derivatives at x = 0
  y1[0] =  0; y1[1] = -1; y1[2] = -2; y1[3] = -3;  // value and derivatives at x = 1

  getHermiteCoeffs3(y0, y1, a);
  getHermiteCoeffsM(y0, y1, a, 3);

  static const int N = 100;
  double x[N], y[4*N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.0);

  // evaluate polynomial and 3 derivatives:
  for(int n = 0; n < N; n++)
    RAPT::rsPolynomial<double>::evaluateWithDerivatives(x[n], a, 7, &y[4*n], 3);
  RAPT::rsArrayTools::deInterleave(y, N, 4);

  plotData(N, x, y);
  //plotData(N, x, y, &y[N]);
  //plotData(N, x, y, &y[N], &y[2*N], &y[3*N]);
}

void rotes::testHermiteTwoPointM()
{
  static const int maxM = 3; // maximum number of derivatives (determines how much memory is allocated)
  static const int M = 3;    // number of derivatives to be actually used
  double a [2*maxM+2];
  double y0[maxM+1];
  double y1[maxM+1];

  y0[0] = +1; y0[1] = +1; y0[2] = +1; y0[3] = +1;  // value and derivatives at x = 0
  //y1[0] =  0; y1[1] = -1; y1[2] = -2; y1[3] = -3;  // value and derivatives at x = 1
  y1[0] =  1; y1[1] =  0; y1[2] =  0; y1[3] = 0;  // value and derivatives at x = 1


  //hermiteTwoPoint3(y0, y1, a);
  getHermiteCoeffsM(y0, y1, a, M);

  static const int N = 400;
  double x[N], y[(M+1)*N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.3);

  // evaluate polynomial and 3 derivatives:
  for(int n = 0; n < N; n++)
    RAPT::rsPolynomial<double>::evaluateWithDerivatives(x[n], a, 7, &y[(M+1)*n], M);
  RAPT::rsArrayTools::deInterleave(y, N, M+1);

  plotData(N, x, y);
  //plotData(N, x, y, &y[N], &y[2*N], &y[3*N]);
}

double* createDeltaImpulse(int length, int impulseIndex)
{
  double *d = new double[length];
  RAPT::rsArrayTools::fillWithZeros(d, length);
  d[impulseIndex] = 1.0;
  return d;
}

rosic::Matrix interpolateOneSidedHermite5M(double *x, int xLength, int oversampling, const int M[5], double shape)
{
  rosic::Matrix y(5, xLength*oversampling);
  for(int i = 0; i < 5; i++)
    upsampleHermiteAsymmetricM(x, xLength, y.m[i], oversampling, M[i], shape);
  return y;
}

rosic::Matrix rotes::createHermiteInterpolatorImpulseResponses(int inLength, int oversampling, const int M[5], double shape)
{
  double *x = createDeltaImpulse(inLength, 1);
  rosic::Matrix y = interpolateOneSidedHermite5M(x, inLength, oversampling, M, shape);
  delete[] x;
  return y;
}




void plotMagnitudeResponsesOf(rosic::Matrix &x, int fftSize, double plotFloor, int plotZoom = 1, double freqScale = 1.0,
                              double ampScale = 1.0)
{
  rassert(x.numRows == 5); // this function is written for matrices with 5 rows - \todo: generalize this

  int length  = x.numColumns;
  int numBins = fftSize/2 + 1;

  rosic::Matrix magH(5, numBins, true);
  rosic::Matrix phsH(5, numBins, true);
  double *xPadded = new double[fftSize];
  RAPT::rsArrayTools::fillWithZeros(xPadded, fftSize);
  for(int m = 0; m < 5; m++)
  {
    RAPT::rsArrayTools::copy(x.m[m], xPadded, length);
    fftMagnitudesAndPhases(xPadded, fftSize, magH.m[m], phsH.m[m], fftSize);
    for(int n = 0; n < numBins; n++)
      magH.m[m][n] = RAPT::rsMax(RAPT::rsAmpToDb(magH.m[m][n] * ampScale * fftSize), plotFloor);
  }

  double *frequencies = new double[numBins];
  RAPT::rsArrayTools::fillWithIndex(frequencies, numBins);
  RAPT::rsArrayTools::scale(frequencies, frequencies, numBins, 2.0*freqScale/fftSize);  // wrong?

  int plotMaxIndex = numBins / plotZoom;

  plotData(plotMaxIndex, frequencies, magH.m[0], magH.m[1], magH.m[2], magH.m[3], magH.m[4]);
  //plotData(plotMaxIndex, frequencies, phsH.m[0], phsH.m[1], phsH.m[2], phsH.m[3], phsH.m[4]);

  delete[] frequencies;
  delete[] xPadded;
}



void rotes::plotOneSidedInterpolatorContinuousResponses(int M[5], double shape)
{
  // parameters:
  static const int inLength     = 16;
  static const int oversampling = 128;

  // create (continuous time) impulse responses:
  double *x = createDeltaImpulse(inLength, 1);
  rosic::Matrix h = interpolateOneSidedHermite5M(x, inLength, oversampling, M, shape);
  double *t = new double[h.numColumns];
  RAPT::rsArrayTools::fillWithIndex(t, h.numColumns);
  RAPT::rsArrayTools::scale(t, t, h.numColumns, 1.0/oversampling);
  RAPT::rsArrayTools::add(t, -1.0, t, h.numColumns);

  plotData(h.numColumns, t, h.m[0], h.m[1], h.m[2], h.m[3], h.m[4]);
  //plotMagnitudeResponsesOf(h, 8192, -120.0, 10, oversampling, 1.0/oversampling);

  delete[] x;
  delete[] t;
}

void rotes::plotOneSidedInterpolatorPolyphaseResponses(int M, double shape, double d[5])
{
  static const int inLength2 = 32;
  double *x = createDeltaImpulse(inLength2, inLength2/2);

  double *t = new double[inLength2];
  RAPT::rsArrayTools::fillWithIndex(t, inLength2);
  rosic::Matrix hd(5, inLength2, true);
  for(int m = 0; m < 5; m++)
  {
    for(int n = M+1; n < inLength2; n++)
      hd.m[m][n] = RAPT::getDelayedSampleAsymmetricHermiteM(d[m], &x[n], M, shape);
  }
  plotData(hd.numColumns, t, hd.m[0], hd.m[1], hd.m[2], hd.m[3], hd.m[4]);
  plotMagnitudeResponsesOf(hd, 8192, -60.0, 1);

  delete[] x;
  delete[] t;
}


void rotes::testAsymmetricPolynomialInterpolatorsOld()
{
  double x0, x1, y0, y1, yd0, yd1, ydd0, ydd1;
  double a5, a4, a3, a2, a1, a0;

  static const int oversampling = 128;
  static const int length       = 4*oversampling+1;
  double x[length], yl[length], yc[length], yq[length];
  RAPT::rsArrayTools::fillWithRangeLinear(x,  length, -1.0, 3.0);
  RAPT::rsArrayTools::fillWithZeros(      yl, length           );  // linear
  RAPT::rsArrayTools::fillWithZeros(      yc, length           );  // cubic
  RAPT::rsArrayTools::fillWithZeros(      yq, length           );  // quintic
  int n;

  // linear:
  x0  = -1.0;  y0  =  0.0;  yd0 = +1.0;
  x1  =  0.0;  y1  = +1.0;  yd1 = +1.0;
  RAPT::fitCubicWithDerivative(x0, x1, y0, y1, yd0, yd1, &a3, &a2, &a1, &a0);
  for(n = 0; n < length/4; n++)
    yl[n] = a3*x[n]*x[n]*x[n] + a2*x[n]*x[n] + a1*x[n] + a0;

  x0  =  0.0;  y0  = +1.0;  yd0 = -1.0;
  x1  =  1.0;  y1  =  0.0;  yd1 = -1.0;
  RAPT::fitCubicWithDerivative(x0, x1, y0, y1, yd0, yd1, &a3, &a2, &a1, &a0);
  for(n = length/4; n < 2*length/4; n++)
    yl[n] = a3*x[n]*x[n]*x[n] + a2*x[n]*x[n] + a1*x[n] + a0;

  // cubic:
  y0  =  0.0;  yd0 =  0.0;
  y1  = +1.0;  yd1 = +1.0;
  double xn;
  RAPT::fitCubicWithDerivativeFixedX(y0, y1, yd0, yd1, &a3, &a2, &a1, &a0);
  for(n = 0; n < length/4; n++)
  {
    xn    = x[n] + 1.0;
    yc[n] = a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }

  y0 = +1.0;  yd0 = +1.0;
  y1 =  0.0;  yd1 = -1.0;
  RAPT::fitCubicWithDerivativeFixedX(y0, y1, yd0, yd1, &a3, &a2, &a1, &a0);
  for(n = length/4; n < 2*length/4; n++)
  {
    xn    = x[n] + 0.0;
    yc[n] = a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }

  y0 = 0.0; yd0 = -1.0;
  y1 = 0.0; yd1 =  0.0;
  RAPT::fitCubicWithDerivativeFixedX(y0, y1, yd0, yd1, &a3, &a2, &a1, &a0);
  for(n = 2*length/4; n < 3*length/4; n++)
  {
    xn    = x[n] - 1.0;
    yc[n] = a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }



  // quintic:
  y0  =  0.0;  yd0 =  0.0;  ydd0 =  0.0;
  y1  = +1.0;  yd1 = +1.0;  ydd1 = +1.0;
  RAPT::fitQuinticWithDerivativesFixedX(y0, y1, yd0, yd1, ydd0, ydd1, &a5, &a4, &a3, &a2, &a1, &a0);
  for(n = 0; n < length/4; n++)
  {
    xn    = x[n] + 1.0;
    yq[n] = a5*xn*xn*xn*xn*xn + a4*xn*xn*xn*xn + a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }


  y0 = +1.0;  yd0 = +1.0; ydd0 = +1.0;
  y1 =  0.0;  yd1 = -1.0; ydd1 = -2.0;
  RAPT::fitQuinticWithDerivativesFixedX(y0, y1, yd0, yd1, ydd0, ydd1, &a5, &a4, &a3, &a2, &a1, &a0);
  for(n = length/4; n < 2*length/4; n++)
  {
    xn    = x[n] + 0.0;
    yq[n] = a5*xn*xn*xn*xn*xn + a4*xn*xn*xn*xn + a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }

  y0 =  0.0;  yd0 = -1.0; ydd0 = -2.0;
  y1 =  0.0;  yd1 =  0.0; ydd1 =  1.0;
  RAPT::fitQuinticWithDerivativesFixedX(y0, y1, yd0, yd1, ydd0, ydd1, &a5, &a4, &a3, &a2, &a1, &a0);
  for(n = 2*length/4; n < 3*length/4; n++)
  {
    xn    = x[n] - 1.0;
    yq[n] = a5*xn*xn*xn*xn*xn + a4*xn*xn*xn*xn + a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }

  y0 =  0.0;  yd0 =  0.0; ydd0 = +1.0;
  y1 =  0.0;  yd1 =  0.0; ydd1 =  0.0;
  RAPT::fitQuinticWithDerivativesFixedX(y0, y1, yd0, yd1, ydd0, ydd1, &a5, &a4, &a3, &a2, &a1, &a0);
  for(n = 3*length/4; n < 4*length/4; n++)
  {
    xn    = x[n] - 2.0;
    yq[n] = a5*xn*xn*xn*xn*xn + a4*xn*xn*xn*xn + a3*xn*xn*xn + a2*xn*xn + a1*xn + a0;
  }






  //Plotter::plotData(length, x, yl, yc, yq);



  //Plotter::pl

  // write function to interpolate arbitrary input signals

  // do fourier transform (pad to length 4096 or sth.) -> investigate spectrum (magnitude and phase)

  static const int fftSize = 8192;
  static const int numBins = fftSize/2 + 1;

  double frequencies[numBins];
  double magL[numBins], phsL[numBins], magC[numBins], phsC[numBins], magQ[numBins], phsQ[numBins];

  RAPT::rsArrayTools::fillWithZeros(magL, numBins);
  RAPT::rsArrayTools::fillWithZeros(phsL, numBins);
  RAPT::rsArrayTools::fillWithZeros(magC, numBins);
  RAPT::rsArrayTools::fillWithZeros(phsC, numBins);
  RAPT::rsArrayTools::fillWithZeros(magQ, numBins);
  RAPT::rsArrayTools::fillWithZeros(phsQ, numBins);
  RAPT::rsArrayTools::fillWithIndex(frequencies, numBins);
  RAPT::rsArrayTools::scale(frequencies, frequencies, numBins, (double) 2*oversampling/fftSize);


  double tmp[fftSize];
  double plotFloor = -120.0;
  RAPT::rsArrayTools::fillWithZeros(tmp, fftSize);

  RAPT::rsArrayTools::copy(yl, tmp, length);
  fftMagnitudesAndPhases(tmp, fftSize, magL, phsL, fftSize);
  for(n=0; n<numBins; n++)
    magL[n] = RAPT::rsMax(RAPT::rsAmpToDb(4 * magL[n] * fftSize/length), plotFloor);

  RAPT::rsArrayTools::copy(yc, tmp, length);
  fftMagnitudesAndPhases(tmp, fftSize, magC, phsC, fftSize);
  for(n=0; n<numBins; n++)
    magC[n] = RAPT::rsMax(RAPT::rsAmpToDb(4 * magC[n] * fftSize/length), plotFloor);

  RAPT::rsArrayTools::copy(yq, tmp, length);
  fftMagnitudesAndPhases(tmp, fftSize, magQ, phsQ, fftSize);
  for(n=0; n<numBins; n++)
    magQ[n] = RAPT::rsMax(RAPT::rsAmpToDb(4 * magQ[n] * fftSize/length), plotFloor);


  plotData(600, frequencies, magL, magC, magQ);


    /*
  fftMagnitudesAndPhases(impulseResponse, length, magnitudes, NULL, length);

  for(n=0; n<length; n++)
    decibels[n] = amp2dB(length*magnitudes[n]);

  Plotter::plotData(length/4, frequencies, decibels);
  */



  int dummy = 0;
}
