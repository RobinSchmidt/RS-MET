
typedef std::complex<double> rsComplexDbl;

template<class T>
void interleaveWithZeros(T *x, T *y, int xLength, int factor)
{
  RAPT::rsArrayTools::fillWithZeros(y, factor*xLength);
  for(int i = 0; i < xLength; i++)
    y[factor*i] = x[i];
}
// move into the RSCore library


bool testSmbFFT()
{
  bool testResult = true;

  static const int N = 16;
  float x[N];    // signal
  float X[2*N];  // computed spectrum
  float T[2*N];  // target spectrum

  // create test signal:
  int n;
  for(n = 0; n < N; n++)
    x[n] = (float) (1.0/(n+1) - 0.25);

  // create target spectrum:
  T[0]  = (float) -0.619271006771007;    // DC component
  T[1]  = (float)  0.0;                  // always 0 (for real signals?), formally sine-part of DC
  T[2]  = (float)  1.345577785360745;
  T[3]  = (float) -0.767331087655140;
  T[4]  = (float)  0.986973606551556;
  T[5]  = (float) -0.573743664928213;
  T[6]  = (float)  0.834508128766654;
  T[7]  = (float) -0.429285375589382;
  T[8]  = (float)  0.754267954267954;
  T[9]  = (float) -0.317261904761905;
  T[10] = (float)  0.708179602278973;
  T[11] = (float) -0.224849965506606;
  T[12] = (float)  0.681402461824512;
  T[13] = (float) -0.144306435490984;
  T[14] = (float)  0.667290039149184;
  T[15] = (float) -0.070587985264671;
  T[16] = (float)  0.662871850371850;    // Nyquist frequency component
  T[17] = (float)  0.0;                  // sine-part of Nyquist frequency component
  T[18] = (float)  0.667290039149184;
  T[19] = (float)  0.070587985264671;
  T[20] = (float)  0.681402461824512;
  T[21] = (float)  0.144306435490984;
  T[22] = (float)  0.708179602278973;
  T[23] = (float)  0.224849965506606;
  T[24] = (float)  0.754267954267954;
  T[25] = (float)  0.317261904761905;
  T[26] = (float)  0.834508128766654;
  T[27] = (float)  0.429285375589382;
  T[28] = (float)  0.986973606551556;
  T[29] = (float)  0.573743664928213;
  T[30] = (float)  1.345577785360745;
  T[31] = (float)  0.767331087655140;

  // run FFT, calculate numerical error and check if below some margin:
  interleaveWithZeros(x, X, N, 2);
  smbFft(X, N, -1);
  double error = rsArrayTools::maxDeviation(X, T, 2*N);
  testResult &= (error < 1.e-6);  // a rather large margin is required

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testRsFFT(std::string &reportString)
{
  std::string testName = "RsFFT";
  bool testResult = true;

  static const int N = 16;
  rsComplexDbl x[N];  // signal
  rsComplexDbl X[N];  // computed spectrum
  rsComplexDbl T[N];  // target spectrum
  rsComplexDbl y[N];  // reconstructed signal

  // create test signal:
  int n;
  for(n = 0; n < N; n++)
  {
    x[n].real(1.0/(n+1)-0.25);
    x[n].imag(1.0/(n+1)-0.10);
  }

  // create target spectrum:
  T[0]  = rsComplexDbl(-0.619271006771007, 1.780728993228993);
  T[1]  = rsComplexDbl( 2.112908873015884, 0.578246697705605);
  T[2]  = rsComplexDbl( 1.560717271479769, 0.413229941623343);
  T[3]  = rsComplexDbl( 1.263793504356036, 0.405222753177271);
  T[4]  = rsComplexDbl( 1.071529859029859, 0.437006049506050);
  T[5]  = rsComplexDbl( 0.933029567785579, 0.483329636772367);
  T[6]  = rsComplexDbl( 0.825708897315496, 0.537096026333528);
  T[7]  = rsComplexDbl( 0.737878024413855, 0.596702053884513);
  T[8]  = rsComplexDbl( 0.662871850371850, 0.662871850371850);
  T[9]  = rsComplexDbl( 0.596702053884513, 0.737878024413855);
  T[10] = rsComplexDbl( 0.537096026333529, 0.825708897315496);
  T[11] = rsComplexDbl( 0.483329636772367, 0.933029567785579);
  T[12] = rsComplexDbl( 0.437006049506050, 1.071529859029859);
  T[13] = rsComplexDbl( 0.405222753177271, 1.263793504356036);
  T[14] = rsComplexDbl( 0.413229941623343, 1.560717271479769);
  T[15] = rsComplexDbl( 0.578246697705605, 2.112908873015884);

  // compute spectrum via FFT:
  rsArrayTools::copy(x, X, N);
  rsFFT(X, N);

  // check maximum deviation between target- and computed spectrum:
  double error = 0.0;
  for(n = 0; n < N; n++)
    error = rsMax(error, abs(T[n]-X[n]));
  testResult &= (error < 1.e-15);

  // compute and check spectrum via DFT:
  rsArrayTools::copy(x, X, N);
  rsDFT(X, N);
  error = 0.0;
  for(n = 0; n < N; n++)
    error = rsMax(error, abs(T[n]-X[n]));
  testResult &= (error < 1.e-14); // needs more tolerance than FFT

  // check a forward/inverse turnaround cycle:
  rsArrayTools::copy(x, y, N);
  //rsFFT( y, N);
  //rsIFFT(y, N);
  rsLinearTransforms::fourierRadix2DIF(   y, N);
  rsLinearTransforms::fourierInvRadix2DIF(y, N);
  error = 0.0;
  for(n = 0; n < N; n++)
    error = rsMax(error, abs(x[n]-y[n]));
  testResult &= (error < 1.e-15);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFourierTransformerRadix2(std::string &reportString)
{
  std::string testName = "FourierTransformerRadix2";
  bool testResult = true;

  static const int N = 128; // FFT size
  rsComplexDbl x[N];        // complex input signal
  rsComplexDbl T[N];        // target spectrum
  rsComplexDbl X[N];        // computed spectrum

  typedef rsFourierTransformerRadix2<double> FT;

  // create test signal:
  int n;
  for(n = 0; n < N; n++)
  {
    x[n].real(1.0/(n+1) - 0.05);
    x[n].imag(1.0/(n+1) - 0.07);
  }

  // create target spectrum (we assume here that rsFFT computes a correct result):
  rsArrayTools::convert(x, T, N);
  rsFFT(T, N);

  // create the rsFourierTransformerRadix2 object, set it up and let it compute the spectrum:
  FT ft;
  ft.setBlockSize(N);
  ft.setDirection(FT::FORWARD);
  ft.setNormalizationMode(FT::NORMALIZE_ON_INVERSE_TRAFO);
  ft.transformComplexBuffer(x, X);

  // check maximum deviation between target- and computed spectrum:
  double error = 0.0;
  for(n = 0; n < N; n++)
    error = rsMax(error, abs(T[n]-X[n]));
  testResult &= (error < 1.e-13); // in MSVC, we can use 1.e-14 -> it's more precise there

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFourierTrafoRadix2(int N) // maybe rename to testComplexFourierTrafoRadix2
{
  bool r = true;
  double tol = 1.e-13; // 1.e-13 works up to N=32
  typedef RAPT::rsArrayTools AR;
  typedef RAPT::rsFourierTransformerRadix2<double> FT;
  std::vector<complex<double>> x(N), X(N);  // signal and spectrum
  std::vector<complex<double>> t(N), T(N);  // target signal and target spectrum


  // fill x and t with the same N random values:
  rsFillWithComplexRandomValues(x, -1.0, 1.0);
  AR::copy(&x[0], &t[0], N);

  // compute target spectrum by using the naive DFT implementation:
  AR::copy(&t[0], &T[0], N); // because DFT works in place
  RAPT::rsDFT(&T[0], N);

  // compute spectrum via RAPT::rsRadix2FFT and compare:
  AR::copy(&x[0], &X[0], N);
  //RAPT::rsRadix2FFT(&X[0], N);
  rsLinearTransforms::fourierRadix2DIF(&X[0], N);
  r &= rsAlmostEqual(T, X, tol);

  // use the rsFourierTransformerRadix2 object:
  FT ft;
  ft.setBlockSize(N);
  ft.setNormalizationMode(FT::NORMALIZE_ON_INVERSE_TRAFO); // is actually the default setting anyway
  ft.setDirection(FT::FORWARD);
  ft.transformComplexBuffer(&x[0], &X[0]);
  r &= rsAlmostEqual(T, X, tol);

  // now test inverse trafos:
  AR::copy(&T[0], &x[0], N); // target spectrum into signal buffer x for in-place iFFT
  //RAPT::rsIFFT(&x[0], N);
  rsLinearTransforms::fourierInvRadix2DIF(&x[0], N);
  r &= rsAlmostEqual(t, x, tol);

  ft.setDirection(FT::INVERSE);
  ft.transformComplexBuffer(&T[0], &x[0]);
  r &= rsAlmostEqual(t, x, tol);

  return r;
}

bool testFourierTrafoArbitrary(int N)
{
  bool r = true;
  double tol = 1.e-12; // 1.e-13 works up to N=32 in msc, gcc needs 1.e-12
  typedef RAPT::rsArrayTools AR;
  typedef RAPT::rsFourierTransformerRadix2<double> FTR2;
  typedef RAPT::rsFourierTransformerBluestein<double> FTB;
  std::vector<complex<double>> x(N), X(N);  // signal and spectrum
  std::vector<complex<double>> t(N), T(N);  // target signal and target spectrum

  // fill x and t with the same N random values:
  rsFillWithComplexRandomValues(x, -1.0, 1.0);
  AR::copy(&x[0], &t[0], N);

  // compute target spectrum by using the naive DFT implementation:
  AR::copy(&t[0], &T[0], N); // because DFT works in place
  RAPT::rsDFT(&T[0], N);

  //// compute spectrum via RAPT::rsBluesteinFFT and compare (not yet implemented):
  //AR::copy(&x[0], &X[0], N);
  //RAPT::rsBluesteinFFT(&X[0], N);
  //r &= rsAlmostEqual(T, X, tol);

  // use the rsFourierTransformerBluestein object:
  FTB ft;
  ft.setBlockSize(N);
  ft.setNormalizationMode(FTR2::NORMALIZE_ON_INVERSE_TRAFO); // is actually the default setting anyway
  ft.setDirection(FTR2::FORWARD);
  ft.transformComplexBuffer(&x[0], &X[0]);
  r &= rsAlmostEqual(T, X, tol);

  // now test inverse trafos:
  //AR::copy(&T[0], &x[0], N); // target spectrum into signal buffer x for in-place iFFT
  //RAPT::rsIFFT(&x[0], N);
  //r &= rsAlmostEqual(t, x, tol);

  ft.setDirection(FTR2::INVERSE);
  ft.transformComplexBuffer(&T[0], &x[0]);
  r &= rsAlmostEqual(t, x, tol);

  ft.transformComplexBufferInPlace(&T[0]);
  r &= rsAlmostEqual(t, T, tol);

  return r;
}

bool testVariousFourierTransforms(std::string &reportString)
{
  std::string testName = "various Fourier transforms";
  bool testResult = true;

  // We create random complex signal buffers of various lengths and transform them to the
  // frequency domain using various implementations and compare the results - they should
  // all be the same. Then we also compare the corresponding inverse FFT implementations.

  // ToDo: maybe test also for single-precision, i.e. templatize the functions called in the loops
  // below:

  int minTrafoSize = 2;  // todo: allow trafo size = 1 in rsFourierTransformerRadix2
  int maxTrafoSize = 32;


  // test radix-2 transforms:
  for(int i = minTrafoSize; i <= maxTrafoSize; i *= 2)
    testResult &= testFourierTrafoRadix2(i);

  // test arbitrary size transforms:
  for(int i = minTrafoSize; i <= maxTrafoSize; i++)
    testResult &= testFourierTrafoArbitrary(i);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testNumberTheoreticTransform()
{
  bool ok = true;

  using ModInt = rsModularIntegerNTT_64;
  using VecU64 = std::vector<rsUint64>;
  using VecI32 = std::vector<rsInt32>;
  using AT     = rsArrayTools;

  // Test the magic numbers:
  static const int numRoots = ModInt::numRoots;


  int maxN  = rsPowInt(2, numRoots);
  ModInt a, b, c;
  ModInt one  = ModInt(1);
  for(int i = 0; i < numRoots; i++)
  {
    int n = rsPowInt(2, i+1);
    a = ModInt(n);
    b = ModInt(ModInt::lengthsInv[i]);
    c = a * b;
    ok &= c == one;

    a = ModInt(ModInt::roots[i]);
    b = ModInt(ModInt::rootsInv[i]);
    c = a * b;
    ok &= c == one;

    int k = i+1;
    c = a;
    for(int j = 1; j <= k; j++) {
      ok &= c != one;
      c *= c;  }
    ok &= c == one;

    //c = a;
    //for(int j = 1; j < n; j++) {
    //  ok &= c != one;              // Root should be primitive, we should not get 1 for any power
    //  c *= a;  }                   // ...less than n (i.e. n/2, n/3, etc.)
    //ok &= c == one;
    // This loop takes long for high n. That's not sruprising since n grows exponentially, but i 
    // think, the repeated squaring loop above (over k, which runs only up to log2(n)) is enough to
    // ensure that it all works. 
  }

  // Test NTT convolution:
  auto testConvNTT = [](int Nx, int Nh, int mask = RAPT::allBits, int seed = 0)
  {
    int Ny = Nx+Nh-1;  // length of output signal 
    VecI32 x(Nx), h(Nh), y(Ny);
    rsNoiseGenerator<unsigned int> ng;
    ng.setSeed(seed);
    for(int i = 0; i < Nx; i++) 
      x[i] = ng.getSampleRaw() & mask;
    for(int i = 0; i < Nh; i++) 
      h[i] = ng.getSampleRaw() & mask;
    AT::convolve(&x[0], Nx, &h[0], Nh, &y[0]);
    VecI32 y2 = rsConvolveNTT(x, h);
    //rsPlotVectors(x, h, y, y2); 
      // The signals are actually periodic with period mask+1! Why? (assuming that mask mask away
      // a couple of leading bits)
    return y2 == y;  
  };
  ok &= testConvNTT(900, 90,    63);
  ok &= testConvNTT(100, 20,    15);
  ok &= testConvNTT(100, 20,    31);
  ok &= testConvNTT(100, 20,    63);
  ok &= testConvNTT(100, 27,    63);  // Ny = 126
  ok &= testConvNTT(100, 28,    63);  // Ny = 127
  ok &= testConvNTT(100, 29,    63);  // Ny = 128
  ok &= testConvNTT(100, 30,    63);  // Ny = 129
  ok &= testConvNTT(100, 20,  2047);  // mask = 2^11 - 1
  ok &= testConvNTT(100, 20,  8191);  // mask = 2^13 - 1
  ok &= testConvNTT(100, 20, 16383);  // mask = 2^14 - 1
  //ok &= testConvNTT(100, 20, 32767);  // mask = 2^15 - 1  -> fails
  //ok &= testConvNTT(100, 20, 65535);  // mask = 2^16 - 1  -> fails

  return ok;

  // ToDo: 
  // -Figure out at which point we get overflow. This may depend on the mask as well as the 
  //  lengths. When overflow occurs, probably both, naive and NTT convolution will produce
  //  total garbage - but probably different garbage. When using this in production, we need to
  //  take care of avoiding overflow, maybe by clipping the inputs and restricting the lengths.
  //  ...ok - done: for Nx = 100, Nh = 20, the limit is mask = 16383
  // -We get overflow even for moderate lengths: (950,95,63) overflows, (900,90,63) still works,
  //  so it seems, the NTT convolution is not practical for audio, where we have to deal with
  //  long impulse responses. We could still chop the data into very short segments, though...maybe
  //  NTT for short segments, FFT for the longer ones - but then, the additional complexity of
  //  having different methods is hard to justify.
  // -Is it possible to avoid overflow or make it less likely by using signed (modular) integers?
  //  Is that even possible? Or do we already implícitly represent negative numbers by their 
  //  equivalence to modulus - x, for a number -x where x is positive?
  // -We should be able to push the limit of failure for direct convolution by converting to int64.
  //  But that doesn't really help with the NTT convolution.

  // Also: measure the performance of complex vs modular arithmetic. write optimized NTT routines 
  // that use reducing to the modulus not after every add/sub but only after a group, if possible. 
  // Maybe that works better with high-radix FFT algos. We should also prefer addition over 
  // subtraction because subtraction is a bit more expensive.
};



bool testCorrelation(std::string &reportString)
{
  std::string testName = "Correlation";
  bool testResult = true;

  static const int N = 3;
  float x[N] = {1, -3, 2};
  float y[N] = {1,  5, 3};
  float rxx[N];
  float rxy[N];
  float ryx[N];
  float ryy[N];

  // test direct computation of auto- and cross-correlations:
  rsCrossCorrelationDirect(x, x, 3, rxx);
  testResult &= (rxx[0] == 14);
  testResult &= (rxx[1] == -9);
  testResult &= (rxx[2] ==  2);
  rsCrossCorrelationDirect(x, y, 3, rxy);
  testResult &= (rxy[0] == -8);
  testResult &= (rxy[1] ==  7);
  testResult &= (rxy[2] ==  2);
  rsCrossCorrelationDirect(y, x, 3, ryx);
  testResult &= (ryx[0] == -8);
  testResult &= (ryx[1] == -4);
  testResult &= (ryx[2] ==  3);
  rsCrossCorrelationDirect(y, y, 3, ryy);
  testResult &= (ryy[0] == 35);
  testResult &= (ryy[1] == 20);
  testResult &= (ryy[2] ==  3);

  // test FFT based computation of auto- and cross-correlations:
  float tol = (float) 1.e-5;  // quite a large tolerance necesarry
  rsCrossCorrelationFFT(x, x, 3, rxx);
  testResult &= fabs(rxx[0] - 14) < tol;
  testResult &= fabs(rxx[1] - -9) < tol;
  testResult &= fabs(rxx[2] -  2) < tol;
  rsCrossCorrelationFFT(x, y, 3, rxy);
  testResult &= fabs(rxy[0] - -8) < tol;
  testResult &= fabs(rxy[1] -  7) < tol;
  testResult &= fabs(rxy[2] -  2) < tol;
  rsCrossCorrelationFFT(y, x, 3, ryx);
  testResult &= fabs(ryx[0] - -8) < tol;
  testResult &= fabs(ryx[1] - -4) < tol;
  testResult &= fabs(ryx[2] -  3) < tol;
  rsCrossCorrelationFFT(y, y, 3, ryy);
  testResult &= fabs(ryy[0] - 35) < tol;
  testResult &= fabs(ryy[1] - 20) < tol;
  testResult &= fabs(ryy[2] -  3) < tol;
  rsAutoCorrelationFFT(x, 3, rxx);
  testResult &= fabs(rxx[0] - 14) < tol;
  testResult &= fabs(rxx[1] - -9) < tol;
  testResult &= fabs(rxx[2] -  2) < tol;
  rsAutoCorrelationFFT(y, 3, ryy);
  testResult &= fabs(ryy[0] - 35) < tol;
  testResult &= fabs(ryy[1] - 20) < tol;
  testResult &= fabs(ryy[2] -  3) < tol;

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFitQuadratic(std::string &reportString)
{
  std::string testName = "FitQuadratic";
  bool testResult = true;

  double x[3] = {2, 3, 5};
  double y[3] = {1, 4, 2};
  double a[3];

  rsPolynomial<double>::fitQuadratic(a, x, y);
  testResult &= (a[0] == -13.000000000000000 );
  testResult &= (a[1] ==   9.6666666666666661);
  testResult &= (a[2] ==  -1.3333333333333333);

  x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;
  rsPolynomial<double>::fitQuadratic(a, x, y);
  testResult &= (a[0] ==  1.0);
  testResult &= (a[1] ==  5.5);
  testResult &= (a[2] == -2.5);

  rsPolynomial<double>::fitQuadratic_0_1_2(a, y);
  testResult &= (a[0] ==  1.0);
  testResult &= (a[1] ==  5.5);
  testResult &= (a[2] == -2.5);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFitQuarticWithDerivatives(std::string &reportString)
{
  std::string testName = "FitQuarticWithDerivatives";
  bool testResult = true;

  double y[3] = {1, 2, 0};
  double s0   =  1;
  double s2   = -2;
  double a[5];

  rsPolynomial<double>::fitQuarticWithDerivatives(a, y, s0, s2);

  double tmp[2];
  rsPolynomial<double>::evaluateWithDerivatives(0.0, a, 4, tmp, 1);
  testResult &= (tmp[0] == 1.0);
  testResult &= (tmp[1] == 1.0);
  rsPolynomial<double>::evaluateWithDerivatives(1.0, a, 4, tmp, 1);
  testResult &= (tmp[0] == 2.0);
  rsPolynomial<double>::evaluateWithDerivatives(2.0, a, 4, tmp, 1);
  testResult &= (tmp[0] ==  0.0);
  testResult &= (tmp[1] == -2.0);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

/*
bool testLinearSystem3x3(std::string &reportString)
{
  std::string testName = "LinearSystem3x3";
  bool testResult = true;

  double x[3];
  double y[3]    = {-4, 15, 11};
  double A[3][3] = {{ 1, 2, -3},
                    {-5, 7,  2},
                    {-2, 2,  3}};

  rsSolveLinearSystem3x3(A, x, y);

  testResult &= (x[0] == 1.0);
  testResult &= (x[1] == 2.0);
  testResult &= (x[2] == 3.0);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
*/

bool testTransforms()
{
  std::string testName = "Transforms";
  std::string dummy;
  bool ok = true;

  // test FFT routines:
  ok &= testSmbFFT();
  ok &= testRsFFT(dummy);
  ok &= testFourierTransformerRadix2(dummy);
  ok &= testVariousFourierTransforms(dummy);
  ok &= testNumberTheoreticTransform();  // still fails, i think due to overflow -> we need a smaller modulus

  return ok;
}
