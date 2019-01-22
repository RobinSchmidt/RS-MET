bool testSpectrogramResynthesis(int blockSize, int hopSize, int signalLength, int fftSize, 
  int windowType = RAPT::rsWindowFunction::HANNING_WINDOW_ZN)
{
  bool r = true;

  int B = blockSize;
  int H = hopSize;
  int N = signalLength;
  int M = fftSize;

  // compute the complex spectrogram of a sequence of random numbers:
  std::vector<double> x = rsRandomVector(N, -1, +1);
  RAPT::rsSpectrogram<double> sp;  // spectrogram processor
  sp.setAnalysisWindowType(windowType);
  sp.setSynthesisWindowType(windowType);
  sp.setBlockSize(B);
  sp.setHopSize(H);
  //sp.setTrafoSize(M);
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(&x[0], N);
  //rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(&x[0], N, &w[0], B, H, 1);
  //int F = s.getNumRows();
  // todo: let the function take an FFT-size parameter instead of a zero-padding factor (maybe)
  // facilitates having an FFT size independent from the block-size

  // resynthesize signal from spectrogram:
  std::vector<double> y = sp.synthesize(s);

  // check, if resynthesized matches original signal:
  double tol = 1.e-13;
  int n;
  std::vector<double> err(N);
  for(n = 0; n < N; n++)         // create error signal
    err[n] = x[n] - y[n];
  for(n = 0; n < N; n++)
    r &= abs(err[n]) <= tol;
  //plotVector(y);
  if(r == false)
    plotVector(err);  // plot error signal, if something goes wrong

  return r;
}

// a subclass of rsSpectrogram, so we can access protected variables from test code
class rsSpectrogramUnitTest : public RAPT::rsSpectrogram<double>
{

public:

  bool testForwardTrafo(int N, double tol)
  {
    std::vector<std::complex<double>> y, x = rsComplexRandomVector(N, -1.0, +1.0); y = x;
    fft(  &y[0], N);  // actual
    rsFFT(&x[0], N);  // target
    return rsAlmostEqual(x, y, tol);
  }

  bool testInverseTrafo(int N, double tol)
  {
    std::vector<std::complex<double>> y, x = rsComplexRandomVector(N, -1.0, +1.0); y = x;
    ifft(  &y[0], N);  // actual
    rsIFFT(&x[0], N);  // target
    return rsAlmostEqual(x, y, tol);
  }

  bool testTrafos(int N, double tol)
  {
    bool r = testForwardTrafo(N, tol);
    r &= testInverseTrafo(N, tol);
    return r;
  }




  bool testTransforms()
  {
    bool r = true;

    r &= testTrafos(  8, 1.e-15);
    r &= testTrafos( 16, 1.e-15);
    r &= testTrafos(128, 1.e-14);

    return r;
  }

  bool runTests()
  {
    bool r = true;
    r &= testTransforms();
    return r;
  }

};


bool spectrogramUnitTest()
{
  bool r = true;      // test result

  // last parameter not yet used in function..

  rsSpectrogramUnitTest tester;
  r &= tester.runTests();



  r &= testSpectrogramResynthesis(8, 4, 20, 8);
  //r &= testSpectrogramResynthesis(7, 3, 20, 8);  // odd size window
  r &= testSpectrogramResynthesis(128, 64, 500, 128); // H = B/2
  r &= testSpectrogramResynthesis(128, 32, 500, 128); // H = B/4

  // todo: try block-sizes that are not a power of two, window functions that don't have natural
  // perfect reconstruction properties (i.e. demodulate the output)
  // test resynthesis with zero-padding

  // let the spectrogram class have a direct setFftSize function, let it use an 
  // rsFourierTransformer object, allow arbitrary FFT sizes



  return r;
}