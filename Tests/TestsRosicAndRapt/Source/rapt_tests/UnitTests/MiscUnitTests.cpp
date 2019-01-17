bool testSpectrogramResynthesis(int blockSize, int hopSize, int signalLength, int fftSize, 
  int windowType = RAPT::rsWindowFunction::HANNING_WINDOW_ZN)
{
  bool r = true;

  int B = blockSize;
  int H = hopSize;
  int N = signalLength;
  int M = fftSize;

  // generate the window function (todo: let the spectrogram class do that by itself):
  std::vector<double> w(B);
  RAPT::rsWindowFunction::createWindow(&w[0], B, windowType, true);

  // generate a sequence of random numbers:
  std::vector<double> x(N);
  int n;
  RAPT::rsNoiseGenerator<double> prng;
  for(n = 0; n < N; n++)
    x[n] = prng.getSample();

  typedef RAPT::rsSpectrogram<double> SP; // spectrogram processor
  SP sp;

  // compute the complex spectrogram:
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(&x[0], N, &w[0], B, H, 1);
  //int F = s.getNumRows();
  // todo: let the function take an FFT-size parameter instead of a zero-padding factor (maybe)
  // facilitates having an FFT size independent from the block-size

  // resynthesize signal from spectrogram:
  std::vector<double> y  = sp.synthesize(s, &w[0], B, H, &w[0]);

  // check, if resynthesized matches original signal:
  double tol = 1.e-13;
  std::vector<double> err(N);
  for(n = 0; n < N; n++)         // create error signal
    err[n] = x[n] - y[n];
  for(n = 0; n < N; n++)
    r &= abs(err[n]) <= tol;
  //plotVector(err);  // may be uncommented to plot error signal, if something is going wrong

  return r;
}


bool spectrogramUnitTest()
{
  bool r = true;      // test result

  // last parameter not yet used in function..


  r &= testSpectrogramResynthesis(8, 4, 20, 8);
  //r &= testSpectrogramResynthesis(7, 3, 20, 8);  // odd size window
  r &= testSpectrogramResynthesis(128, 64, 500, 128); // H = B/2
  r &= testSpectrogramResynthesis(128, 32, 500, 128); // H = B/4

  // todo: try block-sizes that are not a power of two, window functions that don't have natural
  // perfect reconstruction properties (i.e. demodulate the output)

  // let the spectrogram class have a direct setFftSize function, let it use an 
  // rsFourierTransformer object, allow arbitrary FFT sizes



  return r;
}