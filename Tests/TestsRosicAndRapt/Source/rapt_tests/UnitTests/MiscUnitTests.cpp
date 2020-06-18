bool blepUnitTest()
{
  bool r = true;

  rsTableLinBlep<double, double> linBlep;  // rename to lb, mb
  rsTableMinBlep<double, double> minBlep;

  // we test various settings for the kernel-length and table precision - this test is only to rule
  // out access violations (if it doesn't trigger assertions, everything should be fine) but 
  // doesn't look at the produced signals:
  int maxLength = 5;   // maximum kernel length
  int maxPrec   = 5;   // maximum table precision
  for(int length = 0; length <= maxLength; length++) {
    for(int prec = 1; prec <= maxPrec; prec++) {

      linBlep.setTablePrecision(prec);
      minBlep.setTablePrecision(prec);

      linBlep.setLength(length);
      minBlep.setLength(length);

      // delay = 0:
      linBlep.prepareForStep(0.0, 1.0);
      minBlep.prepareForStep(0.0, 1.0);

      // delay = 0.5:
      linBlep.prepareForStep(0.5, 1.0);
      minBlep.prepareForStep(0.5, 1.0);

      // delay = 1:
      linBlep.prepareForStep(1.0, 1.0);
      minBlep.prepareForStep(1.0, 1.0);
    }
  }

  // check linear polyblep:
  rsPolyBlep1<double, double> pb1;
  pb1.reset();
  pb1.prepareForStep(0.0, 1.0);
  r &= pb1.getDelayed()   ==  0.0;
  r &= pb1.getCorrector() == -0.5;
  pb1.reset();
  pb1.prepareForStep(0.25, 1.0);
  r &= pb1.getDelayed()   ==  0.03125;
  r &= pb1.getCorrector() == -0.28125;
  pb1.reset();  
  pb1.prepareForStep(0.5, 1.0);
  r &= pb1.getDelayed()   ==  0.125;
  r &= pb1.getCorrector() == -0.125;
  pb1.reset();  
  pb1.prepareForStep(0.75, 1.0);
  r &= pb1.getDelayed()   ==  0.28125;
  r &= pb1.getCorrector() == -0.03125;
  pb1.reset();  
  pb1.prepareForStep(1.0, 1.0);
  r &= pb1.getDelayed()   ==  0.5;
  r &= pb1.getCorrector() == -0.0;
  // sanity-check, if these are really the values that should be expected - these are the ones that
  // currently come out and look good in an informal anti-aliasing test

  // check cubic polyblep:
  //rsPolyBlep2<double, double> pb2;

  return r;
}


class rsTestSyncPhasor : public rsSyncPhasor<double, rsPolyBlep1<double, double>>
{
public:
  void setState(double newMasterPos, double newMasterInc, double newSlavePos, double newSlaveInc, 
    bool resetBlep = true) 
  { 
    masterPos = newMasterPos;
    masterInc = newMasterInc;
    slavePos  = newSlavePos;
    slaveInc  = newSlaveInc;
    if(resetBlep)
      blep.reset();
  }
};
bool syncUnitTest()
{
  bool r = true;

  // We put an rsSyncPhasor into a well known state (consisting of the positions ind increments of 
  // the master and slave and the state of embedded blep object) and let it produce two samples and 
  // compare the results with the prediction/expectation. Two samples are enough, because we use 
  // only a linear polyblep which corrects only two samples. We do this for different sorts of 
  // states that execute different code paths (for wrap-arounds for master-only, slave-only, 
  // master-then-slave, slave-then-master). We also compare the naive (uncorrected) synced output
  // samples to our expectations

  rsTestSyncPhasor sp;
  double x0, x1;  // naive outputs
  double y0, y1;  // corrected outputs


  // no reset:
  sp.setState(0.3, 0.1, 0.3, 0.01);
  x0 = sp.getSampleNaive();   //  0.31    correct
  y0 = sp.applyBlep(x0);      //  0.0     correct because of delay
  x1 = sp.getSampleNaive();   //  0.32    correct
  y1 = sp.applyBlep(x1);      //  0.31    correct

  // master reset only at d=0.5:
  sp.setState(0.95, 0.1, 0.3, 0.01);
  x0 = sp.getSampleNaive();   //  0.005    correct
  y0 = sp.applyBlep(x0);      // -0.03875
  x1 = sp.getSampleNaive();   //  0.015    correct
  y1 = sp.applyBlep(x1);      //  0.04375

  // slave reset only at d=0.5:
  sp.setState(0.3, 0.1, 0.995, 0.01);
  x0 = sp.getSampleNaive();   //  0.005    correct
  y0 = sp.applyBlep(x0);      // -0.125
  x1 = sp.getSampleNaive();   //  0.015    correct
  y1 = sp.applyBlep(x1);      //  0.13

  // master and slave resets occur simultaneously with sample-delay of d=0.5:
  sp.setState(0.95, 0.1, 0.995, 0.01);
  x0 = sp.getSampleNaive();   //  0.005     correct: master-reset has no effect
  y0 = sp.applyBlep(x0);      // -0.125625
  x1 = sp.getSampleNaive();   //  0.015     correct
  y1 = sp.applyBlep(x1);      //  0.130625

  // master reset first (d=0.5), then slave reset (d=0.4):
  sp.setState(0.95, 0.1, 0.994, 0.01);
  x0 = sp.getSampleNaive();   //  0.004    correct: master-reset has no effect, bcs slave reset overwrites state
  y0 = sp.applyBlep(x0);      // -0.12558
  x1 = sp.getSampleNaive();   //  0.014    correct
  y1 = sp.applyBlep(x1);      //  0.12968
  // what happens in x0 = sp.getSampleNaive() to the state/position of the slave phasor
  // t in samples: -1: previous, 0: current, the (0) in the middle means where it would end up at 
  // time 0 if there would not be a second reset, p is position
  // t: -1       -0.5    (0)      -0.4     0        
  // p:  0.994    0       0.005    0       0.004

  // slave reset first (d=0.6), then master reset (d=0.5):
  sp.setState(0.95, 0.1, 0.996, 0.01);
  x0 = sp.getSampleNaive();   //  0.005    correct: master-rest overrides slave reset
  y0 = sp.applyBlep(x0);      // -0.18075
  x1 = sp.getSampleNaive();   //  0.015    correct
  y1 = sp.applyBlep(x1);      //  0.08575




  // i have written the produced numbers into the comments - figure out, if these match the 
  // expectations (with pencil and paper), if so, turn comments into r &= ... checks
  // check also, master-only and slave-only cases ...and maybe the trivial no-reset case, too


  //return false;  // this test is not yet complete and the phasor still doesn't work correctly
  return r;
}

bool testSpectrogramResynthesis(int blockSize, int hopSize, int signalLength, int fftSize, 
  RAPT::rsWindowFunction::WindowType windowType 
  = RAPT::rsWindowFunction::WindowType::hanningZN)
{
  bool r = true;

  int B = blockSize;
  int H = hopSize;
  int N = signalLength;
  int M = fftSize;

  // compute the complex spectrogram of a sequence of random numbers:
  std::vector<double> x = rsRandomVector(N, -1, +1);
  RAPT::rsSpectrogramProcessor<double> sp;  // spectrogram processor
  sp.setAnalysisWindowType(windowType);
  sp.setSynthesisWindowType(windowType);
  //sp.setTrafoSize(M);
  //sp.setBlockSize(B);
  sp.setBlockAndTrafoSize(B, M);
  sp.setHopSize(H);
  rsMatrix<rsComplexDbl> s = sp.getComplexSpectrogram(&x[0], N);
  int numFrames = s.getNumRows();  
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
  //if(r == false)
  //  plotVector(err);  // plot error signal, if something goes wrong

  //// plot the modulation signal resulting from analysis/resynthesis roundtrip:
  //std::vector<double> mod = sp.getRoundTripModulation(numFrames);
  //plotVector(mod);

  return r;
}

// a subclass of rsSpectrogramProcessor, so we can access protected variables from test code
class rsSpectrogramUnitTest : public RAPT::rsSpectrogramProcessor<double>
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
    //ifft(  &y[0], N); RAPT::rsArrayTools::scale(&y[0], N, 1.0/N); // actual
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
    r &= testTrafos(  2, 1.e-15);
    r &= testTrafos(  4, 1.e-15);
    r &= testTrafos(  8, 1.e-15);
    r &= testTrafos( 16, 1.e-14);
    r &= testTrafos(128, 1.e-13);
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

  rsSpectrogramUnitTest tester;
  r &= tester.runTests();

  // parameters to the function (in order)
  // B: block size
  // H: hop size
  // N: signal length
  // M: FFT size

  // maybe use the parameter order: N, M, B, H (values should decrease, 
  // i.e. N >= M >= B > H...well, maybe not necessarrily N >= M - but that will be the typical
  // case for real world signals)


  // try under which conditions the analysis/resynthesis roundtrip is an identity operation:
  r &= testSpectrogramResynthesis(  8,  4, 20,  8);   // B = M = 2^k, H = B/2   - works
  r &= testSpectrogramResynthesis(  6,  3, 20,  6);   // B = M = even, H = B/2  - works
  r &= testSpectrogramResynthesis( 12,  4, 50, 16);   // B even, M = 2^k > B, H = B/3 - works
  //r &= testSpectrogramResynthesis(7, 3, 20, 8);     // B odd, M=2^k > B - fails!
  //r &= testSpectrogramResynthesis(7, 3, 20, 7);       // B = M odd, H = floor(B/2) - fails spectacularly!!
  //r &= testSpectrogramResynthesis( 13,  4, 50, 16);   // fails
  r &= testSpectrogramResynthesis(128, 64, 500, 128);  // H = B/2
  r &= testSpectrogramResynthesis(128, 32, 500, 128);  // H = B/4 - windows add up to 1.5 and there's a fade-in artifact
  r &= testSpectrogramResynthesis(128, 32, 500, 256);  // M = 2*B -> zero padding by factor 2

  // ...eventually the ideal goal would be that identity resynthesis works for any combination of
  // blockSize, trafoSize >= blockSize, hopSize < blockSize

  // todo: try block-sizes that are not a power of two, window functions that don't have natural
  // perfect reconstruction properties (i.e. demodulate the output)
  // test resynthesis with zero-padding
  // test resynthesis with and without demodulation (without, we should use 
  // window/blocksize/hopsize that doesn't need demodulation)

  // let the spectrogram class have a direct setFftSize function, let it use an 
  // rsFourierTransformer object, allow arbitrary FFT sizes

  // why do we need to scale the output of the ifft? in getComplexSpectrogram, there is already a 
  // scaling by 2 / rsArrayTools::sum(w, B); ...aahh - but it's applied to the STFT matrix *after* the
  // STFT has been computed - so, we should probably use no normalziation

  // how it should work:
  // -on forward FFT, scale spectral values by 1 / sum(window) (or 2 / sum(..) bcs of negative 
  //  freqs)
  // -on inverse FFT, no scaling should be applied

  // how it actually works:
  // -on forward FFT, it is scaled as is should
  // -on inverse FFT, it is scaled by 1/N - why does this work? there must be a hidden inverse
  //  scaling somewhere - maybe in applyDemodulation? -> figure out


  // this fails the transformer object seems to compute a scrambled forward FFT when it should 
  // compute an inverse FFT - but testVariousFourierTransforms passes - so where is the bug?

  return r;
}

bool harmonicAnalyzerUnitTest()
{
  bool r = true;      // test result

  // test filling the FFT buffer (with zero-apdding and shifting for zero-phase at center):

  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;
  using HA  = RAPT::rsHarmonicAnalyzer<double>;

  Vec sig = { 1,2,3,4,5,6,7,8 };  // 8 elements

  // test with zero-padding factor = 4:
  Vec buf(32);
  AT::fillWithNaN(&buf[0], (int) buf.size());
  HA::prepareBuffer(sig, buf);
  Vec target = { 5,6,7,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4 };
  r &= buf == target;

  // test with zero-padding factor = 2:
  buf.resize(16);
  AT::fillWithNaN(&buf[0], (int) buf.size());
  HA::prepareBuffer(sig, buf);
  target = { 5,6,7,8,0,0,0,0,0,0,0,0,1,2,3,4 };
  r &= buf == target;

  // test with zero-padding factor = 1:
  buf.resize(8);
  AT::fillWithNaN(&buf[0], (int) buf.size());
  HA::prepareBuffer(sig, buf);
  target = { 5,6,7,8,1,2,3,4 };
  r &= buf == target;

  return r;
}

bool singleSineModelerUnitTest()
{
  bool r = true;

  using Vec = std::vector<double>;
  using SSM = rsSingleSineModeler<double>;

  int N = 1000; // number of samples
  double tol = 1.e-14;

  // Test to resynthesize white noise - the analysis data may be meaningless in this case, but 
  // identity resynthesis should work nevertheless:
  Vec x = rsRandomVector(N, -1.0, 1.0);  // input signal
  Vec a(N), w(N), p(N), pm(N);           // analysis data
  Vec y(N);                              // resynthesized signal
  Vec err;

  SSM ssm;

  ssm.analyzeAmpAndPhase(&x[0], N, &a[0], &p[0]);
  ssm.synthesizeFromAmpAndPhase(&a[0], &p[0], N, &y[0]);
  r &= rsAreVectorsEqual(x, y, tol);
  err = x-y;
  // todo: synthesize and compare

  //y[0] = 0;
  rsFill(y, 0.0);

  ssm.analyzeAmpFreqAndPhaseMod(&x[0], N, &a[0], &w[0], &pm[0]);
  ssm.synthesizeFromAmpFreqPhaseMod(&a[0], &w[0], &pm[0], N, &y[0]);

  err = x-y;
  // only the first sample is wrong


  // todo: test to analyze a perfect sinewave and see, if the analysis data makes sense....then 
  // maybe make it more difficult by introducing a frequency sweep, amplitude fade, etc....



  return r;
}

bool sineModelingUnitTest()
{
  bool r = true;

  r &= harmonicAnalyzerUnitTest();
  r &= singleSineModelerUnitTest();

  return r;
}
