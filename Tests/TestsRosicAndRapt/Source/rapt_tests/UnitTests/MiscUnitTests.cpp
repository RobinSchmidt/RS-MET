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

template<class TOsc>
void plotSyncOutput(TOsc& osc, int numSamples)
{
  std::vector<double> yNaive(numSamples), yBlep(numSamples);
  for(int n = 0; n < numSamples; n++) {
    yNaive[n] = osc.getSampleNaive();
    yBlep[n]  = osc.applyBlep(yNaive[n]);
  }
  RAPT::rsArray::shift(&yBlep[0], numSamples, -osc.getBlep().getDelay());
  rsPlotVectors(yNaive, yBlep);
}
bool syncUnitTest()
{
  bool r = true;

  // we use PolyBlep1 because it has a blep corrector whose desired state can be most easily 
  // predicted
  typedef rsSyncPhasor<double, rsPolyBlep1<double, double>> SP;
  SP sp;

  // maybe this code should go into an experiment...like
  // testSyncPhasorPlot(masterInc, slaveInc, numSamples)
  // testSnycPhasorWrite ...does the same, but writes result to wavefile

  sp.setMasterIncrement( 31./1024);  //  < 1/32  ..try 30/1024 too
  sp.setSlaveIncrement( 128./1024);  // == 1/8
  sp.reset();
  plotSyncOutput(sp, 3200);
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 4x
  // hmm - actually here even the naive signal looks strange


  sp.setMasterIncrement( 32./1024);  // == 1/32
  sp.setSlaveIncrement( 129./1024);  //  > 1/8
  sp.reset();
  plotSyncOutput(sp, 500);
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 3x
  // that looks pretty good

  // slaveFreq = 4*masterFreq: master and slave resets always occur exactly simultaneously - we 
  // expect the master reset to be inconsequential, i.e. the code for master-reset is called and 
  // the blep is prepared, but always with zero step-amplitude:
  sp.setMasterIncrement( 32./1024);  // == 1/32
  sp.setSlaveIncrement( 128./1024);  // == 1/8
  sp.reset();
  plotSyncOutput(sp, 100);  // just for development - produce an output array and plot it
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 3x

  // we expect the blep-signal to be equal to the naive signal (but one sample delayed) because the
  // slave resets occur exactly at integer sample instants and master resets are supposed to be 
  // inconsequential anyway

  // maybe we should manually walk through the samples that are supposed to be generated and define
  // the desired output signals and then check, if they are actually produced (don't check internal
  // states of the osc and its embedded blep directly - just look at the output signal)

  // hmm...that doesn't look like a pure delay
  // ahh - not we should expect a delay and some sort of lowpass character - that seems to fit with
  // what we see
  // the slave-wraparounds write 0.5 nonzero numbers into corrector of the blep (i think, the delayed
  // sample remains untouched). the master wraparounds are indeed inconsequential

  // maybe we should first unit-test ts rsPolyBlep1 object - why would it write 0.5 into the corrector 
  // when the wraparound occurs exactly at the sample-instant?
  // rsPolyBlep1::prepareForStep does this:
  //   corrector += a * (-d2/2 + d - 1./2);
  // this creates the .5 in the corrector - is this really correct? ..but the informal anti-aliasing 
  // tests with it look good - they are probably correct and what we see is the lowpass character
  // that affects the passband


  sp.setMasterIncrement( 32./1024);  // == 1/32
  sp.setSlaveIncrement( 127./1024);  //  < 1/8
  sp.reset();
  plotSyncOutput(sp, 500);
  // alternately enters the "master-only" branch 1x and "slave-only" branch 3x
  // that looks wrong! i think, it exposes the artifact that i'm battling with
  // ..or..well - the black (naive) looks wrong but the blue (blepped) looks actually nice



  // we need a case that sometimes enters the "master-then-slave" branch - i think, we need a 
  // master increment that causes master-wrap-arounds at non-integers
  // ...maybe 31./1024, 33./1024, 31./1023 or so








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
  RAPT::rsSpectrogram<double> sp;  // spectrogram processor
  sp.setAnalysisWindowType(windowType);
  sp.setSynthesisWindowType(windowType);
  //sp.setTrafoSize(M);
  //sp.setBlockSize(B);
  sp.setBlockAndTrafoSize(B, M);
  sp.setHopSize(H);
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(&x[0], N);
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
    //ifft(  &y[0], N); RAPT::rsArray::scale(&y[0], N, 1.0/N); // actual
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

  // why do we need to scale the output of the ifft? in complexSpectrogram, there is already a 
  // scaling by 2 / rsArray::sum(w, B); ...aahh - but it's applied to the STFT matrix *after* the
  // STFT has been computed - so, we should probably use no nromalziation

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

bool sineModelingUnitTest()
{
  bool r = true;      // test result

  // test filling the FFT buffer (with zero-apdding and shifting for zero-phase at center):

  typedef std::vector<double> Vec;
  Vec sig = { 1,2,3,4,5,6,7,8 };  // 8 elements

  // test with zero-padding factor = 4:
  Vec buf(32);
  RAPT::rsArray::fillWithNaN(&buf[0], (int) buf.size());
  rsHarmonicAnalyzer<double>::prepareBuffer(sig, buf);
  Vec target = { 5,6,7,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4 };
  r &= buf == target;

  // test with zero-padding factor = 2:
  buf.resize(16);
  RAPT::rsArray::fillWithNaN(&buf[0], (int) buf.size());
  rsHarmonicAnalyzer<double>::prepareBuffer(sig, buf);
  target = { 5,6,7,8,0,0,0,0,0,0,0,0,1,2,3,4 };
  r &= buf == target;

  // test with zero-padding factor = 1:
  buf.resize(8);
  RAPT::rsArray::fillWithNaN(&buf[0], (int) buf.size());
  rsHarmonicAnalyzer<double>::prepareBuffer(sig, buf);
  target = { 5,6,7,8,1,2,3,4 };
  r &= buf == target;


  return r;
}