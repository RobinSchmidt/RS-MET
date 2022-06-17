//#include "rosic_FilterTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
//#include "../Shared/Plotting/rosic_Plotter.h"
using namespace rosic;

template<class T>
void getImpulseResponse(T &module, double *h, int N)
{
  module.reset();
  h[0] = module.getSample(1.0);
  for(int n = 1; n < N; n++)
    h[n] = module.getSample(0.0);
}

template<class T>
void plotImpulseResponse(T &module, int N, double fs)
{
  double *t = new double[N];  // time axis
  double *h = new double[N];  // impulse response

  RAPT::rsArrayTools::fillWithRangeLinear(t, N, 0.0, (N-1)/fs);
  getImpulseResponse(module, h, N);
  plotData(N, t, h);

  delete[] t;
  delete[] h;
}

template<class T>
void writeImpulseResponseToFile(const char *path, T &module, int N, int fs, int numBits)
{  
  double *h = new double[N];  
  getImpulseResponse(module, h, N);

  RAPT::rsArrayTools::normalize(h, N, 1.0);

  rosic::writeToMonoWaveFile(path, h, N, fs, numBits);
  delete[] h;
}

void rotes::testLadderFilter()
{
  double fs   = 44100;  // samplerate in Hz
  double fcHi = 4000;   // 1st cutoff frequency
  double fcLo = 250;    // 2nd cutoff frequency 
  double fSaw = 120;    // frequency of the input sawtooth - 120 seems to be a problem-freq

  rosic::LadderFilterOld ladder;
  ladder.setSampleRate(fs);
  ladder.setResonance(0.95);
  //ladder.setResonance(0.0);
  //ladder.setResonance(0.5);

  /*
  static const int N = 500;  
  double t[N];  // time axis
  double yLo[N], yHi[N], yHiLo[N], yLoHi[N];  
    // output signals for low and high cutoff and the two switches lo->hi, hi->lo

  rosic::fillWithIndex(t, N);
  int n;

  // low cutoff output without switch:
  ladder.reset();
  ladder.setCutoff(fcLo, true);
  ladder.getSampleTest(1.0);  // only to se up the internal states, output irrelevant
  for(n = 0; n < N; n++)
    yLo[n] = ladder.getSampleTest(0.0);

  // low cutoff output with switch from high cutoff:
  ladder.reset();
  ladder.setCutoff(fcHi, true);
  ladder.getSampleTest(1.0);   
  ladder.setCutoff(fcLo, true);
  for(n = 0; n < N; n++)
    yHiLo[n] = ladder.getSampleTest(0.0);

  Plotter::plotData(N, t, yLo, yHiLo);
  */

  static const int N = 6000;
  double t[N];  // time axis
  double x[N];  // input signal
  double y[N];  // output signal

  rosic::synthesizeWaveform(x, N, rosic::SAW, fSaw, fs);
  RAPT::rsArrayTools::fillWithIndex(t, N);
  RAPT::rsArrayTools::fillWithZeros(y, N);

  int n;
  ladder.setCutoff(fcHi, true);
  for(n = 0; n < N/3; n++)
    y[n] = ladder.getSampleTest(x[n]);
  ladder.setCutoff(fcLo, true);
  for(n = N/3; n < 2*N/3; n++)
    y[n] = ladder.getSampleTest(x[n]);
  ladder.setCutoff(fcHi, true);
  for(n = 2*N/3; n < N; n++)
    y[n] = ladder.getSampleTest(x[n]);

  plotData(N, t, x, y);
  // middle section looks terrible - filter seems to explode? -> investigate!

  int dummy = 0;
}

void rotes::testModalFilter()
{
  // user parameters:
  static const int N = 20000;   // # samples to plot
  double fs  = 44100;  // samplerate in Hz
  double td  = 0.1;    // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.5;    // amplitude as raw factor
  double pm  = 0.0;    // phase-modulation (for the 2nd implementation)


  //ModalFilter mf;
  RAPT::rsModalFilter<double, double> mf;
  mf.setModalParameters(f, A, td, phs, fs);

  //ModalFilter2 mf2;
  RAPT::rsNonlinearModalFilter<double, double> mf2;
  mf2.setModalParameters(f, A, td, phs, fs);
  mf2.setPhaseModulation(pm);


  double t[N], h[N], h2[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  h[0]  = mf.getSample(1.0);
  h2[0] = mf2.getSample(1.0);
  for(int n = 1; n < N; n++)
  {
    //mf.setModalParameters(f+0.5*f*sin(0.01*n), A, td, phs, fs);
    h[n]  = mf.getSample(0.0);

    //mf2.setModalParameters(f, A, td, phs+n, fs);
    //mf2.setModalParameters(f+0.5*f*sin(0.01*n), A, td, phs, fs);
    h2[n] = mf2.getSample(0.0);
  }

  plotData(N, t, h, h2);
  //plotData(N, t, h2);


  //writeImpulseResponseToFile("d:\\TmpData\\ModalImpulseResponse.wav", mf2, (int)fs, (int)fs, 16);
}

void rotes::testModalFilterWithAttack()
{
  // user parameters:
  double td  = 0.100;  // asymptotic decay time constant in seconds
  double tp  = 0.050;  // time of peak occurence in seconds (must be < td)
  double f   = 440;    // frequency in Hz
  double dt  = +0.0;   // detune for the quickly decaying sinusoid in semitones (adds roughness)
  double phs = 90;     // phase in degrees
  double A   = 1.0;    // amplitude as raw factor
  double fs  = 44100;  // samplerate in Hz

  double df  = RAPT::rsPitchOffsetToFreqFactor(dt);

  //ModalFilterWithAttack mf;
  RAPT::rsModalFilterWithAttack<double, double> mf;
  mf.setModalParameters(f, A, tp, td, phs, fs, df);
  //plotImpulseResponse(mf, (int)fs, (int)fs); // why fs two times?
  plotImpulseResponse(mf, (int)fs, 1.0);
  //writeImpulseResponseToFile("d:\\TmpData\\ModalImpulseResponse.wav", mf, fs, fs, 16);
}

void rotes::testBiquadPhasePlot()
{
  int i;
  static const int numBiquads = 10;
  static const int numBins    = 500;

  // allpass parameters:
  double fs = 44100.0;        // sample-rate
  double fc = 0.5*fs/PI;      // characteristic frequency (such that wc = 1)
  double  q = 1.0/sqrt(2.0);  // quality factor
  q = 5.0;

  // create biquad cascade coefficients:
  double b0[numBiquads], b1[numBiquads], b2[numBiquads], a1[numBiquads], a2[numBiquads];
  //rosic::BiquadDesigner::calculateFirstOrderLowpassCoeffs(b0[0], b1[0], b2[0], a1[0], a2[0], 1.0/fs, fc);
  //rosic::BiquadDesigner::calculateFirstOrderLowpassCoeffsBilinear(b0[0], b1[0], b2[0], a1[0], a2[0], 1.0/fs, fc);
  rosic::BiquadDesigner::calculateCookbookAllpassCoeffs(b0[0], b1[0], b2[0], a1[0], a2[0], 1.0/fs, fc, q);
  a1[0] = -a1[0]; // needed because of different sign conventions in BiquadDesigner and FilterAnalyzer -> fix this
  a2[0] = -a2[0];
  for(i = 1; i < numBiquads; i++)
  {
    b0[i] = b0[0];
    b1[i] = b1[0];
    b2[i] = b2[0];
    a1[i] = a1[0];
    a2[i] = a2[0];
  }

  // create normalized radian frequency axis:
  double w[numBins];
  RAPT::rsArrayTools::fillWithRangeLinear(w, numBins, 0.0, PI);

  // obtain magnitude phase response:
  double m[numBins];
  double p[numBins];
  rosic::rsFilterAnalyzerD::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numBiquads, w, m, numBins, false, false);
  rosic::rsFilterAnalyzerD::getBiquadCascadePhaseResponse(b0, b1, b2, a1, a2, numBiquads, w, p, numBins, false);

  // plot the (magnitude- and) phase response:
  plotData(numBins, w, p);
  //plotData(numBins, w, m, p);
}

void rotes::testFiniteImpulseResponseDesigner()
{
  using FD = rosic::FiniteImpulseResponseDesigner;
  using WD = rosic::WindowDesigner;
  using AT = RAPT::rsArrayTools;
  using Window = WD::windowTypes;
  using Mode   = FD::modes;


  // User parameters:
  static const int length    = 201;      // Filter kernel length
  static const int fftLength = 8192;     // FFT length for plot
  double plotMinDb  = -200;              // Minimum for y-axis in the plot in dB
  double plotMaxDb  =  +10;              // Maximum ...
  double sampleRate = 44100.0;           // Sample rate
  double frequency  = 10000.0;           // Cutoff or center frequency
  Mode   mode       = Mode::BANDPASS;
  Window window     = Window::BLACKMAN;


  // Plot the impulse response:
  double impulseResponse[length];
  FD designer;
  designer.setMode(mode);
  designer.setWindowType(window);
  designer.setFrequency(frequency);
  designer.getImpulseResponse(impulseResponse, length);
  //designer.spectralReversal(impulseResponse, length);   // ?
  double indices[length];
  AT::fillWithIndex(indices, length);
  //plotData(length, indices, impulseResponse);

  // Plot the magnitude and phase response:
  static const int N = fftLength;  // shorter name for convenience
  double freqs[N];
  double mags[N];
  double phases[N];
  AT::fillWithIndex(freqs, N);
  AT::scale(freqs, freqs, N, sampleRate/N);
  fftMagnitudesAndPhases(impulseResponse, length, mags, phases, N);
  //plotData(fftLength/2, frequencies, magnitudes);
  //plotData(fftLength/2, frequencies, phases);




  // Plot more magnitude responses wit different windows:
  double mags1[N];
  double mags2[N];
  double mags3[N];
  double mags4[N];
  double mags5[N];

  FiniteImpulseResponseFilter filter;
  filter.setMode(mode);
  filter.setImpulseResponseLength(length);
  filter.setFrequency(frequency);

  filter.setWindowType(WD::RECTANGULAR);
  filter.getMagnitudeResponse(freqs, mags1, N, true, false);
  AT::clip(mags1, N, plotMinDb, plotMaxDb);

  filter.setWindowType(WD::BLACKMAN);
  filter.getMagnitudeResponse(freqs, mags2, N, true, false);
  AT::clip(mags2, N, plotMinDb, plotMaxDb);

  filter.setWindowType(WD::HAMMING);
  filter.getMagnitudeResponse(freqs, mags3, N, true, false);
  AT::clip(mags3, N, plotMinDb, plotMaxDb);

  filter.setWindowType(WD::HANN);
  filter.getMagnitudeResponse(freqs, mags4, N, true, false);
  AT::clip(mags4, N, plotMinDb, plotMaxDb);

  filter.setWindowType(WD::COSINE_SQUARED); // isn't this the same as Hann? ...not exactly..why?
  filter.getMagnitudeResponse(freqs, mags5, N, true, false);
  AT::clip(mags4, N, plotMinDb, plotMaxDb);

  plotData(N/2, freqs, mags1, mags2, mags3, mags4, mags5); // plot all 5
  plotData(N/2, freqs, mags4, mags5);                      // Hann vs cos^2
  plotData(N/2, freqs, mags4, mags2);                      // Hann vs Blackman
  plotData(N/2, freqs, mags2);                             // Blackman only

  /*
  // create some noise, filter it and write it into a file:
  static const int testLength = 44100;
  double noise[testLength];
  double filteredNoise[testLength+length-1];
  for(int n=0; n<testLength; n++)
    noise[n] = rosic::random(-1.0, 1.0);
  rosic::convolve(noise, testLength, impulseResponse, length, filteredNoise);

  writeToStereoWaveFile("D:/TmpAudio/FilteredNoise.wav", noise, filteredNoise, testLength, 44100, 16);
  */

  // test cosine power window overlap:


  int dummy = 0;
}


bool rotes::testConvolverPartitioned()
{
  bool ok = true;

  using AT = RAPT::rsArrayTools;

  static const int impulseLength  = 5000;
  static const int responseLength = 201;
  static const int resultLength   = impulseLength+responseLength-1;

  double impulse[impulseLength];
  double impulseResponse[responseLength];
  double resultTrue[resultLength];
  double result[resultLength];
  double indices[resultLength];
  AT::fillWithIndex(indices, impulseLength);

  AT::fillWithZeros(impulse, impulseLength);
  impulse[1000] = 1.0;

  // Create an impulse response (a ramp down from 1.5 to 1.0) and convole it with our input impulse
  // via the naive convolution algorithm for reference:
  AT::fillWithRangeLinear(impulseResponse, responseLength, 1.5, 1.0);
  AT::convolve(impulse, impulseLength, impulseResponse, responseLength, resultTrue);

  // Compute the convolution also via the partitioned FFT algorithm:
  ConvolverPartitioned convolver;
  convolver.setImpulseResponse(impulseResponse, responseLength);
  for(int n = 0; n < impulseLength; n++)
    result[n] = convolver.getSample(impulse[n]);
  for(int n = impulseLength; n < resultLength; n++)
    result[n] = convolver.getSample(0.0);

  // Check, if the results match:
  double tol = 1.e-15;  // that's a tight tolerance - it works with msc, though
  double maxDiff = AT::maxDeviation(resultTrue, result, resultLength);
  ok &= maxDiff <= tol;

  // plot:
  //plotData(impulseLength,  indices, impulse);
  //plotData(responseLength, indices, impulseResponse);
  //plotData(resultLength,   indices, resultTrue, result);
  // the result plot looks wrong but the data looks ok in the debugger - maybe a problem with the
  // plotting?

  // ToDo: 
  // -try it with various random inputs and random impulse responses of different lengths, be 
  //  sure to cover special cases (like length = 0,1,2,4,8 (powers of 2), 3,7,15, 5,9,17, ...)
  // -if all works, maybe move the convolver to rapt, test it with float, double, rsFloat64x2,
  //  rsSimdVector<double, 2>, rsSimdVector<float, 2>

  return ok;
}

bool rotes::testFiniteImpulseResponseFilter()
{
  // ToDo: move to unit tests

  bool ok = true;
  using AT  = RAPT::rsArrayTools;
  using Vec = std::vector<double>;

  // Create the filter and retrieve its impulse response h with length L:
  FiniteImpulseResponseFilter filter;
  filter.setMode(FiniteImpulseResponseDesigner::BANDPASS);
  filter.setFrequency(1000.0);
  int L = filter.getKernelLength();
  const double* h = filter.getKernelPointer();
  rsPlotArray(h, L);

  // Create some noise and apply the filter:
  static const int N = 2000;
  int M = N+L+10;
  Vec x(N), y(M);
  for(int n = 0; n < N; n++) {
    x[n] = random(-1.0, 1.0);
    y[n] = filter.getSample(x[n]); }
  for(int m = N; m < M; m++)
    y[m] = filter.getSample(0.0);

  // Create reference signal by naive convolution and compare:
  Vec r(M);
  AT::convolve(&x[0], N, &h[0], L, &r[0]);
  double maxDiff = AT::maxDeviation(&r[0], &y[0], M);
  double tol = 1.e-15; 
  ok &= maxDiff <= tol;
  rsPlotVectors(x, y, r);
  
  //writeToStereoWaveFile("D:/TmpAudio/FilteredNoise.wav", noise, filteredNoise, testLength, 44100, 16);

  // todo: test even and odd lengths

  return ok;
}

  
void rotes::testFilterAnalyzer()
{
  static const int maxOrder = 20;
  double filterFrequency    = 5000;
  double sampleRate         = 44100.0;
  int prototypeOrder        = 3;

  // the filter-designer:
  rosic::rsInfiniteImpulseResponseDesignerD designer;
  designer.setApproximationMethod(rosic::rsPrototypeDesignerD::BUTTERWORTH);
  designer.setMode(rosic::rsInfiniteImpulseResponseDesignerD::LOWPASS);
  designer.setFrequency(filterFrequency);
  designer.setSampleRate(sampleRate);
  designer.setPrototypeOrder(prototypeOrder);

  // biquad coeffs:
  static const int maxNumBiquads = maxOrder/2;   // maxOrder should be even - otherwise we need '+1'
  double b0[maxNumBiquads];
  double b1[maxNumBiquads];
  double b2[maxNumBiquads];
  double a1[maxNumBiquads];
  double a2[maxNumBiquads];
  RAPT::rsArrayTools::fillWithValue(b0, maxNumBiquads, -100.0);
  RAPT::rsArrayTools::fillWithValue(b1, maxNumBiquads, -100.0);
  RAPT::rsArrayTools::fillWithValue(b2, maxNumBiquads, -100.0);
  RAPT::rsArrayTools::fillWithValue(a1, maxNumBiquads, -100.0);
  RAPT::rsArrayTools::fillWithValue(a2, maxNumBiquads, -100.0);

  int numBiquads = designer.getNumBiquadStages();
  designer.getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);

  // frequencies and response-stuff:
  static const int numBins = 1000;
  double  frequencies[numBins];
  double  omegas[numBins];
  double  magnitudes[numBins];
  Complex H[numBins];
  for(int k=0; k<numBins; k++)
  {
    frequencies[k] = 0.5 * k * sampleRate / (numBins-1);
    omegas[k]      = 2.0 * PI * frequencies[k] / sampleRate;
  }

  //FilterAnalyzer::getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numBiquads, omegas, magnitudes, numBins, true);

  rosic::rsFilterAnalyzerD::getBiquadCascadeFrequencyResponse(         b0, b1, b2, a1, a2, numBiquads, omegas, rsCastPointer(H), numBins);
  rosic::rsFilterAnalyzerD::multiplyWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, omegas, rsCastPointer(H), numBins);
  rosic::rsFilterAnalyzerD::addWithBiquadCascadeFrequencyResponse(     b0, b1, b2, a1, a2, numBiquads, omegas, rsCastPointer(H), numBins);

  rosic::rsFilterAnalyzerD::getMagnitudes( rsCastPointer(H), magnitudes, numBins);
  rosic::rsFilterAnalyzerD::convertToDecibels(magnitudes, numBins);

  RAPT::rsArrayTools::clip(magnitudes, numBins, -60.0, 20.0);


  plotData(numBins, frequencies, magnitudes);

  int dummy = 0;
}

void rotes::testBiquadCascade()
{
  double sampleRate = 44100.0;

  rosic::rsInfiniteImpulseResponseDesignerD designer;
  rosic::rsBiquadCascadeDD            biquadCascade;

  designer.setApproximationMethod(rosic::rsPrototypeDesignerD::BUTTERWORTH);
  designer.setMode(rosic::rsInfiniteImpulseResponseDesignerD::BANDPASS);
  designer.setPrototypeOrder(5);
  designer.setSampleRate(sampleRate);
  designer.setFrequency(5000.0);
  designer.setBandwidth(1.0);

  biquadCascade.setNumStages(designer.getNumBiquadStages());
  designer.getBiquadCascadeCoefficients(biquadCascade.getAddressB0(), biquadCascade.getAddressB1(), 
    biquadCascade.getAddressB2(), biquadCascade.getAddressA1(), biquadCascade.getAddressA2() );
  //biquadCascade.turnIntoAllpass();

  // frequencies and response-stuff:
  static const int numBins = 1000;
  double  frequencies[numBins];
  double  omegas[numBins];
  double  magnitudes[numBins];
  double  phases[numBins];
  Complex H[numBins];
  for(int k=0; k<numBins; k++)
  {
    frequencies[k] = 0.5 * k * sampleRate / (numBins-1);
    omegas[k]      = 2.0 * PI * frequencies[k] / sampleRate;
  }

  biquadCascade.getFrequencyResponse(omegas, rsCastPointer(H), numBins);
  rosic::rsFilterAnalyzerD::getMagnitudes(rsCastPointer(H), magnitudes, numBins);
  rosic::rsFilterAnalyzerD::getPhases(    rsCastPointer(H), phases,     numBins);
  //rsFilterAnalyzer::convertToDecibels(magnitudes, numBins);

  biquadCascade.getMagnitudeResponse(omegas, magnitudes, numBins, true, false);
  biquadCascade.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);


  plotData(numBins, frequencies, magnitudes, phases);

  int dummy = 0;
}

void rotes::testCrossover4Way()
{
  using AT = RAPT::rsArrayTools;
  using FA = RAPT::rsFilterAnalyzer<double>;
  using Complex = std::complex<double>;

  rsCrossOver4Way<double, double> crossover;

  double sampleRate        = 44100.0;
  int    slope11           = 24;   // middlemost slope
  int    slope21           = 48;
  int    slope22           = 96;
  double lowCrossoverFreq  = 2000.0;
  double midCrossoverFreq  = 5000.0;
  double highCrossoverFreq = 10000.0;

  double lowClipValue      = -80.0;   // lower magnitude limit in dB for the plots
  crossover.setSampleRate(sampleRate);

  // frequencies and response-stuff:
  static const int numBins = 1000;
  double  frequencies[numBins];
  double  omegas[numBins];
  Complex H_LP_1_1[numBins];
  Complex H_HP_1_1[numBins];
  Complex H_LP_2_1[numBins];
  Complex H_HP_2_1[numBins];
  Complex H_LP_2_2[numBins];
  Complex H_HP_2_2[numBins];
  Complex H_Low[numBins];
  Complex H_LowMid[numBins];
  Complex H_HighMid[numBins];
  Complex H_High[numBins];
  Complex H_Sum[numBins];
  double  magnitudesLP_1_1[numBins];
  double  magnitudesHP_1_1[numBins];
  double  magnitudesLow[numBins];
  double  magnitudesLowMid[numBins];
  double  magnitudesHigh[numBins];
  double  magnitudesHighMid[numBins];
  double  magnitudesSum[numBins];
  for(int k=0; k<numBins; k++)
  {
    frequencies[k] = 0.5 * k * sampleRate / (numBins-1);
    omegas[k]      = 2.0 * PI * frequencies[k] / sampleRate;
  }

  //-----------------------------------------------------------------------------------------------
  // 2 ways:

  crossover.setBandActive(false, 1, 0);
  crossover.setBandActive(false, 1, 1);
  crossover.setCrossoverFrequency(lowCrossoverFreq,  1, 0);
  crossover.setCrossoverFrequency(highCrossoverFreq, 1, 1);
  crossover.setCrossoverFrequency(midCrossoverFreq,  0, 0);

  crossover.setSlope(slope11, 0, 0);
  crossover.setSlope(slope21, 1, 0);
  crossover.setSlope(slope22, 1, 1);

  // Obtain lowpass frequency response:
  crossover.getStage1().getLowpassFrequencyResponse(frequencies, H_LP_1_1, numBins, false);
  FA::getMagnitudes(H_LP_1_1, magnitudesLP_1_1, numBins);
  FA::convertToDecibels(magnitudesLP_1_1, numBins);
  AT::clip(magnitudesLP_1_1, numBins, lowClipValue, 10.0);

  // Obtain highpass frequency response:
  crossover.getStage1().getHighpassFrequencyResponse(frequencies, H_HP_1_1, numBins, false);
  FA::getMagnitudes(H_HP_1_1, magnitudesHP_1_1, numBins);
  FA::convertToDecibels(magnitudesHP_1_1, numBins);
  AT::clip(magnitudesHP_1_1, numBins, lowClipValue, 10.0);

  // obtain sum frequency response:
  AT::add(H_LP_1_1, H_HP_1_1, H_Sum, numBins);
  FA::getMagnitudes(H_Sum, magnitudesSum, numBins);
  FA::convertToDecibels(magnitudesSum, numBins);

  AT::copy(magnitudesLP_1_1, magnitudesLow,  numBins);  // why?
  AT::copy(magnitudesHP_1_1, magnitudesHigh, numBins);  // why?
  plotData(numBins, frequencies, magnitudesLow, magnitudesHigh, magnitudesSum);
  // Strange: the lowpass looks linear even though we have a linear freq-axis (it should look
  // linear in a log-log plot)

  //-----------------------------------------------------------------------------------------------
  // 3 ways (lower band split further):

  crossover.setBandActive(true, 1, 0);

  // obtain frequency response for low band:
  crossover.getStage1().getLowpassFrequencyResponse(   frequencies, H_LP_1_1, numBins, false);
  crossover.getStage2(0).getLowpassFrequencyResponse(frequencies, H_LP_2_1, numBins, false);
  AT::multiply(H_LP_1_1, H_LP_2_1, H_Low, numBins);
  FA::getMagnitudes(H_Low, magnitudesLow, numBins);
  FA::convertToDecibels(magnitudesLow, numBins);
  AT::clip(magnitudesLow, numBins, lowClipValue, 10.0);

  // obtain frequency response for low-mid band:
  crossover.getStage2(0).getHighpassFrequencyResponse(frequencies, H_HP_2_1, numBins, false);
  AT::multiply(H_LP_1_1, H_HP_2_1, H_LowMid, numBins);
  FA::getMagnitudes(H_LowMid, magnitudesLowMid, numBins);
  FA::convertToDecibels(magnitudesLowMid, numBins);
  AT::clip(magnitudesLowMid, numBins, lowClipValue, 10.0);

  // obtain frequency response for high band:
  crossover.getStage1().getHighpassFrequencyResponse(frequencies, H_HP_1_1, numBins, false);
  AT::copy(H_HP_1_1, H_High, numBins);
  crossover.getHighBranchAllpass().getFrequencyResponse(omegas, H_High, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_High, magnitudesHigh, numBins);
  FA::convertToDecibels(magnitudesHigh, numBins);
  AT::clip(magnitudesHigh, numBins, lowClipValue, 10.0);

  // obtain the summed frequency response:
  AT::add(H_Low, H_LowMid, H_Sum, numBins);
  AT::add(H_Sum, H_High,   H_Sum, numBins);
  FA::getMagnitudes(H_Sum, magnitudesSum, numBins);
  FA::convertToDecibels(magnitudesSum, numBins);

  plotData(numBins, frequencies, magnitudesLow, magnitudesLowMid, magnitudesHigh, magnitudesSum);

  //-----------------------------------------------------------------------------------------------
  // 3 ways (upper band split further):

  crossover.setBandActive(false, 1, 0);
  crossover.setBandActive(true,  1, 1);

  // obtain frequency response for low band:
  crossover.getStage1().getLowpassFrequencyResponse(frequencies, H_LP_1_1, numBins, false);
  AT::copy(H_LP_1_1, H_Low, numBins);
  crossover.getLowBranchAllpass().getFrequencyResponse(omegas, H_Low, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_Low, magnitudesLow, numBins);
  FA::convertToDecibels(magnitudesLow, numBins);
  AT::clip(magnitudesLow, numBins, lowClipValue, 10.0);

  // obtain frequency response for high-mid band:
  crossover.getStage1().getHighpassFrequencyResponse(  frequencies, H_HP_1_1, numBins, false);
  crossover.getStage2(1).getLowpassFrequencyResponse(frequencies, H_LP_2_2, numBins, false);
  AT::multiply(H_HP_1_1, H_LP_2_2, H_HighMid, numBins);
  FA::getMagnitudes(H_HighMid, magnitudesHighMid, numBins);
  FA::convertToDecibels(magnitudesHighMid, numBins);
  AT::clip(magnitudesHighMid, numBins, lowClipValue, 10.0);

  // obtain frequency response for high band:
  crossover.getStage2(1).getHighpassFrequencyResponse(frequencies, H_HP_2_2, numBins, false);
  AT::multiply(H_HP_1_1, H_HP_2_2, H_High, numBins);
  FA::getMagnitudes(H_High, magnitudesHigh, numBins);
  FA::convertToDecibels(magnitudesHigh, numBins);
  AT::clip(magnitudesHigh, numBins, lowClipValue, 10.0);

  // obtain the summed frequency response:
  AT::add(H_Low, H_HighMid, H_Sum, numBins);
  AT::add(H_Sum, H_High,    H_Sum, numBins);
  FA::getMagnitudes(H_Sum, magnitudesSum, numBins);
  FA::convertToDecibels(magnitudesSum, numBins);

  plotData(numBins, frequencies, magnitudesLow, magnitudesHighMid, magnitudesHigh, magnitudesSum);

  //-----------------------------------------------------------------------------------------------
  // 4 ways (upper band split further):

  crossover.setBandActive(true, 1, 0);
  crossover.setBandActive(true, 1, 1);

  // obtain frequency response for low band:
  crossover.getStage1().getLowpassFrequencyResponse(   frequencies, H_LP_1_1, numBins, false);
  crossover.getStage2(0).getLowpassFrequencyResponse(frequencies, H_LP_2_1, numBins, false);
  AT::multiply(H_LP_1_1, H_LP_2_1, H_Low, numBins);
  crossover.getLowBranchAllpass().getFrequencyResponse(omegas, H_Low, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_Low, magnitudesLow, numBins);
  FA::convertToDecibels(magnitudesLow, numBins);
  AT::clip(magnitudesLow, numBins, lowClipValue, 10.0);

  // obtain frequency response for low-mid band:
  crossover.getStage2(0).getHighpassFrequencyResponse(frequencies, H_HP_2_1, numBins, false);
  AT::multiply(H_LP_1_1, H_HP_2_1, H_LowMid, numBins);
  crossover.getLowBranchAllpass().getFrequencyResponse(omegas, H_LowMid, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_LowMid, magnitudesLowMid, numBins);
  FA::convertToDecibels(magnitudesLowMid, numBins);
  AT::clip(magnitudesLowMid, numBins, lowClipValue, 10.0);

  // obtain frequency response for high-mid band:
  crossover.getStage1().getHighpassFrequencyResponse(  frequencies, H_HP_1_1, numBins, false);
  crossover.getStage2(1).getLowpassFrequencyResponse(frequencies, H_LP_2_2, numBins, false);
  AT::multiply(H_HP_1_1, H_LP_2_2, H_HighMid, numBins);
  crossover.getHighBranchAllpass().getFrequencyResponse(omegas, H_HighMid, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_HighMid, magnitudesHighMid, numBins);
  FA::convertToDecibels(magnitudesHighMid, numBins);
  AT::clip(magnitudesHighMid, numBins, lowClipValue, 10.0);

  // obtain frequency response for high band:
  crossover.getStage2(1).getHighpassFrequencyResponse(frequencies, H_HP_2_2, numBins, false);
  AT::multiply(H_HP_1_1, H_HP_2_2, H_High, numBins);
  crossover.getHighBranchAllpass().getFrequencyResponse(omegas, H_High, numBins, 
    FA::MULTIPLICATIVE_ACCUMULATION);
  FA::getMagnitudes(H_High, magnitudesHigh, numBins);
  FA::convertToDecibels(magnitudesHigh, numBins);
  AT::clip(magnitudesHigh, numBins, lowClipValue, 10.0);

  // obtain the summed frequency response:
  AT::add(H_Low, H_LowMid,  H_Sum, numBins);
  AT::add(H_Sum, H_HighMid, H_Sum, numBins);
  AT::add(H_Sum, H_High,    H_Sum, numBins);
  FA::getMagnitudes(H_Sum, magnitudesSum, numBins);
  FA::convertToDecibels(magnitudesSum, numBins);

  plotData(numBins, frequencies, magnitudesLow, magnitudesLowMid, magnitudesHighMid, 
    magnitudesHigh, magnitudesSum);
}

void rotes::testCrossover4Way2()
{ 
  using AT = RAPT::rsArrayTools;

  static const int N       = 4096;    // numSamples, numBins
  double sampleRate        = 44100.0;
  int    slope11           = 48;      // middlemost slope
  int    slope21           = 48;
  int    slope22           = 48;
  double lowCrossoverFreq  = 250.0;
  double midCrossoverFreq  = 1000.0;
  double highCrossoverFreq = 6000.0;
  double lowClipValue      = -40.0;   // lower magnitude limit in dB for the plots

  // Set up the crossover:
  rosic::rsCrossOver4WayStereo crossover;
  crossover.setSampleRate(sampleRate);
  crossover.setCrossoverFrequency(lowCrossoverFreq,  1, 0);
  crossover.setCrossoverFrequency(highCrossoverFreq, 1, 1);
  crossover.setCrossoverFrequency(midCrossoverFreq,  0, 0);
  crossover.setSlope(slope11, 0, 0);
  crossover.setSlope(slope21, 1, 0);
  crossover.setSlope(slope22, 1, 1);
  crossover.setBandActive(true, 1, 0);
  crossover.setBandActive(true, 1, 1);

  // Obtain and plot the impulse-responses:
  double indices[N];
  AT::fillWithIndex(indices, N);
  double impulseResponses[8*N];
  double impulseResponseSum[N];
  auto recordSample = [&](int n, double x)
  {
    double tmp[8];
    AT::fillWithValue(tmp, 8, x);
    crossover.processSampleFrameStereo(tmp);
    impulseResponses[0*N+n] = tmp[0];
    impulseResponses[2*N+n] = tmp[2];
    impulseResponses[4*N+n] = tmp[4];
    impulseResponses[6*N+n] = tmp[6];
    impulseResponseSum[n] = tmp[0] + tmp[2] + tmp[4] + tmp[6];
  };
  recordSample(0, 1.0);
  for(int n = 1; n < N; n++)
    recordSample(n, 0.0);
  plotData(N, indices, &impulseResponses[0], &impulseResponses[2*N], 
    &impulseResponses[4*N], &impulseResponses[6*N], impulseResponseSum);
 
  // Compute and plot magnitude responses:
  double frequencies[N];
  double omegas[N];
  double magnitudesLow[N];
  double magnitudesLowMid[N];
  double magnitudesHigh[N];
  double magnitudesHighMid[N];
  double magnitudesSum[N];
  for(int k = 0; k < N; k++)
  {
    frequencies[k] = k * sampleRate / N;
    omegas[k]      = 2.0 * PI * frequencies[k] / sampleRate;
  }
  fftMagnitudesAndPhases(&impulseResponses[0],   N, magnitudesLow,     NULL, N, true);
  fftMagnitudesAndPhases(&impulseResponses[2*N], N, magnitudesLowMid,  NULL, N, true);
  fftMagnitudesAndPhases(&impulseResponses[4*N], N, magnitudesHighMid, NULL, N, true);
  fftMagnitudesAndPhases(&impulseResponses[6*N], N, magnitudesHigh,    NULL, N, true);
  fftMagnitudesAndPhases(impulseResponseSum,     N, magnitudesSum,     NULL, N, true);
  plotData(N/2, frequencies, magnitudesLow, magnitudesLowMid, magnitudesHighMid,
    magnitudesHigh, magnitudesSum);

  // OK - looks good. ToDo: plot in loglog scale
}

void rotes::testCrossoverNewVsOld()
{
  /*
  // just a function to check, if the new, templatized version of the crossover works as expected
  // by comparing it with the old implementation - when  the old implementation will be gone, this 
  // test code may be deleted, too

  double sampleRate        = 44100.0;
  double lowCrossoverFreq  = 250.0;
  double midCrossoverFreq  = 1000.0;
  double highCrossoverFreq = 6000.0;
  int    slope11           = 60;   // middlemost slope
  int    slope21           = 48;
  int    slope22           = 48;
  bool   active10          = true;
  bool   active11          = true;

  rosic::rsCrossOver4Way coOld;
  coOld.setSampleRate(44100);
  coOld.setSampleRate(sampleRate);
  coOld.setCrossoverFrequency(lowCrossoverFreq,  1, 0);
  coOld.setCrossoverFrequency(highCrossoverFreq, 1, 1);
  coOld.setCrossoverFrequency(midCrossoverFreq,  0, 0);
  coOld.setSlope(slope11, 0, 0);
  coOld.setSlope(slope21, 1, 0);
  coOld.setSlope(slope22, 1, 1);
  coOld.setBandActive(active10, 1, 0);
  coOld.setBandActive(active11, 1, 1);
  coOld.setMonoMode(false);

  rosic::rsCrossOver4WayStereo coNew;
  coNew.setSampleRate(44100);
  coNew.setSampleRate(sampleRate);
  coNew.setCrossoverFrequency(lowCrossoverFreq,  1, 0);
  coNew.setCrossoverFrequency(highCrossoverFreq, 1, 1);
  coNew.setCrossoverFrequency(midCrossoverFreq,  0, 0);
  coNew.setSlope(slope11, 0, 0);
  coNew.setSlope(slope21, 1, 0);
  coNew.setSlope(slope22, 1, 1);
  coNew.setBandActive(active10, 1, 0);
  coNew.setBandActive(active11, 1, 1);


  RAPT::rsNoiseGenerator<double> noiseGen;
  static const int N = 10000;  // number of samples
  double yo[8];
  double yn[8];
  double tol = 1.e-12;
  for(int n = 0; n < N; n++)
  {
    yo[0] = yn[0] = noiseGen.getSample();
    yo[1] = yn[1] = noiseGen.getSample();

    coOld.processSampleFrame(yo);
    coNew.processSampleFrameStereo(yn);
    rsAssert(RAPT::rsArrayTools::almostEqual(yo, yn, 8, tol));

    int dummy = 0;
  }
  */
}

void rotes::testSlopeFilter()
{
  double sampleRate = 44100.0;

  rosic::SlopeFilter filter;
  filter.setSampleRate(sampleRate);

  using Vec = std::vector<double>;

  // Settings:
  int N = 501;
  int numStages = 1;  // simulates using multiple stages
  Vec slopes = { -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 };
  //Vec slopes = { -18, -15, -12, -9, -6, -3, 0, +3, 6, 9, 12, 15, 18 };
  //Vec slopes = { -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50 };

  // Create and plot magnidtue responses:
  Vec freqs(N), mags(N), decibels(N);
  RAPT::rsArrayTools::fillWithRangeExponential(&freqs[0], N, 20.0, 20000.0);
  GNUPlotter plt2;
  for(size_t i = 0; i < slopes.size(); i++)
  {
    filter.setSlope(slopes[i]);
    for(int k = 0; k < N; k++)
    {
      mags[k] = filter.getMagnitudeAt(freqs[k]);
      mags[k] = pow(mags[k], numStages);
      decibels[k] = RAPT::rsAmpToDb(mags[k]);
    }
    plt2.addDataArrays(N, &freqs[0], &decibels[0]);
  }
  plt2.setLogScale("x");
  plt2.plot();
  int dummy = 0;


  // Observations:
  // -For smaller settings between -6,..,+6, the approximation looks quite good, in particular
  //  between 100 Hz and 10 kHz. Beyond that, the slope levels off (gets shallower) which might 
  //  actually be desirable - we may not want the boost to continue indefinitely into subsonic 
  //  and supersonic frequencies.
  // -The pivot or hinge frequency where all graphs meet is slightly above 1 kHz (around 1012 Hz). 
  //  I actually expected it to be exactly at 1 kHz. -> figure out, maybe fix
  // -For the more extreme settings like -40 or -50, the slopes seem to be smaller than they should
  //  be when just trying to achieve these slopes with one stage. Maybe the slope of one stage 
  //  should be limited to +-10 dB/oct and higher slopes should be realized by a cascade of 
  //  multiple slope filters (with identical coeffs).
}

void rotes::testPrototypeDesigner()
{
  /*
  static const int N = 3;
  double  coeffs[N+1];  
  Complex roots[N+5]; 
  fillWithZeros(coeffs, N+1);
  makeBesselPolynomial(coeffs, N);
  reverse(coeffs, N+1);
  findPolynomialRoots(coeffs, N, roots);
  */

  std::complex<double> poles[10];
  std::complex<double> zeros[10];
  rosic::rsPrototypeDesignerD designer;
  designer.setApproximationMethod(rosic::rsPrototypeDesignerD::BESSEL);
  designer.getPolesAndZeros(poles, zeros);

  int dummy = 0;
}

void rotes::testLowpassToLowshelf()
{
  // user parameters:
  static const int N = 5; // prototype filter order - above 7, it becomes weird for elliptics, Bessel works up to 25
  double Gp  = 0.8;       // passband ripple for elliptic lowpass prototype
  double Gs  = 0.2;       // stopband gain for elliptic lowpass prototype
  double G0  = 0.2;       // reference gain for shelving filter
  double G   = 1.3;       // gain for shelving filter

  // test:
  G0 = 0.2;
  G  = 1.4;
  //G0 = 1.4;
  //G  = 0.2;

  // analog unit cutoff lowpass prototype design:
  std::complex<double> z[N], p[N]; // arrays of lowpass prototype poles and zeros
  double  k;          // lowpass prototype filter gain
  //rsPrototypeDesigner::getBesselLowpassZerosPolesAndGain(z, p, &k, N);
  rosic::rsPrototypeDesignerD::getEllipticLowpassZerosPolesAndGain(z, p, &k, N, Gp, Gs);
  //Plotter::plotAnalogMagnitudeResponse(z, p, k, N, 0.0, 3.0, 1000);

  // transform lowpass to low-shelving:
  std::complex<double> zs[N], ps[N]; // arrays of lowshelf prototype poles and zeros
  double  ks;           // lowshelf prototype filter gain
  rosic::rsPoleZeroMapperD::sLowpassToLowshelf(z, p, &k, zs, ps, &ks, N, G0, G);

  //Plotter::plotAnalogMagnitudeResponse(zs, ps, ks, N, 0.0, 3.0, 2000);
    // needs to be updated for std::complex<double>
}


void rotes::testBesselPrototypeDesign()
{
  // this needs to be updated

  // user parameters:
  static const int N = 9;  //  Bessel works up to 25 before it gets funky
  double fc  = 1000.0;     // cutoff frequency
  double fl  =  500.0;     // lower bandedge frequency for bandpasses
  double fu  = 1500.0;     // upper bandedge frequency for bandpasses
  double G0  = 2.2;        // reference gain for shelving filter
  double G   = 1.4;        // gain for shelving filter

  // radian cutoff frequencies:
  double wc  = 2*PI*fc;   
  double wl  = 2*PI*fl;   
  double wu  = 2*PI*fu;   

  // todo: use pre-warped radian cutoff frequencies:
  double wcw = wc;        
  double wlw = wl;
  double wuw = wu;

  // just for test:
  wcw = 1000.0;           
  wlw =  500.0;
  wuw = 1500.0;

  // analog low-shelving design:
  std::complex<double> z[N], p[N]; // arrays of lowshelf prototype poles and zeros
  double  k;          // lowshelf prototype filter gain
  //rosic::rsPrototypeDesignerD::getBesselLowShelfZerosPolesAndGain(z, p, &k, N, G, G0);
  //Plotter::plotAnalogMagnitudeResponse(z, p, k, N, 0.0, 3.0, 2000);
  //Plotter::plotAnalogPhaseResponse(z, p, k, N, 0.0, 3.0, 2000);


  std::complex<double> zp[N], pp[N]; // arrays of lowpass prototype poles and zeros
  //double  kp;          // lowpass prototype filter gain
  //rosic::rsPrototypeDesignerD::getBesselLowpassZerosPolesAndGain(zp, pp, &kp, N);
  //Plotter::plotAnalogMagnitudeResponse(zp, pp, kp, N, 0.0, 3.0, 2000);
  //Plotter::plotAnalogPhaseResponse(zp, pp, kp, N, 0.0, 3.0, 2000);


  // s-plane frequency transformations:  

  std::complex<double> zLP[N], pLP[N], zHP[N], pHP[N], zBP[2*N], pBP[2*N], zBR[2*N], pBR[2*N];
  double  kLP, kHP, kBP, kBR; 

  // LP -> LP:
  rosic::rsPoleZeroMapperD::sLowpassToLowpass(z, p, &k, zLP, pLP, &kLP, N, wcw);
  //Plotter::plotAnalogMagnitudeResponse(zLP, pLP, kLP, N, 0.0, 3.0*wcw, 200);

  // LP -> HP:
  rosic::rsPoleZeroMapperD::sLowpassToHighpass(z, p, &k, zHP, pHP, &kHP, N, wcw);
  //Plotter::plotAnalogMagnitudeResponse(zHP, pHP, kHP, N, 0.0, 3.0*wcw, 200);
  //Plotter::plotAnalogPhaseResponse(zHP, pHP, kHP, N, 0.0, 3.0*wcw, 200);

  // LP -> BP:
  rosic::rsPoleZeroMapperD::sLowpassToBandpass(z, p, &k, zBP, pBP, &kBP, N, wlw, wuw);
  //Plotter::plotAnalogMagnitudeResponse(zBP, pBP, kBP, 2*N, 0.0, 3.0*wcw, 200);

  // LP -> BR:
  rosic::rsPoleZeroMapperD::sLowpassToBandreject(z, p, &k, zBR, pBR, &kBR, N, wlw, wuw);
  //Plotter::plotAnalogMagnitudeResponse(zBR, pBR, kBR, 2*N, 0.0, 3.0*wcw, 200);
}


void rotes::testPapoulisPrototypeDesign()
{
  // user parameters:
  static const int N = 6; // 
  double G0  = 0.2;        // reference gain for shelving filter
  double G   = 1.4;        // gain for shelving filter

  // with N = 12, G0 = 1, G >= 85 or so, the shelving numerator zeros come out wrong - there are repeated left-halfplane zeros
  // that shouldn't be there - ill numerical conditioning?
  G0 = 1.0;
  G  = 100.0;

  // analog low-shelving design:
  std::complex<double> z[N], p[N]; // arrays of lowshelf prototype poles and zeros
  //double  k;          // lowshelf prototype filter gain
  //rsPrototypeDesigner::getPapoulisLowpassZerosPolesAndGain(z, p, &k, N);
  //rosic::rsPrototypeDesignerD::getPapoulisLowShelfZerosPolesAndGain(z, p, &k, N, G, G0);
  //Plotter::plotAnalogMagnitudeResponse(z, p, k, N, 0.0, 3.0, 2000);
}

void rotes::testEngineersFilter()
{
  // filter parameters:
  double fs  = 44100.0;       // sample-rate
  double fc  = 500.0;         // characteristic frequency in Hz
  double bw  = 2.0;           // bandwidth in octaves (for bandpass, bandreject and bandstop)
  double g   = 2.0;           // gain for (peak and shelving filters)
  double rp  = 1.0;           // passband ripple in dB
  double rs  = 60.0;          // stopband rejection in dB
  int order  = 7;             // order of the prototype
  int mode   = 1;             // LP, HP, BP, BR, LS, HS, PK - todo: make allpass mode available
  int method = 5;             // approximation method

  // plotting parameters:
  static const int numSamples = 1000;
    
  // allocate and initialize memory for the responses:
  double timeAxis[numSamples];  // time axis in samples
  double stepResp[numSamples];  // step response
  double impResp[numSamples];   // impulse response
  RAPT::rsArrayTools::fillWithRangeLinear(timeAxis, numSamples, 0.0, numSamples-1.0);

  // create and set up the filter:
  rosic::rsEngineersFilterMono filter;
  filter.setSampleRate(fs);
  filter.setFrequency(fc);
  filter.setBandwidth(bw);
  filter.setGain(g);
  filter.setRipple(rp);
  filter.setStopbandRejection(rs);
  filter.setPrototypeOrder(order);  
  filter.setMode(mode);
  filter.setApproximationMethod(method);

  // obtain and plot impulse- and step-response:
  int n;          // sample counter
  filter.reset(); // necesarry? -> shouldn't
  impResp[0] = filter.getSample(1.0);
  for(n = 1; n < numSamples; n++)
    impResp[n] = filter.getSample(0.0);
  filter.reset();
  for(n = 0; n < numSamples; n++)
    stepResp[n] = filter.getSample(1.0);
  plotData(numSamples, timeAxis, impResp, stepResp);

  // todo: maybe plot magnitude, phase, phase-delay, group-delay, poles/zeros, etc. also
}

void rotes::testPoleZeroMapping()
{
  static const int N  = 10;      // prototype filter order
  double           fc = 5000.0;  // cutoff frequency in Hz
  double           fs = 44100.0; // samplerate in Hz

  //Complex z[2*N], p[2*N];        // arrays of poles and zeros
  std::complex<double> z[2*N], p[2*N];        // arrays of poles and zeros
  double  k = 1.0;               // filter gain


  rosic::rsInfiniteImpulseResponseDesignerD iirDesigner;
  iirDesigner.setApproximationMethod(rosic::rsPrototypeDesignerD::ELLIPTIC);
  iirDesigner.setMode(rosic::rsInfiniteImpulseResponseDesignerD::LOWPASS);
  iirDesigner.setPrototypeOrder(N);
  iirDesigner.setSampleRate(fs);
  iirDesigner.setFrequency(fc);
  iirDesigner.setRipple(3.0);
  iirDesigner.setStopbandRejection(60.0);
  iirDesigner.getPolesAndZeros(p, z);

  //double c rosic::PoleZeroMapper::getAllpassWarpCoefficient(

  int dummy = 0;
}

void rotes::highOrderFilterPolesAndZeros()
{
  static const int N = 6;    // prototype filter order
  double fs    = 44100.0;  // samplerate
  double fc    =  1000.0;  // cutoff frequency
  double Ap    =     1.0;  // passband ripple in dB
  double As    =    50.0;  // stopband rejection in dB
  double bw    =     1.0;  // bandwidth in octaves
  int mode   = rosic::rsInfiniteImpulseResponseDesignerD::BANDPASS;
  int method = rosic::rsPrototypeDesignerD::ELLIPTIC;

  // create and set up the filter designer object:
  rosic::rsInfiniteImpulseResponseDesignerD designer;
  designer.setPrototypeOrder(N);
  designer.setSampleRate(fs);
  designer.setFrequency(fc);
  designer.setRipple(Ap);
  designer.setStopbandRejection(As);
  designer.setBandwidth(bw);
  designer.setMode(mode);
  designer.setApproximationMethod(method);

  // compute poles and zeros:
  //Complex z[2*N], p[2*N];        // arrays of poles and zeros
  std::complex<double> z[2*N], p[2*N];        // arrays of poles and zeros
  designer.getPolesAndZeros(p, z);

  int dummy = 0;
}
