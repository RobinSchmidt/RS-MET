#include "PhaseVocoderExperiments.h"
// rename to SineModelExperiments

void phaseRepresentation()
{
  // We investigate the importance of representing the unwrapped phase not merely by a double
  // variable but by a double for the argument and an integer k that tells in which cycle we are.
  // The unwrapped phase is pu = p + 2*pi*k for some integer k where p is the actual argument (for 
  // the sine function) and k gives us the cycle in which we are. Using just a double variable for 
  // pu would lead to precision loss if k is large, so it is advidable to represent the phase by 
  // the pair (p,k) instead of just the single value pu. In this experiment we create a sinusoid 
  // once using the (p,k) representation for the unwrapped phase (which means we give just p as 
  // argument to the sine function) and using pu as argument for the sine function. As the effect 
  // becomes more pronounce for high k, we start at some start value k0.

  // user parameters:
  static const int N = 1000;  // number of samples to generate
  double fs  = 96000;         // samplerate in Hz
  double f   = 20000;         // frequency of the sinusoid
  double t0  = 1000000;       // start time in samples

  // internal parameters:
  double w  = 2*PI*f/fs;                   // normalized radian frequency
  double pu = t0*w;                        // unwrapped start phase
  double p  = fmod(pu, 2*PI);              // start phase argument
  int    k  = rsRoundToInt((pu-p)/(2*PI)); // start cycle index

  // generate sinusoids, once using the wrapped and once the unwrapped phase:
  double x1[N], x2[N];
  for(int n = 0; n < N; n++)
  {
    x1[n] = sin(pu);
    x2[n] = sin(p);
    pu += w;
    p  += w;
    if(p > 2*PI)
    {
      p -= 2*PI;
      k++;
    }
  }
    
  // obtain the difference between the two signals - it represents the error:
  double d[N];
  RAPT::rsArray::subtract(x2, x1, d, N);
  double e = RAPT::rsArray::maxAbs(d, N);
  int dummy = 0;

  // Observations:
  // With a samplerate of 96kHz, a sinusoid of 20kHz and a start-time (in samples) of a million
  // (i.e. around 10 seconds into the signal), the error is of the order of 1.e-7 and thus, 
  // significantly larger than the machine epsilon. Moreover, we didn't even consider the actual
  // phase error but rather the error in the resynthesized signal (i guess, the actual phase error
  // might be larger than that). So, it seems indeed a good idea to represent phase with a wrapped
  // argument and an integer k that gives the cycle.
}


/** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero value in
w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again. That means, 
the window has a period length of N. Such a window is suitable for applications where it is 
important that suitably overlapped windows sum up to a constant, like in the phase-vocoder. */
/*
void rsHanningWindowZN(double *w, int N)
{
  double s = 2*PI/N; // for a window that ends with zero: w[N-1]=0, this would be s=2*PI/(N-1)
  for(int n = 0; n < N; n++)
    w[n] = 0.5*(1-cos(s*n));
}
*/

void grainRoundTrip()
{
  // Check, if doing an STFT and an inverse STFT of a grain leads to an ouput grain that euqals
  // the input grain times the product of the analysis- and synthesis window (later, this may be
  // moved into a unit test)

  static const int N = 16;       // number of samples in x
  static const int B = 16;       // blocksize
  static const int K = B/2;      // number of bins (= FFTSize/2)
  static const int M = 2*K;      // FFT size
  int n0 = B/2;                  // time index of STFT frame

  rsSpectrogramD pv;             // for conveniently calling the static functions
                                       
  // create analysis and synthesis windows:
  double wa[B], ws[B];
  rsWindowFunction::hanningZN(wa, B);
  rsWindowFunction::hanningZN(ws, B);
    // try other windows - actually, the roundtrip work even with random sequences

  // create the test signal:
  double x[N];
  RAPT::rsArray::fillWithValue(x, N, 1.0);

  // obtain short-time spectrum:
  rsComplexDbl X[M];
  pv.setBlockSize(B);
  //pv.setTrafoSize(B);
  pv.setAnalysisWindowType(RAPT::rsWindowFunction::HANNING_WINDOW_ZN);
  //pv.setZeroPaddingFactor(1);
  pv.shortTimeSpectrum(x, N, n0, X);

  // under construction - not yet complete

  int dummy = 0;
}


//-------------------------------------------------------------------------------------------------
// Dolph-Chebychev window generation code from here:
// http://practicalcryptography.com/miscellaneous/machine-learning/implementing-dolph-chebyshev-window/
// not recommended for production use because the complexity is O(N^2) - instead use an iFFT 
// approach
// References:
// [1] Lyons, R., "Understanding Digital Signal Processing", Prentice Hall, 2004.
// [2] Antoniou, A., "Digital Filters", McGraw-Hill, 2000.

// This function computes the chebyshev polyomial T_n(x)
double cheby_poly(int n, double x){
  double res;
  if (fabs(x) <= 1) res = cos(n*acos(x));
  else              res = cosh(n*acosh(x));
  return res;
}
// calculate a chebyshev window of size N, store coeffs in out as in Antoniou
//  -out should be array of size N 
//  -atten is the required sidelobe attenuation (e.g. if you want -60dB atten, use '60')
void cheby_win(double *out, int N, double atten){
  int nn, i;
  double M, n, sum = 0, max=0;
  double tg = pow(10,atten/20);  /* 1/r term [2], 10^gamma [2] */
  double x0 = cosh((1.0/(N-1))*acosh(tg));
  M = (N-1)/2;
  if(N%2==0) M = M + 0.5; /* handle even length windows */
  for(nn=0; nn<(N/2+1); nn++){
    n = nn-M;
    sum = 0;
    for(i=1; i<=M; i++){
      sum += cheby_poly(N-1,x0*cos(PI*i/N))*cos(2.0*n*PI*i/N);
    }
    out[nn] = tg + 2*sum;
    out[N-nn-1] = out[nn];
    if(out[nn]>max)max=out[nn];
  }
  for(nn=0; nn<N; nn++) out[nn] /= max; /* normalise everything */
  return;
}
//-------------------------------------------------------------------------------------------------

void plotWindows()
{
  // Plots the overlapping windows an their sum

  static const int B  = 512;            // blocksize
  static const int H  = B/4;          // hopsize
  static const int N  = 1000;           // number of samples in the test signal


  int J = rsSpectrogramD::getNumFrames(N, H);
                                      
  // create the window function:
  double wa[B], ws[B], w[B];
  rsWindowFunction::hanningZN(wa, B);
  rsWindowFunction::hanningZN(ws, B);
  RAPT::rsArray::multiply(wa, ws, w, B);

    // todo: try different window functions: Hmaming, Blackman, versions with both ends nonzero

  double yw[N];             // sum of windows
  RAPT::rsArray::fillWithZeros(yw, N);
  double t[B];              // time indices for current window

  GNUPlotter plt;
  plt.setRange(-B/2, J*H+B/2, -0.1, 1.1);
  plt.addCommand("set xtics " + to_string(H)); 
  int j;
  for(j = 0; j < J; j++)
    plt.setGraphColor(j+1, "505050");
  plt.setGraphColor(j+1, "000000");
  int n = 0;
  for(j = 0; j < J; j++)
  {
    RAPT::rsArray::fillWithRangeLinear(t, B, n-B/2., n+B/2.-1);
    plt.addDataArrays(B, t, w);
    RAPT::rsArray::addInto(yw, N, w, B, n-B/2);
    n += H;
  }
  double s = rsSpectrogramD::getWindowSum(wa, ws, B, H);
  RAPT::rsArray::scale(yw, N, 1/s);
  plt.addDataArrays(N, yw);


  //int n0; 

  plt.plot();

  // ToDo: extend the analysis further to the left and right such that it starts and ends with the
  // first and last frame that would yield a nonzero spectrogram result. this implies that we get 
  // some frames that would formally have negative time indices - we need to take that into account
  // for plotting. But it makes the spectrogram processing easier - we may assume zero-extension to
  // left and right and periodic extension to below and above (take care if the Nyquist freqeuncy
  // value is repeated or not - depends on even/oddness of B)

  // B = 512, H = B/4, N = 1000 is good for a plot for the manual

  // maybe we can derive a formula for a synthesis window given an analysis window, blocksize and
  // hopsize, such that sum(wa*ws) = const ...but i think, that is exactly what the demodulation
  // function already does
}

void spectrogramSine()
{
  static const int B  = 512;            // blocksize
  static const int H  = B/4;            // hopsize
  static const int N  = 10000;          // number of samples in the test signal
  static const int P  = 1;              // zero-padding factor
  static const int M  = B*P;            // FFT size
  static const int K  = M/2 + 1;        // number of non-redundant bins
  double           fs = 44100;          // samplerate
  double           f  = 5000;           // sinusoid frequency
  int W = RAPT::rsWindowFunction::HANNING_WINDOW_ZN;


  // A hopsize of B/4 will result in a constant when overlapping successive frames, assuming that
  // the window is applied twice (once in the analysis stage and once in the synthesis stage). This
  // is the desired condition for perfect resynthesis without modulation artifcats - i.e. no 
  // explicit demodulation will be necessarry.

  // create the test signal:
  double x[N];
  createSineWave(x, N, f, 1.0, fs);
  RAPT::rsArray::scale(x, N/2, 0.1);   // amplitude switch in the middle of the signal

  // compute the complex spectrogram:
  rsSpectrogramD sp;
  sp.setBlockSize(B);
  sp.setHopSize(H);
  sp.setAnalysisWindowType(W);
  sp.setSynthesisWindowType(W); 
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(x, N);
  int F = s.getNumRows();

  // compute (magnitude) spectrogram and phasogram:
  double **mag, **phs, **dB;
  MatrixTools::rsAllocateMatrix(mag, F, K);
  MatrixTools::rsAllocateMatrix(phs, F, K);
  MatrixTools::rsAllocateMatrix(dB,  F, K);

  int i, j;
  for(i = 0; i < F; i++)
  {
    for(j = 0; j < K; j++)
    {
      mag[i][j] = abs(s(i, j));
      dB[i][j]  = rsMax(rsAmp2dB(mag[i][j]), -50.0);
      //dB[i][j]  = rsMax(rsAmp2dB(mag[i][j]), -250.0);
      phs[i][j] = arg(s(i, j));
    }
  }

  // compute frequency- and time-reassignment matrices:
  // ...

  // plot the magnitude spectrogram (later: with or without reassignment):
  plotSpectrogram(F, K, dB, fs, H);

  // resynthesize and plot signal:
  std::vector<double> y  = sp.synthesize(s);
  plotVector(y);

  // free dynamically memory:
  MatrixTools::rsDeAllocateMatrix(mag, F, K);
  MatrixTools::rsDeAllocateMatrix(phs, F, K);
  MatrixTools::rsDeAllocateMatrix(dB,  F, K);

  // create error signal:
  double err[N];
  for(int n = 0; n < N; n++)
    err[n] = x[n] - y[n];

  // todo: experiment with non-power-of-2 blocksizes, maybe also use odd blocksizes, etc.
  // turn into unit test - use noisy input signals
  int dummy = 0;
}

void sineParameterEstimation()
{
  // A testbed for experimenting with different block sizes, transform sizes and window functions
  // to investigate, how these choices affect the accuracy of the estimation of sinusoidal 
  // parameters. We feed a single cosine wave with known frequency, amplitude and phase into the 
  // analysis and visualize the results.

  typedef RAPT::rsWindowFunction WF;
  typedef SinusoidalAnalyzer<double> SA;

  // signal parameters:
  double sampleRate = 10000;    // sample rate
  double frequency  =  1000;    // signal frequency
  double amplitude  =  1;       // signal amplitude
  double startPhase =  0;       // start phase (of the cosine)
  double length     =  0.1;     // length in seconds

  // analysis parameters:
  double anaTime = length/2; // time instant of analysis (exact only up to sample-rate)
  int blockSize  = 500;
  int trafoSize  = 500;
  //int window     = WF::RECTANGULAR_WINDOW;
  //int window     = WF::TRIANGULAR_WINDOW;
  int window     = WF::HAMMING_WINDOW;
  //int window     = WF::HANNING_WINDOW_ZN;
  //int window     = WF::BLACKMAN_WINDOW;
  //int window     = WF::BLACKMAN_HARRIS;

  // tests:                             // errors (with Hamming window, 1 kHz @ 10 kHz):

  //blockSize = 50; trafoSize =  99;    // f: -1.07, a: -0.003, p: -0.003    ok
  //blockSize = 50; trafoSize = 100;    // f: -1.22, a: -0.001, p: -9.9e-5   good
  //blockSize = 50; trafoSize = 101;    // f: -1.38, a: -0.001, p: 0.003     ok

  //blockSize = 49; trafoSize =  97;    // f: -1.08, a: 0.013,  p: 1.5e-14   good
  blockSize = 49; trafoSize =  98;    // f: -1.08, a: 0.01,   p: -0.641    bad
  //blockSize = 49; trafoSize =  99;    // f: -1.19, a: 0.021,  p: 1.7e-14   good
  //blockSize = 49; trafoSize = 100;    // f: -1.34, a: 0.023,  p: -0.628    bad    p = -2*PI/10

  //blockSize = 500; trafoSize = 500;   // f: -7.7e-8, a: -1.1e-6, p: -9.3e-8
  //blockSize = 250; trafoSize = 500;   // f: -0.049,  a: -1.9e-5, p: -7.4e-7
  //blockSize = 249; trafoSize = 500;   // f: -0.050,  a: 0.005,   p: -0.628    large phase error (actually -2*PI/10)
  //blockSize = 500; trafoSize = 4000;  // highly oversampled spectrum

  // generate cosine wave:
  int numSamples = (int) ceil(length*sampleRate); // number of samples
  double w = 2*PI*frequency/sampleRate;
  std::vector<double> x(numSamples);
  for(int n = 0; n < numSamples; n++)
    x[n] = amplitude * cos(w*n + startPhase);
  //plotVector(x);

  // obtain short time spectrum at the desired time-instant:
  int anaIndex = (int) round(anaTime*sampleRate); // sample index at which to center the analysis window
  rsSpectrogramD sp;
  sp.setAnalysisWindowType(window);
  sp.setTimeOriginAtWindowCenter(true);
  sp.setBlockAndTrafoSize(blockSize, trafoSize);
  std::vector<std::complex<double>> trafoBuffer(trafoSize), complexSpectrum(trafoSize);
  sp.prepareTrafoBuffer(&x[0], numSamples, anaIndex, &trafoBuffer[0]);
  sp.shortTimeSpectrum( &x[0], numSamples, anaIndex, &complexSpectrum[0]);
  std::vector<double> magnitudes(trafoSize), decibels(trafoSize), phases(trafoSize);
  double scaler = sp.getAnalysisScaler();
  for(int k = 0; k < trafoSize; k++) {
    magnitudes[k] = scaler * abs(complexSpectrum[k]);
    decibels[k]   = rsMax(rsAmp2dB(magnitudes[k]), -100.0);
    phases[k]     = arg(complexSpectrum[k]);
  }
  //plotComplexVectorReIm(trafoBuffer);
  //plotVector(magnitudes);
  //plotVector(decibels);
  //plotVector(phases);

  // estimate the sine parameters:
  std::vector<int> peaks = SA::peakIndices(&magnitudes[0], trafoSize/2, 0.1);
  int peakIndex = peaks[0]; // we should find one and only one peak, if everything works as it should
  double peakBin, freqEstimate, ampEstimate, phaseEstimate; 
  SA::spectralMaximumPositionAndValue(&magnitudes[0], peakIndex, &peakBin, &ampEstimate);
  freqEstimate  = peakBin*sampleRate/sp.getFftSize(); // maybe we should have a function sp.getBinFrequency(peakBin)
  phaseEstimate = phases[peakIndex]; // maybe use interpolation later

  // compute errors with respect to true values:
  double targetPhase = startPhase + 2*PI*frequency*anaTime; // actual phase of cosine at anaTime
  targetPhase = RAPT::rsWrapToInterval(targetPhase, -PI, PI);
  double freqError  = frequency - freqEstimate;
  double levelError = rsAmp2dB(amplitude / ampEstimate);
  double phaseError = targetPhase - phaseEstimate;

  // Observations:
  // for blockSize = 249; trafoSize = 500; we get a phase error of exactly -2*PI/10 - this is not a
  // random error - something must be systematically wrong

  // todo: 
  // -try different blockSizes and trafoSizes - in particular, try even/odd numbers for both, etc.
  //  -check the phase-error for odd blocksize and trafoSize = 2*blockSize
  // -plot the analyzed short-time spectrum, a highly oversampled version of it (->zero padding),
  //  the fitted parabola and a mark at the found frequency/amplitude (and maybe another mark at
  //  the correct freq/amp)
  // -maybe somehow also visualize the phase estimation (not yet sure, how)
  GNUPlotter plt;
}

// convenience function:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate)
{
  SinusoidalSynthesizer<double> synth;
  synth.setSampleRate(sampleRate);
  return synth.synthesize(model);
}

// rename to testSinusoidalSynthesis1
void sinusoidalSynthesis1()
{
  double stretch = 1.0;
  double fs = 44100;  // use for writing to wavefile
  //double fs = 5000; // use for plotting the signal

  // create the datapoints for a singole partial:
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  partial.appendDataPoint(ISP(stretch*0.0, 100.0, 0.4,  PI/2));    // time, freq, amp, phase
  partial.appendDataPoint(ISP(stretch*0.4, 100.0, 0.2,  PI/2));
  partial.appendDataPoint(ISP(stretch*0.8, 150.0, 0.8, -PI/2));
  partial.appendDataPoint(ISP(stretch*1.2, 100.0, 0.4,  0.0));
  partial.appendDataPoint(ISP(stretch*1.6, 200.0, 0.2,  PI));
  partial.appendDataPoint(ISP(stretch*2.0, 100.0, 0.8,  PI/2));

  // create a model, add the partial to it and synthesize the sound:
  RAPT::rsSinusoidalModel<double> model;
  model.addPartial(partial);
  SinusoidalSynthesizer<double> synth;
  synth.setSampleRate(fs);
  synth.setCubicAmplitudeInterpolation(true);
  synth.setCubicPhaseInterpolation(true);
  std::vector<double> x = synth.synthesize(model);

  // make a sinusoidal analysis of the sound that we have just created and re-create the sound
  // from the model that results from this analysis:
  // maybe try the following parameters: 
  // maxFreqDeviation df = 100
  // windowType = Hamming -> B = 4, L = -42.7
  // windowSize M >= B * fs / df = 4 * 44100 / 100 = 1764 -> use 1801 for odd size
  // threshold t >= L = -42.7 -> use t = -30
  // i think, the Blackman-Nutall window is superior to the Blackman-Harris window - it has a
  // slightly narrower mainlobe and better sidelobe rejection
  RAPT::rsSinusoidalModel<double> model2;
  SinusoidalAnalyzer<double> sa;
  sa.setWindowType(RAPT::rsWindowFunction::HAMMING_WINDOW);
  sa.setMaxFreqDeltaBase(100);
  sa.setTrafoSize(4096);
  //sa.setBlockSize(256);    // gives total nonsense results
  sa.setBlockSize(512);      // can't really follow the freq-trajectory (does some averaging)
  //sa.setBlockSize(1024);
  //sa.setBlockSize(2048);
  //sa.setBlockSize(1801);  // does not yet work
  //sa.setHopSize(225);
  //sa.setHopSize(256);
  sa.setHopSize(128);
  sa.setRelativeLevelThreshold(-15);
  //sa.setRelativeLevelThreshold(-40);
  model2 = sa.analyze(&x[0], (int)x.size(), fs);            // try to recover model from sound
  //std::vector<double> y = synthesizeSinusoidal(model2, fs); // resynthesize sound from recovered model

  //plotSinusoidalAnalysisResult(sa, &x[0], (int) x.size(), fs); // should plot spectrogram under tracks - but that doesn't work yet
  //plotSineModel(model, fs);
  plotTwoSineModels(model, model2, fs);

  // in some way, the plot of the original model is actually wrong - it does not really represent 
  // what is synthesized because the synthesizer uses spline interpolation for the frequency tracks 
  // and  the plot uses linear interpolation - we should perhaps plot a sort of spline-interpolated
  // model as reference (maybe have a function upsampleSinusoidalModel or something)

  //rosic::writeToMonoWaveFile("SinusoidalSynthesisTest.wav", &x[0], (int)x.size(), (int)fs, 16);

  // Observations:
  // -when the sound is synthesized with linear phase interpolation, the analyzed model shows 
  //  sections of constant frequency (almost, the frequency jitters a bit around the constant) 
}

void sinusoidalSynthesis2()
{
  // a sinusoid with a stable frequency - we want to test fade-in/out
  double fs = 2000; // sample rate

  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  RAPT::rsSinusoidalModel<double> model, model2;

  partial.prependDataPoint(ISP(-0.1, 100.0, 0.5, 0.0));
  partial.appendDataPoint( ISP( 0.0, 100.0, 1.0, 0.0));
  partial.appendDataPoint( ISP( 0.1, 100.0, 0.5, 0.0));
  partial.appendDataPoint( ISP( 0.2, 100.0, 1.0, 0.0));

  model.addPartial(partial);
  std::vector<double> x = synthesizeSinusoidal(model, fs);

  plotVector(x);
}
  // maybe try a single linear sweep


// (move to test tools):
std::vector<double> createSinusoid(size_t N, double frequency, double sampleRate,
  double amplitude = 1.0, double phase = 0.0)
{
  std::vector<double> x(N);
  double w = 2*PI*frequency/sampleRate;
  for(size_t n = 0; n < x.size(); n++)
    x[n] = amplitude * sin(w*n + phase);
  return x;
}
// adds a sinusoid with given parameters into the vector:
void addSinusoid(std::vector<double>& x, double frequency, double sampleRate,
  double amplitude = 1.0, double phase = 0.0)
{
  double w = 2*PI*frequency/sampleRate;
  for(size_t n = 0; n < x.size(); n++)
    x[n] += amplitude * sin(w*n + phase);
}

void sinusoidalAnalysis1()
{
  // test signal parameters:
  //double sampleRate = 48000;  // in Hz
  double sampleRate = 50000;  // in Hz
  double length     = 0.1;    // in seconds -> 10000 samples at 50kHz
  double frequency  = 1000;   // in Hz
  double amplitude  = 3.0;    // raw factor
  double phase      = 0.0;    // in radians

  // analsis parameters:
  int window = RAPT::rsWindowFunction::HAMMING_WINDOW;
  double freqRes = 100; // frequency resolution

  phase = PI/2; // for test
  //frequency = 800;
  //frequency = 808;
  //frequency = 809;
  //frequency = 810;
  frequency = 1000 * GOLDEN_RATIO/2;
  //frequency = 1000 * SQRT2_INV;

  // create signal:
  double period = sampleRate / frequency;         // in samples
  size_t N = (size_t)ceil(length * sampleRate);   // number of samples
  std::vector<double> x = createSinusoid(N, frequency, sampleRate, amplitude, phase);

  // create and set up analyzer:
  SinusoidalAnalyzer<double> sa;
  sa.setWindowType(window);
  sa.setMaxFreqDeltaBase(100);
  sa.setTrafoSize(4096);
  sa.setBlockSize(1024);
  sa.setHopSize(256);
  //sa.setHopSize(512);
  //sa.setRelativeLevelThreshold(-25);  // some spurious tracks will occur (with Hamming window)
  sa.setRelativeLevelThreshold(-15);    // no spurious tracks will occur (with Hamming window)
  sa.setFadeInTime(0.01);    // -> 500 samples @50kHz
  sa.setFadeOutTime(0.01);   // -> 500 samples @50kHz
  sa.setMinimumTrackLength(0.021);  // should be a little above fadeInTime+fadeOutTime
  //sa.setFadeInTime(0.0);
  //sa.setFadeOutTime(0.0);
  //sa.setMinimumTrackLength(0.001); 
  //plotSineModel(sa, &x[0], (int) x.size(), sampleRate);

  int minBlockSize = sa.getRequiredBlockSize(window, freqRes, sampleRate);
  double maxLevelThresh = sa.getRequiredThreshold(window, 0.0);

  // find model for the signal:
  rsSinusoidalModel<double> model = sa.analyze(&x[0], (int)N, sampleRate);

  // ok - aside from spurious tracks at the start/end (transients?) it looks good -> clean up by
  // deleting spurious tracks and "finalize" tracks by applying fade-outs
  // ->figure out, why we get spurious tracks, even at a threshold of -25 dB which is clearly above
  // the sidelobe level of -42 dB of the Hamming window - hmm - there are indeed sidelobes in the 
  // spectrum that are higher than that - why? is the window messed up? the window spectrum looks 
  // kinda strange anyway - plot the window and its spectrum
  // aaahhh - OK - it is because the first few frames have the switch-on - in the very first frame,
  // this is like imposing an additional rectangular window - zeroing out the left half of the
  // window - this creates higher sidelobes. well within the signal, the window spectrum looks as
  // expected - OK - all is good. these are just transient artifacts and we should clean them up
  // by deleting the spurious tracks

  // create and set up a sinusoidal synthesizer object and plot resynthesis result:
  SinusoidalSynthesizer<double> synth;
  synth.setSampleRate(sampleRate);
  plotSineResynthesisResult(model, synth, &x[0], (int) N);

  // there are discontinuities in the fade-in/out sections because the original signal has the hard
  // switch on/off discontinuities in it which the resynthesized signal doesn't have - so the 
  // discontinuities are to be expected

  // When using zero fade-times, the model is too short - i think this is because we use a small
  // hopSize and the spectrogram assumes a hopSize of blocksize/2 to compute the required number of
  // frames - check this - hmm...but it happens also when the hopSize actually is blocksize/2
  // ...it seems that generally, out models are one frame (or more?) to short? -> experiment with
  // this - it should all work also when there is no fade-in/out whatsoever
}

void sinusoidalAnalysis2()
{
  // Two sinusoids at 1000 and 1100 Hz

  // input signal parameters:
  double fs = 44100; // sample rate
  double f1 = 1000;  // 1st frequency
  double f2 = 1107;  // 2nd frequency
  double a1 = 1.0;   // 1st amplitude
  double a2 = 1.0;   // 2nd amplitude
  double p1 = PI/2;  // 1st start phase
  double p2 = PI/2;  // 2nd start phase
  double L  = 0.2;   // length in seconds

  // analysis parameters:
  double freqRes  = f2-f1;   // frequency resolution
  double dBmargin = 20;      // dB margin over sidelobe level
  int padFactor = 2;         // zero padding factor
  int window    = RAPT::rsWindowFunction::HAMMING_WINDOW;


  // create a model and synthesize the sound:

  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial1, partial2;
  RAPT::rsSinusoidalModel<double> model;

  double p1e = p1 + 2*PI*f1*L; // end phase for 1st partial
  partial1.prependDataPoint(ISP(0.0, f1, a1, p1));
  partial1.appendDataPoint( ISP(L,   f1, a1, p1e));

  double p2e = p2 + 2*PI*f2*L; // end phase for 2nd partial
  partial2.prependDataPoint(ISP(0.0, f2, a2, p2));
  partial2.appendDataPoint( ISP(L,   f2, a2, p2e));

  model.addPartial(partial1);
  model.addPartial(partial2);

  std::vector<double> x = synthesizeSinusoidal(model, fs);
  //plotSineModel(model, fs); // we need a margin for the y-axis - otherwise it looks like nothing
  //plotVector(x);

  // analyze the produced sound again and compare the original model and analysis result
  // ...

  // create and set up analyzer and try to recover the model from the sound
  SinusoidalAnalyzer<double> sa;
  int blockSize = sa.getRequiredBlockSize(window, freqRes, fs); // +1 for test
  int hopSize   = blockSize/8; // small hopSize needed, otherwise analyzed model has too short tracks
  //int trafoSize = RAPT::rsNextPowerOfTwo(blockSize);  // or maybe use blockSize itself
  int trafoSize = padFactor*blockSize;
  double levelThresh = sa.getRequiredThreshold(window, dBmargin);
  sa.setWindowType(window);
  sa.setRelativeLevelThreshold(levelThresh);
  sa.setBlockAndTrafoSize(blockSize, trafoSize);
  //sa.setBlockSize(blockSize);
  //sa.setTrafoSize(trafoSize);
  sa.setHopSize(hopSize);
  sa.setMaxFreqDeltaBase(10);  // 10 Hz variation allowed between frames
  sa.setFadeInTime(0.01);
  sa.setFadeOutTime(0.01);
  sa.setMinimumTrackLength(0.021);  // should be a little above fadeInTime+fadeOutTime
  rsSinusoidalModel<double> model2 = sa.analyze(&x[0], (int)x.size(), fs);
  //plotTwoSineModels(model, model2, fs);
  // ...try to use the lowest possible values that give a good result


  // create and set up a sinusoidal synthesizer object and plot resynthesis result:
  SinusoidalSynthesizer<double> synth;
  synth.setSampleRate(fs);
  plotSineResynthesisResult(model2, synth, &x[0], (int)x.size());

  // Observations:
  // -when the hopSize is too small (like blockSize/2), the analyzed tracks are too short
  // -for 1000 and 1107 Hz, there's a largish amplitude error in the estimates (over 0.06) when
  //  using a padFactor of 1, using a larger padFactor, the amplitude error gets smaller but we get
  //  a larger phase error (so large, that the residual is actually louder with zero-padding)
  //  -could this be a bug that messes up the phase when zero padding is used? it seems to be a
  //   a systematic error
  //  -the computed blockSize is even (1649) - when used as is, the phase error appears, when 
  //   adding 1 to get 1650, the phase error disappears - there's clearly something wrong with 
  //   zero-padding when the blockSize is odd

  // maybe allow time varying frequencies and amplitudes - but maybe in next test - ramp up the
  // complexity
}


// make various tests for the sinusoidal analysis of increasing level of difficulty:
// 1: single sinusoid with stable frequency and amplitude
//  -check, if the frequency and ampltude is etsimated accurately
//  -try frequencies that coincide with bin-centers (best case) and those that fall halfway 
//   between bins (worts case) and some intermediate cases
// 2: single sinuosoid with time-variying amplitude
// 3: two sinusoids with stable freq and amp
//  -check, how each influences the analysis of the other - as function of their frequencies
// 4: several sinusoids with various frequencies and start and endpoints
// 5: periodic sounds like triangle, square, saw
// ....

// Notes (move to documentation of SinusoidalAnalyzer):
// when the frequency is systematically under- or overestimated (i.e. estimation errors do not 
// average out to zero during a sinusoidal track), we may still hope that the time-domain 
// resynthesis works well, if the phase-values can compensate that error. if that's not possible
// (it may be especially problematic at high frequencies) - try using a smaller hop size (halving 
// the hop-size should double the frequency at which freq-error compensation by phase still 
// works)

/*

Experiment with different window versions: symmetric (w[0]==w[N-1]), periodic (w[1]==w[N-1])
I think, the periodic one produces a little error in the computed phase (maybe also magnitude) but is
better suited for perfect reconstruction (overlapped windows add up to a constant, for properly chosen
shape and hopsize). Perhaps, the periodic version should be used in a blocked overlap/add FFT 
analysis/resynthesis scheme, but when oscillator bank resynthesis is used, symmetric windows might be 
preferable, due to less error in the analysis data (thes are all hypotheses that need to be checked
experimentally)

We could use symmetric windows (that start and end either with a zero or nonzero sample) and increase
or decrease the hopsize by one sample in order to achieve the desired periodicity. if both ends are 
nonzero, increase, if both are zero, decrease. this could solve the symmetry issue.



Experiment with the following test signals:

synthetic signals:
DC, impulses, sines of various frequencies, mixes of sine waves, sweeping sines (and mixes thereof, 
also simultaneous up- and downward sweeps that cross each other), noise, saw, saw sweeps, chords, 
sines + impulses (to observe time- and frequency localized signals simultaneously), damped sines,
filter sweeps

natural sounds (roughly ascending complexity): 
tonal instrument samples (steady and impulsive), drum/percussion samples, 
solo instrument performances, speaking and singing voice, drumloops, full (mastered) mixdowns of 
different musical genres

Ideas:
In order to match the analysis blocksize to the periodicity of the signal, the signal can be 
autotuned to a fixed period. For example, if the fundamental is around 440 Hz, corresponding to 
roughly 100 samples at 44.1 kHz samplerate, the period could be autotuned to 128 samples and the 
phase-vocoder would use a blocksize of some multiple of 128. Or we could fix the analysis 
blocksize to, say, 512 and tune the period to 512/3 such that each block contains 3 cycles.


Miscelleneous remarks:
-When we reynthesize a signal with the overlap/add procedure, we have to a apply a global gain 
 inverse to summed (overlapped) products of the anylysis and synthesi windows. If these sum up to a
 constant, i think, the gain can be calculated as:
 g = 0;
 for(i = 0; i < B; i++)
   g += wa[i] * ws[i];
 where  B: blocksize, H: hopsize, wa: analysis window, ws: synthesis window. The inverse of this 
 computed g value is the desired gain to be applied

Seperation of sinusoids: For a spectrogram value s(j, k) where j is the frame index and k is the 
bin index, a sinusoid will satisfy some constraints for how the s(j, k) relates to its local
neighbourhood. The horizontal constraints depend of the overlap factor (defined as (B-H)/B) and the
vertical constraints depend on the window-shape and zero-padding factor. These constraints tighten
as the overlap or zero-padding increases, respectively - so any decision function that separates
sinusoids should take these parameters into account. The vertical constraints also tighten with
increasing mainlobe-width of the window.
vertical constraints:
-not a minimum: s(j, k) >= s(j, k-1), s(j, k) >= s(j, k+1), but it's not necessarily a maximum 
 because a sinuosid spreads over several bins
-a sinusoid in the freqeuncy domain looks like the transform of the window, so maybe a 
 pattern-matching algorithm can be used -> cross-correlate STFT spectrum with the window transform, 
 maybe only with the mainlobe, but i think, it's not valid to just crosscorrelate the magnitudes
-2D-crosscorrelate similarity with image-processing kernels that correspond to sinusoids, width 
 depends on overlap, height on zero-padding
-compare forward and backward instantaneous frequency estimate (should be similar for sines)
 -maybe we can also use the reassigned frequency values (and maybe also their time reassignments)
-bandpass-filter horizontal lines with the expected modulation frequency for the respective bin 
 (real and imaginary parts are expected to undulate with a frequency characteristic for each bin)
-maybe to modify the spectrogram values, compute an expected value from the neighbouhood for 
 magnitude and phase and the do a weighted sum of actual and expected value (weights should add
 up to unity)
-for experimentation, we should at various cases. stable sinusoid, sweeping amplitude, sweeping 
 frequency, both is sweeping - we may derive
-maybe make use of the time/frequency reassignment values
-there are (at least) 3 ways to determine an instantaneous frequency: phase-differences of 
 successive frames, frequency-reassignment, parabolic magnitude interpolation - maybe the decision,
 whether a pixel (time/frequency point) belongs to a sinusoid can be based on how much these 
 different values agree.
 -the most efficient would be to compare magnitude-interpolation values to phase-difference values
  because no additional FFTs are needed
 -other ideas for finding a frequency: parabolic interpolation of autocorrelation sequence
-if a spectral peak is considered sinusoidal, we should automatically assume some neighbouring peaks
 to belong to that sinuosid, too - based on the spectral width of the FFT of the analysis window
-deviatograms: idea: if there's a stable sinusoid, we could compute an expected value for each 
 pixel based on neighborhood pixels - the (absolute or squared) difference between the expected and
 actual value is called the "deviation". plotting the deviation for each pixel gives the
 "deviatogram". for example, for a magnitude deviatogram, we could define the expected magnitude 
 value se(j, k) = (s(j-1, k) + s(j+1, k)) / 2. this is the mean of the left and right neighbour's
 magnitude and implies an assumed linear amplitude envelope. for an assumption of an exponential
 envelope, we could use the geometric mean. we could also assume that ascension are linear and 
 decays exponential by using geometric mean if s(j-1) > s(j+1), else use arithmetic mean (should
 be used only for impulsive sounds)


*/