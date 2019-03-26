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
  //static const int K = B/2;      // number of bins (= FFTSize/2)
  //static const int M = 2*K;      // FFT size
  static const int M = 16;       // FFT size
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
  pv.setTrafoSize(B);
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
/*
// This function computes the chebyshev polyomial T_n(x)
double cheby_poly(int n, double x)
{
  double res;
  if (fabs(x) <= 1) res = cos(n*acos(x));
  else              res = cosh(n*acosh(x));
  return res;
}

// calculate a chebyshev window of size N, store coeffs in out as in Antoniou
//  -out should be array of size N 
//  -atten is the required sidelobe attenuation (e.g. if you want -60dB atten, use '60')
void cheby_win(double *out, int N, double atten)
{
  int nn, i;
  double M, n, sum = 0, max=0;
  double tg = pow(10,atten/20);         // 1/r term [2], 10^gamma [2]
  double x0 = cosh((1.0/(N-1))*acosh(tg));
  M = (N-1)/2;
  if(N%2==0) M = M + 0.5;               // handle even length windows 
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
  for(nn=0; nn<N; nn++) out[nn] /= max; // normalise everything
  return;
}
*/
// code moved to Prototypes and can be deleted



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

//-------------------------------------------------------------------------------------------------

void spectrogramSine()
{
  static const int B  = 512;            // blocksize
  static const int H  = B/4;            // hopsize
  static const int N  = 10000;          // number of samples in the test signal
  static const int P  = 2;              // zero-padding factor
  static const int M  = B*P;            // FFT size
  static const int K  = M/2 + 1;        // number of non-redundant bins
  double           fs = 44100;          // samplerate
  double           f  = 5000;           // sinusoid frequency
  int W = RAPT::rsWindowFunction::HANNING_WINDOW_ZN;


  // A hopsize of B/4 will result in a constant when overlapping successive frames, assuming that
  // the window is applied twice (once in the analysis stage and once in the synthesis stage). This
  // is the desired condition for perfect resynthesis without modulation artifcats - i.e. no 
  // explicit demodulation will be necessarry. ..todo: currently, this constant is not unity but 
  // rather 1.5 - so without demodulation, the resynthesized signal will be too loud by factor 1.5
  // compared to the original (this does not happen with demodulation, because the demodulation
  // compensates for that factor as well)

  // Note that at the beginning and end of the signal, the demodulation signal is not constant 
  // everywhere - there are still fade-in and fade-out artifacts du to the fact that at start and 
  // end, we do not yet have the full number of overlapping windows. ..this fade-in and out is also 
  // compensated by the demodulation but perhaps that may lead to artifacts. As a workaround, it is
  // possible to prepend and append a block of zeros (length B should be enough, perhaps even more 
  // than enough) to the signal before analysis and remove these paddings after resynthesis. 

  // create the test signal:
  double x[N];
  createSineWave(x, N, f, 1.0, fs);
  RAPT::rsArray::scale(x, N/2, 0.1);   // amplitude switch in the middle of the signal

  // compute the complex spectrogram:
  rsSpectrogramD sp;
  sp.setBlockAndTrafoSize(B, M);
  sp.setHopSize(H);
  sp.setAnalysisWindowType(W);
  sp.setSynthesisWindowType(W); 
  //sp.setOutputDemodulation(false); // with appropriate settings, demodulation should be superfluous
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(x, N);
  int F = s.getNumRows();    // number of frames

  // compute (magnitude) spectrogram and phasogram:
  double **mag, **phs, **dB;
  MatrixTools::rsAllocateMatrix(mag, F, K);
  MatrixTools::rsAllocateMatrix(phs, F, K);
  MatrixTools::rsAllocateMatrix(dB,  F, K);
  for(int i = 0; i < F; i++) {
    for(int  j = 0; j < K; j++) {
      mag[i][j] = abs(s(i, j));
      dB[i][j]  = rsMax(rsAmp2dB(mag[i][j]), -50.0);
      phs[i][j] = arg(s(i, j));
    }
  }

  // compute frequency- and time-reassignment matrices:
  // ...

  // plot the magnitude spectrogram (later: with or without reassignment):
  //plotSpectrogram(F, K, dB, fs, H);

  // resynthesize and plot signal:
  std::vector<double> y  = sp.synthesize(s);
  //plotVector(y);

  // create error and plot signal:
  std::vector<double> err(N);
  for(int n = 0; n < N; n++)
    err[n] = x[n] - y[n];
  plotVector(err);

  // Observations:
  // -for B = 512, H = B/4, P = 1,2,4 the resynthesis error is of the order of 5.e-16 (~2*eps), 
  //  with P = 3 about 6.e-13 and with P = 5 about 1.e-12 - so when the zero padding factor is a 
  //  power of two, the numerical roundoff properties are better, but more general zero-padding 
  //  factors also work in principle

  // todo: experiment with non-power-of-2 blocksizes, maybe also use odd blocksizes, etc.
  // turn into unit test - use noisy input signals


  // free dynamically memory:
  MatrixTools::rsDeAllocateMatrix(mag, F, K);
  MatrixTools::rsDeAllocateMatrix(phs, F, K);
  MatrixTools::rsDeAllocateMatrix(dB,  F, K);
}

void sineParameterEstimation()
{
  // A testbed for experimenting with different block sizes, transform sizes and window functions
  // to investigate, how these choices affect the accuracy of the estimation of sinusoidal 
  // parameters. We feed a single cosine wave with known frequency, amplitude and phase into the 
  // analysis and visualize the results.

  typedef RAPT::rsWindowFunction WF;
  typedef rsSinusoidalAnalyzer<double> SA;

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
  //blockSize = 49; trafoSize =  97;    // f: -1.08, a: 0.013,  p: 1.5e-14   good
  //blockSize = 49; trafoSize =  98;    // f: -1.08, a: 0.018,  p: 1.2e-14   good
  //blockSize = 49; trafoSize =  99;    // f: -1.19, a: 0.021,  p: 1.7e-14   good
  //blockSize = 49; trafoSize = 100;    // f: -1.34, a: 0.023,  p: 1.6e-14   good 

  //blockSize = 50; trafoSize =  99;    // f: -1.07, a: -0.003, p: -0.003,-0.0005    ok
  //blockSize = 50; trafoSize = 100;    // f: -1.22, a: -0.001, p: -9.9e-5, -0.0005   ok
  blockSize = 50; trafoSize = 101;    // f: -1.38, a: -0.001, p: 0.003, -0.0004     ok

  //blockSize = 500; trafoSize = 500;   // f: -7.7e-8, a: -1.1e-6, p: -9.3e-8
  //blockSize = 250; trafoSize = 500;   // f: -0.049,  a: -1.9e-5, p: -7.4e-7
  //blockSize = 249; trafoSize = 500;   // f: -0.050,  a: 0.005,   p: -0.628    large phase error (actually -2*PI/10)
  //blockSize = 499; trafoSize = 4000;  // highly oversampled spectrum

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
  int peakIndex = peaks[0]; // peak index, we should find one and only one peak, if everything works as it should
  double peakPos, freqEstimate, ampEstimate, phaseEstimate; 
  SA::spectralMaximumPositionAndValue(&magnitudes[0], peakIndex, &peakPos, &ampEstimate);
  freqEstimate  = peakPos*sampleRate/sp.getFftSize(); // maybe we should have a function sp.getBinFrequency(peakBin)

  //phaseEstimate = phases[peakIndex]; // maybe use interpolation later
  phaseEstimate = RAPT::rsArray::interpolatedValueAt(&phases[0], (int) phases.size(), peakPos);

  // compute errors with respect to true values:
  double targetPhase = startPhase + 2*PI*frequency*anaTime; // actual phase of cosine at anaTime
  targetPhase = RAPT::rsWrapToInterval(targetPhase, -PI, PI);
  double freqError  = frequency - freqEstimate;
  double levelError = rsAmp2dB(amplitude / ampEstimate);
  double phaseError = targetPhase - phaseEstimate;

  // Observations:
  // -when the blockSize is odd, we can estimate the phase much more accurately which is consistent
  //  with what Xavier Serra says in his lectures - interesingly, that remains true even if we use
  //  linear interpolation ofr the phase - although that does tpyically improve the accuracy, it 
  //  does so only by a factor 10 - and in some cases (B = 50, M = 101), it even gets worse (by 
  //  factor 5). ...so we should probably go with odd size windows - although the errors are very 
  //  small in any case, so it probably doesn't matter

  // todo: 
  // -plot the analyzed short-time spectrum, a highly oversampled version of it (->zero padding),
  //  the fitted parabola and a mark at the found frequency/amplitude (and maybe another mark at
  //  the correct freq/amp)
  // -maybe somehow also visualize the phase estimation (not yet sure, how)
  GNUPlotter plt;
}


void phaseInterpolation() // rename to sineModelPhaseInterpolation
{
  // Tests various phase interpolation methods of SinusoidalSynthesizer - we create a sinusoidal 
  // partial and let the synthesizer generate the phases and plot the results
  // todo: maybe optionally plot the (numeric) derivative of the phase arrays instead of the phase
  // arrays theselves)

  double fs = 10000;  // sample rate
  //int N = 1000;       // number of samples

  double ts = 0.01; // timescale

  // create data for some not too boring frequency trajectory:
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  //partial.prependDataPoint(ISP(  0*ts, 1000.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  1*ts,  800.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  2*ts, 1200.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  3*ts, 1100.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  4*ts,  700.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  5*ts,  500.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  6*ts,  500.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  7*ts,  600.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  8*ts,  800.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(  9*ts,  900.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP( 10*ts, 1000.0, 1.0, 0.0));

  partial.prependDataPoint(ISP(  0*ts, 1000.0, 1.0, 0.0));
  partial.appendDataPoint( ISP(  1*ts, 1000.0, 1.0, 0.0));
  partial.appendDataPoint( ISP(  2*ts, 1100.0, 1.0, 0.0));
  partial.appendDataPoint( ISP(  3*ts, 1200.0, 1.0, 0.0));
  partial.appendDataPoint( ISP(  4*ts, 1000.0, 1.0, 0.0));
  partial.appendDataPoint( ISP(  5*ts, 1000.0, 1.0, 0.0));

  //partial.prependDataPoint(ISP( 0*ts, 1000.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP( 5*ts, 1050.0, 1.0, 0.0));
  //partial.appendDataPoint( ISP(10*ts, 1100.0, 1.0, 0.0));


  // create and set up the synth and create time-axis at sample-rate:
  rsSinusoidalSynthesizer<double> synth;
  synth.setSampleRate(fs);

  // synthesize and plot the sound:
  RAPT::rsSinusoidalModel<double> model;
  model.addPartial(partial);
  std::vector<double> x = synth.synthesize(model);
  //plotVector(x);



  int N = (int) x.size();
  std::vector<double> td = partial.getTimeArray();
  std::vector<double> t(N);
  for(size_t n = 0; n < N; n++) 
    t[n] = n / fs;// fill time-array



  // let the synth generate the phases:
  std::vector<double> pi = synth.phasesViaTweakedIntegral(partial, td, t);
  std::vector<double> pc = synth.phasesHermite(partial, td, t, false); 
  std::vector<double> pq = synth.phasesHermite(partial, td, t, true); // quintic looks wrong
  RAPT::rsArray::unwrap(&pc[0], N, 2*PI);
  RAPT::rsArray::unwrap(&pq[0], N, 2*PI);

  // array for plotting the phase datapoints:
  std::vector<double> pd = partial.getPhaseArray();
  std::vector<double> fd = partial.getFrequencyArray();
  int M = (int) pd.size();
  pd = rsSinusoidalProcessor<double>::unwrapPhase(td, fd, pd);
  //pd = synth.unwrapPhase(td, fd, pd);
  //RAPT::rsArray::unwrap(&pd[0], M, 2*PI);

  std::vector<double> dp = (0.5/PI) * (pc-pq); 
  // normalized difference between the algorithms - at the datapoints, it must be an integer 
  // corresponding to the k in the formula pu = p + k*2*PI


  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(M, &td[0], &pd[0]);
  plt.addDataArrays(N, &t[0],  &pi[0]);
  plt.addDataArrays(N, &t[0],  &pc[0]);
  plt.addDataArrays(N, &t[0],  &pq[0]);
  plt.addDataArrays(N, &t[0],  &dp[0]);
  //plt.setGraphStyles("points pt 7 ps 1.2", "lines", "lines");
  plt.plot();

  // Observations:
  // -the cubic hermite phase is 2pi behind the the integrated phase from datapoint 4 onwards
  // -the slope of the cubic phase at the datapoints is still wrong (close to 0 everywhere)
}



/*

// move to RAPT::rsArray, maybe make cubic versions
template<class T>
void applyFadeIn(T* x, int N, int numFadeSamples)
{
  int nf = rsMin(numFadeSamples, N);
  for(int n = 0; n < nf; n++) {
    T t = T(n) / T(nf);
    x[n] *= t;
  }
}
template<class T>
void applyFadeOut(T* x, int N, int numFadeSamples)
{
  int nf = rsMin(numFadeSamples, N);
  for(int n = 0; n < nf; n++) {
    T t = T(n) / T(nf);
    x[N-n-1] *= t;
  }
}
template<class T>
void applyFadeInAndOut(T* x, int N, int numFadeSamples)
{
  applyFadeIn( &x[0], N, numFadeSamples);
  applyFadeOut(&x[0], N, numFadeSamples);
}

// convenience function - move to rs_testing:
std::vector<double> synthesizeSinusoidal(
  const RAPT::rsSinusoidalModel<double>& model, double sampleRate, double fadeTime = 0.0)
{
  rsSinusoidalSynthesizer<double> synth;
  synth.setSampleRate(sampleRate);
  //synth.setCubicAmplitudeInterpolation(true);
  std::vector<double> x = synth.synthesize(model);
  if(fadeTime > 0.0)
    applyFadeInAndOut( &x[0], (int) x.size(), int (fadeTime*sampleRate));
  return x;
}
*/



void writeTwoSineModelOutputsToFile(
  const char* fileName,
  const RAPT::rsSinusoidalModel<double>& model1,
  const RAPT::rsSinusoidalModel<double>& model2,
  const rsSinusoidalSynthesizer<double>& synth, bool writeResidualFile = false)
{
  // we need some awkwardness here to get two time-aligned outputs...
  std::vector<double> x = synth.synthesize(model1);
  std::vector<double> x1, x2;
  getPaddedSignals(&x[0], (int) x.size(), model2, synth, x1, x2);
  rosic::writeToStereoWaveFile(fileName, &x1[0], &x2[0], (int) x1.size(), (int)synth.getSampleRate());

  if(writeResidualFile) {
    std::vector<double> r = x2 - x1;
    rosic::writeToMonoWaveFile("Residual.wav", &r[0], (int)r.size(), (int)synth.getSampleRate());
    // todo: incorporate given filename into the resiudal filename - we really should have 
    // wave-writing function that take a std::string instead of const char*
  }
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
  rsSinusoidalSynthesizer<double> synth;
  synth.setSampleRate(fs);
  synth.setCubicAmplitudeInterpolation(true);
  //synth.setCubicPhaseInterpolation(true);
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
  RAPT::rsSinusoidalAnalyzer<double> sa;
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
  double fadeCycles = 5.0;    // number of cycles for fade in

  phase = PI/2; // for test
  //frequency = 800;
  //frequency = 808;
  //frequency = 809;
  //frequency = 810;
  frequency = 1000 * GOLDEN_RATIO/2;
  //frequency = 1000 * SQRT2_INV;

  double fadeTime = fadeCycles/frequency;  // in seconds


  // analsis parameters:
  int window = RAPT::rsWindowFunction::HAMMING_WINDOW;
  //int window = RAPT::rsWindowFunction::BLACKMAN_WINDOW;
  double freqRes = frequency; // frequency resolution
  int zeroPadFactor = 2;

  // create signal:
  double period = sampleRate / frequency;         // in samples
  size_t N = (size_t)ceil(length * sampleRate);   // number of samples
  std::vector<double> x = createSinusoid(N, frequency, sampleRate, amplitude, phase);
  applyFadeInAndOut(&x[0], (int)N, int(sampleRate*fadeTime));

  // create and set up analyzer:
  RAPT::rsSinusoidalAnalyzer<double> sa;
  int blockSize = sa.getRequiredBlockSize(window, freqRes, sampleRate);
  int hopSize = blockSize/2;
  double maxLevelThresh = sa.getRequiredThreshold(window, 0.0);
  sa.setWindowType(window);
  sa.setMaxFreqDeltaBase(20);
  sa.setBlockAndTrafoSize(blockSize, zeroPadFactor*blockSize);
  sa.setHopSize(hopSize);
  //sa.setRelativeLevelThreshold(-25);  // some spurious tracks will occur (with Hamming window)
  sa.setRelativeLevelThreshold(-15);    // no spurious tracks will occur (with Hamming window)
  //sa.setFadeInTime(0.01);    // -> 500 samples @50kHz
  //sa.setFadeOutTime(0.01);   // -> 500 samples @50kHz
  //sa.setMinimumTrackLength(0.021);  // should be a little above fadeInTime+fadeOutTime
  sa.setFadeInTime(0);
  sa.setFadeOutTime(0);
  sa.setMinimumTrackLength(0); 
  //plotSineModel(sa, &x[0], (int) x.size(), sampleRate);



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
  rsSinusoidalSynthesizer<double> synth;
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
  // Two sinusoids at around 1000 and 1100 Hz

  // input signal parameters:
  double fs = 44100;   // sample rate
  //double fs = 10000;   // sample rate
  double f1 =  1011;   // 1st frequency
  double f2 =  1107;   // 2nd frequency
  double a1 = 1.0;     // 1st amplitude
  double a2 = 1.0;     // 2nd amplitude
  double p1 = PI/2;    // 1st start phase
  double p2 = PI/2;    // 2nd start phase
  double L  = 9.2;     // length in seconds - use 0.2 for plot and 9.2 for wavefile generation
  double fade = 0.03;  // fade-in/out time for original signal

  // analysis parameters:
  typedef RAPT::rsWindowFunction WF;
  double freqRes  = f2-f1;   // frequency resolution
  double dBmargin = 20;      // dB margin over sidelobe level
  double blockSizeFactor = 1.0;  // factor by which to to make blockSize longer than necessary
  int zeroPaddingFactor = 4;         // zero padding factor
  //int window = WF::HANNING_WINDOW_ZN;
  //int window = WF::HAMMING_WINDOW;
  int window = WF::BLACKMAN_WINDOW;
  //int window = WF::BLACKMAN_HARRIS;


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

  std::vector<double> x = synthesizeSinusoidal(model, fs, fade);
  //plotSineModel(model, fs); // we need a margin for the y-axis - otherwise it looks like nothing
  //plotVector(x);

  // analyze the produced sound again and compare the original model and analysis result

  // create and set up analyzer and try to recover the model from the sound:
  RAPT::rsSinusoidalAnalyzer<double> sa;
  int blockSize = int(blockSizeFactor * sa.getRequiredBlockSize(window, freqRes, fs, false));
  int hopSize   = blockSize/2;
  int trafoSize = zeroPaddingFactor*blockSize;
  double levelThresh = sa.getRequiredThreshold(window, dBmargin);
  sa.setWindowType(window);
  sa.setRelativeLevelThreshold(levelThresh);
  sa.setBlockAndTrafoSize(blockSize, trafoSize);
  sa.setHopSize(hopSize);
  sa.setMaxFreqDeltaBase(20);  // Hz variation allowed between frames - todo: have a parameter
                               // that is independent from the hopSize
  sa.setFadeInTime(0);
  sa.setFadeOutTime(0);
  sa.setMinimumTrackLength(0);
  rsSinusoidalModel<double> model2 = sa.analyze(&x[0], (int)x.size(), fs);
  //plotTwoSineModels(model, model2, fs);
  // ...try to use the lowest possible values that give a good result

  // create and set up a sinusoidal synthesizer object and plot resynthesis result:
  rsSinusoidalSynthesizer<double> synth;

  typedef rsSinusoidalSynthesizer<double>::PhaseInterpolationMethod PIM;
  //synth.setPhaseInterpolation(PIM::cubicHermite);
  //synth.setPhaseInterpolation(PIM::quinticHermite);
  synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  synth.setSampleRate(fs);
  //plotSineResynthesisResult(model2, synth, &x[0], (int)x.size());

  //plotModelOutputComparison(model, model2, synth);

  // figure out, if there's a bias in the frequency estimates:
  double f1a = model2.getPartial(0).getMeanFreq();
  double f2a = model2.getPartial(1).getMeanFreq();
  double freqBias1 = f1 - f1a;
  double freqBias2 = f2 - f2a;

  // test - synthesize and resynthesize only one partial
  //model.removePartial( 1);
  //model2.removePartial(1);
  writeTwoSineModelOutputsToFile("TestTwoSinesResynthesis.wav", model, model2, synth, true);
  //plotModelOutputComparison(model, model2, synth);


  // Observations:
  // -with amplitudes == 1, using a Blackman window with zero-padding of 4 and resynthesis with
  //  -tweaked-integral/cubic natural: residual starts at -inf, jumps to -90dB then to -84 in the 
  //   16 bit wavefile
  //  -hermite cubic: resdidual between -60 and -63 all the time
  //  -quintic hermite: resdidual between -57 and -60 all the time
  //  ->the winner is the cubic natural in this case, but the residual gets worse (louder) over 
  //    time, whereas it stays constant for the hermite schemes - i think for signals shorter than
  //    20 sec, i guess cubic natural will outperform the others, for longer samples, the error 
  //    accumulation for the cubic natural may make one of the others preferable 
  //   ->can we have an interpolation that combines the advantages of hermite and natural (very 
  //     good accuracy that doesn't get worse over time)?
  // -could it have to do with the frequency estimation bias (the freq is slightly overestimated,
  //  enforcing a too steep slope at the datapoints which will have to be compensated by a 
  //  shallower slope somewhere in between the datapoints) - could it be that the "stiffness" of
  //  the natural spline due to the 2nd derivative constraint counteracts that effect?
  //  -with a Hamming window, there's still more of the sinusoids left in the residual due to 
  //  estimation error - even going up to zeroPadding = 8 doesn't help much with the Hamming
  //  window
  // -for the quintic, it does not seem to make a big difference, if the numeric frequency 
  //  derivative is used for the target 2nd derivative or if it's just set to zero - which is quite
  //  surprising
  // -interestingly, the Hanning window seems to perform better than the Hamming window
  // -with the Hanning window and oversampling/zero-padding of 4, we observe a residual that
  //  increases in volume over time - apparently there's some error-accumulation over time
  //  going on
  // -whether the window is even or odd does not seem to make a big difference
  // -when using a blockFactor = 2, we can estimate the frequencies more accurately, but the 
  //  amplitude envelope gets an additional fade-in/out character
  // -todo: try a Gaussian window - it may improve the frequency estimate due to the fact of having
  //  and exact parabola as transform
  //  -maybe try to fit a function other than the parabola to the mainlobe
  //  -if a chebychev window is used, maybe we can try to just fit the appropriate chebychev window
  //   spectrum (which is known analytically)

  // -idea: 
  //  -the first few frames have lots of zero valued samples which lead to an amplitude error,
  //   for example, in the first frame only half of the samples in the buffer are nonzero
  //  -maybe we should take this situation into account already in the spectrogram processor
  //   by multiplying the analyzed amplitudes by an appropriate compensation factor
  //  -maybe this compensation factor should be number of inital zero-samples divided by the
  //   blockSize - or by the ratio of the sum of used window samples over the sum of all window
  //   samples
  //  -something like 
  //   if(n <= blockSize/2)  // or <?
  //     a *= getLeftBoundaryCompensation(n)
  //  -maybe have functions in rsSpectrogram setEdgeCompensation
  // ...or actually this error occurs only because in the synthesis of the sinusoid we cut it off
  // abruptly - this is not supposed to happen in natural signals - here, the sample values outside
  // the range 0..N-1 actually *are* zero, so this edge-compensation would make sense only for 
  // spectrograms that are cut out from some original signal
  // -generally, we should probably set up parameters in such a way, that the allowed freq-delta
  //  corresponds to the bin-distance of the FFT bins - the gives a heuristic to select the 
  //  hopSize - as a function of the normalized frequency delta (like, in Hz per second)

  // maybe allow time varying frequencies and amplitudes - but maybe in next test - ramp up the
  // complexity
}


void sinusoidalAnalysis3()
{
  // A sinusoidal sweep

  // input signal parameters:
  double fs = 44100;   // sample rate
  double f1 = 1000;    // start frequency
  double f2 = 1010;    // end frequency
  double a1 = 0.25;    // start amplitude
  double a2 = 0.25;    // end amplitude
  double p1 = PI/2;    // start phase
  double L  = 0.3;     // length in seconds


  // analysis parameters:
  typedef RAPT::rsWindowFunction WF;
  //double freqRes  = std::min(f2, f1); // frequency resolution
  double freqRes  = 1000.0;
  double dBmargin = 20;               // dB margin over sidelobe level
  double blockSizeFactor = 1.0;       // factor by which to to make blockSize longer than necessary
  int zeroPaddingFactor = 4;          // zero padding factor
  //int window = WF::HANNING_WINDOW_ZN;
  int window = WF::HAMMING_WINDOW;
  //int window = WF::BLACKMAN_WINDOW;
  //int window = WF::BLACKMAN_HARRIS;


  // create a model and synthesize the sound:
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  RAPT::rsSinusoidalModel<double> model;

  double p2 = 0; // end phase todo: set it up correctly integrate frequency sweep analytically to
                 // figure out what the end value must be
  partial.prependDataPoint(ISP(0.0, f1, a1, p1));
  partial.appendDataPoint( ISP(L,   f2, a2, p2));
  model.addPartial(partial);
  std::vector<double> x = synthesizeSinusoidal(model, fs);
  //plotSineModel(model, fs);
  //plotVector(x);

  // create and set up analyzer and try to recover the model from the sound:
  RAPT::rsSinusoidalAnalyzer<double> sa;
  int blockSize = int(blockSizeFactor * sa.getRequiredBlockSize(window, freqRes, fs, false));
  int hopSize   = blockSize/2;
  int trafoSize = zeroPaddingFactor*blockSize;
  double levelThresh = sa.getRequiredThreshold(window, dBmargin);
  sa.setWindowType(window);
  sa.setRelativeLevelThreshold(levelThresh);
  sa.setBlockAndTrafoSize(blockSize, trafoSize);
  sa.setHopSize(hopSize);
  sa.setMaxFreqDeltaBase(20);
  sa.setFadeInTime(0);
  sa.setFadeOutTime(0);
  sa.setMinimumTrackLength(0);
  sa.setFreqPhaseConsistency(true);
  rsSinusoidalModel<double> model2 = sa.analyze(&x[0], (int)x.size(), fs);
  plotTwoSineModels(model, model2, fs);


  // figure out, if there's a bias in the frequency estimates:
  double fa = model2.getPartial(0).getMeanFreq(); // maybe have a getMeanFreq that takes two time-stamps as parameters
  double freqBias = fa - (f1+f2)/2;

  // create and set up a sinusoidal synthesizer object and plot resynthesis result:
  rsSinusoidalSynthesizer<double> synth;
  typedef rsSinusoidalSynthesizer<double>::PhaseInterpolationMethod PIM;
  //synth.setPhaseInterpolation(PIM::cubicHermite);
  //synth.setPhaseInterpolation(PIM::quinticHermite);
  synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  synth.setSampleRate(fs);
  //plotSineResynthesisResult(model2, synth, &x[0], (int)x.size());
  writeTwoSineModelOutputsToFile("TestSineSweepResynthesis.wav", model, model2, synth, true);
  plotModelOutputComparison(model, model2, synth);


  // Observations:
  // -it seems the measured phase is offset from actual phase by an amount that depends on the
  //  sweeping speed
  // -the resynthesized signal is shotrer than the original one - why?
  // -even for a very small sweep amount like 1000 to 1002 Hz, the resynthesized signal is very 
  //  much out of phase with the original - is the error in the analysis or synthesis?
  //  ->compute correct phase at the datapoint time-instants and compare with the measured
  //  data
  // -maybe plot a phase spectrogram, see video 02 - Sinusoidal model 2, 11:32 - we should take the 
  //  time derivative of the phase spectrum - i stable sinusoid should have a constant value, a 
  //  sweeping sinuosid a constantly increasing color (and then jump back, like a sawtooth wave)
  // -with 1000 and 1100, there's some symmetry - in the first half, blue is consitently leading in
  //  the second half, blue is lagging
  // -figure out, if the phase, measurement, freq-measurement or both are to blame
  // -if the freq-measurement is too inaccurate, maybe try the reassignment apporach for estimating
  //  the frequency - let the user switch between freq-estimation algos:
  //  parabolic, reassignment, maybe phase-delta? 
  // -maybe after the partials have been figured out roughly, refine the estimates by the 
  //  heterodyning phase-vocoder approach and/or g�rtzel algorithm
  // -maybe try a global de-biasing of the freq-estimates (the local re-adjustment turned out to 
  //  not not work very well due to the tendency of producing alternating corrections, but maybe a 
  //  global dibiasing works better - maybe try to compensate the bias by adding something or by 
  //  multiplying by a factor - try which works better
  //  -or: maybe the tendency to alternate can somehow by counteracted algorithmically? maybe by a 
  //   lowpass/averaging scheme? ...or mayb detect alternations and do averaging, if it's detected
  // -figure out, how the spectrum of a (windowed) sine sweep actually looks like - especially
  //  the phase-spectrum
  // -can we somehow compute a freq-slope for each bin and frame just from the data within the 
  //  frame?
  // -todo: make the makeFreqsConsistentWithPhases a user option that can be set from client code
  //  and experiment with that

  // todo: there's alot of code duplication with respect to sinusoidalAnalysis3 - maybe we should 
  // make a class SineModelExperiment which factors out all the common code - we will likely need 
  // many more experiments with different input signals and lots of code will have to be duplicated
  // without such a class
  // maybe it should have functions like
  // setInputSineSweep(f1, f2, a1, a2, L)
  // setInputTwoSines(f1, f2, a1, a2, L)
  // setInputWobblySine()
  // setInputSawtooth()..
  // runExperiment()
}

/*
// move function somewhere, where we can see it also from the tests repo - maybe make a new
// juce module rs_tests (that should include also the plotting stuff, signal generation, etc.):
void testHarmonicResynthesis(const std::string& name, std::vector<double>& input, 
  double fs, bool writeWaveFiles, bool plotResults)
{
  // analyze, resynthesize and create error signal:
  double* x = &input[0];   // pointer to first sample (for convenience)
  int Nx = (int) input.size();
  rsHarmonicAnalyzer<double> analyzer;
  analyzer.setSampleRate(fs);
  analyzer.setSincInterpolationLength(512);
  RAPT::rsSinusoidalModel<double> mdl = analyzer.analyze(x, Nx);
  //plotSineModel(mdl, fs);
  std::vector<double> output = synthesizeSinusoidal(mdl, fs); 
  std::vector<double> error = output-input;
  double* y = &output[0]; int Ny = (int) output.size(); // again, for convenience
  double* e = &error[0];  int Ne = (int) error.size();  // dito

  // write original, resynthesized and error signals to files, if desired:
  if(writeWaveFiles == true) {
    rosic::writeToMonoWaveFile((name + "Original.wav").c_str(),      x, Nx, (int)fs);
    rosic::writeToMonoWaveFile((name + "Resynthesized.wav").c_str(), y, Ny, (int)fs);
    rosic::writeToMonoWaveFile((name + "Error.wav").c_str(),         e, Ne, (int)fs);
  }

  // plot original, resynthesized and error signals, if desired:
  if(plotResults == true) {
    GNUPlotter plt;
    plt.addDataArrays(Nx, x);
    plt.addDataArrays(Ny, y);
    plt.addDataArrays(Ne, e);
    //plt.addDataArrays(Ne-2000, &e[1000]);  // middle part of error
    plt.plot();
  }
}
*/


RAPT::rsSinusoidalPartial<double> phaseAlternatingPartial(int numDataPoints, double timeDelta,
  double freq, double phase1, double phase2)
{
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  partial.prependDataPoint(ISP(0.0, freq, 1.0, phase1));
  int count = 1;
  double phase;
  while(count < numDataPoints)
  {
    if( RAPT::rsIsOdd(count) )
      phase = phase2;
    else
      phase = phase1;
    partial.appendDataPoint(ISP(count*timeDelta, freq, 1.0, phase));
    count++;
  }
  return partial;
}

void phaseFreqConsistency()
{
  // Tests the phase/frequency consistency algorithm for the sinusoidal model. The partial 
  // frequencies are re-adjusted in a way such that their numeric integral happens to exactly hit
  // the measured phase values.

  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  partial = phaseAlternatingPartial(200, 0.01, 1000, 0.0, 0.2*PI);
  double error = partial.getMaxFreqPhaseInconsistency(); // should be 0.2*pi

  RAPT::rsSinusoidalProcessor<double>::makeFreqsConsistentWithPhases(partial);
  //partial.makeFreqsConsistentWithPhases();

  error = partial.getMaxFreqPhaseInconsistency();        // should be zero now

  // plot the re-adjusted freq-trajectory:
  std::vector<double> freqs = partial.getFrequencyArray();
  rsPlotVector(freqs);
  // this looks wrong! the frequency delta between adjacent datapoints increases toward the ends
  // commenting out rsMinSqrDiffWithGivnSum(&f[0], &sum[0], M); in makeFreqsConsistentWithPhases 
  // has the effect that the alternation amplitude increases toward the start section
  // ...could the target-phase computation be wrong? ..hmm..looks ok - actually, it's totally
  // implausible that the result should be the solution to the least-squares-of-differences
  // problem...maybe the minimization procedure is still wrong? we need to make more tests with 
  // that - create examples in sage, solve them analytically and compare with results from the
  // procedure

  // commenting the phase re-adjustment (if(delta > maxDelta) ...) in rsConsistentUnwrappedValue0
  // gives different results - still worng but in different ways - we need to re-consider the 
  // computation of the measurement-consistent target-phase value - i think it is wrong to compute
  // the preliminaryUnwrappedValue based on the old freq-values - we must use updated freq-values
  //...but we don't have any available yet...what to do? maybe impose the condition that some 
  // frequency (maybe somewhere in the middle where the estimates are supposedly better than at the
  // ends) should stay put...but it's simler to implement to let the first datapoint stay put

  // maybe enforcing the freq-data to be consistent with the phase-data is not sucha good idea, 
  // after all? essentially, it removes the original freq-data and the new freq data is redundant
  // with the phase data - the phase data could actually be thrown away, since we now can 
  // reconstruct it by numerical integration of the freq - that means, we have actually destroyed
  // data that could be potentially meaningful ...hmmm

  // try different approaches to freq-refinement: 
  // -take the numeric derivative of the unwrapped phase
  // -use frequency reassigment, if that makes sense
  // -use longer windows and parabolic interpolation also in the case of rsHarmonicAnalyzer
  // -maybe a combination of two or more approaches is best? like 
  //  parabolic interpolation -> numeric derivative -> consistency enforcement

  int dummy = 0;
}


void testHarmonicResynthesis(const std::string& name, double fs, int N, double f0 = 0)
{
  // setup (comment out "doStuff = true", if you don't want stuff to be done):
  bool writeWaveFiles = false, plotResults = false;
  writeWaveFiles = true;
  plotResults    = true;

  std::vector<double> input = createNamedSound(name, fs, N); 
  //std::string name2 = name + std::to_string(f) + "Hz";
  testHarmonicResynthesis(name, input, fs, f0, writeWaveFiles, plotResults);
}

void harmonicAnalysis1()  // rename to harmonicResynthesis
{

  // todo: fix findCosistentPhase (move it into rsSinousoidalPartial and write unit test)
  // or better: to AudioFunctions, the re-activate freq-refinement, then add the new algo

  //testHarmonicResynthesis("Sine_Freq=500_Amp=0.5",      44100, 5000);
  //testHarmonicResynthesis("Cosine_Freq=500_Amp=0.5",    44100, 5000);

  //testHarmonicResynthesis("LowpassSaw_Freq=523.25_kMax=10",  44100, 15000);



  //testHarmonicResynthesis("TremoloSine_Freq=200_Rate=10_Depth=5", 44100, 15000);
  // try a little bit of amplitude modulation - could the spurious high-freq components be due to 
  // that? the varying amplitude of low-freq components somehow gets translated to spurious 
  // high-freqs? maybe a tremolo sine - yes! that seems to be it! when the sine is amp-mdoulated,
  // it appears not as single spike but as spike-with-a-tail in the short-time spectrum



  //testHarmonicResynthesis("TwoSines",   44100, 5000);
  //testHarmonicResynthesis("ModalPluck", 44100, 5000);
  // convert all calls to include the frequency in the string (done), then get rid of the frequency 
  // parameter of the function


  //testHarmonicResynthesis("TwoSines_Freq1=500_Freq2=1000_Amp1=1.0_Amp2=1.0", 44100, 5000, 200);
  // good for testing if windows may spectrally resolve the harmonics


  //testHarmonicResynthesis("TwoSines_Freq1=200_Freq2=6000_Amp1=1.0_Amp2=1.0", 44100, 5000, 200);

  testHarmonicResynthesis("TwoSines_Freq1=200_Freq2=6100_Amp1=1.0_Amp2=1.0", 44100, 5000, 200);
  // 2nd harmonic gets missed somehow


  //testHarmonicResynthesis("TwoSines_Freq1=100_Freq2=10020_Amp1=0.5_Amp2=0.1", 44100, 5000, 100);
  // -let's figure out, if the resynthesized signal gets out-of-phase with respect to the orginal 
  //  at higher partials, when they are not exactly located at a harmonic
  // -try different phase interpolation methods and see how this influences the phase decoherence 
  //  between the datapoints
  //  -it changes the shape of how the decoherence raises but not the qualitative behavior
  //   ...but the quinitc has some additional weirdness
  // -plot instantaneous phase trajectories - there seems to be a discontinuity in the middle 
  //  between the datapoints ...could it choose a different path in the wrappedInterpolation in
  //  both half-waves?
  //  ->resynthesize only the high-freq component and compare to original high-freq component
  // -doesn't find correct cycle-marks if fundamental is not given
  // -removing all other partials (except at 100 and 10020) removes the "phase decoherence" issue
  //  ->it's not actually the inharmonic partial itself, whose phase decoheres, but the influence
  //  of all the other partials -> try multiple cycles per window (2 or 4) with windows that show
  //  rolloff


  //testHarmonicResynthesis("TwoSines_Freq1=200_Freq2=2050_Amp1=0.3_Amp2=0.2", 44100, 5000);
  // hits assert because the zero freq partial has some nonzero phase values - how they do arise
  // during analysis? does the FFT produce them? yes - assert was commented

  // if we refine the estimated freqs by phase derivatives with(Freq1=200_Freq2=2050): 
  // -the fundamental gets a bit wiggly around 200Hz 
  //  -the edges get worse (134.25..224 with or without ends-extrapolation of numeric derivative)
  // -the 2050Hz frequency gets estimated quite well indeed 
  // -all the other partial estimates (which have zero amplitude) get their estimated 
  //  frequencies also shifted up by 50Hz 
  // -and the "corrected" DC component actually doesn't stay at DC 
  //  -todo: fix the phase at zero for the DC component


  // with 200/2050 Hz we can clearly see the buzzing artifact, the residual looks similar if we
  // use  200/1950 - there are four sorts of artifacts: upward jumps, upward spikes, downward jumps
  // and downward spikes that alternate in that order
  // -> plot the amp, freq and phase trajectories - especially the (interpolated) phase is probably
  // interesting
  // maybe also have a look at those partials that whould have zero amplitude - maybe their 
  // contribution messes up the signal? try to resynthesize without them
  // injecting rsPlotVector(rsDifference(p)); into rsSinusoidalSynthesizer<T>::synthesizePartial 
  // shows that the difference between adjacent phase values is not a (roughly) constant function 
  // as it should be but shows spikes ...there must be some error in the instantaneous phase 
  // computation. these spikes seem to occur at the datapoints for the fundametal, for example at 
  // samples 2535, 2756, .. with 200/2025 - for higher partials, more spikes in between appear
  // using rsPlotVector(rsDifference(p)); shows that the phase array is not properly unwrapped
  // hmm..ok - i think, i'm not using unwrapping currently...switching the phase-interpolation 
  // algo letzts the spikes disappear - phase looks good now, buzz is still there
  // ...could it be that all the other FFT bins that do not belong to any proper partial sort of
  // conspire to create that jump at that particular instant? maybe it would help to use 
  // cycle-marks at zero corssing of the original waveform rather than those of the fundamental?
  // ...check, how i did it in the matlab code - no - it uses a filtered signal, too - but a much
  // higher order filter (10 passes of a bidirectional bandpass) :-O - nope - ramping up the 
  // filter order doesn't help
  // but i think, the matlab code measures the phase *at* the cycle-marks, not in between (it
  // doesn't do the shift of the FFT buffer by one half) ...should that matter - i don't know
  // why but maybe

  // adding mdl.keepOnly({0, 9}); to testHarmonicResynthesis (after mdl.removePartial(0)) 
  // removes the buzzing - it seems to be indeed the combined effect of all other overtones'
  // contributions - the conspire to produce edges - how can we avoid this? maybe obtain datapoints
  // also midway between the current ones? this would force the resynthesized phase to be in sync
  // with the original phase at the instants that are currently problematic ...and/or maybe the 
  // phase-based frequency estimate refinement could help against this?

  // i think, i know why the buzz occurs: try resynthesizing without the fundamental - it tries
  // to synthesize segments the 2030Hz component with discontinuities because in the FFT buffer,
  // the wave is indeed cut off discontinuously - the resynthesis tries to model these 
  // discontinuities via all other partials - frequency refinement may indeed plausibly break up
  // the phase-coherence of all theses partials that produces the edge - with readjusted 
  // frequencies, the model also won't be strictly harmonic anymore - which is a good thing since
  // the input sound is in fact inharmonic

  // ...but no - the buzz seems to be indeed due to the combined contribution of all partials that
  // are supposed to be zero, not due to the frequency estimation error of the inharmonic partial 
  // - because refining freq-estimates or not makes not much difference with regard to the buzz, 
  // but removing all other overtones (with near zero amplitude) does make a difference

  // i think, the best thing to try next is indeed to use blocks with 2 or 4 cycles, windowing, 
  // parabolic freq-estimation and (maybe) magnitude thresholding...maybe the number of cycles
  // pre block can be an integer user parameter (powers of two, but maybe that constraint can be
  // relaxed later)...and maybe use the dolph-chebychev window  - it's optimal for that purpose..
  // ...or maybe not? maybe sidelobe-rolloff is actually desirable in this application? i think so,
  // because it will lower the analyzed amplitudes of non-existing partials far from the actual
  // partial - try my new windows with given sidelobe rolloff - 

  //

  // when Freq2=2100, each block has exactly 10.5 cycles of the inharmonic wave in it and we get 
  // only discontinuities in the derivative as opposed to outright signal jumps - this all makes
  // sense now



  //testHarmonicResynthesis("VibratoSine", 44100, 95000);
  // produces a resiudal that looks like a train of (sort of) triangular spikes, amplitude 
  // modulated by the vibrato frequency - maybe try to use smoothing on the freq- and amp 
  // trajectories after interpolating them - i guess, it's an artifact from the interpolation 
  // process. the amp env of the resynthesized signal shows some slight amplitude modulation, too
  // actually, the resynthesized signal has also some overtones which shouldn't be there

  // maybe let the function not take a frequency parameter but instead pass the desired frequency 
  // and maybe other parameters as part of the name - for example TwoSines_Freq1=200_Freq2=2100
  // make functions getParameter(string soundName, string paramName. this is more flexible. maybe
  // have reasonable default values, if no parameter is given



  // todo: test with less high freq rolloff (makes it more difficult)
  // -when the first cycle mark ends up at 0, we get an access violation - we don't treat the case,
  //  where there is no initial partial cycle
  //  -this happens with N=8000, key=64, fs=44100 and createSineWave -> fixed

  // -try with saw-wave, square-wave, etc

  // -try resynthesizing without DC

  // -try a better sinc-interpolator (longer kernel and/or better window) ...maybe also try worse
  //  interpolators to see, if it indeed has to do with the interpolation
  // -effects of sinc interpolator kernel-length on approximate amplitude of residual for 
  //  TwoSines200 (1st number is kernel length, 2nd number residual amplitude): 
  //      8   :  2e-5
  //     15.5 : 20e-5
  //     16   :  6e-5
  //     16.5 :  8e-5
  //     17   : 15e-5
  //     20   :  8e-5
  //     64   :  4e-5
  //    256   :  1e-5
  //    512   :  5e-6
  //   1024   :  2e-6
  //  there seems to be a general trend of the amplitude being inversely proportional to kernel 
  //  length but there are weird oscillations
  //  todo: make a sinc-interpolator test: 
  //  create noise -> upsample -> downsample -> create difference ...the amplitude of the 
  //  difference should be as low as possible, try different downsampling factors, try different
  //  window functions - maybe collect data of error amplitude as function of kernel length and
  //  plot that
  //  todo: implement an algorithm that doesn't interpolate to a fixed block-size but instead uses
  //   a variable block-size and Bluestein FFT. this bypasses interpolation problems, but may 
  //   introduce other problems due to rounding of cycle-marks (but maybe that's less problematic
  //   - we'll see)

  // -try different synthesis settings (phase- and amplitude interpolation scheme, etc.)
  //  it's actually quite plausible that the error signal may be due to the interpolation of 
  //  the amplitude and/or phase (maybe more so due to amplitude)
  //  -this would also explain, why the error is larger in the attack-section - the original 
  //   envelope is farther away from being linear there
  //  -nope: for the pluck sound, there's no visible difference between linear and cubic amplitude
  //   interpolation (at least, in the later section of the sound, around 4000)
  //   ...this seems reasonable because we see an error signal also for steady (constant amplitude) 
  //   sine resynthesis - and when the amplitude is flat, both interpolators return the same value

  // todo:

  // -check artifact at end - could it be that the phase at the last datapoints is computed 
  //  wrong? ..more likely, it's due to the sudden cutting off in the middle of a cycle, this 
  //  discontinuity introduces high-freq content which is randomly messed up in resynthesis unless
  //  we would actuall do an iFFT of the block - only in this case, we can expect the discontinuity
  //  to be faithfully reconstructed - but with an osc bank, amplitude-interpolation, etc. the step
  //  gets "randomly" messed up


  // verify, if the resynthesizer avoid generating freqs above fs/2, i.e. if it guards against
  // aliasing - nope, it doesn't (!)
  
  // we have very high frequencies in the transient when the ampSlope is set to -3, with
  // -6, it looks very good - i think, at -3, it gets the cycle-marks wrong - maybe tweak the
  // cycle-mark-finder parameters and/or make them available to client code

  // todo: create a new project/repo where we test the analysis/resynthesis framework on real world
  // sample data (which should not go into the RS-MET repo)

  // maybe the framework can be used to compress samples - in areas where the partials are well
  // approximated by interpolated tarjectories, we could reduce the density of datapoints
  // ...the compression would be lossy, though

  // maybe make sure that y has the same size as x...maybe wrap this analysis/resynthesis roundtrip
  // into a convenience class

// look at files decompositionSteelGuitar002.m, testHarmonicAnalysis.M
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