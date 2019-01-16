#include "PhaseVocoderExperiments.h"

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
  pv.shortTimeSpectrum(x, N, n0, wa, B, M, X);

  // under construction - not yet complete

  int dummy = 0;
}

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

  // A hopsize of B/4 will result in a constant when overlapping successive frames, assuming that
  // the window is applied twice (once in the analysis stage and once in the synthesis stage). This
  // is the desired condition for perfect resynthesis without modulation artifcats - i.e. no 
  // explicit demodulation will be necessarry.

  // create the window function:
  double w[B];
  rsWindowFunction::hanningZN(w, B); // todo: create also the time-derivative and the 
                                     // time-ramped window for reassignment later

  // create the test signal:
  double x[N];
  createSineWave(x, N, f, 1.0, fs);
  RAPT::rsArray::scale(x, N/2, 0.1);   // amplitude switch in the middle of the signal

  // compute the complex spectrogram:
  rsSpectrogramD sp;
  // sp.setAnalysisBlockSize(B);
  // sp.setAnalysisHopSize(H);
  // sp.setAnalysisWindowType(W); 
  // ...
  rsMatrix<rsComplexDbl> s = sp.complexSpectrogram(x, N, w, B, H, P);
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
      phs[i][j] = arg(s(i, j));
    }
  }

  // compute frequency- and time-reassignment matrices:
  // ...

  // plot the magnitude spectrogram (later: with or without reassignment):
  plotSpectrogram(F, K, dB, fs, H);

  // resynthesize and plot signal:
  // sp.setSynthesisBlockSize(B);
  // sp.setSynthesisHopSize(H);
  // sp.setSynthesisWindowType(W); 
  // ...
  //std::vector<double> y  = rsSpectrogramD::synthesize(s, w, B, H, w);
  std::vector<double> y  = sp.synthesize(s, w, B, H, w);
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

// rename to testSinusoidalSynthesis1
void sinusoidalSynthesis1()
{
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  RAPT::rsSinusoidalModel<double> model, model2;
  //RAPT::rsSinusoidalSynthesizer<double> synthesizer;

  double stretch = 1.0;
  partial.appendDataPoint(ISP(stretch*0.0, 100.0, 0.4, 0.0));    // time, freq, amp, phase
  partial.appendDataPoint(ISP(stretch*0.4, 100.0, 0.2,  PI/2));
  partial.appendDataPoint(ISP(stretch*0.8, 150.0, 0.8, -PI/2));
  partial.appendDataPoint(ISP(stretch*1.2, 100.0, 0.4, 0.0));
  partial.appendDataPoint(ISP(stretch*1.6, 200.0, 0.2,  PI));
  partial.appendDataPoint(ISP(stretch*2.0, 100.0, 0.8, 0.0));

  // other idea: maybe instead of representing amplitude and phase, we can use a complex amplitude?
  // ...but that can actually be left to the synthesis algorithm - we may have one that uses amp 
  // and phase and another one that uses real/imag

  model.addPartial(partial);
  double fs = 44100;
  //double fs = 5000; // for plot
  //std::vector<double> x = synthesizer.synthesize(model, fs);
  std::vector<double> x = synthesizeSinusoidal(model, fs);

  // make a sinusoidal analysis of the sound that we have just created and re-create the sound
  // from the model that results from this analysis:
  // maybe try the following parameters: 
  // maxFreqDeviation df = 100
  // windowType = Hamming -> B = 4, L = -42.7
  // windowSize M >= B * fs / df = 4 * 44100 / 100 = 1764 -> use 1801 for odd size
  // threshold t >= L = -42.7 -> use t = -30
  // i think, the Blackman-Nutall window is superior to the Blackman-Harris window - it has a
  // slightly narrower mainlobe and better sidelobe rejection
  SinusoidalAnalyzer<double> sa;
  sa.setWindowType(RAPT::rsWindowFunction::HAMMING_WINDOW);
  sa.setMaxFreqDeltaBase(100);
  //sa.setBlockSize(1801);  // does not yet work
  //sa.setHopSize(225);
  sa.setBlockSize(2048);
  sa.setHopSize(256);
  sa.setRelativeLevelThreshold(-40);
  sa.setZeroPaddingFactor(2);
  //model2 = sa.analyze(&x[0], (int)x.size(), fs);
  //std::vector<double> y = synthesizeSinusoidal(model2, fs);
  // there are loads of spurious partials - also, we need to wrap the phase into -pi..pi for the 
  // fade-in/out datapoints, we may also need to apply fade-outs to all tracks that are alive until
  // the end to properly finish them (there are spurious partials with length 2 which shouldn't 
  // happen)
  // there's one non-spurious track - but it seems to not represent the actual signal very well
  // -> make tests with simpler signals and implement plotting facilities
  // maybe make a class SineModelPlotter which can plot the sinusoidal trajectories over the 
  // spectrogram and maybe plot also the time-domain waveform
  // - maybe have also a visualization of input and resynthesized signal

  plotSineModel(sa, &x[0], (int) x.size(), fs);

  //rosic::writeToMonoWaveFile("SinusoidalSynthesisTest.wav", &x[0], (int)x.size(), (int)fs, 16);
}

void sinusoidalAnalysis1()
{
  // test signal parameters:
  double sampleRate = 48000;  // in Hz
  double frequency  = 3000;   // in Hz
  double length     = 0.2;    // in seconds
  double startPhase = 0.0;    // in radians


  // create signal:
  double period = sampleRate / frequency;    // in samples
  int N = (int)ceil(length * sampleRate);    // number of samples
  std::vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = sin(n * 2*PI*frequency/sampleRate + startPhase);


  // create and set up analyzer:
  SinusoidalAnalyzer<double> sa;
  sa.setWindowType(RAPT::rsWindowFunction::HAMMING_WINDOW);
  sa.setMaxFreqDeltaBase(100);
  sa.setBlockSize(1024);
  sa.setHopSize(256);
  sa.setZeroPaddingFactor(4);
  sa.setRelativeLevelThreshold(-25);
  sa.setFadeInTime(0.01);
  sa.setFadeOutTime(0.01);
  sa.setMinimumTrackLength(0.021);  // should be a little above fadeInTime+fadeOutTime
  plotSineModel(sa, &x[0], (int) x.size(), sampleRate);

  // find model for the signal:
  rsSinusoidalModel<double> model = sa.analyze(&x[0], N, sampleRate);

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

  // todo: implement zero-phase windowing, arbitrary window-sizes and hop-sizes for spectrogram 
  // -> always check identity resynthesis with unit test

  // todo: resynthesize and create residual

  int dummy = 0;
}


// make various tests for the sinusoidal analysis of increasing level of difficulty:
// 1: single sinusoid with stable frequency and amplitude
//  -check, if the frequency a nd ampltude is etsimated accurately
//  -try frequencies that coincide with bin-centers (best case) and those that fall halfway 
//   between bins (worts case) and some intermediate cases
// 2: single sinuosoid with time-variying amplitude
// 3: two sinusoids with stable freq and amp
//  -check, how each influences the analysis of the other - as function of their frequencies
// 4: several sinusoids with various frequencies and start and endpoints
// 5: periodic sounds like triangle, square, saw
// ....


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