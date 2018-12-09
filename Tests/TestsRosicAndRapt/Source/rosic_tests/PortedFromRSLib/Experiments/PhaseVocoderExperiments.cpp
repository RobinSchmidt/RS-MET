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

  rsPhaseVocoderD pv;             // for conveniently calling the static functions
                                       
  // create analysis and synthesis windows:
  double wa[B], ws[B];
  pv.hanningWindowZN(wa, B);
  pv.hanningWindowZN(ws, B);
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


  int J = rsPhaseVocoderD::getNumFrames(N, H);
                                      
  // create the window function:
  double wa[B], ws[B], w[B];
  rsPhaseVocoderD::hanningWindowZN(wa, B);
  rsPhaseVocoderD::hanningWindowZN(ws, B);
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
  double s = rsPhaseVocoderD::getWindowSum(wa, ws, B, H);
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
  static const int N  = 10000;           // number of samples in the test signal
  static const int P  = 1;              // zero-padding factor
  static const int M  = B*P;            // FFT size
  static const int K  = M/2 + 1;        // number of non-redundant bins
  double           fs = 44100;          // samplerate
  double           f  = 5000;            // sinusoid frequency

  // A hopsize of B/4 will result in a constant when overlapping successive frames, assuming that
  // the window is applied twice (once in the analysis stage and once in the synthesis stage). This
  // is the desired condition for perfect resynthesis without modulation artifcats - i.e. no 
  // explicit demodulation will be necessarry.

  // create the window function:
  double w[B];
  rsPhaseVocoderD::hanningWindowZN(w, B); // todo: create also the time-derivative and the 
                                         // time-ramped window for reassignment later

  // create the test signal:
  double x[N];
  createSineWave(x, N, f, 1.0, fs);
  RAPT::rsArray::scale(x, N/2, 0.1);   // amplitude switch in the middle of the signal

  // compute the complex spectrogram:
  rsMatrix<rsComplexDbl> s = rsPhaseVocoderD::complexSpectrogram(x, N, w, B, H, P);
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
  std::vector<double> y  = rsPhaseVocoderD::synthesize(s, w, B, H, w);
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
  int dummy = 0;
}


template<class T>
void synthesizePartial(const rsSinusoidalPartial<T>& partial, T* x, int numSamples, 
  T sampleRate)
{
  // figure out number of samples to produce:
  int nStart = (int) floor(sampleRate * partial.getStartTime());
  int nEnd   = (int) ceil( sampleRate * partial.getEndTime());
  nStart = rsClip(nStart, 0, numSamples-1);
  nEnd   = rsClip(nEnd,   0, numSamples-1);
  int N = nEnd - nStart;

  // create arrays for non-interpolated instantaneous parameter data:
  size_t M = partial.getNumDataPoints();
  std::vector<T> td(M), fd(M), ad(M), wpd(M);
  for(size_t m = 0; m < M; m++) {
    rsInstantaneousSineParams<T> dp = partial.getDataPoint(m);
    td[m]  = dp.getTime();          // time data
    fd[m]  = dp.getFrequency();     // frequency data
    ad[m]  = dp.getAmplitude();     // amplitude data
    wpd[m] = dp.getWrappedPhase();  // wrapped phase data
  }

  // obtain preliminary uwrapped phase data points by numerically integrating the frequency:
  std::vector<T> upd(M);
  rsNumericIntegral(&td[0], &fd[0], &upd[0], (int)M, wpd[0]);
  upd = 2*PI*upd; // convert from "number of cycles passed" to radians

  // incorporate the target phase values into the unwrapped phase:
  bool accumulatePhaseDeltas = true;  // make user parameter - experiment, which is better
  for(size_t m = 0; m < M; m++)
  {
    T wp1 = fmod(upd[m], 2*PI); // 0..2*pi
    //T wp2 = wpd[m] + PI;        // 0..2*pi
    T wp2 = fmod(wpd[m], 2*PI); // 0..2*pi
    T d   = wp2-wp1;            // -2*pi..2*pi, delta between target and integrated frequency
    if(d < 0) d += 2*PI;        // 0..2*PI
    if(d > PI)                  // choose adjustment direction of smaller phase difference
      d = 2*PI - d;             // in 0...pi ...check, if this formula is correct 
    upd[m] += d;                // re-adjust final unwrapped phase
    if(accumulatePhaseDeltas)
      for(size_t k = m+1; k < M; k++) // re-adjustment at m should also affect m+1, m+2, ...
        upd[m] += d;
  }
  // maybe the unwrapped phase computation should be factored out into a function 

  // interpolate the amplitude and unwrapped phase data to sample-rate:
  std::vector<T> t(N), f(N), a(N), p(N); // interpolated instantaneous data
  for(size_t n = 0; n < N; n++)          // fill time-array
    t[n] = (nStart + n) / sampleRate;    // ...optimize

  // interpolate phase:
  //rsInterpolateLinear(&td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);
  rsNaturalCubicSpline(&td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);

  // interpolate amplitude:
  //rsNaturalCubicSpline(&td[0], &ad[0],  (int)M, &t[0], &a[0], (int)N);
  rsInterpolateLinear(&td[0], &ad[0],  (int)M, &t[0], &a[0], (int)N);

  // it seems cubic interpolation for the phase and linear for the amplitude is most suitable,
  // although, for the amplitude, we may also use cubic - but linear for the phase leads to aduible
  // artifacts (sort of clicks) at the segement junctions
  // -maybe linear interpolation of frequency with subsequent integration would work (due to the
  //  smoothing effect of integration) - but that would make it difficult to incorporate the 
  //  target phases (i think, we would have to produce an interpolated phase-delta array and add 
  //  that to an interpolated preliminary unwrapped-phase array)
  // -maybe using cubic interpolation for frequency, than integrating and then adding an 
  //  interpolated phase-delta array could give and even smoother freq-trajectory? maybe try it
  // -or maybe use higher order numeric integration on the non-interpolated freq-data?
  // -maybe the synthesizer should have "presets" for most useful combinations of synthesis 
  //  parameters


  // maybe the user should be able to select the interpolation method and maybe a bidirectional 
  // smoothing filter (this should be set separately for amplitude and phase)

  // synthesize the sinusoid and add it to what's already there:
  std::vector<T> s(N); // needed here only for plotting, remove for production code
  for(size_t n = 0; n < N; n++)
    s[n] = x[nStart+n] += a[n] * sin(p[n]);

  GNUPlotter plt;
  //plt.addDataArrays(M, &td[0], &fd[0]);
  //plt.addDataArrays((int)M, &td[0], &upd[0]);
  //plt.addDataArrays((int)N, &t[0],  &p[0]); // interpolated phase
  //plt.addDataArrays((int)N, &t[0],  &a[0]);   // interpolated amplitude
  plt.addDataArrays((int)N, &t[0],  &s[0]);   // produced sinusoid
  plt.plot();
}

template<class T>
std::vector<T> synthesizeSinusoidal(const rsSinusoidalModel<T>& model, T sampleRate)
{
  int N = (int) ceil(sampleRate * model.getEndTime()); // number of samples
  std::vector<T> x(N);
  for(size_t i = 0; i < model.getNumPartials(); i++)
    synthesizePartial(model.getPartial(i), &x[0], N, sampleRate);
  return x;
}

void sinusoidalModel1()
{
  typedef RAPT::rsInstantaneousSineParams<double> ISP;
  RAPT::rsSinusoidalPartial<double> partial;
  RAPT::rsSinusoidalModel<double> model;
  //RAPT::rsSinusoidalSynthesizer<double> synthesizer;

  double stretch = 1.0;
  partial.appendDataPoint(ISP(stretch*0.0, 100.0, 0.4, 0.0)); // time, freq, amp, phase
  partial.appendDataPoint(ISP(stretch*0.4, 100.0, 0.2,  PI/2));
  partial.appendDataPoint(ISP(stretch*0.8, 150.0, 0.8, -PI/2));
  partial.appendDataPoint(ISP(stretch*1.2, 100.0, 0.4, 0.0));
  partial.appendDataPoint(ISP(stretch*1.6, 200.0, 0.2,  PI));
  partial.appendDataPoint(ISP(stretch*2.0, 100.0, 0.8, 0.0));

  //// test with all phases zero:
  //partial.appendDataPoint(ISP(stretch*0.0, 100.0, 0.4, 0.0)); // time, freq, amp, phase
  //partial.appendDataPoint(ISP(stretch*0.4, 100.0, 0.2, 0.0));
  //partial.appendDataPoint(ISP(stretch*0.8, 150.0, 0.8, 0.0));
  //partial.appendDataPoint(ISP(stretch*1.2, 100.0, 0.4, 0.0));
  //partial.appendDataPoint(ISP(stretch*1.6, 200.0, 0.2, 0.0));
  //partial.appendDataPoint(ISP(stretch*2.0, 100.0, 0.8, 0.0));

  // it seems, the phase is still wrong and there's one sample too little in the resulting signal

  // cycles[m] = cycles[m-1] + 0.5*(freq[m]+freq[m-1]) * (time[m]-time[m-1])
  // ...hmm...i really think this should not be part of the data structure - it should be computed 
  // during synthesis

  // other idea: maybe instead of representing amplitude and phase, we can use a complex amplitude?
  // ...but that can actually be left to the synthesis algorithm - we may have one that uses amp 
  // and phase and another one that uses real/imag

  model.addPartial(partial);
  //double fs = 44100;
  double fs = 4000; // for plot
  //std::vector<double> x = synthesizer.synthesize(model, fs);
  std::vector<double> x = synthesizeSinusoidal(model, fs);


  rosic::writeToMonoWaveFile("SinusoidalSynthesisTest.wav", &x[0], (int)x.size(), (int)fs, 16);
  int dummy = 0;
}


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