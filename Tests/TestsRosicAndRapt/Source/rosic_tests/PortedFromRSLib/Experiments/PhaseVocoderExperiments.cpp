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


  int J = rsSpectrogramD::getNumFrames(N, H);
                                      
  // create the window function:
  double wa[B], ws[B], w[B];
  rsSpectrogramD::hanningWindowZN(wa, B);
  rsSpectrogramD::hanningWindowZN(ws, B);
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
  rsSpectrogramD::hanningWindowZN(w, B); // todo: create also the time-derivative and the 
                                          // time-ramped window for reassignment later

  // create the test signal:
  double x[N];
  createSineWave(x, N, f, 1.0, fs);
  RAPT::rsArray::scale(x, N/2, 0.1);   // amplitude switch in the middle of the signal

  // compute the complex spectrogram:
  rsMatrix<rsComplexDbl> s = rsSpectrogramD::complexSpectrogram(x, N, w, B, H, P);
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
  std::vector<double> y  = rsSpectrogramD::synthesize(s, w, B, H, w);
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
  int nEnd   = (int) ceil( sampleRate * partial.getEndTime()) + 1;
  nStart = rsClip(nStart, 0, numSamples);
  nEnd   = rsClip(nEnd,   0, numSamples);
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
  bool accumulatePhaseDeltas = false;  // make user parameter - experiment, which is better - note
  // that accumulation gives rise to an O(M^2) complexity of the algorithm
  for(size_t m = 0; m < M; m++)
  {
    T wp1 = fmod(upd[m], 2*PI); // 0..2*pi
    T wp2 = fmod(wpd[m], 2*PI); // 0..2*pi
    T d   = wp2-wp1;            // -2*pi..2*pi, delta between target phase and integrated frequency
    if(d < 0) d += 2*PI;        // 0..2*PI
    if(d > PI)                  // choose adjustment direction of smaller phase difference
      d -= 2*PI;                // -pi..pi
    upd[m] += d;                // re-adjust final unwrapped phase
    if(accumulatePhaseDeltas)
      for(size_t k = m+1; k < M; k++) // re-adjustment at m should also affect m+1, m+2, ...
        upd[k] += d;
  }
  // maybe the unwrapped phase computation should be factored out into a function
  // maybe provide an alternative implementation that uses the measured (unwrapped) phases y and 
  // the measured instantaneous frequencies as y'' in a Hermite interpolation scheme (this is how
  // it's described in the literature). to unwrap the phase of datapoint m, take the phase of m-1
  // and add the integral over t[m-1]..t[m] of the average frequency (f[m-1]+f[m])/2 and choose
  // p[m] + 2*pi*k as unwrapped where k is chosen such that p[m] is closest to value obtained from
  // integration

  // interpolate the amplitude and unwrapped phase data to sample-rate:
  std::vector<T> t(N), f(N), a(N), p(N); // interpolated instantaneous data
  for(size_t n = 0; n < N; n++)          // fill time-array
    t[n] = (nStart + n) / sampleRate;    // ...optimize

  // interpolate phase:
  //rsInterpolateLinear(&td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);
  rsNaturalCubicSpline(&td[0], &upd[0], (int)M, &t[0], &p[0], (int)N);

  // interpolate amplitude:
  rsNaturalCubicSpline(&td[0], &ad[0],  (int)M, &t[0], &a[0], (int)N);
  //rsInterpolateLinear(&td[0], &ad[0],  (int)M, &t[0], &a[0], (int)N);

  // it seems cubic interpolation for the phase and linear for the amplitude is most suitable,
  // although, for the amplitude, we may also use cubic - but linear for the phase leads to audible
  // artifacts (sort of clicks) at the segement junctions
  // -maybe linear interpolation of frequency with subsequent integration would work (due to the
  //  smoothing effect of integration) - but that would make it difficult to incorporate the 
  //  target phases (i think, we would have to produce an interpolated phase-delta array and add 
  //  that to an interpolated preliminary unwrapped-phase array)
  // -maybe using cubic interpolation for frequency, than integrating and then adding an 
  //  interpolated phase-delta array could give and even smoother freq-trajectory? maybe try it
  // -or maybe use higher order numeric integration on the non-interpolated freq-data?
  // -but when smoothing comes into play later, linear interpolation of the phase might be not so
  //  bad as it is without smoothing
  // -however - when the data comes from an analysis, it will be typically much denser in time than
  //  in this test (like one datapoint per cycle), so the difference between the interpolation
  //  methods may be not that important anymore - but for sparse data, it is crucial
  // -maybe the synthesizer should have "presets" for most useful combinations of synthesis 
  //  parameters


  // maybe the user should be able to select the interpolation method and maybe a bidirectional 
  // smoothing filter (this should be set separately for amplitude and phase)

  // synthesize the sinusoid and add it to what's already there:
  std::vector<T> s(N); // needed here only for plotting, remove for production code
  for(size_t n = 0; n < N; n++)
    s[n] = x[nStart+n] += a[n] * cos(p[n]);
  // we use the cosine (not the sine) because that's what's used in the literature - probably 
  // because it's consistent with representing real sinusoids as the real part of complex
  // sinusoids (using the sine, we would have to take the imaginary part instead)


  GNUPlotter plt;
  //plt.addDataArrays(M, &td[0], &fd[0]);
  //plt.addDataArrays((int)M, &td[0], &upd[0]);
  //plt.addDataArrays((int)N, &t[0],  &p[0]); // interpolated phase
  //plt.addDataArrays((int)N, &t[0],  &a[0]);   // interpolated amplitude
  //plt.addDataArrays((int)N, &t[0],  &s[0]);   // produced sinusoid
  //plt.plot();
}

template<class T>
std::vector<T> synthesizeSinusoidal(const rsSinusoidalModel<T>& model, T sampleRate)
{
  int N = (int) ceil(sampleRate * model.getEndTime()) + 1; // number of samples
  std::vector<T> x(N);
  for(size_t i = 0; i < model.getNumPartials(); i++)
    synthesizePartial(model.getPartial(i), &x[0], N, sampleRate);
  return x;
}

// move to library:
template<class T>
rsMatrix<T> matrixMagnitudes(const rsMatrix<std::complex<T>>& A)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<T> mags(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      mags(i, j) = abs(A(i, j));
  return mags;
}
template<class T>
rsMatrix<T> matrixPhases(const rsMatrix<std::complex<T>>& A)
{
  int N = A.getNumRows();
  int M = A.getNumColumns();
  rsMatrix<T> phases(N, M);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < M; j++)
      phases(i, j) = arg(A(i, j));
  return phases;
}
// maybe factor out common code...maybe something like applyMatrixFunction with different
// input and output types for the template parameter

template<class T>
std::vector<int> peakIndices(T* x, int N, T threshToMax = 0)
{
  T max = RAPT::rsArray::maxValue(x, N);
  std::vector<int> peaks;
  for(int i = 1; i < N-1; i++)
    if(x[i] > x[i-1] && x[i] > x[i+1] && x[i] > threshToMax*max)
      peaks.push_back(i);
  return peaks;
} 
// move to library, maybe have additional criteria like a threshold with respect to the rms, minimum
// distance between the peaks, etc.

// move to library:
template<class T>
void fitQuadratic_m1_0_1(T *a, T *y)  // x = -1,0,+1
{
  a[0] = y[1];
  a[1] = 0.5*(y[2]-y[0]);
  a[2] = y[2] - a[0] - a[1];
}
template<class T>
T quadraticExtremumPosition(T *a)
{
  return T(-0.5) * a[1]/a[2];
}
template<class T> 
void spectralMaximumPositionAndValue(T *x, int k, T* pos, T* val)
{
  // coeffs of parabolic interpolant:
  T lowAmp = 0.0000001; // -140 dB - to prevent log-of-zero
  T a[3], y[3];
  y[0] = rsAmpToDbWithCheck(x[k-1], lowAmp);   // left
  y[1] = rsAmpToDbWithCheck(x[k],   lowAmp);   // mid
  y[2] = rsAmpToDbWithCheck(x[k+1], lowAmp);   // right
  fitQuadratic_m1_0_1(a, y);

  // find maximum position and evaluate parabola there:
  T d  = quadraticExtremumPosition(a);
  *pos = k + d;
  *val = rsDbToAmp(a[0] + (a[1] + a[2]*d)*d);
}

// interpolate between two values that are supposed to wrapping around and alway be in xMin..xMax, 
// for example -1..1, 0..1, 0..2*pi, -pi..pi, etc. t is the interpolation parameter between 0..1 
// where such that x = (1-t)*x0 + t*x1 = x0 + t*(x1-x0)
template<class T> 
T rsInterpolateWrapped(T x0, T x1, T t, T xMin, T xMax)
{
  T r  = xMax - xMin;  // range
  T du = x1 + r - x0;  // upper difference
  T dm = x1     - x0;  // middle difference
  T dl = x1 - r - x0;  // lower difference
  T au = rsAbs(du);
  T am = rsAbs(dm);
  T al = rsAbs(dl);
  T x;

  // add an appropriate fraction the delta that has the smallest absolute difference to x0:
  if(au < am && au < al)
    x = x0 + t*du;
  else if(am < al)
    x = x0 + t*dm;
  else
    x = x0 + t*dl;

  // re-wrap result into allowed interval:
  if(x > xMax)
    x -= r;
  else if(x < xMin)
    x += r;
  return x;
}
// move to library and check if it is correct (maybe by a unit test?)...and if it can be simplified



template<class T>
size_t findBestMatch(T freq, std::vector<RAPT::rsSinusoidalPartial<double>>& tracks,
  T maxFreqDeviation, const std::vector<bool>& trackContinued)
{
  T dfMin = RS_INF(T);  
  size_t bestIndex = 0;
  for(size_t i = 0; i < tracks.size(); i++) {
    T trackFreq = tracks[i].getEndFreq();
    T df = rsAbs(freq - trackFreq);
    if(df < dfMin && trackContinued[i] == false) { // look only in those that are not already continued
      dfMin = df;
      bestIndex = i;
    }
  }
  if(dfMin < maxFreqDeviation)
    return bestIndex;
  else
    return tracks.size();
}

// this implements the peak continuation step - for all current spectral peaks in newPeakData, find 
// a corresponding continuation partner among the activeTracks - 3 situations have to be handled:
// -when a partner is found, continue the track
// -when no partner is found, create a new track ("birth")
// -all active tracks that have not been used in this continuation are killed (i.e. moved to 
//  finishedTracks
template<class T>
void continuePartialTracks(
  std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeakData,
  std::vector<RAPT::rsSinusoidalPartial<T>>& activeTracks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& finishedTracks,
  T maxFreqDeviation, T frameTimeDelta, int direction) // additionally needed information
{

  // initializations:
  typedef std::pair<size_t, size_t> IndexPair;
  std::vector<IndexPair> continuationPairs;  // 1st: activeTracks, 2nd: newPeakData
  std::vector<size_t> killTrackIndices;      // index into activeTracks
  std::vector<size_t> birthPeakIndices;      // index into newPeakData
  size_t pkIdx;                              // peak index
  size_t trkIdx;                             // track index
  size_t numActiveTracks = activeTracks.size();
  std::vector<bool> trackContinued(numActiveTracks);
  for(trkIdx = 0; trkIdx < numActiveTracks; trkIdx++)
    trackContinued[trkIdx] = false;

  // Figure out which tracks have to be continued (and if so, which new peak data should be 
  // appended), which tracks have to be discontinued ("death") and for which peaks a new track
  // has to be created ("birth"):

  // loop over the new peaks to figure out birthes and continuations:
  for(pkIdx = 0; pkIdx < newPeakData.size(); pkIdx++) {

    trkIdx = findBestMatch(newPeakData[pkIdx].freq, activeTracks, 
      maxFreqDeviation, trackContinued); // looks only in those that are not already continued

    if(trkIdx == activeTracks.size())    // no match found
      birthPeakIndices.push_back(pkIdx);
    else {
      continuationPairs.push_back(IndexPair(trkIdx, pkIdx));
      trackContinued[trkIdx] = true;
    }
  }
  // Actually, it would be better to run the loop over the active tracks and find a match for
  // each track among the new peaks. This would allow to search for a frequency that is not exactly 
  // the last value in the respective track but one that is obtained by extrapolation of the 
  // frequency trajectory so far (linear prediction or polynomial extrapolation). Also, the search 
  // could take advantage of the fact that the peak array is sorted by frequency (to use binary 
  // instead of linear search). This would imply that we would have to work with a boolean 
  // "peakUsed" array instead of "trackContinued". Maybe keep both variants in a prototype 
  // implementation.

  // loop over the "trackContinued" array to figure out deaths:
  for(trkIdx = 0; trkIdx < numActiveTracks; trkIdx++) {
    if(trackContinued[trkIdx] == false)
      killTrackIndices.push_back(trkIdx);
  }

  // We have figured out the desired continuations, deaths and birthes. Now, we actually do them:

  size_t i;

  // continue matched tracks with new peaks:
  for(i = 0; i < continuationPairs.size(); i++) {
    trkIdx = continuationPairs[i].first;
    pkIdx  = continuationPairs[i].second;
    activeTracks[trkIdx].appendDataPoint(newPeakData[pkIdx]);
  }

  // kill discontinued tracks (where no matching peak was found for a track):
  if(killTrackIndices.size() > 0) {
    for(i = killTrackIndices.size()-1; i >= 0; i--) {
      trkIdx = killTrackIndices[i];
      rsAppend(finishedTracks, activeTracks[trkIdx]);
    }

    // todo: append an additional datapoint with zero amplitude to the killed track for a smooth
    // fade out (needs to take into account frameTimeDelta and direction to figure out the time for
    // the datapoint...actually just their product - maybe passed as one parameter)...but maybe the
    // fade-in/out may be shorter than one hopSize?

    rsRemove(activeTracks, trkIdx);
  }

  // create new tracks (where no matching track was found for a peak):
  for(i = 0; i < birthPeakIndices.size(); i++)
  {
    pkIdx = birthPeakIndices[i];
    RAPT::rsInstantaneousSineParams<T> newData = newPeakData[pkIdx];
    RAPT::rsSinusoidalPartial<T> newTrack;
    newTrack.appendDataPoint(newData);

    // todo: append an additional datapoint with zero amplitude for smooth fade-in

    activeTracks.push_back(newTrack);
  }
}

template<class T>
rsSinusoidalModel<T> analyzeSinusoidal(T* sampleData, int numSamples, T sampleRate)
{
  // -maybe pre-process the input signal by flattening the pitch and make the period coincide with
  //  a fraction of the analysis frame size 
  // -we should keep the time warping map around and when done with the analysis, use it to 
  //  re-map the time-instants in the model
  // -we should also keep the corresponding (instantaneous) frequency modifier around and apply
  //  the inverse to the freq data in the model
  //  ...but later - first, let's see how far we get without such a pre-processing and try to make 
  //  the algorithm work well under this (suboptimal) condition, too

  // make function parameters:
  int blockSize = 4096;
  int hopSize = 2048;
  int zeroPaddingFactor = 2;

  //// test:
  //blockSize = 1024;
  //hopSize = 32;
  //zeroPaddingFactor = 1;

  // Obtain a spectrogram:
  rsSpectrogram<T> sp;  // maybe it should be called rsSpectrogramProcessor because
                        // the spectrogram itself it the STFT matrix
  sp.setBlockSize(blockSize);
  sp.setHopSize(hopSize);
  sp.setZeroPaddingFactor(zeroPaddingFactor);
  //size_t numBins = sp.getNumNonRedundantBins();
  rsMatrix<std::complex<T>> stft = sp.complexSpectrogram(sampleData, numSamples);
  rsMatrix<T> mag = matrixMagnitudes(stft);
  rsMatrix<T> phs = matrixPhases(stft);



  // Initializations:
  typedef RAPT::rsInstantaneousSineParams<double> InstParams;
  typedef RAPT::rsSinusoidalPartial<double> Partial;
  //std::vector<Partial> partials;     // array of active partials (empty at first)
  int numFrames  = stft.getNumRows();
  int numBins    = stft.getNumColumns();
  int frameStep  = +1;  // +1: scan through frames forward, -1: backward
  int firstFrame = 0; 
  int lastFrame  = numFrames-1;
  if(frameStep == -1)
    rsSwap(firstFrame, lastFrame);
  int frameIndex  = firstFrame;
  double binDelta   = sampleRate / sp.getFftSize();
  double frameDelta = sp.getHopSize() / sampleRate;


  // ...maybe plot the spectrogram here...
  //plotPhasogram(numFrames, numBins, phs.getDataPointer(), sampleRate, sp.getHopSize());


  std::vector<Partial> activeTracks, finishedTracks;
  std::vector<InstParams> instPeakParams;

  // maybe we need separate arrays for activePartials and finishedPartials

  //bool plotShortTimeSpectra = true;

  // needed for plotting only:
  std::vector<T> freqs(numBins);
  for(int i = 0; i < numBins; i++)
    freqs[i] = i * sampleRate / sp.getFftSize();

  // loop over the frames:
  while(frameIndex != lastFrame) {

    //T time = frameIndex * sp.getHopSize() / sampleRate;
    T time = frameIndex * frameDelta;
    T* pMag = mag.getRowPointer(frameIndex);
    T* pPhs = phs.getRowPointer(frameIndex);
    std::complex<T>* pCmp = stft.getRowPointer(frameIndex);  // pointer to complex short-time spectrum

    //plotData(numBins/64, &freqs[0], pMag); // for development
    //plotData(numBins/64, &freqs[0], pPhs);

    // find spectral peaks:
    T peakThresh = rsDbToAmp(-30.0); // should be somewhere above the sidelobe level
    std::vector<int> peaks = peakIndices(pMag, numBins, peakThresh);

    // determine exact peak frequencies, amplitudes (exact phases are left for later):
    instPeakParams.resize(peaks.size());
    for(size_t i = 0; i < peaks.size(); i++) {
      T peakBin, peakAmp;
      spectralMaximumPositionAndValue(pMag, peaks[i], &peakBin, &peakAmp);
      T peakFreq = peakBin*sampleRate/sp.getFftSize();
      T peakPhase = 0; // preliminary - see comment below function for ideas 
      instPeakParams[i].time  = time;
      instPeakParams[i].freq  = peakFreq;
      instPeakParams[i].gain  = peakAmp;
      instPeakParams[i].phase = peakPhase;
    }

    // peak continuation, birth or death:
    double maxFreqDelta = 2*binDelta; // replace factor 2 by adjustable parameter
    continuePartialTracks(instPeakParams, activeTracks, finishedTracks, 
      maxFreqDelta, frameDelta, frameStep);

    frameIndex += frameStep;
  }

  rsSinusoidalModel<T> model;
  // model.addPartials(finishedTracks);
  // model.addPartials(activeTracks);

  // if(frameStep == -1) time-reverse all partials (well, reverse the arrays such that the time 
  // runs forward)

  return model;

  // algorithm:


  // -scan through this spectrogram from right to left - i.e. start at the end
  // -for each frame m, do:
  //  -find the spectral peaks in the frame
  //  -for each peak k in the frame, do:
  //   -figure out its exact frequency, amplitude and phase (and maybe time, if re-assignment is 
  //    used)
  //   -try to find a match in the set of active partials:
  //    -if a match is found, prepend the new datapoint to that matched partial
  //    -if no match is found, start a new partial at frame m+1 with zero amplitude and prepend
  //     the new datapoint (time m)
  //   -if there are active partials that have not been used up by this matching procedure, prepend
  //    a datapoint with amplitude zero at time m into them - they end now
  // 
  // -whenever a partial is created at time m, it will fade in from amplitude 0 at m-1 and whenever
  //  one is killed (i.e. not continued from m+1 to m), it will fade out to reach zero at time m

  // -maybe after this is all done, we may go over the data a second time to see if we can merge
  //  and/or remove partials


  // -find peak frequency by (parabolic) interpolation (or maybe cubic?)
  // -evaluate complex amplitude at peak-freq
  // -compute magnitude phase from interpolated complex amplitude and store in one of the partials

  // -start scanning forward to the spectrogram after the transient has passed, later extend the
  //  partials towards the start by scanning backward from the start position
  // -maybe later find (multiple) transients via the onset detector


  // references:
  // http://www.cerlsoundgroup.org/Loris/
  // https://ccrma.stanford.edu/~juan/ATS_manual.html
  // http://clam-project.org/


}
/*
      // how do we best compute the instantaneous phase - linear interpolation?
      int binInt  = peaks[i];
      T   binFrac = peakBin-binInt;
      if(binFrac < 0) {
        binInt  -= 1;
        binFrac  = 1-binFrac;
      }
      //T amp2 = (1-binFrac)*pMag[binInt] + binFrac*pMag[binInt+1]; // another way to compute amplitude - vary bad
      //std::complex<T> binCmpVal = (1-binFrac)*pCmp[binInt] + binFrac*pCmp[binInt+1]; 
      //T binMag = abs(binCmpVal);
      //T binPhs = arg(binCmpVal);
      // complex value at peak-bin
      T phs0 = pPhs[binInt];   // just for debug
      T phs1 = pPhs[binInt+1]; // dito
      peakPhase = rsInterpolateWrapped(pPhs[binInt], pPhs[binInt+1], binFrac, -PI, PI);
      // hmm...it doesn't seem to make sense to interpolate the phase between bind like that - it 
      // does not seem to be a smooth function - it jumps erratically between bins - maybe we need 
      // a completely different way to obtain the instantaneus phase...maybe by comparing values in
      // different frames? or maybe by correlating the function with a complex sine that has 
      // exactly the frequency? maybe the Goertzel algo can be used?:
      // https://en.wikipedia.org/wiki/Goertzel_algorithm
      // https://www.embedded.com/design/configurable-systems/4006427/A-DSP-algorithm-for-frequency-analysis
      // https://www.st.com/content/ccc/resource/technical/document/design_tip/group0/20/06/95/0b/c3/8d/4a/7b/DM00446805/files/DM00446805.pdf/jcr:content/translations/en.DM00446805.pdf
      // https://stackoverflow.com/questions/13499852/scipy-fourier-transform-of-a-few-selected-frequencies
      // https://stackoverflow.com/questions/13499852/scipy-fourier-transform-of-a-few-selected-frequencies
      // https://stackoverflow.com/questions/11579367/implementation-of-goertzel-algorithm-in-c
      //...if we go down that route, we may also use the obtained amplitude
      // for refinig our amplitude estimate...ok - for now, just leave the instantaneous phase 
      // measurement at zero

      https://ccrma.stanford.edu/~jos/sasp/Phase_Interpolation_Peak.html
*/


// rename to testSinusoidalSynthesis1
void sinusoidalModel1()
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
  model2 = analyzeSinusoidal(&x[0], (int)x.size(), fs);
  std::vector<double> y = synthesizeSinusoidal(model2, fs);

  rosic::writeToMonoWaveFile("SinusoidalSynthesisTest.wav", &x[0], (int)x.size(), (int)fs, 16);
  int dummy = 0;
}

void sinusoidalAnalysis1()
{
  double sampleRate = 48000;  // in Hz
  double frequency  = 100;    // in Hz
  double length     = 0.2;    // in seconds
  double startPhase = 0.0;    // in radians
  int    blockSize  = 1024;   // a bit larger than 2 cycles (1 is 480 samples with 100Hz@48kHz)
  int    hopSize    = 256;
  int    zeroPad    = 1;


  //// test
  //frequency = 5000;
  //length = 0.05; 

  // create signal:
  double period = sampleRate / frequency;    // in samples
  int N = (int)ceil(length * sampleRate);    // number of samples
  std::vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = sin(n * 2*PI*frequency/sampleRate + startPhase);

  // find model for the signal:
  rsSinusoidalModel<double> model = analyzeSinusoidal(&x[0], N, sampleRate);
    // we need to pass the desired blockSize, hopSize and zeroPadding factor to the function

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