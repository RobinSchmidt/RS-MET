#include "ResamplingExperiments.h"

void fadeOut()
{
  static const int N = 500;  // number of samples to plot
  double t[N], x[N], y[N];
  createTimeAxis(N, t, 1.0);
  rsFillWithValue( x, N, 1.0);
  rsCopyBuffer(x, y, N);
  rsFadeOut(y, 100, 300);
  plotData(N, t, x, y);

  //Plotter plt; plt.plotData(N, t, x, y);
}

void resampler()
{
  // As test signal, we create a pulse-wave at frequency "f" at samplerate "fs" with a duty cycle 
  // of "d" and amplitude "A". We create "xN" samples of that signal and the resample that signal 
  // with a ratio of "r".
  static const int xN = 21000; // number of samples in input signal
  static const int yN = xN;    // maximum number of samples in output signal
  double f  = 1000.0;          // frequency of test signal
  double fs = 44100.0;         // sample rate
  double d  = 0.4;             // duty cycle of pulse-wave
  double A  = 0.5;             // amplitude
  double r  = 2.1;             // resampling ratio
  double x[xN];                // input signal
  double y[yN];                // output signal

  // create the test signal and write it into a wavefile:
  synthesizePulseWave(x, xN, f, d, fs, 0.0, true);
  rsScale(x, xN, A);
  writeToMonoWaveFile("ResamplingInput.wav",  x, xN, (int) fs, 16);

  // resample it with various settings and write the results into wavefiles:
  rsResampler::transposeLinear(x, xN, y, yN, r);
  writeToMonoWaveFile("ResamplingOutputLinear.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 2, true);
  writeToMonoWaveFile("ResamplingOutputSinc2.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 4, true);
  writeToMonoWaveFile("ResamplingOutputSinc4.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 8, true);
  writeToMonoWaveFile("ResamplingOutputSinc8.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 16, true);
  writeToMonoWaveFile("ResamplingOutputSinc16.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 64, true);
  writeToMonoWaveFile("ResamplingOutputSinc64.wav",  y, yN, (int) fs, 16);

  rsResampler::transposeSinc(x, xN, y, yN, r, 512, true);
  writeToMonoWaveFile("ResamplingOutputSinc512.wav",  y, yN, (int) fs, 16);

  // actually, the simple cosin (non-squared) window seems to giev better results than
  // squared - do the comparison with a longer signal, so we can see the spectrum at higher 
  // resolution

  /*
  // plot signal:
  double t[xN];
  createTimeAxis(xN, t, fs);
  plotData(xN/50, t, x, y);
  */
}


void sincResamplerAliasing()
{
// todo: test sinc interpolation aliasing systematically - for example like this: 
// -create a sinewave at 11025 Hz (fs = 44100)
// -make a sweep from 1.0-5.0 (transpose factors)
// -where the sweep reaches 2, the output signal crosses the Nyquist frequency
// -plot the output signal against the sweep-factors (maybe dB scaled)
// -around 2, we should see a rapid drop of the output signal amplitude (because any signal beyond
//  2 is pure aliasing)
// -if the drop begins long before 2, we have too much passband attenuation
// -if it doesn't drop below a threshold, we have not enough stopband rejection
// -if it doesn't drop fast enough, we have a too wide transition band
// -passband ripple and -attenuation and stopband rejection depend on window-shape
// -transition width depends on window length (and window shape - but we can always compensate
//  the window shape effect by using a longer length)
// so, we should do this test with different window-shapes and lengths and find the sweet spot.
// maybe we should also try the Blackman-Nutall window because it has good stopband rejection 
// which is the most important feature for supressing the aliasing components

// this ist still flawed. we should really implement teh TimeWarp class first


  static const int xN = 2000;  // number of input samples
  static const int yN = 8000;  // maximum number of output samples

  double f  = 10000.0;         // input sinusoid frequency
  double fs = 44100.0;         // sample rate
  double A  = 0.5;             // input amplitude

  double *x = new double[xN];  // input signal
  double *y = new double[yN];  // output signal
  double *w = new double[yN];  // warping map

  // create the input signal and write it into a wavefile:
  synthesizeWaveform(x, xN, SINE, f, fs, 0.0, true);
  rsScale(x, xN, A);
  //writeToMonoWaveFile("SincResamplerAliasingInput.wav",  x, xN, (int) fs, 16);

  // create warping map for linear frequency sweep:
  double rL = 1.0;             // lower transposition ratio (at start)
  double rU = 5.0;             // upper transposition ratio (at end)
  rsFillWithRangeLinear(w, yN, rL, rU);
  rsCumulativeSum(w, w, yN); 



  // fixed sinc length of 64:
  rsTimeWarper::timeWarpSinc(x, xN, y, w, yN, 64, 1.0, true);
  //writeToMonoWaveFile("SincResamplerAliasing64.wav",  y, yN, (int) fs, 16);

  // recover frequency sweep data from time-warp map (as x-axis for plot)
  rsDifference(w, yN);


  //plotData(yN, w, y);


  delete[] x;
  delete[] y;
  delete[] w;
}


void sincResamplerModulation()
{
  static const int N = 2000;  // number of samples
  int    sincLength  = 10;    // length of the sinc function
  double ratio       = 2.95;  // resampling ratio
  bool   demodulate  = true;  // switch instantaneous DC normalization on/off

  // select window-function:
  FunctionPointer3DoublesToDouble wnd;
  wnd = rsRaisedCosineWindow;
  //wnd = rsExactBlackmanWindow;

  // select window-parameter (only relevant, if wnd = rsRaisedCosineWindow is used):
  double wp;
  //wp = 0.0;    // Hann window
  //wp = 0.08;   // Hamming window
  wp = 1.0;    // rectangular window

  double t[N];
  createTimeAxis(N,t, 1.0);

  // resample a DC signal:
  int n;
  double xDC[N];                // DC input signal 
  double yDC[N];                // resampled DC output signal
  rsFillWithValue(xDC, N, 1.0);
  rsFillWithZeros(yDC, N);
  double tr = 0.0;
  for(n = 0; n < N; n++)
  {
    yDC[n] = signalValueViaSincAt(xDC, N, tr, sincLength, ratio, wnd, wp, demodulate);
    tr += ratio;
  }

  // resample a sinusoid:
  double xSin[N], ySin[N];
  double w = 2*PI/400;
  for(n = 0; n < N; n++)
    xSin[n] = sin(w*n);
  tr = 0.0;
  for(n = 0; n < N; n++)
  {
    ySin[n] = signalValueViaSincAt(xSin, N, tr, sincLength, ratio, wnd, wp, demodulate);
    tr += ratio;
  }

  //plotData(N, t, xDC, yDC);
  plotData(N, t, xDC, yDC, xSin, ySin);
  // Observations:
  // When "demodulate" is set to false, the output shwos a modulation, the period of which 
  // (measured in samples) is given by 1/min(f, 1-f) where f is the fractional part of the 
  // resampling ratio [verify this formula]. This modulation artifact is most pronounced for the 
  // rectangular window, it gets better with Hann and Hamming window and still better with the 
  // Blackman window. Increasing the sinc-length has also the tendency to lessen the modulation, 
  // but it doesn't seem to decrease monotonically (ratio = 0.95, 20 is better than 10 but 40 is 
  // again worse but not as bad as 10). Additionally to the modulation, we make an overall 
  // amplitude error - the signal tends to gets more quiet for higher ratios. Both of these 
  // problems are solved by setting "demodulate" to true. This will lead to normalization by
  // the sum of the tap-weights, such that each instantaneous interpolation filter has a DC gain
  // of unity.

  int dummy = 0;
}


void sincResamplerPassbandRipple()
{
  static const int N = 10000;  // number of samples
  int    sincLength  = 30;    // length of the sinc function
  double ratio       = 0.71;  // resampling ratio
  bool   demodulate  = true;  // switch instantaneous DC normalization on/off

  // select window-function:
  FunctionPointer3DoublesToDouble wnd;
  //wnd = rsRaisedCosineWindow;
  wnd = rsExactBlackmanWindow;

  // select window-parameter (only relevant, if wnd = rsRaisedCosineWindow is used):
  double wp;
  wp = 0.0;    // Hann window
  //wp = 0.08;   // Hamming window
  //wp = 1.0;    // rectangular window

  double t[N];
  createTimeAxis(N,t, 1.0);

  // create a sine sweep:
  double x[N], y[N];
  double w[N];
  double pMax = 20;  // maximum period for sweep
  double pMin = 5;   // minimum period for sweep
  rsFillWithRangeLinear(w, N, 1/pMax, 1/pMin); // sweep 
  int n;
  for(n = 0; n < N; n++)
    x[n] = sin(w[n]*n);


  // resample the sine-sweep:
  double tr = 0.0;
  for(n = 0; n < N; n++)
  {
    y[n] = signalValueViaSincAt(x, N, tr, sincLength, ratio, wnd, wp, demodulate);
    tr += ratio;
  }

  //plotData(N, t, x, y);
  plotData(rsMin((int)floor(N/ratio),N), t, y);
  // Observations:
  // when using a rectangular window, we see a slight modulation of the amplitude of the output due
  // to the passband ripple of the interpolation filter. This artifact becomes invisible already 
  // with the Hann window. Theoretically, Hamming and even more so Blackman should give even better
  // results.


  int dummy = 0;
}


double windowedSincWeightSum(double sincLength, double stretch, double tf, 
  FunctionPointer3DoublesToDouble windowFunction = rsRaisedCosineWindow, 
  double windowParameter = 0.0)
{
  int L = (int) floor(sincLength);
  double ws = 0.0;      // sum of all weights
  double w;             // weight for a sample
  double ts;            // time instant with shift
  int mMin = -(L/2-1);  // minimum shift
  int mMax = +(L/2);    // maximum shift
  for(int m = mMin; m <= mMax; m++)
  {
    ts = m-tf;
    w = rsNormalizedSinc(ts/stretch) * windowFunction(ts, sincLength, windowParameter);
    ws += w;
  }
  return ws; 
  // ws approximately equals stretch so we return the "normalized" value instead:

  //return ws/stretch;
  
  // the approximation gets worse when the sincLength gets shorter and/or the stretch-value gets 
  // larger - this ration ripples around unity and at some points falls off - pretty much like
  // a lowpass-filter with passband ripple. the "passband" (i.e. the range of strecth-factors 
  // for which this ration wiggles around unity before it falls off) is smaller for shorter 
  // sinc-lengths. that means shorter sinc-lengths allow for less stretch before amplitude error 
  // becomes severe

  // these were just for some experimentation:
  //return 2*ws/sincLength;
  //return stretch*ws/sincLength;
  //return sincLength*ws/stretch;
}
void sincResamplerSumOfTapWeights()
{
  // we plot the sum of the weights for different (fixed) sinc lengths as function of the stretch 
  // factor:
  static const int N = 1000; // number of values to plot
  double sMax = 100.0;        // maximum stretch factor
  double tf   = 0.2;         // fractional part of time instant

  // select window-function:
  FunctionPointer3DoublesToDouble wnd;
  //wnd = rsRaisedCosineWindow;
  wnd = rsExactBlackmanWindow;

  // select window-parameter (only relevant, if wnd = rsRaisedCosineWindow is used):
  double wp;
  //wp = 0.0;    // Hann window
  //wp = 0.08;   // Hamming window
  wp = 1.0;    // rectangular window

  // array of stretch factors and resulting weight sums for different sinc lengths 
  // (8, 16, ...)
  double s[N], ws8[N], ws16[N], ws32[N], ws64[N], ws128[N];
  rsFillWithRangeLinear(s, N, 1.0, sMax);
  for(int i = 0; i < N; i++)
  {
    ws8[i]   = windowedSincWeightSum( 8,  s[i], tf, wnd, wp); // works up to s=2.38
    ws16[i]  = windowedSincWeightSum(16,  s[i], tf, wnd, wp); // up to s=4.78
    ws32[i]  = windowedSincWeightSum(32,  s[i], tf, wnd, wp); // 9.55
    ws64[i]  = windowedSincWeightSum(64,  s[i], tf, wnd, wp); // 19.0
    ws128[i] = windowedSincWeightSum(128, s[i], tf, wnd, wp); // 38.0
  }


  // plot the weight-sums as function of the stretch factor for various sinc-lengths:
  plotData(N, s, ws8, ws16, ws32, ws64, ws128);
  // Observations:
  // For low stretch factors, the weight-sum approximately equals the stretch factor but wiggles   
  // around it (is sometimes smaller, sometimes larger). As the stretch factor gets higher, the 
  // weight-sum saturates at the half of the sinc-length. 
  
  // plot the ratio of the weight-sum over the stretch-factor: ws[i]/s[i]
  rsDivide(ws8,   s, ws8,   N);
  rsDivide(ws16,  s, ws16,  N);
  rsDivide(ws32,  s, ws32,  N);
  rsDivide(ws64,  s, ws64,  N);
  rsDivide(ws128, s, ws128, N);
  plotData(N, s, ws8, ws16, ws32, ws64, ws128);
  // Observations:
  // As expected, the ratio first wiggles around unity until some cutoff point beyond which it 
  // falls off. The cutoff point depends on the sinc length and is higher for longer lengths. This
  // means, longer lengths support larger stretch amounts before they start to attenuate the output 
  // signal, so longer is better. Also, longer windows reduce the amount of the wiggle for a given 
  // fixed stretch factor (i.e. when looking at the wiggles at the left end of the plot) - the 
  // maximum deviation, more to the right of the plot where shorter lengths are already dropped 
  // off tends to grow with longer sinc-lengths (but not so much). The overall amount of the 
  // wiggles can be reduced by using an appropriate window function - it gets better in that order:
  // rectangular < Hann < Hamming < Blackman. The ripple of the Blackman window is only one tenth 
  // of the ripple of the Hamming window, but it falls off earlier than Hamming, so it should be 
  // used with slightly longer lengths (about 25% longer)
  
  // we measure the maximum of the wiggles of the 64-sample long sinc and express it in dB:
  double d64 = rsAmp2dB(rsMaxValue(ws64, N));
    // rect: 1.43, Hann: 0.11, Hamming: 0.034, Blackman: 0.0029829146655588447,
    // exact Blackman: 0.00084361968017353120 -> that's what we should use for the interpolator

  int dummy = 0;
}

void timeWarp()
{
  int    xN = 400000;          // number of samples in input signal
  int    yN = 88200;           // maximum number of samples in output signal
  double f  = 1000.0;          // frequency of test signal
  double fs = 44100.0;         // sample rate
  double d  = 0.4;             // duty cycle of pulse-wave
  double A  = 0.5;             // amplitude

  double *x = new double[xN];  // input signal
  double *y = new double[yN];  // output signal
  double *w = new double[yN];  // warping map

  // create the input signal and write it into a wavefile:
  synthesizePulseWave(x, xN, f, d, fs, 0.0, true);
  rsScale(x, xN, A);
  writeToMonoWaveFile("TimeWarpInput.wav",  x, xN, (int) fs, 16);

  // create warping map for linear frequency sweep:
  double fMin = 0.25;  // transposition minimum
  double fMax = 8.0;   // transposition maximum
  rsFillWithRangeLinear(w, yN, fMin, fMax);
  rsCumulativeSum(w, w, yN); 
    // make this function take an output buffer as second argument (which may or may not be equal
    // to the input buffer), write a function that does the inverse (first difference)
   
  /*
  // plot warping map:
  double *t = new double[yN];
  createTimeAxis(yN, t, 1.0);
  plotData(yN, t, w);
  delete[] t;
  */

  // without anti-aliasing, using a fixed sinc length of 64:
  rsTimeWarper::timeWarpSinc(x, xN, y, w, yN, 64, 1.0, false);
  writeToMonoWaveFile("TimeWarpOutputNoAntiAliasSinc64.wav",  y, yN, (int) fs, 16);

  // with anti-aliasing, using a fixed sinc length of 16:
  rsTimeWarper::timeWarpSinc(x, xN, y, w, yN, 16, 1.0, true);
  writeToMonoWaveFile("TimeWarpOutputSinc16.wav",  y, yN, (int) fs, 16);

  // with anti-aliasing, using a variable sinc length between 16 and 4.0*16=64:
  rsTimeWarper::timeWarpSinc(x, xN, y, w, yN, 16, 4.0, true);
  writeToMonoWaveFile("TimeWarpOutputSinc16-64.wav",  y, yN, (int) fs, 16);

  // with anti-aliasing, using a fixed sinc length of 64:
  rsTimeWarper::timeWarpSinc(x, xN, y, w, yN, 64, 1.0, true);
  writeToMonoWaveFile("TimeWarpOutputSinc64.wav",  y, yN, (int) fs, 16);

  delete[] x;
  delete[] y;
  delete[] w;
}

/*
std::vector<double> removePitchModulation(std::vector<double> input, double targetFrequency)
{
  double fs = 44100.0; // shoould be irrelevant
  
  input.

  //rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(x, fm, xN, fs, 20.0, 
  //  5000.0)
}
*/

// The function takes a reference to a pointer for the output signal y and will allocate memory for
// y and the return-value its number of samples in the y-array. The caller is responsible for 
// deleting the y-array eventually:
int removePitchModulation(double *x, int N, double *&y,  double targetFrequency)
{
  double fs = 44100.0;       // samplerate - irrelevant but required as argument
  double *f = new double[N]; // instantaneous frequency

  // Measure instantaneous frequency:
  rsInstantaneousFundamentalEstimator::
    measureInstantaneousFundamental(x, f, N, fs, 20.0, 5000.0);

  // Compute desired readout speed for each sample which is given by the ratio of the desired 
  // instantaneous frequency and the actual instantaneous frequency of the input signal (re-use the
  // f-array for this)
  for(int n = 0; n < N; n++)
    f[n] = targetFrequency / f[n];

  // Compute required length for output, allocate memory and compute output:
  int Ny = rsTimeWarper::getPitchModulatedLength(f, N);
  y = new double[Ny]; 
  rsTimeWarper::applyPitchModulation(x, f, N, y, 16.0, 4.0, true); 

  // Clean up and return the length of the output:
  delete[] f;
  return Ny;
}

void pitchDemodulation()
{
  // As test input, we create a sinewave that sweeps linearly from some start frequency to some 
  // end frequency and also has some vibrato in it. This swept and vibrato'ed sine wave will then  
  // have its pitch-modulation removed by time-varying resampling.

  int    xN = 88200;    // number of samples in input signal
  double fs = 44100.0;  // sample rate
  double f1 = 430.0;    // start value for input signal frequency
  double f2 = 470.0;    // end value for input signal frequency
  double vf =   7.0;    // vibrato frequency in Hz
  double vd =  40.0;    // vibrato depth in Hz (as frequency deviation)
  double fy = 440.0;    // target value for output frequency
  double a  =   0.5;    // signal amplitude

  double *f  = new double[xN]; // instantaneous frequency of input signal
  double *fm = new double[xN]; // measured instantaneous frequency
  double *rl = new double[xN]; // reliability of measurements
  double *x  = new double[xN]; // input signal
  double *r  = new double[xN]; // readout-ratios



  // create the array with instantaneous frequencies of input signal:
  rsFillWithRangeLinear(f, xN, f1, f2);
  int n;
  for(n = 0; n < xN; n++)
    f[n] += 0.5 * vd * sin((2*PI*vf/fs)*n);

  // create input signal:
  createSineWave(x, xN, f, a, fs);
  writeToMonoWaveFile("PitchDemodulationInput.wav",  x, xN, (int) fs, 16);

  // compute desired readout speed for each sample which is given by the ratio of the desired 
  // instantaneous frequency and the actual instantaneous frequency of the input signal:
  for(n = 0; n < xN; n++)
    r[n] = fy / f[n];

  // compute the output signal:
  int yN = rsTimeWarper::getPitchModulatedLength(r, xN);
  double *y = new double[yN]; 
  rsTimeWarper::applyPitchModulation(x, r, xN, y, 16.0, 4.0, true);  
  writeToMonoWaveFile("PitchDemodulationOutputKnownF0.wav",  y, yN, (int) fs, 16);


  // measure the instantaneous frequency and do the same thing with measured values:
  rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(x, fm, xN, fs, 20.0, 
    5000.0, rl);
  for(n = 0; n < xN; n++)
    r[n] = fy / fm[n];  
  int yN2 = rsTimeWarper::getPitchModulatedLength(r, xN);
  double *y2 = new double[yN2]; 
  rsTimeWarper::applyPitchModulation(x, r, xN, y2, 16.0, 4.0, true);  
  writeToMonoWaveFile("PitchDemodulationOutputMeasuredF0.wav",  y2, yN2, (int) fs, 16);

  // now, the same, but using the convenience function removePitchModulation
  double *y3 = nullptr;
  int yN3 = removePitchModulation(x, xN, y3, fy);
  writeToMonoWaveFile("PitchDemodulationOutputMeasuredF0-2.wav", y3, yN3, (int) fs, 16);
  delete[] y3;


  plotData(xN, 0, 1/fs, f, fm); // plot actual and measured fundamental
  //plotData(xN, 0, 1/fs, rl);    // plot reliability

  /*
  double *t = new double[xN]; // time axis for plot
  //createTimeAxis(xN, t, fs);
  createTimeAxis(xN, t, 1);
  //plotData(5000, t, wi, w);
  //plotData(xN, t, fx);
  //plotData(xN, t, x);
  delete[] t;
  */

  // todo: instead of just removing existing pitch modulation, we could generalize this approach 
  // to apply a new pitch envelope which may or may not be constant 

  delete[] f;
  delete[] fm;
  delete[] x;
  delete[] y;
  delete[] y2;
  delete[] r;
}

void phaseLockedCrossfade()
{
  // We create two sinusoids with different frequencies and different vibrato rates and apply a 
  // crossfade between them where, during the crossfade, the sinusoids are resampled in a way so
  // as to ensure a phase-match between them at all time-instants. 

  static const int N1 = 20000, N2 = 15000; // length of 1st and 2nd signal
  double fs  = 8000;                       // samplerate
  double fc1 = 140, fc2 = 260;             // center frequencies
  double a1  = 1.0, a2  = 0.5;             // amplitudes
  double fv1 = 6,   fv2 = 6;               // vibrato frequencies
  double vd1 = 80,  vd2 = 80;              // vibrato depths
  double x1[N1], x2[N2];                   // 1st and 2nd signal
  double f1[N1], f2[N2];                   // instantaneous frequencies of x1, x2
  int n;                                   // sample index

  // create instantaneous frequency arrays:
  for(n = 0; n < N1; n++)
    f1[n] = fc1 + 0.5 * vd1 * sin(n * 2*PI*fv1/fs);
  for(n = 0; n < N2; n++)
    f2[n] = fc2 + 0.5 * vd2 * sin(n * 2*PI*fv2/fs);

  // create signals:
  createSineWave(x1, N1, f1, a1, fs);
  createSineWave(x2, N2, f2, a2, fs);

  // target frequency for flattened signals:
  double ft = 0.5 * (fc1 + fc2); 

  // create and set up the rsPhaseLockedCrossfader object:
  rsPhaseLockedCrossfader plc;
  plc.setInputs(x1, f1, N1, x2, f2, N2, ft);

  // retrieve flattened signals:
  vector<double> x1f = plc.getFlattenedSignal1();
  vector<double> x2f = plc.getFlattenedSignal2();

  // Having obtained the pitch-flattened signals, it is now possible to set up the crossfade start
  // and end as well as an overall time-shift of the flattened x2f with respect to the flattened 
  // x1f. All values are in samples but not necessarily integers:
  double start = 2000;       // crossfade start in flattened signal x1
  double end   = 12000;      // crossfade end in flattened signal x1
  double shift = 1000;       // time-shift of flattened x2 with respect to flattened x1

  // Having chosen our crossfade section in the flattened time-domain, we set it up in the 
  // crossfader object and retrieve the crossfaded signal:
  plc.setFlattenedCrossfade(start, end, shift);
  vector<double> y = plc.getOutput();

  // write wavefiles:
  writeToMonoWaveFile("PhaseLockCrossfadeInput1.wav", x1, N1, (int) fs, 16);
  writeToMonoWaveFile("PhaseLockCrossfadeInput2.wav", x2, N2, (int) fs, 16);
  writeToMonoWaveFile("PhaseLockCrossfadeOutput.wav", &y[0], y.size(), (int) fs, 16);

  // plot crossfaded signal:
  GNUPlotter plt;
  //plt.addDataArrays(x1f.size(), &x1f[0]);
  //plt.addDataArrays(x2f.size(), &x2f[0]);
  plt.addDataArrays(y.size(), &y[0]);
  //plt.addDataArrays(pmc.t1.size(), &pmc.t1[0]);
  //plt.addDataArrays(pmc.t2.size(), &pmc.t2[0]);

  plt.plot();

  // Observations:
  // When both signals have different vibrato-speeds, the vibrato in the output during the
  // crossfade is messed up. It would be desirable to see the vibrato speed just sweep between
  // both input speeds during the crossfade - i don't know, if that's possible, though. It's
  // all a question of the exact formula used in computeReadoutTimes - at the moment, this is
  // double txw = (1-t)*tx1 + t*tx2; 
}

//template<class T>
//void rsCrossfade(T *x1, size_t N1, T *x2, size_t N2, T *y, size_t start, size_t end, 
//  size_t shift = 0)
//{
//  rsAssert(end >= start);
//  rsAssert(end < N1);
//  rsAssert(shift <= start);
//  rsAssert(end+shift < N2);
//  size_t Ny = N2+shift;
//}

// Applies a linear crossfade to x1, x2. The start, end and shift values are all with respect to 
// the time axis of x1. The shift parameter can be used to shift x2 with respect to x1.
vector<double> linearCrossfade(vector<double> x1, vector<double> x2, int start, int end, 
  int shift)
{
  int N1 = x1.size();                    // length of x1
  int N2 = x2.size();                    // length of x2
  int Ny = N2+shift;                     // length of output
  vector<double> y(Ny);
  int n;
  for(n = 0; n < start; n++)             // leading section
    y[n] = x1[n]; 
  double scl = 1.0 / (end-start+1);
  for(n = start; n < end+1; n++)         // crossfade section
  {
    double c = scl * (n-start+1);
    y[n] = (1-c)*x1[n] + c*x2[n-shift];
  }
  for(n = end+1; n < Ny; n++)            // trailing section
  {
    int n2 = n - shift;
    y[n] = x2[n-shift];
  }
  return y;
}

// Given two speed arrays which are assumed to be applicable to two pitch-flattened signals y1, y2
// in order to unflatten them (giving back the original non-flat signals x1, x2), this function 
// modifies the speed arrays in a way so as to when applied to y1, y2, you get two signals z1, z2
// which can be crossfaded and will be phase-locked during this crossfade. You also have to pass
// the crossfade parameters (in terms of the flattened time-axis of y1). All parameters are 
// pointers because they all may potentially be modified. That means, your crossfade start, end and
// shift values will also be updated.
void crossfadeUnflatteningSpeeds(vector<double>* speeds1, vector<double> *speeds2, int *start, 
  int *end, int *shift)
{
  vector<double> s = linearCrossfade(*speeds1, *speeds2, *start, *end, *shift);
  rsVariableSpeedPlayer vsp;
  vsp.setInputAndSpeed(nullptr, &s[0], s.size());
  for(int n = 0; n < speeds1->size(); n++)
    (*speeds1)[n] = 1 / (vsp.warpTime(n+1) - vsp.warpTime(n));
  for(int n = 0; n < speeds2->size(); n++)
    (*speeds2)[n] = 1 / (vsp.warpTime(n+1+*shift) - vsp.warpTime(n+*shift));
  *start = (int) vsp.warpTime(*start);
  *end   = (int) vsp.warpTime(*end);
  *shift = (int) vsp.warpTime(*shift);
  if(*shift > *start)
    *shift = *start;  // shift should be <= start
}

void phaseLockedCrossfade2()
{
  static const int N1 = 2000, N2 = 1500;   // length of 1st and 2nd signal
  double fs  = 1000;                       // samplerate
  double fc1 = 20,  fc2 = 30;              // center frequencies
  double a1  = 1.0, a2  = 0.5;             // amplitudes
  double fv1 = 4,   fv2 = 6;               // vibrato frequencies
  double vd1 = 10,  vd2 = 15;              // vibrato depths
  vector<double> x1(N1), x2(N2);           // 1st and 2nd signal
  vector<double> f1(N1), f2(N2);           // instantaneous frequencies of x1, x2
  int n;                                   // sample index                  

  // create instantaneous frequency arrays:
  for(n = 0; n < N1; n++)
    f1[n] = fc1 + 0.5 * vd1 * sin(n * 2*PI*fv1/fs);
  for(n = 0; n < N2; n++)
    f2[n] = fc2 + 0.5 * vd2 * sin(n * 2*PI*fv2/fs);

  // create signals:
  createSineWave(&x1[0], N1, &f1[0], a1, fs);
  createSineWave(&x2[0], N2, &f2[0], a2, fs);

  // target frequency for flattened signals:
  double ft = 0.5 * (fc1 + fc2); 

  // create read-speed arrays from frequency arrays:
  vector<double> s1(N1), s2(N2);
  for(n = 0; n < N1; n++)
    s1[n] = ft / f1[n];
  for(n = 0; n < N2; n++)
    s2[n] = ft / f2[n];

  // create pitch-flattened signals:
  rsVariableSpeedPlayer vsp1, vsp2;
  vsp1.setInputAndSpeed(&x1[0], &s1[0], N1);
  vsp2.setInputAndSpeed(&x2[0], &s2[0], N2);
  vector<double> y1 = vsp1.getOutput();
  vector<double> y2 = vsp2.getOutput();

  // create the inverse speed-arrays (which turn y1, y2 back into x1, x2):
  vector<double> s1i = rsVariableSpeedPlayer::invertSpeeds(s1);
  vector<double> s2i = rsVariableSpeedPlayer::invertSpeeds(s2);

  // create flattened crossfaded signal y:
  int start = 400;   // crossfade start in y1
  int end   = 1000;  // crossfade end in y1
  int shift = 200;   // time-shift of y2 with respect to y1
   
  // create flat crossfaded signal
  vector<double> y = linearCrossfade(y1, y2, start, end, shift);          

  // unflatten the crossfaded signal:
  vector<double> si = linearCrossfade(s1i, s2i, start, end, shift);
  vector<double> yu = rsVariableSpeedPlayer::applyPlaybackSpeed(y, si);

  // obtain unflattening arrays to be applied to y1, y2 separately, pre-crossfade, apply them 
  // (giving z1, z2) and then do the crossfade post-unflattening.
  crossfadeUnflatteningSpeeds(&s1i, &s2i, &start, &end, &shift); // call modifies all arguments
  vector<double> z1 = rsVariableSpeedPlayer::applyPlaybackSpeed(y1, s1i);
  vector<double> z2 = rsVariableSpeedPlayer::applyPlaybackSpeed(y2, s2i);
  vector<double> z  = linearCrossfade(z1, z2, start, end, shift);

  //// write wavefiles:
  //writeToMonoWaveFile("PhaseLockCrossfadeInput1.wav",     &x1[0], N1, (int) fs, 16);
  //writeToMonoWaveFile("PhaseLockCrossfadeInput2.wav",     &x2[0], N2, (int) fs, 16);
  writeToMonoWaveFile("PhaseLockCrossfadeOutputPost.wav", &yu[0], yu.size(), (int) fs, 16);
  writeToMonoWaveFile("PhaseLockCrossfadeOutputPre.wav",  &z[0],  z.size(),  (int) fs, 16);

  // plot:
  GNUPlotter plt;

  // flattened signals:
  //plt.addDataArrays((int)y1.size(), &y1[0]);
  //plt.addDataArrays((int)y2.size(), &y2[0]);

  //// flattened crossfaded signal:
  //plt.addDataArrays((int)y.size(), &y[0]);

  //// unflattening speed arrays:
  //plt.addDataArrays((int)s1i.size(), &s1i[0]);
  //plt.addDataArrays((int)s2i.size(), &s2i[0]);
  //plt.addDataArrays((int)si.size(),  &si[0]);

  // unflattened crossfade output:
  plt.addDataArrays((int)yu.size(), &yu[0]);
  plt.addDataArrays((int)z.size(),  &z[0]);
  plt.plot();
  // The two arrays yu and z should be roughly equal, but possibly not exactly due to the 
  // rounding step that is necessary
}

//void phaseLockedCrossfade2()
//{
//  static const int N1 = 2000, N2 = 1500;   // length of 1st and 2nd signal
//  double fs  = 1000;                       // samplerate
//  double fc1 = 20,  fc2 = 30;              // center frequencies
//  double a1  = 1.0, a2  = 0.5;             // amplitudes
//  double fv1 = 4,   fv2 = 6;               // vibrato frequencies
//  double vd1 = 10,  vd2 = 15;              // vibrato depths
//
//  vector<double> x1(N1), x2(N2);           // 1st and 2nd signal
//  vector<double> f1(N1), f2(N2);           // instantaneous frequencies of x1, x2
//  int n;                                   // sample index                  
//                                           
//  // create instantaneous frequency arrays:
//  for(n = 0; n < N1; n++)
//    f1[n] = fc1 + 0.5 * vd1 * sin(n * 2*PI*fv1/fs);
//  for(n = 0; n < N2; n++)
//    f2[n] = fc2 + 0.5 * vd2 * sin(n * 2*PI*fv2/fs);
//
//  // create signals:
//  createSineWave(&x1[0], N1, &f1[0], a1, fs);
//  createSineWave(&x2[0], N2, &f2[0], a2, fs);
//
//  // target frequency for flattened signals:
//  double ft = 0.5 * (fc1 + fc2); 
//
//  // create read-speed arrays from frequency arrays:
//  vector<double> s1(N1), s2(N2);
//  for(n = 0; n < N1; n++)
//    s1[n] = ft / f1[n];
//  for(n = 0; n < N2; n++)
//    s2[n] = ft / f2[n];
//
//  // create pitch-flattened signals:
//  rsVariableSpeedPlayer vsp1, vsp2;
//  vsp1.setInputAndSpeed(&x1[0], &s1[0], N1);
//  vsp2.setInputAndSpeed(&x2[0], &s2[0], N2);
//  vector<double> y1 = vsp1.getOutput();
//  vector<double> y2 = vsp2.getOutput();
//
//  // create the inverse speed-arrays (which turn y1, y2 back into x1, x2):
//  vector<double> s1i = rsVariableSpeedPlayer::invertSpeeds(s1);
//  vector<double> s2i = rsVariableSpeedPlayer::invertSpeeds(s2);
//
//  // create flattened crossfaded signal y:
//  int start = 400;   // crossfade start in y1
//  int end   = 1000;  // crossfade end in y1
//  int shift = 200;   // time-shift of y2 with respect to y1
//
//  // create flat crossfaded signal
//  vector<double> y = linearCrossfade(y1, y2, start, end, shift);          
//
//  // unflatten the crossfaded signal:
//  vector<double> si = linearCrossfade(s1i, s2i, start, end, shift);      // unflattening speed array
//  vector<double> yu = rsVariableSpeedPlayer::applyPlaybackSpeed(y, si);  // unflattened crossfade signal
//
//  // From the unflattening speed array si (to be applied to the crossfade output y), generate two 
//  // separate unflattening arrays si1, si2, to be applied to y1, y2 pre-crossfade. This also 
//  // requires readjustment of the crossfade parameters (start, end, shift):
//  rsVariableSpeedPlayer vsp;
//  vsp.setInputAndSpeed(nullptr, &si[0], si.size());
//  double start2 = vsp.warpTime(start);
//  double end2   = vsp.warpTime(end);
//  double shift2 = vsp.warpTime(shift);
//  vector<double> si1 = s1i, si2 = s2i;
//  for(n = 0; n < si1.size(); n++)
//    si1[n] = 1 / (vsp.warpTime(n+1) - vsp.warpTime(n));
//  for(n = 0; n < si2.size(); n++)
//    si2[n] = 1 / (vsp.warpTime(n+1+shift) - vsp.warpTime(n+shift));
//
//  // apply the two unflattening arrays to the flattened signals:
//  vector<double> z1 = rsVariableSpeedPlayer::applyPlaybackSpeed(y1, si1);
//  vector<double> z2 = rsVariableSpeedPlayer::applyPlaybackSpeed(y2, si2);
//
//  // obtain an unflattened crossfaded signal z by applying a crossfade to z1, z2:
//  vector<double> z = linearCrossfade(z1, z2, (int)start2, (int)end2, (int)shift2);
//
//  int dummy = 0;
//
//  //// write wavefiles:
//  //writeToMonoWaveFile("PhaseLockCrossfadeInput1.wav", &x1[0], N1, (int) fs, 16);
//  //writeToMonoWaveFile("PhaseLockCrossfadeInput2.wav", &x2[0], N2, (int) fs, 16);
//  //writeToMonoWaveFile("PhaseLockCrossfadeOutput.wav", &yu[0], yu.size(), (int) fs, 16);
//
//  // plot
//  GNUPlotter plt;
//
//  // flattened signals:
//  //plt.addDataArrays((int)y1.size(), &y1[0]);
//  //plt.addDataArrays((int)y2.size(), &y2[0]);
//
//  //// flattened crossfaded signal:
//  //plt.addDataArrays((int)y.size(), &y[0]);
//
//  //// unflattening speed arrays:
//  //plt.addDataArrays((int)s1i.size(), &s1i[0]);
//  //plt.addDataArrays((int)s2i.size(), &s2i[0]);
//  ////plt.addDataArrays((int)si.size(),  &si[0]);
//  //plt.addDataArrays((int)si1.size(), &si1[0]);
//  //plt.addDataArrays((int)si2.size(), &si2[0]);
//
//
//  //// unflattening warping maps:
//  //plt.addDataArrays((int)m1.size(), &m1[0]);
//  //plt.addDataArrays((int)m2.size(), &m2[0]);
//  //plt.addDataArrays((int)m.size(),  &m[0]);
//
//  // unflattened crossfade output:
//  plt.addDataArrays((int)yu.size(), &yu[0]);
//  plt.addDataArrays((int)z.size(), &z[0]);
//
//  plt.plot();
//
//  // Observations:
//  // yu is the unflattened output of the crossfade of the two flattened signals y1, y2
//  // z is the crossfaded output of z1, z2 where z1, and z2 are separately unflattened versions
//  // of y1, y2
//  // The two arrays yu and z should be roughly equal, but possibly not exactly due to the 
//  // rounding step that is necessary
//
//  dummy = 0;
//}

//void phaseLockedCrossfade2()
//{
//  static const int N1 = 20000, N2 = 15000; // length of 1st and 2nd signal
//  double fs  = 8000;                       // samplerate
//
//
//  //double fc1 = 200, fc2 = 200;             // center frequencies
//  //double a1  = 1.0, a2  = 0.5;             // amplitudes
//  //double fv1 = 8,   fv2 = 8;               // vibrato frequencies
//  //double vd1 = 80,  vd2 = 40;              // vibrato depths
//
//  double fc1 = 150, fc2 = 250;             // center frequencies
//  double a1  = 1.0, a2  = 0.5;             // amplitudes
//  double fv1 = 8,   fv2 = 6;               // vibrato frequencies
//  double vd1 = 80,  vd2 = 60;              // vibrato depths
//
//  vector<double> x1(N1), x2(N2);           // 1st and 2nd signal
//  vector<double> f1(N1), f2(N2);           // instantaneous frequencies of x1, x2
//  int n;                                   // sample index
//
//  // create instantaneous frequency arrays:
//  for(n = 0; n < N1; n++)
//    f1[n] = fc1 + 0.5 * vd1 * sin(n * 2*PI*fv1/fs);
//  for(n = 0; n < N2; n++)
//    f2[n] = fc2 + 0.5 * vd2 * sin(n * 2*PI*fv2/fs);
//
//  // create signals:
//  createSineWave(&x1[0], N1, &f1[0], a1, fs);
//  createSineWave(&x2[0], N2, &f2[0], a2, fs);
//
//  // target frequency for flattened signals:
//  double ft = 0.5 * (fc1 + fc2); 
//
//  // up to here, the code is (almost) duplicated from phaseLockedCrossfade
//  // maybe factor out a function...
//
//  // create read-speed arrays from frequency arrays:
//  vector<double> s1(N1), s2(N2);
//  for(n = 0; n < N1; n++)
//    s1[n] = ft / f1[n];
//  for(n = 0; n < N2; n++)
//    s2[n] = ft / f2[n];
//
//  // create pitch-flattened signals:
//  rsVariableSpeedPlayer vsp1, vsp2;
//  vsp1.setInputAndSpeed(&x1[0], &s1[0], N1);
//  vsp2.setInputAndSpeed(&x2[0], &s2[0], N2);
//  vector<double> y1 = vsp1.getOutput();
//  vector<double> y2 = vsp2.getOutput();
//
//  // unflatten the flattened signals again (to check, if roundtrip works):
//  vector<double> s1i = rsVariableSpeedPlayer::invertSpeeds(s1);
//  vector<double> s2i = rsVariableSpeedPlayer::invertSpeeds(s2);
//  vector<double> xx1 = rsVariableSpeedPlayer::applyPlaybackSpeed(y1, s1i);
//  vector<double> xx2 = rsVariableSpeedPlayer::applyPlaybackSpeed(y2, s2i);
//
//  // create flattened crossfaded signal y:
//  int start = 4000;                  // crossfade start in y1
//  int end   = 12000;                 // crossfade end in y1
//  int shift = 2000;                  // time-shift of y2 with respect to y1
//
//  // test:
//  //start = 0;
//  shift = 0;
//
//  vector<double> y = crossfade(y1,  y2, start, end, shift); // flat crossfade signal
//
//  // unflatten the crossfaded signal:
//  vector<double> si = crossfade(s1i, s2i, start, end, shift);            // unflattening speed array
//  vector<double> yu = rsVariableSpeedPlayer::applyPlaybackSpeed(y, si);  // unflattened crossfade signal
//
//  // obtain two separate unflattening speed arrays from si:
//  vector<double> si1 = s1i, si2 = s2i;
//  for(n = start; n <= end; n++)
//    si1[n] = si[n];
//  for(n = start; n <= end; n++)
//    si2[n-shift] = si[n];            // verify the shift
//
//  // apply the two unflattening arrays to the flattened signals:
//  vector<double> z1 = rsVariableSpeedPlayer::applyPlaybackSpeed(y1, si1);
//  vector<double> z2 = rsVariableSpeedPlayer::applyPlaybackSpeed(y2, si2);
//
//  rsVariableSpeedPlayer vsp;
//  vsp.setInputAndSpeed(nullptr, &si[0], si.size());
//  int start2 = vsp.warpTime(start);
//  int end2   = vsp.warpTime(end);
//  int shift2 = vsp.warpTime(shift);
//  //int start2 = vsp.unwarpTime(start);
//  //int end2   = vsp.unwarpTime(end);
//  //int shift2 = vsp.unwarpTime(shift);
//    // check, if these truncations (from double to int) may lead to problems...
//
//
//  // obtain an unflattened crossfaded signal zu by applying a crossfade to z1, z2:
//  //vector<double> zu = crossfade(z1, z2, start, end, shift);
//  vector<double> zu = crossfade(z1, z2, start2, end2, shift2);
//    // (maybe, we will have to modify the start, end, shift values)
//    // zu should be (almost) equal to yu, a little error may occur due to rounding of the warped
//    // start, end, shift values ...or maybe only the shift has to be modified by the 
//    // time-warping map corresponding to si1?
//    // ...try with shift 0 - nope - makes no difference - i think, we have to modify all 3 values
//    // possibly according to the warping map corresponding to...
//
//
//  //// create the unflattening warping map for y:
//  //vector<double> m1 = vsp1.getTimeWarpMapYX();  // same length as y1
//  //vector<double> m2 = vsp2.getTimeWarpMapYX();  // same length as y2
//  //vector<double> m(Ny);
//  //for(n = 0; n < start; n++)       // leading section
//  //  m[n] = m1[n]; 
//  //for(n = start; n < end; n++)     // crossfade section
//  //{
//  //  double c = (n-start) / (end-start-1.0);
//  //  m[n] = (1-c)*m1[n] + c*m2[n-shift];
//  //}
//  //for(n = end; n < Ny; n++)        // trailing section
//  //  m[n] = m2[n-shift];
//
//
//  //// write wavefiles:
//  //writeToMonoWaveFile("PhaseLockCrossfadeInput1.wav", &x1[0], N1, (int) fs, 16);
//  //writeToMonoWaveFile("PhaseLockCrossfadeInput2.wav", &x2[0], N2, (int) fs, 16);
//  //writeToMonoWaveFile("PhaseLockCrossfadeOutput.wav", &yu[0], yu.size(), (int) fs, 16);
//
//
//  // plot
//  GNUPlotter plt;
//
//  // flattened signals:
//  //plt.addDataArrays((int)y1.size(), &y1[0]);
//  //plt.addDataArrays((int)y2.size(), &y2[0]);
//
//  //// original and reconstructed input 1:
//  //plt.addDataArrays((int)x1.size(),  &x1[0]);
//  //plt.addDataArrays((int)xx1.size(), &xx1[0]);
//
//  //// flattened crossfaded signal:
//  //plt.addDataArrays((int)y.size(), &y[0]);
//
//  //// unflattening speed arrays:
//  ////plt.addDataArrays((int)s1i.size(), &s1i[0]);
//  ////plt.addDataArrays((int)s2i.size(), &s2i[0]);
//  ////plt.addDataArrays((int)si.size(),  &si[0]);
//  //plt.addDataArrays((int)si1.size(), &si1[0]);
//  //plt.addDataArrays((int)si2.size(), &si2[0]);
//
//
//  //// unflattening warping maps:
//  //plt.addDataArrays((int)m1.size(), &m1[0]);
//  //plt.addDataArrays((int)m2.size(), &m2[0]);
//  //plt.addDataArrays((int)m.size(),  &m[0]);
//
//  // unflattened corssfade output:
//  plt.addDataArrays((int)yu.size(), &yu[0]);
//  plt.addDataArrays((int)zu.size(), &zu[0]);
//
//  plt.plot();
//
//  // Observations: 
//  // xx1, xx2 are each 1 sample shorter than x1, x2 respectively - the warp/unwarp roundtrip 
//  // looses the last sample
//  // at the start and end of xx1, xx2, there's some ripple due to the sinc-interpolator
//  //
//
//  int dummy = 0;
//}


void sineShift()
{
  // input signal parameters:
  static const int N = 3000;  // number of samples
  double fs = 44100.0;        // sample rate
  double f  = 100.0;          // sinusoid frequency
  double a  = 0.5;            // amplitude
  double p  = 2.5;            // start-phase

  // user parameters:
  int    n0 = 1000;           // sample-instant, where target phase should be effective
  double p0 = PI/2;           // target phase at n0


  // create input signal:
  double x[N];
  int n;
  double w = 2*PI*f/fs;         // (true) normalized radian frequency
  for(n = 0; n < N; n++)
    x[n] = a * sin(w*n + p);

  // compute desired shift amount:
  double dn = rsSineShiftAmount(x, N, n0, p0);

  // create shifted signal using a rounded integer shift value:
  double yi[N];
  rsCopyBuffer(x, yi, N);
  rsShift(yi, N, (int) rsRound(dn));

  // create shifted signal using the exact shift value and sinc-interpolation:
  double y[N];
  rsCopyBuffer(x, y, N);

  rsResampler::shiftSinc(y, y, N, dn, 64.0);

  // plot:
  plotData(N, 0, 1, x, yi, y); // input and output signals with integer and noninteger shift


  int dummy = 0;
}

void sineShift2()
{
  // user parameters:
  static const int numSamples   = 3000;  // number of samples
  static const int numHarmonics = 5;     // number of harmonics
  double fs = 44100.0;                   // sample rate
  double f0 = 100.0;                     // fundamental frequency
  double f[numHarmonics];                // harmonic frequencies
  double a[numHarmonics];                // harmonic amplitudes (randomized)
  double p[numHarmonics];                // harmonic start-phases (randomized)
  double aMin = 0.1;                     // minimum amplitude for a harmonic
  double aMax = 0.2;                     // maximum amplitude for a harmonic
  double p0   = 1.57;                    // target phase (0: odd symmetry, PI/2: even symmetry)
  int    n0   = 1000;                    // sample, where we want to see target phase pt
  int    randomSeed = 2;                 // seed for random number generator

  // allocate 2D-array for the harmonics:
  double x[numHarmonics][numSamples];    // input harmonics
  double y[numHarmonics][numSamples];    // output harmonics

  // fill arrays of frequencies, amplitudes, start-phases and create input harmonics (with 
  // randomized amplitudes and initial phases)
  rsRandomUniform(0.0, 1.0, randomSeed);
  int h, n;                              // loop indices for harmonics and samples
  for(h = 0; h < numHarmonics; h++)
  {
    f[h] = f0 * (h+1);
    a[h] = rsRandomUniform(aMin, aMax);
    p[h] = rsRandomUniform(0.0, 2*PI); 
    for(n = 0; n < numSamples; n++)
    {
      double w = 2*PI*f[h]/fs;            // normalized radian frequency of current sinusoid
      x[h][n]  = a[h] * sin(w*n + p[h]);
    }
  }

  // shift-input signals and store results in y:
  for(int h = 0; h < numHarmonics; h++)
  {
    double shiftAmount = rsSineShiftAmount(x[h], numSamples, n0, p0); 
    rsResampler::shiftSinc(x[h], y[h], numSamples, shiftAmount, 64.0); 
  }

  // mix the (original and shifted) harmonics together and write the results into wavefiles:
  double xMix[numSamples], yMix[numSamples];
  rsFillWithZeros(xMix, numSamples);
  rsFillWithZeros(yMix, numSamples);
  for(h = 0; h < numHarmonics; h++)
  {
    rsAdd(xMix, x[h], xMix, numSamples);
    rsAdd(yMix, y[h], yMix, numSamples);
  }
  writeToMonoWaveFile("SineShiftInputMix.wav",  xMix, numSamples, (int) fs, 16);
  writeToMonoWaveFile("SineShiftOutputMix.wav", yMix, numSamples, (int) fs, 16);

  // plot:
  //plotData(numSamples, 0, 1, x[0], x[1], x[2], x[3], x[4]); // input harmonics
  //plotData(numSamples, 0, 1, y[0], y[1], y[2], y[3], y[4]); // output harmonics
  //plotData(numSamples, 0, 1, xMix, yMix);
}

void pitchDetectWithSilence()
{
  int N = 10000;
  double fs = 44100;

  // not yet implemented


  int dummy = 0;
}



// x: input, y: output, N: number of samples, f0: fundamental frequency (in Hz), bw: bandwidth for 
// harmonic splitter (in Hz), fs: samplerate (in Hz), n0: adjustment sample index, 
// tp: target-phase, nh: number of harmonics:
void rsHarmonicPhaseAdjust(double *x, double *y, int N, double f0, double bw, double fs, int n0, 
  double tp, int nh)
{
  if( nh == -1 )
    nh = rsFloorInt(0.5*fs/f0);   // all harmonics
  double f;                       // frequency of current harmonic
  double shift;                   // amount of shift for the harmonic in samples
  double *tmp = new double[N];    // memory for one harmonic

  // extract, shift and accumulate one harmonic at a time:
  rsFillWithZeros(y, N);
  for(int h = 0; h < nh; h++)
  {
    f = (h+1)*f0;
    rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, tmp, N, f, bw, fs, 10);
    shift = rsSineShiftAmount(tmp, N, n0, tp, 2*PI*f/fs);
    rsResampler::shiftSinc(tmp, tmp, N, shift, 64);
    rsAdd(y, tmp, y, N);
  }

  delete[] tmp;
}

void pitchDetectA3()
{
  //const char *path = "../../TestInputs/Sustain`n=A3`dyn=loud`register=1.wav";
  const char *path = "../../TestInputs/SustainA3Cut.wav";

  int N, fs;
  double* x = readMonoWaveFile(path, N, fs);
  double* f = new double[N];

  rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(
    x, f, N, fs, 20, 5000, nullptr);

  GNUPlotter plt;
  plt.addDataArrays(N, f);
  plt.plot();

  delete[] x;
  delete[] f;
}

void phaseLockSaxophone()
{
  const char *path = "../../TestInputs/autotuned_saxaphone_D#2_01.wav";
    // the sample is actually at D#3 - just a different convention?

  int numChannels;
  int N;
  int n0;          // alignment position
  double f0;       // fundamental frequency
  double bw;
  int fs;

  double **px = readFromWaveFile(path, numChannels, N, fs);
  double *x = new double[N];
  rsCopyBuffer(px[0], x, N);

  n0 = rsMaxAbsIndex(x, N); // corresponds to sAlignPos = FindHighestAmplitudeInSignal(data[0]);
                            // it finds n0 = 1872
  f0 = rsPitchToFreq(51);   // frequency for D#3 (MIDI-key 51)
  bw = 0.7*f0;

  double *y = new double[N];  // output signal
  rsHarmonicPhaseAdjust(x, y, N, f0, bw, fs, n0, PI/2, -1);

  writeToMonoWaveFile("PhaseLockedSaxophone.wav",  y, N, fs, 16);

  //plotData(10000, 0, 1.0/fs, x, y);
  //plotData(10000, 0, 1.0/fs, y);


  // cleanup:
  for(int i = 0; i < numChannels; i++)
    delete px[i];
  delete[] px;
  delete[] x;
  delete[] y;
}


void phaseLockSaxophone2()
{
  const char *path = "../../TestInputs/autotuned_saxaphone_D#2_01.wav";
    // the sample is actually at D#3 - just a different convention?

  int numChannels;
  int N;
  int nh = 5;       // number of harmonics
  int n0;           // alignment position
  double f0;        // fundamental frequency
  double f;
  double bw;
  double tp = PI/2; // target phase
  double shift;
  int fs;
  int i;

  double **px = readFromWaveFile(path, numChannels, N, fs);
  double *x = new double[N];
  rsCopyBuffer(px[0], x, N);

  n0 = rsMaxAbsIndex(x, N); // corresponds to sAlignPos = FindHighestAmplitudeInSignal(data[0]);
                            // it finds n0 = 1872
  f0 = rsPitchToFreq(51);   // frequency for D#3 (MIDI-key 51)
  bw = 0.7*f0;

  double **y = new double*[nh];
  for(i = 0; i < nh; i++)
    y[i] = new double[N];

  for(int h = 0; h < nh; h++)
  {
    f = (h+1)*f0;
    rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y[h], N, f, bw, fs, 10);
    shift = rsSineShiftAmount(y[h], N, n0, tp, 2*PI*f/fs);
    rsResampler::shiftSinc(y[h], y[h], N, shift, 64);
  }
  
  plotData(3000, 0, 1.0, y[0], y[1], y[2], y[3], y[4]);

  // cleanup:
  for(i = 0; i < numChannels; i++)
    delete px[i];
  delete[] px;
  delete[] x;
  for(i = 0; i < nh; i++)
    delete[] y[i];
  delete[] y;
}

void autoTuneHorn()
{
  const char *path = "../../TestInputs/Horn_Muted_F1.wav";
  int numChannels;
  int N;
  int fs;
  double **px = readFromWaveFile(path, numChannels, N, fs);
  double *x = new double[N];
  rsCopyBuffer(px[0], x, N);

  double f0 = 87.0; // target frequency
  int i, n;

  double *f, *rl, *r;
  f   = new double[N];
  rl  = new double[N];
  r   = new double[N];

  // measure the instantaneous frequency:
  rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(x, f, N, fs, 20.0, 
    5000.0, rl);

  // tune the original signal:
  for(n = 0; n < N; n++)
    r[n] = f0 / f[n];  
  int Ny = rsTimeWarper::getPitchModulatedLength(r, N);
  double *y = new double[Ny]; 
  rsTimeWarper::applyPitchModulation(x, r, N, y, 16.0, 4.0, true);  
  writeToMonoWaveFile("HornAutotuned.wav",  y, Ny, (int) fs, 16);

  // filter out fundamental from original:
  double bw = 0.7*f0;
  double *xf = new double[N];
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, xf, N, f0, bw, fs, 10);
  //writeToMonoWaveFile("HornFundamental.wav",  xf, N, (int) fs, 16);

  // write original and extracted fundamental into file:
  writeToStereoWaveFile("HornAndFundamental.wav", x, xf, N, fs, 16);


  // cleanup:
  for(i = 0; i < numChannels; i++)
    delete px[i];
  delete[] px;
  delete[] x;
  delete[] f;
  delete[] rl;
  delete[] r;
  delete[] y;
  delete[] xf;
}

void autoTuneHorn2()
{
  int numChannels;
  int N;
  int fs;
  int n;

  // read full horn sample:
  double **px = readFromWaveFile("../../TestInputs/Horn_F1.wav", numChannels, N, fs);
  double *xFull = new double[N];
  rsCopyBuffer(px[0], xFull, N);
  for(n = 0; n < numChannels; n++)
    delete px[n];
  delete[] px;

  // read chunk of 1st harmonic of horn sample:
  px = readFromWaveFile("../../TestInputs/Horn_F1_H1_Chunk.wav", numChannels, N, fs);
  double *xH1 = new double[N];
  rsCopyBuffer(px[0], xH1, N);
  for(n = 0; n < numChannels; n++)
    delete px[n];
  delete[] px;

   // todo: use readMonoWaveFile - get rid of that code duplication

  // measure the instantaneous frequency using the chunk - we use only the nonzero part of the 
  // signal which is the range from nStart to nEnd, the length of our frequency array is 
  // nEnd-nStart+1:
  int nStart = 130801;        // this is where the non-silence part starts
  int nEnd   = N-1;           // this is where the non-silence part ends
  int L      = nEnd-nStart+1; // length of the non-silence part
  double *f, *rl, *r;
  f   = new double[L];
  rl  = new double[L];
  rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(xH1+nStart, f, L, fs, 20.0, 
    5000.0, rl);

  // create array of readout speed factors:
  double f0     = 87.0;         // target frequency
  double thresh = 0.85;         // reliability threshold for tuning
  r = new double[N];
  int m;
  for(n = 0; n < nStart; n++)   // use unit speed for the part before our chunk
    r[n] = 1.0;
  for(n = n; n <= nEnd; n++)    // now we are into the part where the analysis signal is nonzero
  {
    m = n - nStart;
    if(rl[m] > thresh)
      r[n] = f0/f[m];           // use the detected frequency only if reliability is good enough
    else
      r[n] = 1.0;
  }
  for(n = n; n < N; n++)        // use unit speed for the part after our chunk
    r[n] = 1.0;

  // tune the full signal using the readout-speed array that we have just created:
  int Ny = rsTimeWarper::getPitchModulatedLength(r, N);
  double *y = new double[Ny]; 
  rsTimeWarper::applyPitchModulation(xFull, r, N, y, 16.0, 4.0, true);  
  writeToMonoWaveFile("HornAutotunedViaH1Chunk.wav",  y, Ny, (int) fs, 16);

  delete[] xFull;
  delete[] xH1;
  delete[] f;
  delete[] rl;
  delete[] r;
  delete[] y;
}

void sylophoneCycleMarks()
{
  int N, fs;
  char* path = "../../TestInputs/Sylophone/Sustain`n=A2`rr=3-shortChunk.wav";
  double* x = readMonoWaveFile(path, N, fs);

  // find cycle marks by different algorithms:
  rsCycleMarkFinder cmf(fs, 20, 5000);
  vector<double> cm1, cm2;
  cmf.setAlgorithm(rsCycleMarkFinder::F0_ZERO_CROSSINGS); 
  cm1 = cmf.findCycleMarks(&x[0], N);
  cmf.setAlgorithm(rsCycleMarkFinder::CYCLE_CORRELATION); 
  cm2 = cmf.findCycleMarks(&x[0], N);

  // plot signal and cycle marks:
  int Nz = cm1.size();       // # cycle marks
  vector<double> cmy(Nz);    // y values for plotting (all zero)
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(Nz, &cm1[0], &cmy[0]);
  plt.addDataArrays(Nz, &cm2[0], &cmy[0]);
  plt.setGraphStyles("lines", "points", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
  delete[] x;
}

void autoTuneSylophone()
{
  char* path = "../../TestInputs/Sylophone/Sustain`n=A2`rr=3.wav";

  int N, fs;
  double* x = readMonoWaveFile(path, N, fs);

  //int mode = rsInstantaneousFundamentalEstimator::F0_ZERO_CROSSINGS;
  //int mode = rsInstantaneousFundamentalEstimator::CYCLE_CORRELATION;

  std::vector<double> f1(N), f2(N);
  rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(x, &f1[0], N, fs, 20, 5000,
    nullptr, 0);
  //rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(x, &f2[0], N, fs, 20, 5000,
  //  nullptr, 1);

  // plot:
  int start  = 75000;
  int length = 35000;
  GNUPlotter plt;
  //plt.addDataArrays(length, &x[start]);
  plt.addDataArrays(length, &f1[start]);
  //plt.addDataArrays(length, &f2[start]);
  plt.setPixelSize(800, 200);
  plt.setRange(0, length, 216.0, 217.0);
  plt.plot();
  delete[] x;
}

void bestMatchShift()
{
  int numChannels;
  int N1, N2;        // lengths of the 2 signals
  double *x1, *x2;   // the 2 signals
  int fs;            // samplerate
  int n;
            
  // read signals:
  x1 = readMonoWaveFile("../../TestInputs/autocorrelate_test-001.wav", N1, fs);
  x2 = readMonoWaveFile("../../TestInputs/autocorrelate_test-002.wav", N2, fs);

  double shift = rsGetShiftForBestMatch(x1, x2, rsMin(N1, N2));

  delete[] x1;
  delete[] x2;
}
