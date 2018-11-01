#include "PartialExtractionExperiments.h"

/*
idea for analysis algorithm for modal modeling:
-(maybe) pre-process by seprating excitation from resonator via linear prediction
-take a big chunk from the sample and FFT it (maybe N = 65536) with suitable window, possibly, the
 window should be only one-sided (right) to preserve information in the transient
-find spectral peaks:
 -estimate their frequencies and amplitudes by quadratic interpolation on dB values
 -maybe sort them according to ascending frequency (peak picking will probably sort with descending 
  magnitudes - but that might not be so bad either)
 -remember the globally maximal magnitude: mMax
-isolate and analyze one mode at a time:
 -apply a partial isolation filter for the given mode frequency (see below for details on this 
  filter)
 -refine the estimated mode frequency fM by counting periods (upward zero crossings) of the 
  filtered signal and divide by the time-interval
  -the first and last zero crossings should be determined with subsample precision via cubic 
   interpolation to get a more accurate (noninteger) time-interval
 -since our mode isolation filter may have nonunity gain at fM due to passband ripple, compute 
  the actual gain at fM and multiply signal with reciprocal to ensure unity gain
 -apply 1st order allpass with 90� frequency adjusted to fM to obtain a quadrature signal
  (maybe in reverse -> experiment, if this gives an advantage)
 -extract the instantaneous envelope as sqrt(y*y + q*q) where y is the isolated mode and q is the
  quadrature component
 -optionally smooth the instantaneous envelope with a 1st order bi-directional lowpass adjusted to
  a cutoff frequency of fM/s, where s is some smoothing parameter, maybe s = 4
 -to this extracted instantaneous mode envelope, we fit our model envelope (least-squares fit of
  a (constrained) sum of exponentials)
 -the phase of the mode can be obtained by looking at the first zero-crossing (or at some other
  zero crossing later)

---------------------------------------------------------------------------------------------------

The Partial Isolation Filter

The main challenge for this filter is to isolate one partial in the frequency domain without 
smearing its amplitude envelope too much in the time domain. Filters with very high Q and/or great 
steepness will be able to isolate the partial better in the frequency domain at the expense of 
introducing more time smearing due to longer ringing. In order to preserve the phase of the 
partial, the filter should be applied bi-directionally (...which should have the side-effect of 
symmetrizing the time-smearing of the partial's amplitude envelope (?) - but this is not really 
what is observed - but maybe this is because the envelope itslef is asymmetrical)

The filter may be applied multiple times. With biquad bandpass filters, it has been found that the
time-smearing can be reduced by using multiple passes with lower Q instead of a single pass with 
high Q.

Maybe the time smearing (shifting of the envelope peak to the right) can be explained in terms of
group delay - but this should be zero in the bidirectional case (or not?). Or - more plausibly - 
the shifting of the peak can be related to the asymmetry of the input signal's envelope. Perhaps 
its useful to consider the convolution of an asymmetric signal like our the enveloped partial with 
our symmetric, zero-centered impulse response of the bidirectional filter. It's intuitive that the 
output is again asymmetric, skewed to where the input signal has more "weight", which in our case 
is the right side.

Idea for elliptic partial isolation filter
 -design an elliptic bandpass filter with its lower cutoff between the mode frequency fM and its 
  left neighbour and upper cutoff between fM and the right neighbour
  -the desired stopband rejection should be chosen such as to suppress all other modes to a certain 
   amount (like 40dB maybe), we need to take into account ratio of the mode amplitude mM and the 
   global maximum mMax
  -we may tolerate some passpand ripple too, since we will compensate for nonunity gain at the 
   mode frequency later
  -the exact cutoff points and the transition width should be chosen such that the neighbouring 
   modes are already in the stopband, ideally we would want to start the stopbands at the 
   neighbouring partials -> see, if elliptic filters can eb designed in terms of their stopband
   frequencies
  -maybe halfway between the frequencies could be used as cutoff points which may be defined as 
   arithmetic or geometric mean - expermient
  -there's a tradeoff to be made between the filter's ability to isolate the mode in the frequency 
   domain and the filters impulse response length - steeper filters will ring longer and therefore
   smear the isolated mode's envelope in the time domain - therefore, choose reasonable stopband
   rejections - experimentation required
 -filter the signal with the obtained mode-isolation filter bi-directionally (to preserve phase)
  -the bi-directional filtering will square the magnitude response of the filter - that should be
   taken into account in the design step before
 -more ideally, we would like to design the elliptic filter in terms of its 2 stopband frequencies 
  which should be taken as lower and upper neighbour mode frequencies, the transtion witdh should
  ensure that fM falls into the passband (with some margin because we don't know the frequency 
  exactly, so we need an error tolerance)


Experiments:

partialExtractionBell():
I tried to run a complex-phasor based implementation of a two-pole one-zero resonator over the 
signal to directly get the instantaneous envelope as the radius of the phasor - this was not very 
promising.

partialExtractionTriple():
To assess the performance of the partial isolation filter, we create a signal containing 3 
frequencies (for example: fL=900, fM=1000, fU=1150 Hz) all with unit amplitude and some kind 
attack/decay envelope and try to recover the partial at the middle frequency. To make it difficult 
for the filter, the middle partial should have a short attack time. We then look at the isolated 
partial's envelope to see how well the peak of the measured envelope matches the actual envelope of 
the partial. Typically, the attack will be smeared in the sense that the peak's location is 
right-shifted and its height is decreased. The smaller this effect, the better.

partialExtractionSample():
Feeds some real-world samples into the partial extraction filter and plots the resulting extracted
partial along with its amplitude envelope.


To assess the performance when dealing with a real world signal, we may do time-domain subtraction
of the isolated partial from the whole sound and see how much the partial gets reduced in the 
spectrum and how little the other partial are affected. Because the transient is important for us,
we should use either a one-sided (right-sided) window or no window at all for the FFT.



Experimental results:

partialExtractionTriple():

The measured envelope generally shows two kinds of artifacts compared to the true envelope: 
(1) time-shifting of the peak to the right (with resulting lowered peak height - the height is 
given by the partial's actual amplitude there) and (2) modulation at the difference frequencies 
fM-fL and fU-fM [check this]

Biquad bandpass
-multi-pass, low-Q filters smear envelope less compared to single-pass high-Q filter
 -the multi-pass application leads to a repeated convolution of the filter's impulse response with
  itself, such that the amplitude envelope of the impulse response approaches a Gaussian shape 
  [verify this]. The Gaussian shape has the minimum product of time-width and frequency-width and
  might be optimal therefore and rather independent from the filter design method [verify this].
  The more passes we run, the less filter-Q do we need in order to obtain a specified supression
  of other partials - which is kinda trivial. However, running more passes with lower-Q (adjusted 
  such that the supression in dB of a neighbour partial is always the same, say -50dB), has the
  effect of contracting the impulse response to a narrower time interval. the neighbour partial
  supression seems a relevant parameter since it determines how much our measured envelope will be
  contaminated by modulation artifacts.
-zero-padding the signal before multi-pass, bi-directional filtering (to take into account the 
 filter's ringing, avoiding edge effects) leads to small improvement for the envelope's peak 
 location (in the test with the partial-triple at 1kHz, improvement was around 10%)
-using unidirectional instead of bidirectional filtering just seems to shift the whole signal left
 (in case of backward running) or right (in case of forward running) but does not really seem
 to help to locate the peak more exactly
-in order to get an envelope that is not contaminated too much with modulation artifacts from other
 partials, it seems to be necessarry to supress all other partials to a level at least -50dB below 
 the partial under investigation. we should probably not supress them much further because that 
 would sacrifice too much time-resolution (the steeper, the more smear)


Experiments to do:
-instead of prepending zeros, try to prepend a reversed (and possibly negated) copy of the signal 
 to put more "weight" to the left side
-for the lowest mode, it might not be the best choice to pass 0 as lower partial neighbour 
 frequency - investigate that with the piano_E2 sample, or create an example sound
-investigate, if a deterministic function can be found that maps the real peak location/height
 to the measured values - this function depends on the number of passes and the filter shape, but 
 once we have settled for a particular filter, we can measure this function and apply its inverse
 to our peak measurements, it may also depend on the decay time of the partial -> investigate
 ->the actual height might possibly be found by extrapolating the decay-envelope to the actual 
   location
-try Butterworth/Papoulis/Elliptic/... bandpasses with lower and upper cutoff adjusted somewhere 
 between lower (or upper) neighbour partial and the middle partial to be isolated (maybe a weighted
 geometric mean could be used)
 -check, how increasing the order behaves compared to increasing the number of passes
-try a combination of elliptic filtering for supressing direct neighbours and Butterworth or 
 Papoulis filtering for supressing more distant partials
-try a combination of bandpass filtering with notch-filtering at other prominent partials
-maybe for each extracted partial, we could do time-domain subtraction before extracting the next
 partial. this may ease the work of subsequent filter. disadvantage: we can't extract just one
 partial indpendently from the others (but that might not matter much)



Extensions for Future Work:

To extract partials with time-varying frequency (assuming its instantaneous frequency to be known
at all samples) we need a time-varying filter, that tracks the instantaneous frequency. Now, time 
variation of a filter's parameters may introduce artifacts, the nature of which depends on the 
implementation structure of the filter. It is desirable to have some idea about what these 
artifacts are. Maybe a filter based on a spiraling complex phasor might be best suited for time 
variation: when we change the frequency, it will just continue to rotate with the new angular 
increment and when we change the damping, it will just continue to rotate but with a different
multiplier for the radius. ...but maybe somehow we have to account for phase response.

We need to find out, if this zero-phase property of bi-directional filtering does (in some sense) 
continue to hold for our time-varying filter.

To investigate behavior under time-variation, we could feed an impulse and at some point within
the filter's impulse response, change a parameter. If we change the damping, we should just see 
that the exponential envelope continues with a different time-constant. If we change the frequency,
we should just see that change in frequency. If we change the (start)phase, we should see a 
discontinuity in the waveform, reflecting a phase jump.

Experiments to do:
-tune filter to 1000Hz with a decay time-constant of 0.2s, feed an impulse
 (1) at 0.2s, change the decay to 0.1s
 (2) at 0.2s, change the decay to 0.4s
 (3) at 0.2s, change the frequency to 500Hz
 (3) at 0.2s, change the frequency to 2000Hz
 -look, if we see what we want to see
 -repeat these experiments with feeding sines of 500, 1000 and 2000 Hz
 -maybe investigate what a change of the phase does


---------------------------------------------------------------------------------------------------

Amplitude Envelope Extraction

Having isolated a partial with frequency fM, we estimate its amplitude envelope by obtaining a 90�
phase-shifted ("quadrature") version of it by allpass filtering and then use sqrt(x*x+q*q) as 
instantaneous envelope where x is the partial signal and q is the quadrature component. 


Experiments to do:
-try to apply the allpass in reverse direction, see if this makes a difference
-

---------------------------------------------------------------------------------------------------

Resynthesis:

-we may decide to neglect modes if they are determined to be inaudible because of masking effects


*/


// to check the analysis algorithm: synthesize a sound with a couple of modes and try to retrieve
// the parameters, for example 110, 130, 170, 190, 230 Hz (difficult)
// 100, 250, 320, 540, 720, 830







// todo: rename to freeArray2D, templatize, move to RSLib:
void freeSampleData(double **sampleData, int numChannels)
{
  for(int i = 0; i < numChannels; i++)
    delete[] sampleData[i];
  delete sampleData;
}

// convenience function that can be used to plot up to 5 different (audio) signals in one plot, N is the 
// number of samples, fs is the sample rate - the time axis will be scaled in seconds
// - move to Common
void plotSignals(int N, double fs, double *y1, double *y2 = nullptr, double *y3 = nullptr, 
  double *y4 = nullptr, double *y5 = nullptr)
{
  double *t = new double[N]; 
  createTimeAxis(N, t, fs);          
  plotData(N, t, y1, y2, y3, y4, y5);
  delete[] t;
}


void biDirectionalFilter()
{
  // tests the rsBiDirectionalFilter::applyConstPeakBandpassBwInHz function by applying the 
  // bandpass to an impulse and plotting the impulse response or magnitude response

  //static const int N = 8192;  // number of samples (power of two for FFT)
  static const int N = RS_POW2_13;  // number of samples (power of two for FFT)
  double fs = 44100;          // samplerate
  double fc = 1000;           // center frequency
  double bw = 200;            // bandwidth in Hz
  //double gc = RS_SQRT2_INV;   // bandedge gain
  double gc = 0.5;            // bandedge gain
  int order = 3;              // (prototype) order of single pass filter

  bool lowpass = true;      // if true, we apply a lowpass, otherwise a bandpass

  // create input signal with centered impulse and filter it (with different values for the number
  // of passes):
  double x[N], y1[N], y2[N], y4[N], y8[N], y16[N];
  RAPT::rsArray::fillWithZeros(x, N);
  x[N/2] = 1.0;
  //x[0] = 1.0;  // for padding test

  if( lowpass == true )
  {
    rsBiDirectionalFilter::applyButterworthLowpass(x, y1,   N, fc, fs, order, 1,  gc);
    rsBiDirectionalFilter::applyButterworthLowpass(x, y2,   N, fc, fs, order, 2,  gc);
    rsBiDirectionalFilter::applyButterworthLowpass(x, y4,   N, fc, fs, order, 4,  gc);
    rsBiDirectionalFilter::applyButterworthLowpass(x, y8,   N, fc, fs, order, 8,  gc);
    rsBiDirectionalFilter::applyButterworthLowpass(x, y16,  N, fc, fs, order, 16, gc);

    /*
    rsBiDirectionalFilter::applyLowpass(x, y1,  N, fc, fs, 1);
    rsBiDirectionalFilter::applyLowpass(x, y2,  N, fc, fs, 2);
    rsBiDirectionalFilter::applyLowpass(x, y4,  N, fc, fs, 4);
    rsBiDirectionalFilter::applyLowpass(x, y8,  N, fc, fs, 8);
    rsBiDirectionalFilter::applyLowpass(x, y16, N, fc, fs, 16);
    */
  }
  else
  {
    rsBiDirectionalFilter::applyButterworthBandpassBwInHz(x, y1,  N, fc, bw, fs, order, 1,  gc);
    rsBiDirectionalFilter::applyButterworthBandpassBwInHz(x, y2,  N, fc, bw, fs, order, 2,  gc);
    rsBiDirectionalFilter::applyButterworthBandpassBwInHz(x, y4,  N, fc, bw, fs, order, 4,  gc);
    rsBiDirectionalFilter::applyButterworthBandpassBwInHz(x, y8,  N, fc, bw, fs, order, 8,  gc);
    rsBiDirectionalFilter::applyButterworthBandpassBwInHz(x, y16, N, fc, bw, fs, order, 16, gc);

    // old:
    //rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y1,  N, fc, bw, fs, 1);
    //rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y2,  N, fc, bw, fs, 2);
    //rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y4,  N, fc, bw, fs, 4);
    //rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y8,  N, fc, bw, fs, 8);
    //rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y16, N, fc, bw, fs, 16);
  }

  // compute magnitude responses:
  double m1[N/2], m2[N/2], m4[N/2], m8[N/2], m16[N/2];
  rsMagnitudeAndPhase(y1,  N, m1);
  rsMagnitudeAndPhase(y2,  N, m2);
  rsMagnitudeAndPhase(y4,  N, m4);
  rsMagnitudeAndPhase(y8,  N, m8);
  rsMagnitudeAndPhase(y16, N, m16);
  /*
  for(int n = 0; n < N; n++)
  {
    m1[n]  = rsAmp2dB(m1[n]);
    m2[n]  = rsAmp2dB(m2[n]);
    m4[n]  = rsAmp2dB(m4[n]);
    m8[n]  = rsAmp2dB(m8[n]);
    m16[n] = rsAmp2dB(m16[n]);
  }
  */

  // plot: 
  //plotData(N, 0.0, 1/fs, y1, y2, y4, y8, y16); // impulse responses


  double freqs[N];
  RAPT::rsArray::fillWithRangeLinear(freqs, N, 0.0, fs);
  //plotDataLogX(N/2, freqs, m1, m2, m4, m8, m16);
  //plotData(N/2, 0.0, fs/N, m1, m2, m4, m8, m16); // magnitude responses
  plotData(N/16, 0.0, fs/N, m1, m2, m4, m8, m16); // magnitude responses - shouldn't it be
                                                  // fs/(N+1)?

  // normalize and write output files:
  double s = 1 / RAPT::rsArray::maxAbs(y1, N);
  RAPT::rsArray::scale(y1,  N, s);
  RAPT::rsArray::scale(y2,  N, s);
  RAPT::rsArray::scale(y4,  N, s);
  RAPT::rsArray::scale(y8,  N, s);
  RAPT::rsArray::scale(y16, N, s);
  writeToMonoWaveFile("BiDirectionalFilterPasses1.wav",    y1, N, (int) fs, 16);
  writeToMonoWaveFile("BiDirectionalFilterPasses2.wav",    y2, N, (int) fs, 16);
  writeToMonoWaveFile("BiDirectionalFilterPasses4.wav",    y4, N, (int) fs, 16);
  writeToMonoWaveFile("BiDirectionalFilterPasses8.wav",    y8, N, (int) fs, 16);
  writeToMonoWaveFile("BiDirectionalFilterPasses16.wav",  y16, N, (int) fs, 16);


  // Observations
  // when the center frequency goes up, the magnitude-squared curves do not cross exactly at the 
  // 0.5 point anymore. this is because the bandwidth-scaling does not take into account bilinear
  // transform frequency warping (its exact only for analog prototype filters).
  // It doesn't seem to make much sense to apply Butterworth filters with orders higher than 1 
  // multiple times because the ringing seems to get worse with more passes. Only for filters of 
  // order 1, the multipass application reduces the ringing - maybe this related to the bandwidth
  // scaling? ...This seems to depend on the bandedge gain - for 1/sqrt(2), it's is true, but for
  // 0.5 it seems to be advantageous to go higher with order

  int dummy = 0;
}

void envelopeDeBeating()
{
  // We create two attack/decay sinusoids with frequencies close to each other such that the 
  // resulting amplitude envelope will show beating effects at the difference frequency. Then, we
  // try to remove the amplitude modulation from the envelope.

  double fs = 44100; // sample rate
  double f1 = 140;   // frequency 1
  double f2 = 150;   //           2
  double a1 = 0.5;   // amplitude 1
  double a2 = 0.5;   //           2
  double d1 = 0.2;   // decay time 1
  double d2 = 0.3;   //            2
  double fc = 10;   // smoothing lowpass cutoff
  int np = 8;        // smoothing lowpass order/number of passes
  int N = 30000;     // number of samples


  // If amplitudes and decay-times of both sinusoids are the same, the beating is most extreme

  // set up a modal filter bank to produce the output:
  RAPT::rsModalFilterBank<double, double> mfb;
  mfb.setSampleRate(fs);
  mfb.setReferenceAttack(0.2);
  mfb.setReferenceFrequency(f1);
  mfb.setModalParameters({ 1, f2/f1 }, { a1, a2 }, { d1, d2 }, { d1, d2 }, { 0, 0 });


  typedef std::vector<double> Vec;

  // create signal and estimate envelope:
  Vec x(N), env(N);
  getImpulseResponse(mfb, &x[0], N);
  //rsSineEnvelopeViaQuadrature(&x[0], &env[0], N, (f1+f2)/2, fs, 2.0);
   // todo: maybe use instantaneous envelope algorithm


  RAPT::rsEnvelopeExtractor<double>::sineEnvelopeWithDeBeating(&x[0], N, &env[0]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &env[0]);
  //plt.addDataArrays((int)envTime.size(),  &envTime[0],  &envValue[0]);
  //plt.addDataArrays((int)envTime2.size(), &envTime2[0], &envValue2[0]);
  plt.plot();

  // Observations:
  // f1 = 140, f2 = 150, a1 = 0.5, a2 = 0.5, d1 = 0.2, d2 = 0.3:
  //  -natural spline without smoothing, startAtZero=1, endAtZero=0, gives reasonable result
  // f1 = 140, f2 = 145, a1 = 0.5, a2 = 0.5, d1 = 0.2, d2 = 0.3:
  //  -natural spline shows strong overshoot
  //  -linear looks reasonable but should probably use smoothing
  //  -
 

  // de-beat:
  // ideas:
  // -connect peaks (by lines or splines)
  // -filter out beating frequencies (requires to know/estimate them first)
  // -use a variation of the "True-Envelope" algorithm applied to the time-domain signal
  // -actually, what we want is an envelope of an envelope - maybe we can appyl the same algorithm
  //  twice

  // Idea for general envelope extraction:
  // -take absolute value
  // -pick peaks - but only those that satisfy additional constraints - for example:
  //  -ratio of peak vs local average must be above some threshold (local average can be obtained
  //   by bidirectional lowpass)
  //  -distance between two neighbouring peaks must be above some threshold
  // -connect peaks by line or spline segments
  //  -maybe the lines/splines should be applied to the log of the envelope and after that, the
  //   the log should be undone via exp (or - more generally - use some monotonic nonliinear 
  //   transformation and its inverse)..maybe take y = log(1+x) and x = exp(y)-1 to avoid
  //   log-of-zero problems

  // perhaps the cubic interpolation is no good idea due to overshoot - instead use linear and 
  // apply a bidirectional Bessel or Gaussian filter to the result

  // the log/exp stuff is actually not good because the attack-phase should be inverted exp in case
  // of exp....so perhaps better to just stick to simple linear
}

void sineRecreation()
{
  // Creates a sine wave with an amplitude envelope and tries to retrieve that envelope and 
  // recreate sine from the envelope (with possibly different frequency and startphase)

  static const int N = 44100;   // number of samples
  double fx     = 100.0;        // input sine frequency
  double fy     = 120.0;        // output sine frequency
  double fs     = 44100.0;      // samplerate
  double px     = PI/2;         // start phase for x
  double py     = PI/4;         // start phase for y
  double smooth = 1.5;          // smoothing for envelope
  int n;                        // sample index

  double *x   = new double[N];  // input signal
  double *xq  = new double[N];  // quadrature part of input signal
  double *et  = new double[N];  // true envelope of x
  double *em  = new double[N];  // measured envelope of x
  double *em2 = new double[N];  // measured envelope of x (2nd algorithm)
  double *y   = new double[N];  // output signal

  // create the amplitude envelope:
  rsBreakpointModulatorD bm;
  bm.setSampleRate(fs);
  bm.initialize();
  bm.modifyBreakpoint(0, 0.0, 0.0, rsModBreakpoint<double>::ANALOG, 0.5);
  bm.modifyBreakpoint(1, 1.0, 0.0, rsModBreakpoint<double>::ANALOG, 0.5);
  bm.insertBreakpoint(0.1, 1.0, rsModBreakpoint<double>::ANALOG, 0.5);
  bm.insertBreakpoint(0.5, 0.5, rsModBreakpoint<double>::ANALOG, 0.5);
  bm.noteOn();
  for(n = 0; n < N; n++)
    et[n] = bm.getSample();

  // create the eneveloped sine:
  double w = 2*PI*fx/fs;
  double p = px;
  for(n = 0; n < N; n++)
  {
    x[n] = et[n] * sin(p);
    p += w;
  }

  // obtain quardrature part and envelope:
  rsSineQuadraturePart(x, xq, N, fx, fs, true);
  rsSineEnvelopeViaQuadrature(x, em,  N, fx, fs, smooth);
  rsSineEnvelopeViaAmpFormula(x, em2, N, fx, fs, smooth);  // other algorithm

  // recreate the sine (with possibly different frequency and phase):
  rsRecreateSine(x, y, N, fx, fy, fs, py, smooth);

  writeToMonoWaveFile("InputSine.wav",  x, N, (int) fs, 16);
  writeToMonoWaveFile("OutputSine.wav", y, N, (int) fs, 16);

  // plot:
  //plotData(N, 0, 1/fs, x, et);  // signal and true envelope
  //plotData(N, 0, 1/fs, x, xq);  // signal and quadrature part
  //plotData(N, 0, 1/fs, x, xq, em);  // signal, quadrature part and measured envelope
  //plotData(N, 0, 1/fs, x, em);  // signal and measured envelope
  //plotData(N, 0, 1/fs, et, em, em2);  // true and measured envelope
  plotData(N, 0, 1/fs, x, y);  // input and output sine

  // Observations:
  // with smooth 0.0, we see a little modulation in the measured envelope, with smooth = 1.0, it's
  // supressed reasonably but still visible, with smooth = 2.0, the envelope is very smooth indeed but 
  // might be already oversmoothed. Maybe 1.5 is the sweet spot
  // The algorithm based on the formula seems to give better results than the one based on the
  // quadrature signal. The output of the algo using the quadrature signal seems to be a bit 
  // time-shifted (to the left)

  delete[] x;
  delete[] xq;
  delete[] et;
  delete[] em;
  delete[] em2;
  delete[] y;
}

void sineWithPhaseCatchUp()
{
  static const int N = 5000;    // number of samples
  double f  = 100.0;            // sine frequency
  double fs = 44100.0;          // samplerate
  double p0 = -PI/2;            // phase at sample 0
  double pk = PI/2;             // target phase at sample index k
  int    k  = 1000;             // sample index where pk should be reached
  double as = 0.5;              // amplitude at start
  double ae = 0.3;              // amplitude at end
  int    d  = 0;                // direction: -1: down, +1: up, 0: automatic

  double a[N], x[N];            // amplitude envelope and signals

  // create amplitude envelope (and copy to x):
  RAPT::rsArray::fillWithRangeLinear(a, N, as, ae);
  RAPT::rsArray::copyBuffer(a, x, N);

  // create the sine with envelope - we pass x as envelope to see if it also works when x and a
  // point to the same array:
  rsEnvelopedPhaseCatchSine(x, N, f, fs, p0, pk, k, x, d);

  // write to file:
  writeToMonoWaveFile("SineCatchUp.wav", x, N, (int) fs, 16);

  // plot:
  //plotData(N, 0, 1, x, a);


  // test: create a non-swept sinusoid with the frequency f, obtain its envelope and recreate a 
  // sine with the same envelope but different frequency fy with sweep
  double fy = 170.0;
  double y1[N], y2[N];
  rsEnvelopedSine(y1, N, f, fs, p0, a); 
  rsRecreateSineWithPhaseCatch(y1, y2, N, f, fy, fs, p0, pk, k, 0.0, 0);
  plotData(N, 0, 1, x, y2, a);
}


void partialExtractionTriple()
{
  // create and mix 3 partials at 900, 1000, 1150 Hz, try to retrieve the middle partial at 1000Hz
  // via multipass bidirectional bandpass filtering

  double fs = 44100;          // sample rate
  int    N  = (int)fs;        // number of samples

  // lower partial parameters:
  double fL  = 900;    // frequency in Hz
  double aL  =  0.3;   // amplitude as raw factor
  double pL  =  0.0;   // phase in radians (really radians? - check this)
  double atL =  0.01;  // attack in seconds
  double dcL =  0.2;   // decay in seconds

  // middle partial parameters:
  double fM  = 1000;   // frequency in Hz
  double aM  =  0.3;   // amplitude as raw factor
  double pM  =  0.0;   // phase in radians (really radians? - check this)
  //double atM =  0.025; // attack in seconds
  //double atM =  0.02; // attack in seconds
  //double atM =  0.015; // attack in seconds
  //double atM =  0.005; // attack in seconds
  double atM =  0.001; // attack in seconds
  double dcM =  0.15;  // decay in seconds

  // upper partial parameters:
  double fU  = 1150;   // frequency in Hz
  double aU  =  0.3;   // amplitude as raw factor
  double pU  =  0.0;   // phase in radians (really radians? - check this)
  double atU =  0.02;  // attack in seconds
  double dcU =  0.1;   // decay in seconds

  // filter parameters:
  double bw = 70;      // filter bandwidth in Hz
  int    np = 10;       // number of passes
  double gp = 0.5;     // passband edge gain


  // memory allocation:
  double *xL = new double[N]; // lower mode signal
  double *xM = new double[N]; // middle mode signal
  double *xU = new double[N]; // upper mode signal
  double *x  = new double[N]; // all 3 modes mixed - our signal
  double *yM = new double[N]; // middle mode isolated back from mix
  double *ye = new double[N]; // measured envelope of middle mode

  // synthesize and mix the 3 partials:
  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(fL, aL, atL, dcL, pL, fs);
  getImpulseResponse(mf, xL, N);
  mf.setModalParameters(fM, aM, atM, dcM, pM, fs);
  getImpulseResponse(mf, xM, N);
  mf.setModalParameters(fU, aU, atU, dcU, pU, fs);
  getImpulseResponse(mf, xU, N);
  RAPT::rsArray::add(xL, xM, x, N);
  RAPT::rsArray::add(x,  xU, x, N);

  // retrieve the middle partial via multipass bidirectional bandpassing
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, yM, N, fM, bw, fs, np, gp);

  // write all 3 partials seperately, the mix of all 3 partials and the retrieved middle partial
  // into wave files:
  writeToMonoWaveFile("InputLowerPartial.wav",   xL, N, (int) fs, 16);
  writeToMonoWaveFile("InputMiddlePartial.wav",  xM, N, (int) fs, 16);
  writeToMonoWaveFile("InputUpperPartial.wav",   xU, N, (int) fs, 16);
  writeToMonoWaveFile("InputPartialMix.wav",     x,  N, (int) fs, 16);
  writeToMonoWaveFile("OutputMiddlePartial.wav", yM, N, (int) fs, 16);

  // get envelope of recovered middle mode:
  rsSineEnvelopeViaQuadrature(yM, ye, N, fM, fs, 0.0);

  // get peak location and height of envelope (later, maybe use quadratic interpolation to find it 
  // with subsample precision - this is actually overkill, but however)
  int    nPeak = RAPT::rsArray::maxIndex(ye, N);
  double tPeak = nPeak / fs;
  double aPeak = ye[nPeak];



  //plotSignals(N/2, fs, xL, xM, xU); // plot 3 modes seperately
  //plotSignals(N/2, fs, x);          // plot mix of modes
  //plotSignals(N/4, fs, yM, ye);       // plot isolated mode with measured envelope
  plotSignals(N/4, fs, xM, yM, ye);     // plot original and recovered middle mode and measured 
                                        // envelope


  // Observations:
  // when the mode is not seperated cleanly in the frequency domain because of too little filter 
  // narrowness/steepness, we see modulations in the envelope with frequencies fM-fL and fU-fM (check
  // the frequencies)

  // with np = 10:
  // bw = 100: the partial is not cleanly separated, there's a lot of modulation visible in the 
  // envelope
  // bw = 10: the partial is cleanly separated, but the envelope envelope is smeared very much
  // the sweet spot seems to be somewhere between 50 and 80 Hz
  // generally, maybe it's best to use bandwidths of k * df with k = 0.5...0.8 and df: distance
  // to closest neighbour partial in Hz

  // ideas for improving the envelope mesurement: 
  // the "smoothing" filter that is applied to the envelope could have notches at the modulation 
  // frequencies. more generally, we could perform frequency analysis on the envelope and determine 
  // modulation frequencies and place notches according to the analysis

  // when modulation is present in the envelope, we might want to choose the first local peak 
  // instead of the global one - maybe depending on the amplitude difference between first and 
  // highest - if the difference is small, choose the earlier. this effect is exemplified  with an 
  // attack of 0.025s - the global maximum is found at 0.029 but there's a local maximum at 0.21 
  // before - hmm - actually there's a local minimum at 0.025. maybe we should choose a higher Q 
  // and/or greater smoothing to get rid of the modulation entirely

  // cleanup:
  delete[] xL;
  delete[] xM;
  delete[] xU;
  delete[] x;
  delete[] yM;
  delete[] ye;
}




// kinda obsolete in this context, but maybe keep it for inspiration of other things
void partialExtractionBell()
{
  // read sample:
  const char *samplePath = "D:/Dissertation/Audio/bell_light_1_short.wav";
  int numChannels;         // (some) mode frequencies: 980, 1320, 1800, 2350
  int N;                   // number of sample-frames
  int fs;                  // sample rate
  double **sampleData = readFromWaveFile(samplePath, numChannels, N, fs);

  // allocate internal buffers:
  double *x  = new double[N]; RAPT::rsArray::copyBuffer(sampleData[0], x, N); // input signal
  double *t  = new double[N]; createTimeAxis(N, t, fs);          // time axis for plot
  double *yr = new double[N];                                    // real part of filtered signal
  double *yi = new double[N];                                    // imaginary part
  double *ye = new double[N];                                    // instantaneous envelope

  // set up filter:
  rsNonlinearModalFilterDD mf;
  double f = 980.0;
  double Q = 4.0;  // ad hoc
  mf.setModalParameters(f, 1.0/Q, Q/f, 0.0, fs); // check tau ?= Q/f -> look up definition of Q, tau
                                                 // amplitude scaling of 1/Q ad hoc

  // calculate output:
  rsComplexDbl z;
  for(int n = 0; n < N; n++)
  {
    z     = mf.getComplexSample(0.5*rsComplexDbl(x[n], x[n])); // maybe, we should use 1/sqrt(2) instead of 0.5
    yr[n] = z.real();
    yi[n] = z.imag();
    ye[n] = abs(z);
  }

  plotData(N/10, t, x, yr, yi, ye);
  //plotData(N/10, t, yr, yi, ye);
  //plotData(N/10, t, x, yr, yi);

  // observations:
  // the attack portion of the envelope looks indeed exponential, but that might be an artifact of
  // the envelope analysis filter - the higher we choose the Q, the longer the attack phase gets 
  // and the more "clean" the exponential attack looks, Q values of order 10 seem to isolate the 
  // 1st mode reasonably (can this be generalized to higher modes? maybe not - modes tend to get 
  // denser at higher frequencies on a logarithmic frequency axis), but such high Q smears the 
  // envelope wayyy tooo much, also, the envelope should probably be smoothed (maybe a 
  // bi-directional lowpass could be suitable)

  // maybe, we should isolate the mode first by bi-directional bandpass filtering, maybe applied
  // multiple times. after that, we could use the "instantaneous envelope" approach and/or use an
  // envelope follower

  // maybe we could design a "mode isolation filter" that has zeros at all modes except the one, we
  // are going to isolate (maybe, it should have a pole there to compensate for amplitude loss)
  // or maybe not just zeros but a pole/zero pair that creates a narrow notch for supressing modes

  // or maybe we could calculate the mode envelope by a kind of overlapped STFT, but just for one 
  // frequency bin - i.e. a cross-correlation with a windowed complex sinusoid, maybe of M periods
  // length (1 <= M <= 10, for example), with adjustable overlap (also measured in mode-periods)
  // that would be a kind of wavelet/Gabor transform with adaptive frequencies

  // clean up:
  freeSampleData(sampleData, numChannels);
  delete[] x;
  delete[] t;
  delete[] yr;
  delete[] yi;
  delete[] ye;
}




// later, we should use rsFilterBiDirectional - but not before having a unit test in place for it,
// and this very function can be used inside the unit-test later for producing a reference output
template <class T>
void filterBiDirectional(T *x, int xLength, T *y, int yLength, T *b, int bOrder, T *a, 
   int aOrder, int numRingOutSamples = 0, bool backwardFirst = false)
{
  if( backwardFirst )
    rsReverse(x, xLength);

  if( numRingOutSamples == 0 )
  {
    rsFilter(x, xLength, y, yLength, b, bOrder, a, aOrder);
    rsReverse(y, yLength);
    rsFilter(y, yLength, y, yLength, b, bOrder, a, aOrder);
    rsReverse(y, yLength);
  }
  else
  {
    if( numRingOutSamples == -1 )
    {
      // find a suitable number of samples by looking at the impulse-response of the filter - maybe
      // all state variables should be below some threshold (maybe the machine-epsilon?)
      rsError("not yet implemented");
    }
    int zLength = xLength + rsMax(yLength-xLength, numRingOutSamples);
    T *z = new T[zLength];
    filterBiDirectional(x, xLength, z, zLength, b, bOrder, a, aOrder, 0);
    rsCopyBuffer(z, y, yLength);
    delete[] z;
  }

  if( backwardFirst )
  {
    rsReverse(x, xLength);
    rsReverse(y, yLength);
  }
}

// Given a signal x of length that is a mix of partials, this function tries to isolate the 
// partial at frequency fM by means of (bi-directional) filtering and writes the resulting signal
// into y. it returns an improved estimate of the mode frequency. parameters: 
// x:    input signal
// y:    output signal
// N:    signal length
// fL:   frequency of lower neighbour partial
// fM:   frequency of (middle) partial to be isolated
// fU:   frequency of upper neighbour partial
// p:    multiplier for the filter'S Q factor, the higher, the better the frequency domain 
//       seperation at the expense of more time-domain smearing, reasonable: 2...10
// np:   number of (bi-directional) passes of the filter
double isolatePartialWithBiquad(double *x, double *y, int N, double fL, double fM, double fU, 
  double fs, double p, int np)
{
  // very simple and preliminary algorithm using just a RBJ biquad bandpass:

  double Q = p * fM/(fU-fL); 
    // filter Q - tweak factor in front (let's call the factor "p")
    // observations when passing a single partial at 1kHz, unit amplitude, attack: 0.015, 
    // decay: 0.15:
    // at p = 8: the maximum shifts to t = 0.03 with an amplitude a = 0.85
    // p =  4: t = 0.022, a= 0.93
    // p = 16: t = 0.042, a= 0.72
    // ...maybe, we can model this shifting of peak-location and -height as function of p by 
    // collecting more such data and fitting a function to it, the actual values can then be found 
    // by inverting this function
    // i thought that maybe doing the backward-pass first, with appropriate "pre"-padding to give 
    // the filter enough time to ring out, might shift the peak back into correct position but this 
    // doesn't seem to work, unfortunately
    // it makes definitely sense to use multiple passes with a somewhat lower Q compared to a 
    // high-Q single pass, when q is scaled like 1/numPasses, numPasses = 10 gives good results 
    // with this test signal
    // ...move this comment to the experimental results at the bottom of the file


  /*
  // apply multipass bi-directional bandpass filtering at the mode-frequency:
  double a[3], b[3]; // filter coeffs
  rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ
    (b[0], b[1], b[2], a[1], a[2], 1.0/fs, fM, Q);
  rsNegate(a, a, 3);
  a[0] = 1.0;
  filterBiDirectional(x, N, y, N, b, 2, a, 2);
  for(int i = 1; i < np; i++)
    filterBiDirectional(y, N, y, N, b, 2, a, 2);
  */


  /*
  // test: unidirectional (forward or backward) filtering, see, if this is better for retrieving 
  // attack time - no, not really:
  double a[3], b[3]; // filter coeffs
  rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ
    (b[0], b[1], b[2], a[1], a[2], 1.0/fs, fM, Q);
  rsNegate(a, a, 3);
  a[0] = 1.0;
  rsCopyBuffer(x, y, N);
  //rsReverse(y, N);
  for(int i = 0; i < 2*np; i++) // 2*np to have same total number of passes as bidirectional
    rsFilter(y, N, y, N, b, 2, a, 2);
  //rsReverse(y, N);
  */

  // apply multipass bi-directional bandpass filtering at the mode-frequency with padding to take 
  // filter ringing into account:
  int Np = 10000;  // ad-hoc - should be long enough to hold the filter's ringing
  Np = 0;          // for test
  int M  = N+2*Np; // length with pre- and post-padding
  double *tmp = new double[M];
  RAPT::rsArray::fillWithZeros(tmp, Np);
  RAPT::rsArray::copyBuffer(x, &tmp[Np], N);
  RAPT::rsArray::fillWithZeros(&tmp[Np+N], Np);
  double a[3], b[3]; // filter coeffs
  rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ
    (b[0], b[1], b[2], a[1], a[2], 1.0/fs, fM, Q);
  RAPT::rsArray::negate(a, a, 3);
  a[0] = 1.0;
  for(int i = 0; i < np; i++)
  {
    //filterBiDirectional(tmp, M, tmp, M, b, 2, a, 2);
    RAPT::rsArray::filterBiDirectional(tmp, M, tmp, M, b, 2, a, 2);
  }
  RAPT::rsArray::copyBuffer(&tmp[Np], y, N);
  delete[] tmp;

  /*
  // apply multipass bi-directional bandpass filtering at the mode-frequency on a symmetrized
  // signal - symmetrizing doesn't seem to help to locate the envelope peak more accurately, but 
  // maybe something is wrong with the implementation:
  int M  = 2*N; // length with symmetriczed signal
  double *tmp = new double[M];
  rsCopyBuffer(x, &tmp[0], N);
  rsCopyBuffer(x, &tmp[N], N);
  rsReverse(tmp, N);
  rsNegate(tmp, tmp, N);  // optional - gives signal odd symmetry
  double a[3], b[3]; // filter coeffs
  rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ
    (b[0], b[1], b[2], a[1], a[2], 1.0/fs, fM, Q);
  rsNegate(a, a, 3);
  a[0] = 1.0;
  for(int i = 0; i < np; i++)
    filterBiDirectional(tmp, M, tmp, M, b, 2, a, 2);
  rsCopyBuffer(&tmp[N], y, N);
  delete[] tmp;
  */


  fM = fM; // preliminary - todo: refine mode frequency estimate fM by counting zero-crossings

  // re-adjust gain of filter output signal (factor out):
  //rsNegate(a, a, 3); // analysing function's sign convention is different from filtering routine - fix this!
                       // hmm - nope - doesn't seem to be the case
  double w = 2*PI*fM/fs;
  double g = biquadMagnitudeAt(b[0], b[1], b[2], a[1], a[2], 2*PI*fM/fs);
  g = pow(g, 2*np);  // bidirectional: g^2, multipass: g^np -> g^(2*np)
  RAPT::rsArray::scale(y, N, 1.0/g);
    // we are getting HUUUUGE gain factors here, maybe, we should design a constant peak-gain 
    // bandpass instead of constant skirt gain to avoid this


  // compute relative gains for lower and upper neighbour (just during development to find out, how 
  // much attenuation we need for other modes):
  double gL = biquadMagnitudeAt(b[0], b[1], b[2], a[1], a[2], 2*PI*fL/fs);
  gL = pow(gL, 2*np);
  double gU = biquadMagnitudeAt(b[0], b[1], b[2], a[1], a[2], 2*PI*fU/fs);
  gU = pow(gU, 2*np);
  double dBL = rsAmp2dB(gL/g);
  double dBU = rsAmp2dB(gU/g);

  return fM;
}

double isolatePartialWithElliptic(double *x, double *y, int N, double fL, double fM, double fU, 
  double fs, double aM, double aMax, double p)
{
  // later, we should use aM and aMax to determine the required stopband rejection: if they are
  // expressed in decibels (as dbM and dbMax), we should add dbMax-dbM to the required rejection

  return fM;
}



// Given a signal x of length that is a mix of partials, this function tries to isolate the 
// partial at frequency fM by means of (bi-directional) filtering and writes the resulting signal
// into y. it returns an improved estimate of the mode frequency. parameters: 
// x:    input signal
// y:    output signal
// N:    signal length
// fL:   frequency of lower neighbour partial
// fM:   frequency of (middle) partial to be isolated
// fU:   frequency of upper neighbour partial
// aM:   amplitude of partial to be extracted
// aMax: amplitude of loudest partial in the signal
// p:    a time/frequency resolution tradeoff parameter: the higher, the better the frequency 
//       domain seperation at the expense of more time-domain smearing
// np:   number of (bi-directional) passes of the filter
double isolatePartial(double *x, double *y, int N, double fL, double fM, double fU, double fs, 
  double aM, double aMax, double p, int np)
{
  return isolatePartialWithBiquad(x, y, N, fL, fM, fU, fs, p, np);
    // p scales the filter's Q - observations when passing a single partial at 
    // 1kHz, unit amplitude, attack: 0.015, decay: 0.15:
    // at p = 8: the maximum shifts to t = 0.03 with an amplitude a = 0.85
    // p =  4: t = 0.022, a= 0.93
    // p = 16: t = 0.042, a= 0.72
    // ...maybe, we can model this shifting of peak-location and -height as function of p by 
    // collecting more such data and fitting a function to it, the actual values can then be found 
    // by inverting this function
}


// for meaning of parameters, see documentation of isolatePartial
void plotModalEnvelope(const char *samplePath, double fL, double fM, double fU, double aM, 
  double aMax, double p, int np, int numSamplesToPlot)
{
  // things worth to try: use bandpasses with constant bandwidth in Hz, or adjust lower and upper 
  // cutoff frequency somewhere between the mode of interest and its neighbours, the bandpass 
  // should have unity gain at the center frequency, maybe higher order filters (butterworth, 
  // papoulis, elliptic etc. could be used) if the bandpass has non-unity gain at the modal 
  // frequency, calculate the actual gain and apply an appropriate amplification to the output 
  // signal

  int numChannels; 
  int N;         // number of sample-frames
  int fs;        // sample rate
  double **sampleData = readFromWaveFile(samplePath, numChannels, N, fs);

  // allocate internal buffers:
  double *x  = new double[N]; RAPT::rsArray::copyBuffer(sampleData[0], x, N); // input signal
  double *y  = new double[N];                                    // filtered signal
  double *ye = new double[N];                                    // envelope
  double *yn = new double[N];                                    // negated envelope

  // extract partial and envelope:
  fM = isolatePartial(x, y, N, fL, fM, fU, fs, aM, aMax, p, np);
  rsSineEnvelopeViaQuadrature(y, ye, N, fM, (double)fs, 4.0);

  // plot the extracted mode with its envelope:
  RAPT::rsArray::negate(ye, yn, N);
  numSamplesToPlot = rsMin(N, numSamplesToPlot);
  plotSignals(numSamplesToPlot, fs, y, ye, yn);

  // clean up:
  freeSampleData(sampleData, numChannels);
  delete[] x;
  delete[] y;
  delete[] ye;
  delete[] yn;
}



void partialExtractionViaBiquadTriple()
{
  // create 3 partials at 900, 1000, 1150 Hz, try to retrieve the 1000Hz mode

  double fs = 44100;    // sample rate
  int    N  = (int)fs;  // number of samples

  // memory allocation:
  double *xL = new double[N]; // lower mode signal
  double *xM = new double[N]; // middle mode signal
  double *xU = new double[N]; // upper mode signal
  double *x  = new double[N]; // all 3 modes mixed - our signal
  double *yM = new double[N]; // middle mode isolated back from mix
  double *ye = new double[N]; // measured envelope of middle mode

  // lower neighbour parameters:
  double fL  = 900;    // frequency in Hz
  double aL  =  1.0;   // amplitude as raw factor
  double pL  =  0.0;   // phase in radians (really radians? - check this)
  double atL =  0.01;  // attack in seconds
  double dcL =  0.2;   // decay in seconds

  // middle mode parameters:
  double fM  = 1000;   // frequency in Hz
  double aM  =  1.0;   // amplitude as raw factor
  double pM  =  0.0;   // phase in radians (really radians? - check this)
  //double atM =  0.025; // attack in seconds
  //double atM =  0.02; // attack in seconds
  //double atM =  0.015; // attack in seconds
  //double atM =  0.005; // attack in seconds
  double atM =  0.001; // attack in seconds
  double dcM =  0.15;  // decay in seconds

  // upper mode parameters:
  double fU  = 1150;   // frequency in Hz
  double aU  =  1.0;   // amplitude as raw factor
  double pU  =  0.0;   // phase in radians (really radians? - check this)
  double atU =  0.02;  // attack in seconds
  double dcU =  0.1;   // decay in seconds



  // synthesize and mix the 3 modes
  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(fL, aL, atL, dcL, pL, fs);
  getImpulseResponse(mf, xL, N);
  mf.setModalParameters(fM, aM, atM, dcM, pM, fs);
  getImpulseResponse(mf, xM, N);
  mf.setModalParameters(fU, aU, atU, dcU, pU, fs);
  getImpulseResponse(mf, xU, N);
  RAPT::rsArray::add(xL, xM, x, N);
  RAPT::rsArray::add(x,  xU, x, N);


  // just for inspecting the impulse-reponse during development - replace our x-signal with a 
  // unit-impulse at 0.1s, we need to uncomment plotSignals(N/4, fs, yM, ye); below
  //rsFillWithZeros(x, N);
  //x[4410] = 1.0;



  // try to recover middle mode - we have adjusted the Q-scaler p for each number of passes such 
  // that the supression of neighbouring partials is 50dB which seems to be the best compromise 
  // between frequency selectivity and time resolution. this amount of supression reduces the 
  // modulation in the envelope just enough to avoid local maxima:


  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs, 21.0,     1);
  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  4.84,    2);
  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  2.12,    4);
  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  1.212,   8);
  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  0.777,  16);
  //isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  0.5245, 32);
  isolatePartialWithBiquad(x, yM, N, fL, fM, fU, fs,  0.3622, 64);


  // we should find an algorithm that determines the required Q (or p) such that we obtain a 
  // specified supression at some specified frequency

  // for the mode extraction with the biquad, a p-value of 0.6 with 20 passes seems to be the point
  // after which increasing the number of passes further (and reducing p further) does not give
  // much improvement - a plateau seems to be reached. 
  // \todo: investigate, if this set of parameters is also optimal for other mode frequencies
  // like 4900, 5000, 5150

  // as test, feed just the middle mode into the filter, see how it gets smeared:
  //isolatePartialWithBiquad(xM, yM, N, fL, fM, fU, fs, 0.6, 20);


  // get envelope of recovered middle mode:
  rsSineEnvelopeViaQuadrature(yM, ye, N, fM, fs, 0.0);

  // get peak location and height of envelope (later, maybe use quadratic interpolation to find it 
  // with subsample precision - this is actually overkill, but however)
  int    nPeak = RAPT::rsArray::maxIndex(ye, N);
  double tPeak = nPeak / fs;
  double aPeak = ye[nPeak];

   // for fM = 1000, attack = 0.001:
   // passes: 16: tp = 0.013877551020408163, ap = 0.90231040500525828
   // passes: 32: tp = 0.013378684807256236, ap = 0.90712670808762219 -> 5% better
   // passes: 64: tp = 0.012879818594104309, ap = 0.90930771999937432 -> 5% better
   // ...we still get gains but driving the number of passes ever higher suggests that the 
   // bidirectional IIR approach might be replaced by using a gaussian FIR instead (the more passes
   // we do, the more the impulse response approximates a sinusoid with gaussian envelope), or use
   // a higher order gaussian IIR filter to start with -> maybe the gaussian shape is approximated
   // more rapidly (with fewer passes) when the single-pass filter itself approximates a gaussian 
   // shape


  //plotSignals(N/2, fs, xL, xM, xU); // plot 3 modes seperately
  //plotSignals(N/2, fs, x);          // plot mix of modes
  //plotSignals(N/4, fs, yM, ye);       // plot isolated mode with measured envelope
  plotSignals(N/4, fs, xM, yM, ye);     // plot original and recovered middle mode and measured 
                                        // envelope

  // observations:
  // when the mode is not seperated cleanly in the frequency domain because of too little
  // filter steepness, we see modulations in the envelope with frequencies fM-fL and fU-fM (check
  // the frequencies) - this suggests that the "smoothing" filter that is applied to the envelope
  // should have notches there. more generally, we could perfom frequency analysis on the envelope
  // and determine modulation frequencies and place notches according to the analysis

  // when modulation is present in the envelope, we might want to choose the first local peak 
  // instead of the global one - maybe depending on the amplitude difference between first and 
  // highest - if the difference is small, choose the earlier. this effect is exemplified  with an 
  // attack of 0.025s - the global maximum is found at 0.029 but there's a local maximum at 0.21 
  // before - hmm - actually there's a local minimum at 0.025. maybe we should choose a higher Q 
  // and/or greater smoothing to get rid of the modulation entirely

  // cleanup:
  delete[] xL;
  delete[] xM;
  delete[] xU;
  delete[] x;
  delete[] yM;
  delete[] ye;
}

void partialExtractionSample()
{
  int N = 22050; // number of samples to plot

  // Beerglas 
  // mode frequencies: 655, 1310, 1630, 2530, 3100
  // the values passed as aM and aMax are not yet correct (set to unity)

  N = 11025;
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassRubber.wav",    0,  655, 1310, 1, 1, 8, 1, N); // 1st
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassRubber.wav",    0,  655, 1310, 1, 1, 0.6, 20, N); // 1st
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassRubber.wav", 1310, 1630, 2530, 1, 1, 8, 1, N); // 3rd - has amplitude modulation with decreasing frequency
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassRubber.wav", 1310, 1630, 2530, 1, 1, 0.6, 20, N); // 3rd
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassTeaspoon.wav", 0, 655, 1310, 1, 1, 8, 1, N); // 1st
  //plotModalEnvelope("D:/Dissertation/Audio/BeerglassTeaspoon.wav", 0, 655, 1310, 1, 1, 0.6, 20, N);

  // observation: the envelope of the 1st mode looks more like having a fast early decay and slower
  // late decay - might be modeled by adding 3 exponentials, 1 for late decay, 1 for early decay 
  // and one for attack, like ae*exp(-t/te) + al*exp(-t/tl) - (ae+al)*exp(-t/ta), 
  // where te: tau early, tl: tau late, ta: tau attack, we should have tl > te > ta
  // parallel connection of 3 instead of 2 filters

  // when exited with the rubber, the 1st mode reaches its maximum after 0.012s whereas when 
  // excited with the teaspoon, the maximum is reached after 0.017s - so a soft exictation lead to 
  // faster buildup of a low-frequency mode - maybe in an hard excitation, higher modes are excited 
  // first and over time energy is transfered into the lower modes via nonlinear effects?

  //-----------------------------------------------------------------------------------------------

  // Piano E2:
  // mode frequencies: 82, 164, 246, 328, 410, ...
  // the values passed as aM and aMax are not yet correct (set to unity)

  N = 22050;
  //plotModalEnvelope("D:/Dissertation/Audio/piano_E2.wav",   0,  82, 164, 1, 1, 0.6, 20, N); // 1st
  plotModalEnvelope("D:/Dissertation/Audio/piano_E2.wav",  82, 164, 246, 1, 1, 0.6, 20, N); // 2nd - envelope fits model well
  //plotModalEnvelope("D:/Dissertation/Audio/piano_E2.wav", 164, 246, 328, 1, 1, 0.6, 20, N); // 3rd
  //plotModalEnvelope("D:/Dissertation/Audio/piano_E2.wav", 246, 328, 410, 1, 1, 0.6, 20, N); // 4th
}




