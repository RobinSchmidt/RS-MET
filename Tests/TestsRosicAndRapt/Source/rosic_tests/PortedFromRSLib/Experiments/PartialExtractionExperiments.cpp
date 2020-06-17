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
 -apply 1st order allpass with 90° frequency adjusted to fM to obtain a quadrature signal
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

Having isolated a partial with frequency fM, we estimate its amplitude envelope by obtaining a 90°
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
  RAPT::rsArrayTools::fillWithZeros(x, N);
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
  RAPT::rsArrayTools::fillWithRangeLinear(freqs, N, 0.0, fs);
  //plotDataLogX(N/2, freqs, m1, m2, m4, m8, m16);
  //plotData(N/2, 0.0, fs/N, m1, m2, m4, m8, m16); // magnitude responses
  plotData(N/16, 0.0, fs/N, m1, m2, m4, m8, m16); // magnitude responses - shouldn't it be
                                                  // fs/(N+1)?

  // normalize and write output files:
  double s = 1 / RAPT::rsArrayTools::maxAbs(y1, N);
  RAPT::rsArrayTools::scale(y1,  N, s);
  RAPT::rsArrayTools::scale(y2,  N, s);
  RAPT::rsArrayTools::scale(y4,  N, s);
  RAPT::rsArrayTools::scale(y8,  N, s);
  RAPT::rsArrayTools::scale(y16, N, s);
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

void beatingSines()
{
  // We investigate the beating effects that occur when two sinusoids with similar frequencies
  // are played simultaneously.

  // experiment parameters:
  int N = 1000;           // number of samples
  double tMin   =  0.0;   // start time
  double tMax   = 40.0;   // end time
  double freq1  =  0.95;  // frequencies of first...
  double freq2  =  1.05;  // ...and second sine (in cycles per time unit)
  //double amp1   =  1.0;   // amplitudes of first...
  //double amp2   =  1.0;   // ...and second sine (as raw multiplication factor)
  double phase1 =  0.0;   // phases of first...
  double phase2 =  0.0;   // ...and second sine (in degrees)

  // algorithm parameters:
  double w1 = 2*PI*freq1;
  double w2 = 2*PI*freq2;
  double p1 = rsDegreeToRadiant(phase1);
  double p2 = rsDegreeToRadiant(phase2);

  // compute carrier and modulator parameters:
  double wc = (w1+w2)/2;
  double wm = (w1-w2)/2;
  double pc = (p1+p2)/2;
  double pm = (p1-p2)/2 + PI/2;


  // synthesize signals:
  std::vector<double> t(N), s1(N), s2(N), sc(N), sm(N), sum(N), prod(N); // sc,sm: carrier,modulator
  RAPT::rsArrayTools::fillWithRangeLinear(&t[0], N, tMin, tMax);     // time axis
  for(int n = 0; n < N; n++) {
    //s1[n]  = amp1 * sin(w1*t[n] + p1);  // i couldn't find formulas for re-expressing the sum of 
    //s2[n]  = amp2 * sin(w2*t[n] + p2);  // sines as product for arbitrary amplitudes - so use 1
    s1[n]  = sin(w1*t[n] + p1);
    s2[n]  = sin(w2*t[n] + p2);
    sum[n] = s1[n] + s2[n];

    // compute carrier and modulator and their product::
    sc[n]   =   sin(wc*t[n] + pc);
    sm[n]   = 2*sin(wm*t[n] + pm);
    prod[n] = sc[n] * sm[n];
  }

  // now we re-create the signal by means of the strictly positive rectified-sine amplitude
  // function and a time-varying phase function - as the sinusoidal modeling would do:
  std::vector<double> a(N), p(N), ip(N), rs(N), db(N);
  // amplitude- and phase functions, instantaneous phase and recreated signal
  for(int n = 0; n < N; n++) {
    a[n] = fabs(sm[n]);                  // amplitude function is rectified modulator
    p[n] = pc;                           // phase function is...
    if(sm[n] < 0) p[n] += PI;            // ...a square wave
    ip[n] = wc*t[n] + p[n];              // instantaneous phase
    rs[n] = a[n] * sin(ip[n]);           // re-synthesized signal
    db[n] =        sin(ip[n]);           // "de-beated" version
  }


  // plotting:
  GNUPlotter plt;

  // verification - plot pairs of signals that should be equal to one another:
  //plt.addDataArrays(N, &t[0], &sum[0], &prod[0]);         // sum and product - should be equal
  //plt.addDataArrays(N, &t[0], &sum[0], &rs[0]); // sum and recreated sum - should be equal

  // results of adding two components or ring-modulating carrier and modulator:
  //plt.addDataArrays(N, &t[0], &s1[0], &s2[0], &sum[0]);   // component sines and sum
  //plt.addDataArrays(N, &t[0], &sc[0], &sm[0], &prod[0]);  // carrier, modulator and product

  // this is the plot, we actually want:
  //plt.addDataArrays(N, &t[0], &s1[0], &s2[0], &sum[0], &sc[0], &sm[0]);
  // the green signal can be thought of in two ways: 
  // 1: sum of black and blue (how it's actually done)
  // 2: product of red (carrier) and purple (modulator)

  // amplitude- and phase-function used to re-create the final output as-is together with the
  // de-trende intsantaneous phase and re-created output:
  rsDeTrender<double> dtr;
  dtr.removeTrendAndOffset((int)ip.size(), &t[0], &ip[0], &ip[0]);
  //plt.addDataArrays(N, &t[0], &a[0], &ip[0], &rs[0]); // only detrended inst phase
  //plt.addDataArrays(N, &t[0], &a[0], &p[0], &ip[0], &rs[0]);
  plt.addDataArrays(N, &t[0], &a[0], &db[0], &rs[0]);  // de-beated


  plt.setPixelSize(1600, 400);
  plt.plot();

  // todo: figure out the phase-function for the carrier that results from interpreting the sign as
  // part of the carrier - the analysis gives a sawtooth function...which may make sense - we have 
  // the general linear upward trend that any phase-function has together with these jumps
  // ..so maybe it's pc(t) = wc*t - sum_of_delta_functions ...but i think, in the analysis, there
  // was a downward saw...hmmm

  // When thinking of the signal as the product of a carrier and a modulator signal, we may further
  // decompose the modulator into an absolute value (the amplitude) and a sign. We want to enforce 
  // a positive amplitude, so we interpret the sign as belonging to the carrier - it switches 
  // between positive and negative at the zero crossings of the modulator. A negative sign, in 
  // turn, can be interpreted as a 180° phase discontinuity. So, we may interpret beating as a 
  // combination of ring-modulation with a rectified sine-wave combined with phase-modulation with
  // a square wave. Both effects are tuned with respect to each other so as to yield a smooth 
  // overall signal - but if we remove one, the amplitude ring-modulation, say - the other one that
  // remains will make the signal non-smooth. to re-smooth it, we may lowpass filter the 
  // phase-modulation or remove it altogether.

  // If we just lowpass it, a sort of sinusoidal freq/phase 
  // modulation would remain. ...but wait - is the phase modulator really a rectangle wave - in the
  // harmonic analysis framework, it looks more like sawtooth shape. In any case, it has harsh
  // discontinuities that we need to get rid of.

  // When both start-phases are zero, then whenever the modulator crosses zero, the carrier goes 
  // through a zero as well - the carrier always has an upward zero crossing whereas the modulator 
  // has alternatingly and upward and downward zero crossing. Bute when the 2nd sine has a phase of 
  // 180°, the modulator zeros coincide with carrier peaks (all maxima at +1 - but these become 
  // minima at -1, when the phase is -180)

  // I was trying to find a carrier/modulator-like decomposition for the general case, when the two
  // sines each have their own amplitude - not successfully yet:
  // sage:
  //  var("t w1 w2 a1 a2 p1 p2")
  //  f(t) = a1 * sin(w1*t + p1) + a2 * sin(w2*t + p2)
  //  f.simplify_trig() # expands it in terms of sines and cosines (angles go into amplitudes)
  // a1*cos(p1)*sin(t*w1) + a2*cos(p2)*sin(t*w2) + a1*cos(t*w1)*sin(p1) + a2*cos(t*w2)*sin(p2)
  // a1 * (cos(p1)*sin(t*w1) + cos(t*w1)*sin(p1)) + a2 * (cos(t*w2)*sin(p2) + cos(p2)*sin(t*w2))
  // f.expand_trig(), f.reduce_trig()
  // https://brilliant.org/wiki/sum-to-product-trigonometric-identities/

  // How about trying it via a complex exponential:
  // e^(ia) + e^(ib) = e^(-i(a+b)) * ( e^(ia + i(a+b)) + e^(ib + i(a+b)) )
  // ...hmm...but that seems a dead end
}

void envelopeDeBeating()
{
  // We create two attack/decay sinusoids with frequencies close to each other such that the 
  // resulting amplitude envelope will show beating effects at the difference frequency. Then, we
  // try to remove the amplitude modulation from the envelope.

  double fs = 44100; // sample rate
  double f1 = 147;   // frequency 1
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


  //rsEnvelopeExtractor<T>

  typedef RAPT::rsEnvelopeExtractor<double> EE;
  typedef RAPT::rsInterpolatingFunction<double, double> IF;
  EE envExtractor;
  envExtractor.setInterpolationMode(IF::CUBIC_NATURAL);
  //envExtractor.setInterpolationMode(IF::LINEAR);
  //envExtractor.setSampleRate(fs);
  //envExtractor.setSmoothing(20.0, 4);
  envExtractor.setStartMode(EE::FREE_END);
  envExtractor.setEndMode(  EE::FREE_END);
  //envExtractor.setStartMode(EE::ZERO_END);
  //envExtractor.setEndMode(  EE::ZERO_END);
  //envExtractor.setStartMode(EE::EXTRAPOLATE_END);
  //envExtractor.setEndMode(  EE::EXTRAPOLATE_END);
  envExtractor.sineEnvelopeWithDeBeating(&x[0], N, &env[0]);


  //RAPT::rsEnvelopeExtractor<double>::sineEnvelopeWithDeBeating(&x[0], N, &env[0]);

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
  // With smooth 0.0, we see a little modulation in the measured envelope, with smooth = 1.0, it's
  // supressed reasonably but still visible, with smooth = 2.0, the envelope is very smooth indeed 
  // but might be already oversmoothed. Maybe 1.5 is the sweet spot.
  // The algorithm based on the formula seems to give better results than the one based on the
  // quadrature signal. The output of the algo using the quadrature signal seems to be a bit 
  // time-shifted (to the left).

  delete[] x;
  delete[] xq;
  delete[] et;
  delete[] em;
  delete[] em2;
  delete[] y;
}


template<class T>
void rsSineAmplitudeAndPhase(T yL, T y0, T yR, T w, T* a, T* p)
{
  T sw = sin(w);
  T cw = cos(w);

  T pR = atan2( y0*sw, yR - y0*cw);       // phase estimate at y0 using right neighbour yR

  T pL = atan2(-y0*sw, yL - y0*cw) + PI;  // phase estimate at y0 using left neighbour yL
  // get rid of the addition of pi by rotating the input to atan2 by 180° (negate both inputs)
  pL = atan2(y0*sw, y0*cw - yL);


  int dummy = 0;
}
// when the phase has been computed, to decide, which of the 3 values we want to use to compute the
// amplitude, we want to use the one with the largest absolute value of: 
//   sL = sin(p-w), s0 = sin(p), sR = sin(p+w) 
// and then do:
//   a = yL/sL or a = y0/s0 or a = yR/sR. 
// where we should use the value, where the argument for the sine is farthest away from a multiple
// of pi (i think, only -pi, 0, pi need to be considered). hmm...but is perfect resynthesisi 
// guaranteed in all 3 cases or only in the case of using s0?



// like rsSineAmplitudeAndPhase in rapt, but uses left neighbour of y0 instead of right neighbour
template<class T>
void rsSineAmplitudeAndPhaseL(T y0, T yL, T w, T* a, T* p)
{
  T sw = sin(w);
  T cw = cos(w);

  *p = atan2(-y0*sw, yL-y0*cw) + PI;
  *a = y0 / sin(*p);

  int dummy = 0;
}
// ...not yet finished...needs a switch to avoid div-by-zero

bool testSineAmpAndPhaseEstimation2()
{
  bool r = true;
  double tol = 1.e-13;


  double a = 0.3;  // actual amplitude
  double p = 1.0;  // actual phase
  double w = 0.5;  // normalized radian freq


  // the three sample values:
  double yL = a*sin(p - w);
  double y0 = a*sin(p);
  double yR = a*sin(p + w);

  // measurements:
  double aR, pR, aL, pL;
  rsSineAmplitudeAndPhase(y0, yR, w, &aR, &pR); // old library function
  r &= rsIsCloseTo(aR, a, tol);
  r &= rsIsCloseTo(pR, p, tol);
  rsSineAmplitudeAndPhaseL(y0, yL, w, &aL, &pL);
  r &= rsIsCloseTo(aL, a, tol);
  r &= rsIsCloseTo(pL, p, tol);




  return r;
}

void testSineAmpAndPhaseEstimation()
{
  double y0 = 0.2;
  double y1 = 0.3;
  double w  = 0.5;
  double a, p;

  rsSineAmplitudeAndPhase(y0, y1, w, &a, &p);

  // try new formula - it's simpler but has a division in the phase-computation:
  double sw = sin(w);
  double cw = cos(w);
  double p2 = atan2(y0*sw, y1 - y0*cw);   // ...that should fix it
  double a2 = y0 / sin(p2);               // what if p2 == 0?

  // try left-looking formula - should give the same amplitude but the phase should be incremented
  // by w:
  double p3 = atan2(-y1*sw, y0-y1*cw);
  double a3 = y1 / sin(p3);
  // hmm..we get a negative amplitude - the absolute value is correct, though - hmm..so let's flip
  // phase and amplitude:
  p3 += PI;
  a3  = -a3;
  // i think, we could get the same effect by rotating the input to atan2 by 180°

  double w2 = p3-p2; // should reconstruct w

  // ok - it works in this case - but we should really test many more cases with different values
  // for y0, y1, w - maybe this should become a unit test, like this...
  double tol = 1.e-13;
  bool r = true;
  r &= rsIsCloseTo(p2+w, p3, tol);
  r &= rsIsCloseTo(a2,   a3, tol);
  r &= rsIsCloseTo(w,    w2, tol);

  int dummy = 0;
}
// in practice, use both formulas, compute an average of the phases and compute the amplitude for
// exact resynthesis

// Computes the median of 3 values - maybe move to library (near rsMin/rsMax):
template<class T>
T median(T x1, T x2, T x3)
{
  if(x1 >= x2 && x1 >= x3) return rsMax(x2, x3);  // x1 is greatest
  if(x2 >= x1 && x2 >= x3) return rsMax(x1, x3);  // x2 is greatest
  if(x3 >= x1 && x3 >= x2) return rsMax(x1, x2);  // x3 is greatest
  rsError("we should always take one of the branches above");
  return 0;
}
// test with all permutations of 1,2,3 and 1,2,2

// Computes the error function - see SineParameters.txt
template<class T>
T rsSineParametersError(T yLL, T yL, T y0, T yR, T yRR, T a, T p, T w)
{
  T sLL = sin(p-2*w);
  T sL  = sin(p-w);
  T s0  = sin(p);
  T sR  = sin(p+w);
  T sRR = sin(p+2*w);

  T eLL = yLL - a*sLL;
  T eL  = yL  - a*sL;
  T e0  = y0  - a*s0;
  T eR  = yR  - a*sR;
  T eRR = yRR - a*sRR;

  return eLL*eLL + eL*eL + e0*e0 + eR*eR + eRR*eRR;
}
// todo: minimize this function numerically with respect to a,p,w - maybe an analytic solution is
// possible, but it involves solving a messy nonlinear system of equations...
// maybe define a function that computes also the gradient - it can re-use the sLL,.. stuff - but we 
// will need the cosines, too

// numerically optimizes a,p,w so as to minimize the error function given above - returns the 
// number of function evaluations
template<class T>
int rsOptimizeSineParameters(T yLL, T yL, T y0, T yR, T yRR, T* a, T* p, T* w)
{
  double params[3];
  params[0] = *a;
  params[1] = *p;
  params[2] = *w;
  std::function<double(double*)> f = [&](double* params)->double
  { 
    double err = rsSineParametersError(yLL, yL, y0, yR, yRR, params[0], params[1], params[2]);
    return err;
  };
  double tol = 1.e-12;
  double h[3] = { pow(2.0,-15), pow(2.0,-15), pow(2.0,-15) };

  int evals = minimizePartialParabolic(f, params, 3, h, tol); 
  // todo: use better algo - maybe gradient descent with analytic gradient?

  *a = params[0];
  *p = params[1];
  *w = params[2];

  return evals;
}


// can be delete AFTER checking, if all relevant comments have been included into the 
// documentation of rsSineParameterEstimator<T>::sigToAmpsViaPeaks

template<class T>
void rsAmpEnvelope(const T* x, int N, T* a)
{
  rsSineParameterEstimator<T>::sigToAmpsViaPeaks(x, N, a);
  return;


  // todo: take a shadowing-time parameter and use a peak-shadower

  // Algo:
  // -obtain shadowed version of abs(x)
  // -find peaks
  //  -maybe refine their values by using the maximum through a parabola
  // -connect them by linear interpolation

  bool parabolicHeight = true; 
  bool parabolicTime   = false; // makes sense only, if parabolicHeight is also true
  // Env looks smoother with parabolicHeight. parabolicTime may lead to the envelope-estimate 
  // undershooting the signal - so use with care. I think, we should use parabolicHeight but not
  // parabolicTime, if the goal is to estimate the instantaneous phase, too


  for(int n = 0; n < N; n++)
    a[n] = rsAbs(x[n]);  // todo: apply shadower here (shadows are casted only rightward)

  //rsPlotArrays(N, x, a);

  int nL = 0,     nR;  // index of current left and right peak
  T   tL = T(nL), tR;  // position of current left and right peak
  T   aL = a[0],  aR;  // amplitude of current left and right peak
  for(int n = 1; n < N-1; n++)
  {
    if(a[n] >= a[n-1] && a[n] >= a[n+1])
    {
      // there's a peak at a[n]...
      nR = n;
      tR = T(nR);
      aR = a[n];


      if(parabolicHeight)
      {
        // maybe factor out into a function rsParabolicExtremumValue(T* x, int n)
        using Poly = rsPolynomial<T>;
        T c[3]; Poly::fitQuadratic_m1_0_1(c, &a[n-1]);  // c = polynomial coeffs of parabola
        if(c[2] != 0)  // TODO: use a tolerance
        {
          T dt = Poly::quadraticExtremumPosition(c); // time offset of peak between -1..+1
          aR   = Poly::evaluate(dt, c, 2);           // height of peak
          if(parabolicTime)                          // we may or may not use the time offset..
            tR += dt;                                // ..in the linear interpolation loop below
        }
        int dummy = 0;

        // quadraticExtremumPosition computes c[1]/c[2], so the tolerance should be based on the 
        // ratio |c[1]| and |c[2]| - if abs(c[2]) < (small * c1), skip the step

        // todo: modify aR to be the peak of the parabola going through a[n-1], a[n], a[n+1],
        // take care to handle linear sections and plateaus correctly (the parabola becomes
        // degenerate in such cases and the formula will give a division by zero). we do not not 
        // bother to do anything about the location of the peak - it will actually be located at a 
        // subsample position around n but here, we have no way to represent that - but adjusting
        // the height may be beneficial anyway - we assume the amplitude to be approximately 
        // constant between n-1...n+1 ...or wait: we actually *can* use non-integers xL, xR instead
        // of nL, nR...try it!
      }


      for(int i = nL; i < nR; i++)
      {
        a[i] = rsLinToLin(T(i), tL, tR, aL, aR); // optimize!
        // todo: test, if 2nd version is really better ...looks strange at sample 504 - undershoots
        // actual value there - try large flat peaks next to small peaks or maybe flat peaks near
        // sharp peaks of equal height
        // for the phase-formula p[n] = asin(x[n]/a[n]), we need to ensure that a[n] >= |x[n]|, so 
        // such undershoots should be avoided
      }
      //rsPlotArrays(N, x, a);


      // update for next iteration:
      nL = nR;
      tL = tR;
      aL = aR;
    }

  }

  nR = N-1;
  tR = T(nR);
  aR = a[nR];
  for(int i = nL; i < nR; i++)
    a[i] = rsLinToLin(T(i), tL, tR, aL, aR); // optimize!

  //rsPlotArrays(N, x, a);

  // maybe optionally smooth the result by a bidirectional filter

  //int dummy = 0;
}
// In the case of using amplitude as primary (i.e. first estimated) variable, we need to look in 
// the neighbourhood of peaks for other peaks. When frequency is the primary variable, we need to
// look for other zero crossing in the neighborhood of zero-crossings (at least, when the 
// zero-crossing based freq estimation is used). When we use local information only, like just
// the two neighbouring samples of each sample, we have the "most localized" estimation algo.
// ...move all this stuff into a class rsSineParameterEstimator...this function may be called
// signalToAmp and we should also have signalToFreq1, signalToFreq2, signalAndFreqToAmpAndPhase, 
// etc.



template<class T>
void repairPhase(T* p, int N)
{
  for(int n = 0; n < N; n++)
  {
    if(rsIsNaN(p[n])) // what about inf? can this occur? i don't think so
    {
      p[n] = 0.0;   // preliminary
      rsError("We should probably prevent that from happening");
    }
  }

  int dummy = 0;
}
// we must also consider the possibility of a stretch of NaNs - when we see a NaN, we must check, 
// how many NaNs follow and then fill the whole section with "better" data - but these changes 
// will break identity resynthesis - a NaN is actually a sign that the amplitude estimate is wrong
// ...maybe we should repair this too
// encounters far mor NaNs, if we use
//  bool parabolicTime = true  in   rsSineParameterEstimator<T>::connectPeaks
// but even if we don't use it - some NaNs may still happen



template<class T>
void unreflectPhase(T* p, int N)
{
  std::vector<double> tmp = toVector(p, N); // for plotting
  T pi2 = 0.5*PI;
  for(int n = 1; n < N; n++)
  {
    if(p[n] >= 0)
    {
      if(p[n] < p[n-1] || (p[n] > p[n-1] && p[n-1] >= pi2)  )
      {
        p[n] = PI - p[n];
        // this looks ok
      }
    }
    else  // p[n] < 0
    {
      bool cond1 = p[n] < p[n-1];  // 
      bool cond2 = p[n] > p[n-1] && p[n-1] <= -pi2;


      if(cond1 || cond2 /* || cond3*/ )
      {
        p[n] = -PI - p[n];  
        // this looks wrong - test with sine input - maybe we nee a third condition for the
        // wrap-around PI
      }

    }
  }
  rsPlotArrays(N, &tmp[0], p);
  // maybe we should use conditions not based on previous values of the phase but on whether x
  // as ascending or descending - if it's descending, we are either in the range pi/2...pi or
  // -pi...-pi/2

}

template<class T>
void unreflectPhase2(const T* x, T* p, int N)
{
  std::vector<double> tmp = toVector(p, N); // for plotting
  T pi2 = 0.5*PI;
  for(int n = 1; n < N; n++)
  {
    if(x[n] >= 0)
    {
      //if(x[n] < x[n-1] && p[n] < p[n-1]) // why the 2nd codition?
      if(x[n] < x[n-1]) 
        p[n] = PI - p[n];
        // x is positive and going down ->  we are past pi/2 -> phase should continue to go up
    }
    else
    {
      if(x[n] < x[n-1])
        p[n] = -PI - p[n];  
        // x is negative and going down -> we are before -pi/2 -> phase should go up
    }
  }
  //rsPlotArrays(N, x, &tmp[0], p);

}

template<class T>
void sigAndAmpToPhase(const T* x, const T* a, int N, T* p)
{
  for(int n = 0; n < N; n++)
    p[n] = asin(x[n] / a[n]);
  //repairPhase(p, N);
  //unreflectPhase(p, N);
  unreflectPhase2(x, p, N);

  //rsPlotArrays(N, p);




  // make a function unreflectPhase - or 2 - one using phase only and one using also the signal

  int dummy = 0;

  // todo: catch div-by-zero and asin of arguments > 1
  // or maybe let the infs and nans just happen and then pass over the p-array again and replace 
  // the invalid sections with numbers resulting from linearly interopolating between surrounding
  // valid values
}
// the measured instantaneous phase seems to be useless - it tends to oscillate back and forth 
// instead of keeping running forward - maybe we need to take the last value into account and 
// reflect ...maybe if p[n] < p[n-1]..well - asin can only return values in the rang -pi/2...pi/2
// but we need values from -pi...pi

template<class T>
void phaseToFreq(const T* p, int N, T* w)
{


  rsAssert(p != w);
  using AT = rsArrayTools;

  //AT::copy(p, w, N);
  //rsPlotArrays(N, w);

  //AT::add(p, 2*PI, w, N);    // w = p + 2*PI
  //rsPlotArrays(N, w);

  //rsPlotArrays(N, p);

  for(int n = 0; n < N; n++)
    w[n] = rsWrapToInterval(p[n], 0.0, 2*PI);
  //rsPlotArrays(N, w);
  // test, if now one of the simpler functions above works, too

  //for(int n = 1; n < N; n++)
  //  w[n] = rsConsistentUnwrappedValue0(w[n], w[n-1], 2*PI); // does not work
  //rsPlotArrays(N, w);

  AT::unwrap(w, N, 2*PI);  // look at code comment there - optimize!
  //rsPlotArrays(N, w);


  AT::difference(w, N);
  rsPlotArrays(N, w, p);
  // has bipolar spikes - could be repaired with 3-point median - but let's first figure out, why
  // they occur - try with a sine input signal - the do occur at the phase jaggies - are these
  // jaggies an artifact - maybe of the unreflecting algo? or are they a legit signal feature? 
  // figure out - maybe resynthesize signal with unreflected phase and see, if we get identity
  // resynthesis - if we do, it doesn't seem to be an artifact

  // has no effect - why? maybe it works only if the start is zero?


}



void testSineParameterEstimation()
{
  testSineAmpAndPhaseEstimation();
  testSineAmpAndPhaseEstimation2();
  // make it a unit test

  using Vec = std::vector<double>;
  Vec x = Vec({10,11,5,0,1,0});  // shows undershooting problem with parabolicTime
  Vec a = x;
  rsAmpEnvelope(&x[0], (int) x.size(), &a[0]); // try to use it in place with a = x
  //rsPlotVectors(x, a);

  x = Vec({10,11,0,0,10,0}); rsAmpEnvelope(&x[0], (int) x.size(), &a[0]); // not interesting
  //rsPlotVectors(x, a);
}

void sineRecreationBandpassNoise()
{
  testSineParameterEstimation();

  // under construction

  // We create a bandpass noise and try to re-create it with a frequency-, phase-, and amplitude 
  // modulated sinusoid

  // once with instantaneous freq-estimation and once with known instantaneous frequency

  int N  = 2000;
  int fs = 44100;
  double f1  = 5000;    // center frequency of input sine at start
  double f2  = 2000;    // center frequency of input sine at end
  double bw1 = f1/20;  // bandwidth in Hz at start
  double bw2 = f2/20;  // bandwidth in Hz at end
  double amp = 0.25;   // maximum amplitude of noise (the normalization level)
  int numPasses = 2;


  using Vec = std::vector<double>;
  using Flt = rsStateVariableFilter<double, double>;

  // generate the sweeping bandpass noise:
  Vec x(N), fa(N);
  rsNoiseGenerator<double> ng;
  Flt flt;
  flt.setSampleRate(fs);
  flt.setMode(Flt::BANDPASS_PEAK);
  int n;
  for(n = 0; n < N; n++) {
    fa[n] = rsLinToLin(double(n), 0.0, N-1.0, f1, f2);  // actual instantaneous center freq
    x[n]  = ng.getSample(); }
  for(int m = 0; m < numPasses; m++) {
    flt.reset();
    for(n = 0; n < N; n++) {
      double bw = rsLinToLin(double(n), 0.0, N-1.0, bw1, bw2); // instantaneous bandwidth
      double bwOct = rsBandwidthConverter::absoluteBandwidthToOctaves(bw, fa[n]);
      flt.setFrequency(fa[n]);
      flt.setBandwidth(bwOct);
      x[n] = flt.getSample(x[n]);
    }
  }
  rsArrayTools::normalize(&x[0], N, amp, true);
  // maybe (optionally) use two equal filters in series



  using SPE = rsSineParameterEstimator<double>;

  // measure instantaneous frequency (with algo 1 - sine recursion formula):
  Vec fm1(N);
  SPE::sigToOmegasViaFormula(&x[0], N, &fm1[0]); 
  fm1 = (fs/(2*PI)) * fm1;

  // Create cleaned up version via 3-point median filter:
  Vec fm1c(N); 
  for(n = 1; n < N-1; n++)
    fm1c[n] = median(fm1[n-1], fm1[n], fm1[n+1]);

  // Measure instantaneous frequency (with algo 2):
  Vec fm2(N);
  for(n = 0; n < N; n++)
    fm2[n] = rsSineFrequencyAt(&x[0], N, n, false) * (fs/(2*PI));

  // Create a median-filtered version of that also:
  Vec fm2c(N);  
  for(n = 1; n < N-1; n++)
    fm2c[n] = median(fm2[n-1], fm2[n], fm2[n+1]);
  // first an last value look wrong - for the moment, just repeat 2nd and 2nd-to-last:
  fm2c[0]   = fm2c[1];
  fm2c[N-1] = fm2c[N-2];
  // ...until we find a better solution...actually, the median-filtering does not change the 
  // estimates very much anyway - at least when numPasses = 2 - ok, with a single pass, the 
  // difference is more obvious


  // For what follows, we need an array of instantaneous frequencies (or omegas):
  Vec f;
  //f = fa;     // use prefectly correct data
  //f = fm1;  // use data from measurement with algo 1
  //f = fm1c;
  //f = fm2;  // use data from measurement with algo 2
  f = fm2c;   // use cleaned up data from measurement with algo 2
  Vec w = (2*PI/fs) * f;

  // Interestingly, for perfect resynthesis, the content of the f- and w-arrays is actually 
  // irrelevant - it could be anything, even noise (->test this). Of course, if it would be 
  // meaningless data, the content of the p- and a-arrays (measured below in the next step) 
  // would be equally meaningless - their meaning would then be only to compensate for the 
  // nonsense of the f-array....
  // w = rsRandomVector(N, -10.0, 10.0); // test -> yep, resynthesis is still perfect


  // measure instantaneous phase and amplitude:
  Vec p(N), a(N);  // maybe use a1, p1
  for(n = 0; n < N-1; n++)
    rsSineAmplitudeAndPhase(x[n], x[n+1], w[n], &a[n], &p[n]);
  // todo: use a symmetric estimation - looking forard and backward and using an average
  // what about the last sample? should we use extrapolation? or is there a similar formula that 
  // uses the current and the previous sample rather than the current and the next?

  // optimize amp-, freq-, and phase-measurements jointly using numeric optimization:
  Vec ao = a, wo = w, po = p;
  //for(n = 2; n < N-2; n++)
  //  rsOptimizeSineParameters(x[n-2], x[n-1], x[n], x[n+1], x[n+2], &ao[n], &po[n], &wo[n]);
  Vec fo = (fs/(2*PI)) * wo;
  // ahh - but with the optimized parameters, we may not get exact resynthesis - if we want exact
  // resynthesis, we should only optimize w and compute a,p as above
  // whoa - the optimization fails when using only one pass


  // re-create the bandpass noise by a freq-, phase- and amp-modulated sine:
  Vec y(N);   // recreated signal 1 - rename to y1
  for(n = 0; n < N; n++)
    y[n] = a[n] * sin(p[n]);

  // Now we want to use the w-array for synthesis, too:
  Vec wi(N);  // integrated w
  Vec pm(N);  // modified p
  rsArrayTools::cumulativeSum(&w[0], &wi[0], N);  // maybe try trapezoidal integration instead
  for(n = 0; n < N; n++)
    pm[n] = rsWrapToInterval(p[n]-wi[n], -PI, PI);
 
  // actual resynthesis
  Vec z(N);   // recreated signal 2 - rename to y2
  for(int n = 0; n < N; n++)
    z[n] = a[n] * sin(wi[n] + pm[n]);


  // Use algo that estimates the amp-envelope first and bases everything else on that:
  Vec a3(N), p3(N), w3(N);
  rsAmpEnvelope(&x[0], N, &a3[0]);
  sigAndAmpToPhase(&x[0], &a3[0], N, &p3[0]);  // may produce nans - fixed
  phaseToFreq(&p3[0], N, &w3[0]);

  rsPlotArrays(N, &x[0], &a3[0]);
  rsPlotVectors(x, a3, p3);






  //rsPlotVectors(fa, fm1, fm1c); // actual and estimated instantaneous freq
  //rsPlotVectors(fa-fm1, 5000.0*x);  // estimation error together with signal for reference

  //rsPlotVectors(test2, fo);

  //rsPlotVectors(fa, fm1, fm1c);
  //rsPlotVectors(fa, fm2, fm2c);


  //rsPlotVectors(fa, fm2c, fo);
  //rsPlotVectors(fa-fm2, 1000.0*x);
  //rsPlotVectors(fa-fm2);

  //rsPlotVectors(fa, fm1c, fo);
  //rsPlotVectors(fa, fm1_2, fo); // this looks close to the "optimal" local approximation



  //rsPlotVectors(fa, fm1, fm2);

  //rsPlotVectors(a, p);

  Vec err1 = x-y;
  Vec err2 = x-z;
  //rsPlotVectors(err1, err2);

  //rsPlotVectors(x, y, x-y, a);  // last sample wrong
  //rsPlotVectors(x, z, x-z, a);  // dito


  //writeToMonoWaveFile("BandpassNoiseOriginal.wav",  &x[0], N, (int) fs, 16);
  //writeToMonoWaveFile("BandpassNoiseRecreated.wav", &y[0], N, (int) fs, 16);


  // Observations:
  // -With numPasses = 1, we get a very erratic (raw) estimate of the instantaneous frequency with 
  //  the whole range of values from zero up to the Nyquist freq - it gets better with more passes
  //  which is what we should expect. It could make sense to set up an upper bound for the estimate
  //  and/or apply a smoothing lowpass to the data.
  // -The frequency estimation errors of algo 1 do not really look like white noise but more like 
  //  spikes around the true frequency. The maxima of the estimation error seem to be at the 
  //  zero-crossings of the pseudo-sine. Maybe the problem is more ill-conditioned around 
  //  zero-crossings - that would make sense. Maybe around zero-crossings, we should distrust the
  //  estimator and tend to just hold the previous value. Actually, the error spikes seem to be one
  //  sample before the zero-crossings - maybe we should use:
  //     wn = rsSineFrequency(x[n-1], x[n], x[n+1]);
  //  for a more symmetric formula? Maybe we get 1 sample latency with the way we are doing it now?
  //  ...done - yep - look at sample 2485
  // -We should perhaps have a reliability measure based on the ratio of the middle sample and the 
  //  average of left and right sample...this should be between 0 and 1 and if it's 1, we use the 
  //  value as is and if it's zero, we use the estimate from the previous sample and if it's in 
  //  between we use a weighted sum...or something
  // -the freq-estimation errors are largest in region where the overall amplitude is low - for 
  //  both algorithms
  // -we could also estimate the freq from the zero-crossings and interpolate
  // -I wonder, why we never see negatiev estimates to the instantaneous freq - is the formula such
  //  that this can't occur? If so - why?
  // -The numerically optimized frequency trajectory closely resembles the cleaned up algo-1 
  //  trajectory, even when the algo-2 output was for the initial guess (!) - with algo-2 
  //  (zero-crossings), the freq-trajectory is smoother (which does not necessarily mean better - 
  //  maybe the noise *should* be in the freq-trajectory and it's oversmoothed with algo-2?). The 
  //  conclusion is that the algo-1 indeed captures the local frequency better - which is not 
  //  surprising since it's based on local information. So, maybe after all, the 3-sample 
  //  estimation with some clean-up is not such a bad idea - the numerical optimization chooses
  //  similar frequencies (and is actually even more erratic). I think, fm1_2 is closest to the 
  //  local optimimum - but fm2c is closer to the underlying actual frequency - so both are better
  //  by different measures. maybe use algo1 with some additional smoothing of the freq (before)
  //  computing a and p.
  // -algo1 has a very distinctive zizag pattern - in each half-cycle of the input, the error
  //  zigzags through one full cycle - could it be related to the ratio of first and second 
  //  derivative?
  // -An entirely different algo could estimate the amplitude first (take abs, connect peaks, etc.)
  //  and compute instantaneous phase by y[n] = a[n] * sin(p[n]) -> p[n] = asin(y[n]/a[n]). The 
  //  instantaneous phase could then be unwrapped and split into instantaneous frequency and
  //  a phase-modulation component: p[n] = sum_k=0^n(w[k]) + pm[k] (maybe the upper limit of the 
  //  sum should be n-1?). This split could be done by a lowpass filter.
  //  -the amplitude seems easier to estimate and may be more "fundamental"
  //  -when we plot the current amplitude estimates, they look kinda jaggy
  //
  // If we use for resynthesis:
  //   y[n] = a[n] * sin(w*n + p[n]);
  // y is totally different from x - why? there seems to be strong aliasing going on - there's a 
  // very audible sweep in the opposite direction. When we use a fixed frequency, the recreated
  // signal is an octave higher when the phase term is included - without the phase-term, the
  // frequency is correct. But when we leave the phase-term out, the sweep is recreated wrongly.
  // OK - when using *only* the phase-term, we get perfect resynthesis. ...hmmm....maybe we should
  // modify the instantaneos phase measurements to take into account the instantaneous frequency
  // measurements? maybe by subtracting the integrated freq-measurements? The perfect resynthesis 
  // works also, if we use the measured frequency in the analysis of amplitude and phase - in fact,
  // we can use *any* frequency value - as long as we take the formual
  //   y[n] = a[n] * sin(p[n])
  // for resynthesis (i.e. without using the frequency information in resynthesis), we get perfect
  // resyntesis. 
  // ToDo: figure out, how the freq-information can be included into the resynthesis without 
  // breaking perfect resynthesis. We want to use the first formula for resynthesis and have to 
  // modify the phase data in such a way that we can. ...wait no: the first formula is wrong! 
  // Instead of using w[n]*n (where w[n] = 2*PI*f[n]/fs) we should use the integral of w[n] up to n


  //rsPlotVectors(x);
  //rsPlotVectors(x, y1, y2);
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
  RAPT::rsArrayTools::fillWithRangeLinear(a, N, as, ae);
  RAPT::rsArrayTools::copy(a, x, N);

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
  RAPT::rsArrayTools::add(xL, xM, x, N);
  RAPT::rsArrayTools::add(x,  xU, x, N);

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
  int    nPeak = RAPT::rsArrayTools::maxIndex(ye, N);
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
  double *x  = new double[N]; RAPT::rsArrayTools::copy(sampleData[0], x, N); // input signal
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
  RAPT::rsArrayTools::fillWithZeros(tmp, Np);
  RAPT::rsArrayTools::copy(x, &tmp[Np], N);
  RAPT::rsArrayTools::fillWithZeros(&tmp[Np+N], Np);
  double a[3], b[3]; // filter coeffs
  rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ
    (b[0], b[1], b[2], a[1], a[2], 1.0/fs, fM, Q);
  RAPT::rsArrayTools::negate(a, a, 3);
  a[0] = 1.0;
  for(int i = 0; i < np; i++)
  {
    //filterBiDirectional(tmp, M, tmp, M, b, 2, a, 2);
    RAPT::rsArrayTools::filterBiDirectional(tmp, M, tmp, M, b, 2, a, 2);
  }
  RAPT::rsArrayTools::copy(&tmp[Np], y, N);
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
  RAPT::rsArrayTools::scale(y, N, 1.0/g);
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
  double *x  = new double[N]; RAPT::rsArrayTools::copy(sampleData[0], x, N); // input signal
  double *y  = new double[N];                                    // filtered signal
  double *ye = new double[N];                                    // envelope
  double *yn = new double[N];                                    // negated envelope

  // extract partial and envelope:
  fM = isolatePartial(x, y, N, fL, fM, fU, fs, aM, aMax, p, np);
  rsSineEnvelopeViaQuadrature(y, ye, N, fM, (double)fs, 4.0);

  // plot the extracted mode with its envelope:
  RAPT::rsArrayTools::negate(ye, yn, N);
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
  RAPT::rsArrayTools::add(xL, xM, x, N);
  RAPT::rsArrayTools::add(x,  xU, x, N);


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
  int    nPeak = RAPT::rsArrayTools::maxIndex(ye, N);
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




