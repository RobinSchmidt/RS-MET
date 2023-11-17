#include "AnalysisExperiments.h"



void autoCorrelation()
{
  static const int bufferSize = 512;
  double sampleRate = 44100;
  double frequency  = 3024.5;
  double startPhase = 0.0;


  int    k = 8;
  double d = 0.5;

  frequency  = sampleRate/(k+d);

  //frequency  = sampleRate/10;    // exact peak at l=10,20,...
  //frequency  = sampleRate/10.1;
  //frequency  = sampleRate/10.2;
  //frequency  = sampleRate/10.5;
  //frequency  = sampleRate/11;  // exact peak at l=11,22,...

  //frequency  = sampleRate/10.5;
  //frequency  = sampleRate/10.25;
  //startPhase = 0.0*PI;
  //startPhase = PI/8;
  //startPhase = PI/4;
  //startPhase = PI/3;
  //startPhase = PI/2;
  //startPhase = 1.0*PI;

  double t[bufferSize];
  //createTimeAxis(bufferSize, t, sampleRate);
  RAPT::rsArrayTools::fillWithIndex(t, bufferSize);
  double buffer[bufferSize];
  int n;
  for(n = 0; n < bufferSize; n++)
    buffer[n] = sin(n*2*PI*frequency/sampleRate + startPhase);

  // obtain sample autocorrelation function (ACF):
  double acf[bufferSize];
  rsAutoCorrelationFFT(buffer, bufferSize, acf);
  RAPT::rsArrayTools::scale(acf, bufferSize, 1.0/acf[0]);

  // obtain unbiased ACF:
  double uacf[bufferSize];
  rsRemoveCorrelationBias(acf, bufferSize, uacf);

  //for(n = 0; n < bufferSize; n++)
  //  uacf[n] = acf[n]/(bufferSize-n);
  //normalize(uacf, bufferSize, 1.0);


  /*
  // obtain a smoothed acf by a 3-point moving average (MA) filter:
  double sacf[bufferSize];
  copy(acf, sacf, bufferSize);
  double h[3] = {1, 1, 0};
  sacf[0] = h[1]*acf[0] + h[2]*acf[1];
  for(n = 1; n <= bufferSize-2; n++)
    sacf[n] = h[0]*acf[n-1]+h[1]*acf[n]+h[2]*acf[n+1];
  sacf[bufferSize-1] = h[1]*acf[bufferSize-1] + h[2]*acf[bufferSize-2];
  scale(sacf, bufferSize, 1.0/sum(h, 3));
    // A kernel of h[-1] = 1, h[0] = 1, h[+1] = 0 seems to be good to avoid the situation where the
    // 2nd peak becomes higher than the 1st for sinusoidal signals at f=fs/(k+1/2) (k integer).
    // The estimate of the exact peak location shifts by 0.5 due to this filter - this must be
    // taken into account when the peak location is computed.
    */

  // Apply a 2-point moving average (MA) filter to the ACF (without the scaling by 1/2). This
  // avoids octave errors when a sinusoid of f = fs/(k+1/2) (k integer) is applied. Later we will
  // need to take accout of this by moving the exact location of the maximum by -0.5:
  double sacf[bufferSize];
  RAPT::rsArrayTools::copy(acf, sacf, bufferSize);
  /*
  sacf[0] = 2*acf[0];
  for(int n = 1; n < bufferSize; n++)
    sacf[n] = acf[n] + acf[n-1];
  scale(sacf, bufferSize, 1.0/2.0); // the missing scaling of the MA
  */

  // find maximum peak index in the smoothed ACF:
  int maxIndex = rsFindHighestPeakIndex(sacf, bufferSize);

  // compute exact value of the index by parabolic interpolation:
  double a[3];
  RAPT::rsPolynomial<double>::fitQuadratic_0_1_2(a, &sacf[maxIndex-1]);
  double offset = 0.0;
  if( a[2] != 0.0 )
    offset = -0.5*a[1]/a[2];
  double exactIndex = (maxIndex-1)+offset;
  //exactIndex -= 0.5; // -0.5 takes the MA filter into account

  //double cOpt = acf[2*k+1] / acf[k];
  double cOpt = uacf[2*k+1] / uacf[k];

  double c;
  //c = double (bufferSize) / (double) (bufferSize*(1+k));
  //double ks = (double) k / (double) bufferSize;
  c  = 1.0+1.0/(k*k);
  c *= c;


  //plotData(50, t, buffer, acf);
  plotData(bufferSize, t, acf, uacf);
  //plotData(40, t, buffer, acf, sacf);
  //int dummy = 0;
}


// Find the deviation between the peak of an autocorrelation sequence for an input sinusoid
// of frequency f = fs/k and f = fs / (k+1/2). In the former case, there's an exact peak at
// lag k, in the latter case, there are two equal values at lag k and k+1 which imply a peak
// halfway in between. We want to find out the difference (or ratio) between these cases.
void autocorrelationPeakVariation()
{
  static const int N      = 512;
  static const int kMin   = 3;
  static const int kMax   = 50;
  static const int kRange = kMax-kMin;

  double t[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  double x1[N];
  double x2[N];

  double d[kRange+1];
  double r[kRange+1];

  // the ratio depends on k, so we create a loop over k in the hope to "see" the function:
  int k;
  for(k = kMin; k <= kMax; k++)
  {
    // create sinusoid x1 with acf-peak at lag k:
    double fs = 44100;
    double f  = fs/k;
    int n;
    for(n = 0; n < N; n++)
      x1[n] = sin(n*2*PI*f/fs);

    // create sinusoid x2 with acf-peak at lag k+1/2:
    f = fs/(k+0.5);
    for(n = 0; n < N; n++)
      x2[n] = sin(n*2*PI*f/fs);

    // obtain (scaled, unbiased) sample autocorrelation of x1 and x2
    double acf1[N];
    rsAutoCorrelationFFT(x1, N, acf1);
    rsRemoveCorrelationBias(acf1, N, acf1);
    RAPT::rsArrayTools::scale(acf1, N, 1.0/acf1[0]);
    double acf2[N];
    rsAutoCorrelationFFT(x2, N, acf2);
    rsRemoveCorrelationBias(acf2, N, acf2);
    RAPT::rsArrayTools::scale(acf2, N, 1.0/acf2[0]);

    d[k-kMin] = acf1[k]-acf2[k]; // difference
    r[k-kMin] = acf1[k]/acf2[k]; // ratio

    //Plotter::plotData(3*k, t, x1, x2);
    //Plotter::plotData(5*k, t, acf1, acf2);
  }

  // plot the difference and ratio as function of k:
  double kValues[kRange+1];
  RAPT::rsArrayTools::fillWithRangeLinear(kValues, kRange+1, (double)kMin, (double)kMax);
  //Plotter::plotData(kRange, kValues, d, rLog);
  //Plotter::plotData(kRange, kValues, r);
  //Plotter::plotData(kRange, kValues, rLog);

  // model the dependence of the difference d on the peak-location k with a rational function:
  double dModel[kRange+1];
  for(k = kMin; k <= kMax; k++)
    dModel[k-kMin] = 4.5/(3.0+k*k);

  // plot difference and its model as function of k:
  plotData(kRange, kValues, d, dModel);
}

void autoCorrelationPitchDetector()
{
  // parameters:
  static const int blockSize  = 256;   // size of the successive input blocks
  //static const int blockSize  = 290;   // size of the successive input blocks
  static const int numBlocks  = 500;   // number of blocks to process
  static const int bufferSize = 1024;  // internal buffersize of the pitch-detector
  //static const int bufferSize = 512;
  //static const int bufferSize = 2048;


  double sampleRate = 44100;
  double frequency  = 160.0; // frequency of the sinuosid


  frequency = 3520.0;
  //frequency = 3620.0;

  // comments about what happens are at blockSize = 256, bufferSize = 1024:
  //frequency = 3000.0;
  //frequency = 3024.0;      // works well
  //frequency = 3024.5;      // jumps back and forth between correct and sub-octave
  //frequency = 3025.0;      // jumps into sub-octave
  //frequency = 3030.0;    // jumps into sub-octave
  //frequency = 3040.0;
  //frequency = sampleRate / 14.5;  // 3041.379...
  //frequency = 3056.0;    // jumps into sub-octave
  //frequency = 3057.0;    // jumps back and forth between correct and sub-octave
  //frequency = 3058.0;    // jumps back and forth between correct and sub-octave
  //frequency = 3059.0;    // works well
    // with a smaller bufferSize (512) it doesn't jump into the sub-octave

  //frequency = sampleRate/11;  // no bias

  //frequency = 4400.0;    // bias about +0.33, with sqrt about -0.21
  //frequency = 4410.0;    // == fs/10 -> no bias
  //frequency = 4420.0;    // bias about -0.33 Hz, with sqrt about +0.21
    // speculative conclusion: the estimate leans towards the frequencies whicha ere a divider of
    // the sample-rate, the parabolic fit to the autocorrelation function has a maximum that leans
    // towards then inside

  //frequency = sampleRate/9;  // no bias
  //frequency = sampleRate / 9.7;

  rsAutoCorrelationPitchDetector<double> pd;
  pd.setBufferSize(bufferSize);
  pd.setUpdateInterval(blockSize);
  pd.setSampleRate(sampleRate);
  pd.setMaxFundamental(8000.0);

  double x[blockSize];  // temporarily holds a block of input data

  double t[numBlocks];  // time-axis (in units of the input blocksize)
  double ft[numBlocks]; // true frequency
  double fe[numBlocks]; // estimated frequency
  double fm[numBlocks]; // running mean of estimated frequency
  double fs[numBlocks]; // exponentially smoothed estimate
  RAPT::rsArrayTools::fillWithIndex(t,  numBlocks);
  RAPT::rsArrayTools::fillWithValue(ft, numBlocks, frequency);

  double w    = 2*PI*frequency/sampleRate;
  double fSum = 0.0;
  int b, n;
  for(b = 0; b < numBlocks; b++)
  {
    for(n = 0; n < blockSize; n++)
      x[n] = sin( w * (b*blockSize+n) );
    fe[b]  = pd.processBlock(x, blockSize);
    fSum  += fe[b];
    fm[b]  = fSum / (b+1);
  }

  // obtain an exponentially smoothed estimated frequency:
  double y     = ft[0];
  double coeff = 0.95;
  for(b = 0; b < numBlocks; b++)
  {
    y     = coeff*y + (1-coeff)*fe[b];
    fs[b] = y;
  }

  // get the maximum error between true and estimated frequency:
  //double maxError = maxDeviation(ft, fe, numBlocks);
  double maxError = RAPT::rsArrayTools::maxDeviation(&ft[10], &fe[10], numBlocks-10);
  double bias     = RAPT::rsArrayTools::mean(&fs[100], numBlocks-100) - ft[0];

  // plot the true frequency, the instantaneous estimate and the running mean of the current
  // estimate:
  plotData(numBlocks, t, ft, fe, fs);
  //plotData(numBlocks, t, ft, fe);
}


void autoCorrelationPitchDetectorOffline()
{
  // We create a linear sine-sweep with known instantaneous frequency at each sample and try to
  // retrieve that instantaneous frequency by measurement as accurate as possible

  int    N  = 2000;   // number of samples in input signal
  double fs = 44100.0; // sample rate
  double f1 = 440.0;   // start value for input signal frequency
  double f2 = 440.0;   // end value for input signal frequency
  double a  =   0.5;   // signal amplitude

  double *f  = new double[N]; // instantaneous frequency of input signal
  double *x  = new double[N]; // input signal
  double *fm = new double[N]; // measured instantaneous frequency

  // create the array with instantaneous frequencies of input signal:
  RAPT::rsArrayTools::fillWithRangeLinear(f, N, f1, f2);

  // create input signal:
  createSineWave(x, N, f, a, fs);

  // measure the instantaneous frequency (preliminary):
  RAPT::rsArrayTools::fillWithZeros(fm, N);
  int hopSize = 1;
  int n;
  for(n = 0; n < N; n += hopSize)
    fm[n] = rsInstantaneousFundamentalEstimator<double>::estimateFundamentalAt(x, N, n, fs, 400.0, 5000.0);

  // OK - this is way off in terms of accuracy - we need the zero-crossing approach

  double *t = new double[N]; // time axis for plot
  createTimeAxis(N, t, 1);
  plotData(N, t, f, fm);
  //plotData(N, t, x);
  delete[] t;

  delete[] f;
  delete[] x;
  delete[] fm;
}

void crossCorrelationBestMatch()
{
  // we create two sinusoids at the same frequency with different start-phases and find the amount
  // of shift required to align them for the best match

  static const int N1 = 2000;  // length of signal 1
  static const int N2 = 2100;  // length of signal 2
  double fs = 44100.0;         // samplerate
  double f  = 300.0;           // frequency of the sinusoids
  double p1 = 1.0;             // start-phase of signal 1
  double p2 = 4.0;             // start-phase of signal 2
  bool deBias = false;         // switch to use unbiased sample cross-correlation function

  // create the input signals:
  double x1[N1], x2[N2];
  double w = 2*PI*f/fs;
  int n;
  for(n = 0; n < N1; n++)
    x1[n] = sin(p1 + n*w);
  for(n = 0; n < N2; n++)
    x2[n] = sin(p2 + n*w);

  // find the maximum of the cross-correlation sequence with subsample precision:
  int N = rsMin(N1, N2);
  double lag     = rsGetShiftForBestMatch(x1, x2, N, deBias);
  double lagTrue = (p2-p1)/w;
  double error   = (lag-lagTrue)/lagTrue;

  // shift the 2nd signal so as to match the 1st:
  double y2[N2];
  rsResampler<double, double>::shiftSinc(x2, y2, N2, lag, 64);

  // plot:
  plotData(N, 0, 1, x1, y2);     // original 1st and shifted 2nd signal

  // Observations: It turns out that the de-biasing may increase the correlation at later lags by
  // amounts that let them be higher than the actually desired maximum which occurs earlier in the
  // sequence. On the other hand, not de-biasing makes the shift values inaccurate.
  // Generally, the precision increases with increasing signal length and decreases with decreasing
  // frequency. Put another way: lower frequency require longer lengths for a given amount of
  // precision hwich is intuitive

  int dummy = 0;
}


// maybe move to RSMath experiments:
void combineFFTs()
{
  static const int N = 8;
  rsComplexDbl x[N], X[N], X1[N/2], X2[N/2], XX[N];

  // create test signal:
  rsRandomUniform(-1.0, 1.0, 0);  // to init the seed
  for(int n = 0; n < N; n++)
  {
    x[n].real(rsRandomUniform(-1.0, 1.0));
    x[n].imag(rsRandomUniform(-1.0, 1.0));
  }

  // compute an N-point FFT:
  RAPT::rsArrayTools::copy(x, X, N);
  rsFFT(X, N);

  // compute two N/2-point FFTs (of the 1st and 2nd half of the buffer):
  RAPT::rsArrayTools::copy(x, X1, N/2);
  rsFFT(X1, N/2);
  RAPT::rsArrayTools::copy(&x[N/2], X2, N/2);
  rsFFT(X2, N/2);

  // combine the two N/2 point FFTs into one N-point FFT:
  rsComplexDbl wN2 = exp( -rsComplexDbl(0, 1)*rsComplexDbl(PI/N) );
  rsComplexDbl w;
  for(int k = 0; k < N/2; k++)
  {
    w = pow(wN2, rsComplexDbl(k+1));

    XX[2*k]   = X1[k] + X2[k];
    XX[2*k+1] = w*XX[2*k];       // ...the odd-numbered bins are still wrong
                                 // i think, it doesn't work out anyway
  }
}

void envelopeFollower()
{
  double fs = 2000;
  double f  = 100.5;
  int N = 4000;
  int nRelease = 3*N/4;

  // create enveloped sawtooth wave as input signal:
  rsBreakpointModulatorD bm;
  bm.setSampleRate(fs);
  std::vector<double> x(N), e(N);
  createWaveform(&x[0], N, 1, f, fs, 0.0, true); // use anti-aliasing to see Gibbs ripple effects
  //createWaveform(&x[0], N, 2, f, fs, 0.0, true); // use anti-aliasing
  int n;
  bm.noteOn(false, 64, 64);
  for(n = 0; n < nRelease; n++) {
    e[n]  = bm.getSample();
    x[n] *= e[n];;
  }
  bm.noteOff();
  for(n = nRelease; n < N; n++) {
    e[n]  = bm.getSample();
    x[n] *= e[n];;
  }

  // new - do the same thing but now with the new convenience class:
  rsEnvelopeFollower2<double> ef2;
  ef2.setSampleRate(fs);
  ef2.setInputFrequency(f);
  //std::vector<double> e1(N);
  std::vector<double> y1(N), y2(N), y3(N), y4(N);
  for(n = 0; n < N; n++) {
    y1[n] = ef2.getSamplePreFilteredAbs(x[n]);
    y2[n] = ef2.getSampleSlewLimited(y1[n]);
    y3[n] = ef2.getSampleMinMaxSmoothed(y2[n]);
    y4[n] = ef2.getSamplePostFiltered(y3[n]);
    // y4 is the same as if we would just have called y4[n] = ef2.getSample(x[n]) instead of the
    // four calls above
  }
  // todo: try different settings for attack/release, see if the affects the delay, etc. 

  int smoothingLength = (int) ceil(fs/f);  // one cycle - todo: inquire from ef2 object
  int shiftAmount = 3*smoothingLength/2;   // factor 3/2 ad hoc
  RAPT::rsArrayTools::leftShift(&y3[0], N, shiftAmount);
  RAPT::rsArrayTools::leftShift(&y4[0], N, shiftAmount+5);
  // todo: figure out, if these shifts are generally applicable or specific to the settings and/or
  // input signal


  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &e[0]);   // true envelope for reference


  // new - obtained by the convenience class:
  //plt.addDataArrays(N, &y1[0]);
  //plt.addDataArrays(N, &y2[0]);
  //plt.addDataArrays(N, &y3[0]);
  plt.addDataArrays(N, &y4[0]);  

  plt.setPixelSize(1200, 400);
  plt.plot();

  // Observations: 
  // -the anti-aliasing makes the maximum excursion of each cycle different, resulting in an 
  //  undesired modulation of the detected envelope
  // -maybe try instantaneous envelope detecto based on a 90° filter
  //  -hmm - the raw instantaneous envelope is really bad - but maybe applying a lowpass or
  //   slew rate limiter to that could work?
  // -maybe we should suppress the Gibbs-ripples by lowpassing the input signal befor going into
  //  the envelope detector...maybe Butterworth or Bessel filter
  // -maybe try to use a power p at the input and its inverse power 1/p at the output of the 
  //  detector if p < 1, the peaks get less influence - does not seem to be useful
  // -maybe use a post-filter, too (Bessel)
  // -i think, instead of the current "hold" implementation, we need to apply a sort of 
  //  moving-maximum filter to the output of the slewrate limiter (which then doesn't need hold)
  //  -how would we implement that efficiently. i.e. without searching for a maximum in a circular
  //   buffer each sample? maybe on-the-fly decimation?

  // -OK, the factor 1.23 works also for signals that don't show Gibbs ripples - good!
  // -but: it doesn't work for other waveforms - for the sine, we get too large envelope values
  //  -presumable, because the sine-wave is not really affected by the pre-filter
  //
  // ToDo:
  // -For the basic envelope follower, try to achieve different curve shapes. Currently, we have 
  //  the typical exponential shape of the 1-pole filter based on a update equation like
  //  y[n] = c1 * y[n-1]. Maybe try y[n] = max(0, y[n-1] - c0) for linear decay down to zero or
  //  y[n] = max(0, c0 + c1*y[n-1] + c2*y[n-1]*y[n-1])  where c1 is the exponential decay, c0 is a
  //  negative number that is responsible for linear decay and c2 is for changing the shape to make
  //  it hopefully smoother. Or maybe c1 should somehow depend on y[n-1]: when y[n-1] is big, c1 
  //  should be small, I think. This should smooth the upper section such that it doesn't drop off
  //  with a corner. Or what about y[n] = max(0, c0 + c1*y[n-1] + c2*sqrt(y[n-1]))? 

}


void instantaneousFrequency()
{
  // we create a linear frequency sweep of a sinusoidal waveform and measure the instantaneous
  // frequency with rsSineFrequencyAt

  static const int N = 10000;  // number of samples
  double fs = 44100;           // sample rate
  double a  = 0.5;             // amplitude (maybe use a time-varying amplitude later)
  double f1 = 2000;            // start frequency
  double f2 = 1000;            // end frequency

  double smalll = 1.e-8;

  double f[N];                 // true instantaneous frequencies
  double fm1[N];               // measured instantaneous frequencies with core algorithm
  double fm2[N];               // ... with robustified algorithm
  double fm3[N];               // ... with robustified, refined algorithm
  double e1[N], e2[N], e3[N];  // relative measurement errors
  double x[N];                 // signal

  // create frequency array and signal:
  RAPT::rsArrayTools::fillWithRangeLinear(f, N, f1, f2);
  createSineWave(x, N, f, a, fs);

  // measure instantaneous frequency and obtain error:
  for(int n = 0; n < N; n++)
  {
    fm1[n] = rsSineFrequencyAtCore(x, N, n, smalll); // normalized radian frequency
    fm2[n] = rsSineFrequencyAt(x, N, n, false);
    fm3[n] = rsSineFrequencyAt(x, N, n, true);

    fm1[n] = fs*fm1[n] / (2*PI);                    // frequency in Hz
    fm2[n] = fs*fm2[n] / (2*PI);
    fm3[n] = fs*fm3[n] / (2*PI);

    e1[n]  = (f[n] - fm1[n]) / f[n];
    e2[n]  = (f[n] - fm2[n]) / f[n];
    e3[n]  = (f[n] - fm3[n]) / f[n];
  }

  double maxError = RAPT::rsArrayTools::maxAbs(e3, N);
    // c = 0.76: 3.6285047484070877e-005
    // c = 0.75: 3.4832150523579422e-005
    // c = 0.72: 3.0473459643206141e-005

  double meanError = RAPT::rsArrayTools::mean(e3, N);
    // c = 0.66: -1.2387131128786400e-006
    // c = 0.666: -2.4751832800983252e-007

  //plotData(N, 0, 1, f, fm1); // true and measured frequency with core algorithm
  //plotData(N, 0, 1, f, fm2); // true and measured frequency with robust algorithm
  //plotData(N, 0, 1, x, e1, e2);    // signal and measurement errors
  //plotData(N, 0, 1, x, e2);    // signal and measurement error of robust algorithm
  plotData(N, 0, 1, e2, e3);    // measurement error of robust algorithm

  // Observations:
  // with core algorithm (fm1):
  // When the sinusoid is not stable (i.e. it sweeps), there are errors in the measurement. These
  // errors are largest at the zero crossings of the signal (and there, they may become excessive
  // indeed) but are small at the peaks (positive or negative)
  // solution: for any sample, find the peaks before and above that sample, measure there
  // and take the average, this is the fm2 meausrement. In this case, the measurements are musch
  // more stable over time but there is a little bias: the freqeuncies are slightly underestimated
  // in case of an upward sweep and slightly overestimated in case of a downward sweep.

  // more ideas:
  // (we may use an arithmetic or geometric mean - the choice will determine
  // whether we will be more exact in case of linear or exponential sweeps - maybe we can even
  // use a kind of mixed mean (geometric/arithmetic) with a user parameter to adjust the mix)
  // the average should be a weighted average, weighted by the distances to the peaks (maybe
  // we can even use subsample-precision peak locations)

  int dummy = 0;
}

void instantaneousPhase()
{
  static const int N = 10000;  // number of samples
  double fs = 44100;           // sample rate
  double a1 = 0.5;             // start amplitude
  double a2 = 0.5;            // end amplitude
  double f1 = 1000;            // start frequency
  double f2 = 2000;            // end frequency
  double p0 = 0.0;            // start phase

  double w[N];                 // (true) instantaneous normalized radian frequencies
  double wm[N];                // measured values
  double we[N];                // relative measurement error
  double a[N];                 // (true) instantaneous amplitudes
  //double am[N];                // measured values
  //double ae[N];                // relative measurement error
  double p[N];                 // (true) instantaneous phases
  double pm1[N];               // measured values with algorithm based on euqations
  double pe1[N];               // relative measurement error
  double pm2[N];               // measured values with algorithm based on zeros crossings
  double pe2[N];               // relative measurement error
  double x[N];                 // signal

  int precision = 2;           // precision parameter for zero-crossing based measurement

  // fill arrays with true instantaneous (normalized radian) frequencies, phases and amplitudes:
  int n;
  RAPT::rsArrayTools::fillWithRangeLinear(w, N, f1, f2);
  for(n = 0; n < N; n++)
    w[n] *= 2*PI/fs;
  RAPT::rsArrayTools::fillWithRangeLinear(a, N, a1, a2);
  p[0] = p0;
  for(n = 1; n < N; n++)
    p[n] = rsWrapToInterval(p[n-1] + w[n], -PI, PI);

  // create signal:
  for(n = 0; n < N; n++)
    x[n] = a[n] * sin(p[n]);

  // measure instantaneous parameters:
  for(n = 0; n < N; n++)
  {
    wm[n] = rsSineFrequencyAt(x, N, n, true);
    pm1[n] = rsSinePhaseAt(x, N, n, wm[n]);
    pm2[n] = rsSinePhaseAtViaZeros(x, N, n, precision);
    //am[n] = rsSineAmplitudeAt(x, N, n, wm[n]); //we don't have this function yet
  }

  // compute errors
  for(n = 0; n < N; n++)
  {
    we[n] = (w[n] - wm[n]) / w[n];  // relative error
    pe1[n] = (p[n] - pm1[n]);       // absolute error
    pe2[n] = (p[n] - pm2[n]);
  }


  //plotData(N, 0, 1/fs, x);
  //plotData(N, 0, 1, w, wm);  // true and measured normalized radian frequency
  //plotData(N, 0, 1, p, pm1, pm2);  // true and measured instantaneous phase
  //plotData(N, 0, 1, pe1, pe2);       // phase error
  plotData(N-3, 0, 1, pe1, pe2);   // last 3 are outliers - need to be cosidered as special case later


  int dummy = 0;
}

void maxShortTimeRMS()
{
  int N = 100000;
  std::vector<double> x = createSineWave(N, 1000.0, 48000.0);
  double maxRms = RAPT::getMaxShortTimeRMS(&x[0], N, 4800); // 48 samples = 1 cycle, rms = 1/sqrt(2)
  int dummy = 0;
}

void arrayRMS()
{
  int N = 100000;
  std::vector<double> x = createSineWave(N, 1000.0, 48000.0);
  double rms = RAPT::rsArrayTools::rootMeanSquare(&x[0], N); // 48 samples = 1 cycle, rms = 1/sqrt(2)
  int dummy = 0;
}

void peakFinder()
{
  // Create a sinuosid and find its minima and maxima with subsample precision

  int N = 50;
  int oversampling = 40;
  int precision = 4;
  double w = 1.5;  // omega

  // Create the test signal:
  int No = N*oversampling;
  using Vec = std::vector<double>;
  Vec t(N), to(No);
  Vec x(N), xo(No);
  for(int n = 0; n < N; n++)  {
    t[n] = n;
    x[n] = sin(w*t[n]); }
  for(int n = 0; n < No; n++) {
    to[n] = n;
    to[n] /= oversampling;
    xo[n] = sin(w*to[n]);  }

  // Find the minima and maxima:
  using PF = rsPeakFinder<double>;
  using AT = rsArrayTools;
  Vec peakPositions, peakHeights;
  double pos, height;
  for(int n = 1; n < N-1; n++) {
    if( AT::isPeakOrValley(&x[0], n) ) {
      PF::exactPeakPositionAndHeight(&x[0], N, n, precision, &pos, &height);
      peakPositions.push_back(pos);
      peakHeights.push_back(height); }}

  // Plot underlying (pseudo) continuous signal, the sampled signal and the peaks that were 
  // obtained from the sampled signal. They should be in the positions of the actual peaks of the
  // underlying continuous signal:
  GNUPlotter plt;
  plt.addDataArrays(N,  &t[ 0], &x[ 0]);
  plt.addDataArrays(No, &to[0], &xo[0]);
  plt.addDataArrays((int)peakHeights.size(), &peakPositions[0], &peakHeights[0]);
  plt.setGraphStyles("lines", "lines", "points");
  plt.setPixelSize(1200, 400);
  plt.plot();

  // Observations:
  // -With precision = 1 (parabolic) an w = 1, the height error is around 2%, with precision = 0
  //  around 10% - the eestimation error indeed decreases with increasing precision - it works!
}

// convenience function to make the zero-crossing finding work for plain arrays (as required for
// plotting)
void upwardZeroCrossings(double *x, int N, double *z, int Nz, int p)
{
  std::vector<double> zv = rsZeroCrossingFinder::upwardCrossings(x, N, p);
  for(int n = 0; n < rsMin(Nz, (int)zv.size()); n++)
    z[n] = zv[n];
}
void zeroCrossingFinder()
{
  // We create a linear sine-sweep for which the true time-instants of the zero crossings can be
  // computed analytically and compare the measured values with these true values.

  static const int N  = 10000; // number of samples
  static const int Nh = 5;     // number of harmonics
  double fs = 44100;           // samplerate in Hz
  double f1 = 2000.0;          // signal frequency at start
  double f2 = 4000.0;          // signal frequency at end
  double p0 = 1.0;             // start phase
  double A[Nh];                // amplitudes of harmonics
  double *x = new double[N];   // input signal

  // fill array with harmonic amplitudes:
  A[0] = 1.0;
  A[1] = 0.0;
  A[2] = 0.3;
  A[3] = 0.0;
  A[4] = 0.0;

  // normalized radian frequencies:
  double w1 = 2*PI*f1/fs;
  double w2 = 2*PI*f2/fs;

  // parameters of linear function w(t) = a*t + b
  double b  = w1;              // offset
  double a  = (w2-w1)/(N-1);   // slope

  // create signal:
  int n;
  for(n = 0; n < N; n++)
    x[n] = sineSum(a*n*n + b*n + p0, A, Nh);

  // allcate memory for true and measured zero crossings:
  int Nz = rsZeroCrossingFinder::numUpwardCrossings(x, N);
  double *zt = new double[Nz]; // true zero crossing values
  double *z0 = new double[Nz]; // measurement with precision 0 (linear)
  double *z1 = new double[Nz]; // measurement with precision 1 (cubic)
  double *z2 = new double[Nz]; // measurement with precision 2 (quintic)
  double *z3 = new double[Nz]; // measurement with precision 3 (septic)

  // compute true zero-crossings:
  double p = b/a;
  double q;
  int k = 0;
  n = 0;
  while( k < Nz )
  {
    //q = (p0-k*PI)/a;  // finds all zero crossings
    q = (p0-2*k*PI)/a; // finds only upward zero crossings
    if( q < 0.0 )
    {
      zt[n] = -p/2 + sqrt(p*p/4 - q);
      n++;
    }
    k++;
  }
  Nz = n; // ...why?

  // do the measurements:
  upwardZeroCrossings(x, N, z0, Nz, 0);
  upwardZeroCrossings(x, N, z1, Nz, 1);
  upwardZeroCrossings(x, N, z2, Nz, 2);
  upwardZeroCrossings(x, N, z3, Nz, 3);

  // compute measurement errors
  double *e0 = new double[Nz]; // error for precision 0 (linear)
  double *e1 = new double[Nz]; // error for precision 1 (cubic)
  double *e2 = new double[Nz]; // error for precision 2 (quintic)
  double *e3 = new double[Nz]; // error for precision 3 (heptic)
  RAPT::rsArrayTools::subtract(z0, zt, e0, Nz);
  RAPT::rsArrayTools::subtract(z1, zt, e1, Nz);
  RAPT::rsArrayTools::subtract(z2, zt, e2, Nz);
  RAPT::rsArrayTools::subtract(z3, zt, e3, Nz);

  // find maximum errors:
  double eMax0 = RAPT::rsArrayTools::maxAbs(e0, Nz);
  double eMax1 = RAPT::rsArrayTools::maxAbs(e1, Nz);
  double eMax2 = RAPT::rsArrayTools::maxAbs(e2, Nz);
  double eMax3 = RAPT::rsArrayTools::maxAbs(e3, Nz);

  // compute error ratios - these are the precision improvement factors by which choosing a higher
  // precision value actually affects the precision
  double e10 = eMax0/eMax1;
  double e21 = eMax1/eMax2;
  double e32 = eMax2/eMax3;

  // plot:
  //plotData(N, 0, 1, x);
  plotData(Nz, 0, 1, e0, e1, e2, e3);

  // Observations:
  // For a purely sinusoidal sweep (A[0]=1, A[n]=0 for n>0), the error increases linearly with
  // frequency. Using higher precision values leads to a frequency dependent precision improvement.
  // For a sweep between 200 and 400 Hz, the improvement from linear to cubic e10 is around 559,
  // the improvent from cubic to quintic e21 is around 516 and the improvement from quintic to
  // heptic is around 149. For a sweep from 2000 to 4000 Hz: e10=7.0, e21=6.6, e32=6.5 - so the
  // precision improves more slowly at higher frequencies. Generally, it can be said that
  // incrementing the precision setting by 1 will typically increase the actual precision by some
  // factor which will be higher for lower frequencies. It's around of the order of 5 for
  // frequencies around 4kHz and a few hundreds for frequencies around 400Hz. I tend to think, for
  // typical input frequencies around a few hundreds to a few thousands of Hertz, the quintic
  // (p=2) precision setting should be adequate.
  // If there are overtones present, the precision is generally much worse and also does increase
  // more slowly when going higher with the polynomial order. Also, the error increase with respect to
  // frequency looks not like a linear function anymore but more curved. If only a single overtone
  // is present, the error increases with the amplitude of that overtone.

  delete[] x;
  delete[] zt;
  delete[] z0;
  delete[] z1;
  delete[] z2;
  delete[] e0;
  delete[] e1;
  delete[] e2;
}

void zeroCrossingFinder2()
{
  // We just create a small array of -1, +1 values and see, if it finds the correct locations
  // (maybe, instead this could be put into a unit test and maybe use a somewhat longer signal with
  // more zero crossings in the signal (so that we don't just have some at the borders where the
  // detector is forced use a linear interpolant only)

  int p = 3; // precision

  static const int N = 10;
  double x[N];
  RAPT::rsArrayTools::fillWithValue(x,       N/2, +1.0);
  RAPT::rsArrayTools::fillWithValue(&x[N/2], N/2, -1.0);
  x[0]   = -1.0;
  x[N-1] = +1.0;

  // find zero crossings:
  std::vector<double> z = rsZeroCrossingFinder::upwardCrossings(x, N, p);
   // z should be: 0.5, 8.5
  int dummy = 0;
}

void zeroCrossingFinder3()
{
  // Demonstrates the edge-case where we have sequences of many zeros in the signal.

  // user parameters:
  static const int N  = 20000;  // number of samples
  double fs = 44100;            // samplerate in Hz
  double f  = 100.0;            // signal frequency
  double thresh = 0.5;          // threshold below which values are set to zero
  int prec = 1;                 // precision

  vector<double> x = sineAndDeacyingInharmonic(N, f, fs, 0);
  for(int n = 0; n < N; n++) {
    if(fabs(x[n]) < thresh)
      x[n] = 0;
  }
  std::vector<double> z = rsZeroCrossingFinder::upwardCrossings(&x[0], N, prec);
  vector<double> zy(z.size()); // all zeros

  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays((int)z.size(), &z[0], &zy[0]);
  plt.setGraphStyles("lines", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();

  // Observations:
  // -the detected zero-marks are at the starts of the zero-sequences
  // ..we may want them in the centers...
}

// todo: create a contrived test signal consisting of 2 inharmonic sines and use:
// if(abs(x[n] < thresh)
//   x[n] = 0;
// i.e. clamp small values to zero - see, if the zero markers are placed intuitively - try both
// conventions x[n-1] < 0 && x[n] >= 0 or x[n-1] <= 0 && x[n] > 0

/*
std::vector<double> twoTonesAndDecayingDc(int N, double f, double fs, double overtoneRatio,
  double overtoneAmplitude, double dcAmount, double dcDecay)
{
  double w = 2*PI*f/fs;
  vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = sin(w*n) + overtoneAmplitude * sin(overtoneRatio*w*n) + dcAmount * exp(-(n/fs)/dcDecay);
  return x;
}
*/

void cycleMarkFinder()
{
  // user parameters:
  static const int N  = 50000;  // number of samples
  double fs = 44100;            // samplerate in Hz
  double f  = 100.0;           // signal frequency
  double corrLength = 1.0;      // length of correlation (in terms of cycles)
  //fs = 44000;                  // test: cycle exactly 44 samples long
  //fs = 44300;
  //fs = 44500;                  // test: cycle exactly 44.5 samples long
  //fs = 44700;
  //fs = 45000;

  double period = fs/f;

  // create test input signal:
  vector<double> x;
  //x = createSineWave(N, f, fs);   // sine wave at frequency f
  //x = sineAndDeacyingInharmonic(N, f, fs, 0);  // sine at f + sine at f*(1+sqrt(5))
  //x = sineAndDeacyingInharmonic(N, f, fs, 20); // sine at f + decaying sine at f*(1+sqrt(5))
  //x = cycleMarkTestSignal(N, f, fs, 20, 5);
  //x = twoSinesAndDecayingDc(N, f, fs, 15, 0.1, 0.9, 0.2);
    // with 0.2 for the overtone amplitude, the cycle-mark finder thinks, the fundamental
    // is at the overtone - bad!
  //x = twoSinesAndDecayingDc(N, f, fs, 9.9, 0.2, 0.0, 0.2);
  x = sawAndSquare(N, fs, f, 1.0, 10.1*f, 0.2, false);



  // find cycle marks by different algorithms:
  typedef rsCycleMarkFinder<double> CMF;
  CMF cmf(fs, 20, 5000);
  vector<double> cm1, cm2;
  cmf.setRelativeCorrelationLength(corrLength);
  //cmf.setRelativeCorrelationHighpassFreq(0.5);
  cmf.setSubSampleApproximationPrecision(0);  // linear
  cmf.setAlgorithm(CMF::F0_ZERO_CROSSINGS);
  cm1 = cmf.findCycleMarks(&x[0], N);
  //cmf.setAlgorithm(CMF::WINDOWED_CORRELATION);
  cmf.setAlgorithm(CMF::CYCLE_CORRELATION);
  //cmf.setAlgorithm(CMF::ZERO_CROSSINGS);
  //cmf.setAlgorithm(CMF::CORRELATED_ZERO);
  cm2 = cmf.findCycleMarks(&x[0], N);

  vector<double> deltas(cm1.size());
  //RAPT::rsArrayTools::subtract(&cm1[0], &cm2[0], &deltas[0], (int)cm1.size());

  // compute average distances between the cycle-marks - the correct value would be fs/f for a
  // periodic input
  double dAvg1, dAvg2;
  dAvg1 = RAPT::rsArrayTools::meanDifference(&cm1[0], (int)cm1.size());
  dAvg2 = RAPT::rsArrayTools::meanDifference(&cm2[0], (int)cm2.size());
  // maybe compute the variance of the difference, too - maybe with a steady inharmonic input, this
  // value may say something about the stability of the pitch estimate over time in the presence
  // of inharmonicity - maybe a good cycle-mark finder should give a low variance?
  // maybe compute the min/max values of the differences and the maximum absolute deviation from
  // the true/desired value, maybe make a data structure, rsCycleMarkQualityMeasures (maybe one
  // that can be returned from rsCycleMarkFinder itself (using as input an array of cycle marks and
  // a period length)


  // plot signal and cycle marks:
  int Nz = (int) std::max(cm1.size(), cm2.size()); // # cycle marks
  vector<double> cmy(Nz);    // y values for plotting (all zero)
  GNUPlotter plt;

  // there seems to be a small bias for the cross-correlation approach leading of a small drift
  // of the cycle marks (unless the period length happens to be integer). seems like the length
  // is estimated systematically too short (verify this)

  // check, if there are any (asymmetric) biases in the cross-correlation function...but maybe
  // even a sysmmetric bias could have a bad effect?

  // hmm . currently, it seems that increasing the corrLength has a negative effect on the quality
  // of the marks (3 is worse than 1, for example) - investigate, why - seems, it goes away when
  // increasing the length (something about the mea-compuatation?)

  // todo: the CYCLE_CORRELATION algorithm may fail due to the refined marks overrunning
  // the unrefined ones
  // -find a test scenario where this happens
  // -fix it by not considering the f0-zero-crossings as preliminary estimates but starting
  //  from scratch
  // -begin in the middle -> set a cycle mark there
  // -look for the next cycle mark (left and right) near the position that would be predicted
  //  by periodicity (find f0 estimate via autocorrelation), then refine that position by
  //  cyclic cross correlation
  // -do the same thing for the next neighbour cycle (for prediction, we can now use the length
  //  of the previous cycle)
  // -this method has some ambiguity to it - we could have started somewhere els near the center
  //  ...maybe it makes sense to shift the start-position such that the first cycle mark has the
  //  same distance to sample 0 as the distance between the first and second...or something
  //  ...but maybe close to the beginning, we should not set cycle-marks at all because this is
  //  the transient and very probably aperiodic...or maybe carry along a reliability for each
  //  cycle mark given by the cross-correlation value

  // more ideas:
  // -use a nonlinear (monotonic) function on the autocorrelation values before fitting a parabola
  // -use a quartic instead of quadratic parabola
  // -use step-size - don't jump one cycle forward but M and put M-1 markers in between at equal
  //  distances
  // -when it's all done, check, how well the first and last cycle correlate. is there an offset?
  //  if so, can this be compensated by multiplying all cycle-lengths by a factor?

  //plt.addDataArrays(deltas.size(), &deltas[0]); // test
  //plt.plot();

  // 1: blue crosses, 2: green stars

  plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays((int)cm1.size(), &cm1[0], &cmy[0]);
  plt.addDataArrays((int)cm2.size(), &cm2[0], &cmy[0]);
  plt.setGraphStyles("lines", "points", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
}
void cycleMarkErrors()
{
  // user parameters:
  static const int N  = 40000;   // number of samples
  double fs = 44100;             // samplerate in Hz
  double minPeriod  = 40;        // minimum signal period in samples
  double maxPeriod  = 50;
  double corrLength = 1.0;       // length of correlation (in terms of cycles)
  //int numPeriods    = 41;        // number of signal periodicities between min and max
  //int numPeriods = 101;
  //int numPeriods = 201;
  int numPeriods = 401;

  // maybe have a minPeriod and maxPeriod, for example 99..101 and a stepsize and check for various
  // periods in between (99.0, 99.1, 99.2, ..., 100.9, 101.0) and plot the errors as function
  // of the period


  typedef rsCycleMarkFinder<double> CMF;
  CMF cmf(fs, 20, 5000);
  cmf.setRelativeCorrelationLength(corrLength);
  cmf.setSubSampleApproximationPrecision(2); // 0: linear, 1: cubic, ...
  vector<double> x;
  vector<double> cm1, cm2;
  vector<double> periods, meanErrors1, meanErrors2, maxErrors1, maxErrors2, minErrors1, minErrors2,
    maxAbsErrors1, maxAbsErrors2;
  for(int i = 0; i < numPeriods; i++)
  {
    // create test input signal:
    double period = minPeriod + i*(maxPeriod-minPeriod)/(numPeriods-1);
    double f = fs/period;           // signal frequency
    x = createSineWave(N, f, fs);   // sine wave at frequency f

    // find cycle marks by different algorithms:
    cmf.setAlgorithm(CMF::F0_ZERO_CROSSINGS);
    cm1 = cmf.findCycleMarks(&x[0], N);
    //cmf.setAlgorithm(CMF::CYCLE_CORRELATION);
    cmf.setAlgorithm(CMF::CYCLE_CORRELATION);
    cm2 = cmf.findCycleMarks(&x[0], N);

    // get errors:
    CMF::ErrorMeasures errors1, errors2;
    errors1 = cmf.getErrorMeasures(cm1, period);
    errors2 = cmf.getErrorMeasures(cm2, period);

    // add measured errors to data arrays for plotting
    periods.push_back(period);
    meanErrors1.push_back(errors1.mean);
    meanErrors2.push_back(errors2.mean);
    maxAbsErrors1.push_back(errors1.maxAbs);
    maxAbsErrors2.push_back(errors2.maxAbs);
  }

  GNUPlotter plt;
  int M = (int)periods.size();
  plt.addDataArrays(M, &periods[0], &meanErrors1[0]);
  plt.addDataArrays(M, &periods[0], &meanErrors2[0]);
  //plt.addDataArrays(M, &periods[0], &maxAbsErrors1[0]);
  //plt.addDataArrays(M, &periods[0], &maxAbsErrors2[0]);
  plt.plot();

  // Observations:
  // -the correlation algorithm seems to have a bias towards overestimating the period length
  //  when the correlation length is 1, it seems to go down with longer lengths, at 3, it's
  //  disapperaed
  //  -try different window functions
}

void applyBellFilter(double *x, double *y, int N, double f, double fs, double b, double g)
{
  // create and set up bell filter:
  rsStateVariableFilterDD flt; // we may use a biquad later as well (c/p legacy code to RSLib)
  flt.setMode(rsStateVariableFilterDD::BELL);
  flt.setFrequency(f);
  flt.setSampleRate(fs);
  flt.setBandwidth(b);
  flt.setGain(rsDB2amp(g));
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(x[n]);
}
void zeroCrossingPitchDetector()
{
  // we create a test signal with a given frequency (pulse-wave with 2 formants) and try to
  // measure the instantaneous frequency at each sample with different algorithms

  static const int N = 30000;  // number of samples
  double fs  = 44100;          // samplerate in Hz
  double f   = 100.0;          // signal frequency
  double *x  = new double[N];  // signal
  double *f1 = new double[N];  // measured instantaneous frequency by realtime algo
  double *f2 = new double[N];  // ...realtime algo on preprocessed signal
  double *f3 = new double[N];  // new algorithm - interpolate before convert
  int n;

  // synthesize input signal - a pulse-wave with 2 formants:
  synthesizePulseWave(x, N, f, 0.25, fs, 0.0, false);
  //synthesizePulseWave(x, N, f, 0.5, fs, 0.0, true);
  applyBellFilter(x, x, N,  750, fs, 0.2, 15);
  applyBellFilter(x, x, N, 1500, fs, 0.2, 15);
  RAPT::rsArrayTools::scale(x, N, 0.25);
  //writeToMonoWaveFile("PitchDetectorInput.wav",  x, N, (int) fs, 16);

  // get initial estimate of fundamental by using an autocorrelation based algorithm at the center
  // of the input signal:
  double fi = rsInstantaneousFundamentalEstimatorD::estimateFundamentalAt(x, N, N/2, fs, 20.0,
    5000.0);

  // detect the instantaneous frequency with the realtime algorithm, intialized with our initial
  // estimate - this gives our f1 array:
  rsZeroCrossingPitchDetectorD pd;
  pd.reset(fi);
  for(n = 0; n < N; n++)
    f1[n] = pd.estimateFundamentalFrequency(x[n]);

  // apply bidirectional bandpass filter to the signal, tuned to the estimated fundamental:
  double *y = new double[N];

  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y, N, fi, fi, fs, 3);
    // bw and np chosen by a bit of trial and error - experiment more to find sweet spot

  // use realtime algo again, but this time on the filtered signal:
  pd.reset(fi);
  for(n = 0; n < N; n++)
    f2[n] = pd.estimateFundamentalFrequency(y[n]);

  // use the new nonrealtime algo:
  rsInstantaneousFundamentalEstimatorD::measureInstantaneousFundamental(x, f3, N, fs, 20.0, 5000.0);


  double *ft = new double[N];  // target value
  RAPT::rsArrayTools::fillWithValue(ft, N, f);
  //plotData(N, 0, 1/fs, x); // plot original signal
  //plotData(N, 0, 1/fs, y); // plot filtered signal
  //plotData(N, 0, 1/fs, x, y); // plot original and filtered signal
  plotData(N, 0, 1/fs, ft, f1, f2, f3);  // plot actual frequency and measurements
  //plotData(N, 0, 1/fs, ft, f3);  // plot frequency and estimate
  delete[] ft;

  // Observations:

  // f1 - the raw frequency etsimate from the realtime algo:
  // depending on the input frequency, we may see a periodic error in the raw
  // frequency estimate f1 obtained by the realtime algorithm - for example, if f=200, the estimate
  // oscillates between 199.94 and 200.06 in a squarewave-like fashion, around the correct value of
  // 200. For f=220 or f=250, the pattern is more complicated. It seems to be a good idea to apply
  // a moving-average filter adapted to the length of the period of this oscillation - but how do
  // we determine this period? i think, it is a multiple of the signal's period - but what
  // multiple? Or, maybe we should just apply (non-ringing) lowpass smoothing filter to the raw
  // pitch-detector output.

  // f2 - the filtered input signal passed to the realtime algo:
  // for bw=0.3, np=5, f=220: the oscillations are greatly reduced, but we have problems at the
  // edges of the signal - possibly due to not letting the bandpass filter ring out appropriately
  // in the bidirectional filtering process

  // f3: at f=100 Hz, it oscillates around the true value and is actually worse that the realtime
  // algorithm - there must be some bug - or maybe not?
  // Ahhh! i think, the two halfcycles have different lengths. maybe, we should only use the upward
  // zero crossings instead of both. The oscillation gets worse, when the pulsewave becomes more
  // asymmetric

  // Notes on the offline algorithm:
  // There are "edge-effects", i.e. the measured pitch is (sometimes) slightly off at the very
  // start and end of the signal. This could be traced to a slight mistuning of the bandpass filter
  // by which the signal is filtered first. It is tuned to a preliminary estimate of the
  // fundamental given by an autocorrelation algorithm taken in the middle of the signal. If this
  // preliminary estimate is a bit off, the final measurement will be off in the same direction at
  // the edges. This is because at the edges, we see more of the transient response of the bandpass
  // and less of the actual signal content near the bandpass'es center frequency. A possible way to
  // counteract would be to apply a second pass of the whole algorithm, but this time with a
  // refined bandpass center frequency taken as the average of the actual measured pitch (or some
  // middle portion thereof, in order to not take the poorly measured edges into account in the
  // averaging).

  // Further ideas:
  // In case of a missing (or weak) fundamental, we could apply a steep (elliptic?) bandpass to
  // isolate the 2nd and 3rd harmonic, square this signal (this creates the fundamental as
  // intermodulation product) and then use this algorithm on the resulting signal.


  delete[] x;
  delete[] f1;
  delete[] f2;
  delete[] f3;
  delete[] ft;
  delete[] y;
}

void zeroCrossingPitchDetectorTwoTones()
{
  // Creates 2 sinusoids with different frequencies and tries to detect the pitch of the sum of
  // those signals.

  static const int N = 22050;  // number of samples

  double fs = 44100.0;  // sample rate
  double f1 = 95;       // frequency of 1st sinusoid
  double f2 = 105;      // frequency of 2nd sinusoid
  double a1 = 1.0;      // amplitude of 1st sinusoid
  double a2 = 1.0;      // amplitude of 2nd sinusoid
  double p1 = 0.0;      // phase of 1st sinusoid
  double p2 = 0.0;      // phase of 2nd sinusoid

  double x[N], f[N], r[N]; // input signal, detected frequencies and reliabilities

  // create the signal:
  double w1 = 2*PI*f1/fs;
  double w2 = 2*PI*f2/fs;
  for(int n = 0; n < N; n++)
    x[n] = a1*sin(w1*n+p1) + a2*sin(w2*n+p2);


  // todo: refactor the function, such that we can use the zero-crossing based pitch detection
  // also on the unfiltered signal.

  rsInstantaneousFundamentalEstimatorD::measureInstantaneousFundamental(
    x, f, N, fs, 20.0, 5000.0, r);

  // plot:
  //plotData(N, 0.0, 1.0/fs, x);  // input signal
  plotData(N, 0.0, 1.0/fs, f);  // measured frequency
  //plotData(N, 0.0, 1.0/fs, r);  // reliability


  int dummy = 0;
}









/*
template<class T>
class rsPeakTrailer // or rsPeakTrailDragger or just rsTrailDragger
{

public:


  void setDecaySamples(T d, T r = 0.5) // or maybe default value for r should be 1/e
  {
    c = pow(r, T(1)/d);
    // d-th root of r -> after d successive multiplications, we reach r, if we init with 1
  }


  T getSample(T x)
  {
    y *= c;    // or maybe try y = (1-c)*x + c*y; 
    if(x > y)
      y = x;
    return y;
  }


  T getSample2(T x)  
  {
    T out = y;
    if(x > y)
      y = x;
    y *= c;
    return out;
  }
  // with artificial delay for trying to avoid having two overlapping trails exactly at the spikes

  // T getSample(T x, T dt) ...for usage with non-equidistant data

  void reset() { y = T(0); }

protected:

  T c = T(0);
  T y = T(0);

};
// this is kinda nice - each peak carries an exponentially decaying trail of influence - minor 
// peaks below this trail will be dismissed as irrelevant
// moved to class RAPT::rsPeakMasker
*/

template<class T>
void ropeway(const T* x, int N, T halfTime, T* y, T* w) // w is workspace
{
  RAPT::rsPeakMasker<T> pm;
  pm.setDecaySamples(halfTime, T(0.5));
  int n;
  for(n = 0;   n <  N;  n++) w[n] = pm.getSample(x[n]);  pm.reset();   // forward pass
  for(n = N-1; n >= 0;  n--) y[n] = pm.getSample(x[n]);                // backward pass
  // todo: use pm.applyForward(x, N); pm.applyBackward(x, N);
  //rsPlotArrays(N, x, w, y);

  for(n = 0;   n <  N;  n++) y[n] = rsMax(w[n], y[n]);               // maximum of both passes - sorta works
  //for(n = 0;   n <  N;  n++) y[n] = T(0.5) * (w[n]+y[n]);            // average - nope!
  //for(n = 0;   n <  N;  n++) y[n] = w[n]+y[n]-x[n];                  // nope
  //for(n = 0;   n <  N;  n++) y[n] = w[n]+y[n];
  //for(n = 0;   n <  N;  n++) y[n] = rsMax(x[n],T(0.5) * (w[n]+y[n]));
  //for(n = 0;   n <  N;  n++) y[n] = rsMax(x[n],w[n]+y[n]);

  // average doesn't work, max sort of works but is not smooth between spikes - try bidiretional
  // serial passes with averaging between forward-first and backward first
  // ...the non-smoothness is actually not a problem in peak-finding applications


  // bidirectional parallel - todo: maybe try a serial bidirectional application of the filter
  // ...but maybe average between forward-first and backward-first pass - due to the nonlinearity,
  // the results of both may not be the same...or are they? -> figure out!
}

// averaging does not work because we multiply by 0.5 - the sum would be more more appropriate, 
// except at positions of the spikes themselves
// what if the traildragger had a 1-sample delay to respond to the spikes and at the end, we take
// the maximum of sum and original signal?

// i'd really like to have an algorithm that produces proper catenary shaped trails, but averaging 
// forward/backward exponetial trails doesn't really work and taking the max of forward/backward 
// gives non-smooth junctions between the peaks - however, for a peak-picker, we do not care about
// this - so max of forward and backward pass is perhaps most reasonable

// this works only correctly, if x[n] >= 0 for all n - maybe we should fix that by using an offset
// before and after....


template<class T>
void ropeway2(const T* x, int N, T halfTime, T* y, T* w) // w is workspace
{
  //rsPeakTrailer<T> pt;
  RAPT::rsPeakTrailDragger<T> pt;
  pt.setDecaySamples(halfTime, T(0.5));
  int n;

  // forward first:
  for(n = 0;   n <  N;  n++) w[n] = pt.getSample(x[n]);   // forward pass
  for(n = N-1; n >= 0;  n--) w[n] = pt.getSample(w[n]);   // backward pass

  // backward first:
  pt.reset();
  for(n = N-1; n >= 0;  n--) y[n] = pt.getSample(x[n]);   // backward pass
  for(n = 0;   n <  N;  n++) y[n] = pt.getSample(y[n]);   // forward pass

  //rsPlotArrays(N, x, w, y);

  // average of forward-first and backward-first:

  for(n = 0;   n <  N;  n++) y[n] = T(0.5) * (w[n]+y[n]);
  //rsPlotArrays(N, x, y);
}
// result look the same like the other version

template<class T>
void ropeway3(const T* x, int N, T* y, int numPasses)
{
  rsArrayTools::copy(x, y, N);
  for(int i = 0; i < numPasses; i++) {
    rsArrayTools::movingAverage3pt(y, N, &y[0]);
    rsArrayTools::maxElementWise(x, &y[0], N, &y[0]);
  }
}

template<class T>
std::vector<T> ropeway(const std::vector<T>& x, T halfTime)
{
  int N = (int) x.size();
  std::vector<T> y(N), w(N);
  ropeway(&x[0], N, halfTime, &y[0], &w[0]);
  //ropeway2(&x[0], N, halfTime, &y[0], &w[0]);
  return y;
}

void ropewayAlgo()
{
  // Simple demo of the ropeway algo on some hand-created artificial data ...tbc...

  typedef std::vector<double> Vec;

  double halfTime = 5;

  Vec x(101);
  x[20] = 0.5;
  x[30] = 1.0;
  x[40] = 0.2; 
  x[50] = 0.2;  //x[50] = 1;
  x[60] = 0.6;
  x[70] = 1.0; 
  x[80] = 0.4;
  x[90] = 0.1;
  Vec y = ropeway(x, halfTime);

  rsPlotVectors(x, y);
}

bool testPeakPicker()  // move to unit tests
{
  bool result = true;

  typedef std::vector<double> VecD;
  typedef std::vector<int>    VecI;

  rsPeakPicker<double> pp;
  VecD x;
  VecI p;

  // check against one neighbor to each side (that's the default setting):
  x = { 1,2,3,4,3,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({3});
  x = { 1,2,4,4,3,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({2,3});
  x = { 1,2,4,4,4,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({2,3,4});
  x = { 4,2,4,4,4,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,2,3,4});
  x = { 4,2,4,4,4,2,4 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,2,3,4,6});
  x = { 4,4,4,4,4,4,4 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,1,2,3,4,5,6});

  // check against two neighbours to each side:
  pp.setNumNeighbors(2);
  x = { 1,2,3,4,3,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({3});
  x = { 1,2,4,4,3,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({2,3});
  x = { 1,2,4,4,4,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({2,3,4});
  x = { 4,2,4,4,4,2,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,2,3,4});
  x = { 4,2,4,4,4,2,4 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,2,3,4,6});
  x = { 4,4,4,4,4,4,4 }; p = pp.getPeakCandidates(x); result &= p == VecI({0,1,2,3,4,5,6});
  x = { 1,2,3,4,3,4,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({3,5});
  x = { 1,2,3,4,3,5,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({5});
  x = { 1,4,3,4,3,5,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({1,5});
  x = { 1,4,3,4,3,4,1 }; p = pp.getPeakCandidates(x); result &= p == VecI({1,3,5});

  // todo: test with an asymmetric setting like 1 left, two right neighbors


  // how should we handle plateaus? maybe for a plateau, consider the start of the plateau as peak?
  // ...or maybe the center....but then what about plateaus of even length? maybe the idices should
  // be real numbers, then we could just use some x.5 value for the peak in case of even-length 
  // plateaus

  // maybe use two different comparator for functions for left and right neighbours, like > for 
  // left and >= for right - with this, we could control how plateaus are treated: using
  // >,>= would pick the 1st value of the plateau as peak, >=,> would pick the last, >,> would pick
  // none and >=,>= would pick all - maybe if one wants the center of a plateau, we should 
  // preliminarily pick all and the post-process

  // how should we treat the boundaries/edge-cases? maybe in the same way as interior points but 
  // the loop over the left neighbours is shortened at the left boundary

  // I think, the best is to treat plateau-values *all* as peak values and to include values at the 
  // borders, when they are greater than their existing neighbors. This is what we want for 
  // amplitude envelopes and also spectral envelopes.

  // the peak-picker could also be useful for spectral envelope estimation...maybe is this case, 
  // the smoothing width whould be a function of frequency - like proportional to frequency?

  return result;
}




// 3-point moving average, handles edges by taking one-sided 2-point average there - to be applied 
// iteratively to figure out, how many iterations of smoothing a peak survives and use that as 
// measure for the peak's relevance
template<class T>
void rsSmooth(const T* x, int N, T* y)
{
  rsAssert(x != y, "algorithm not suitable for in-place processing");
  // maybe make an even stronger assertion that x and y do not overlap - implement a function
  // int rsArrayTools::getOverlap(T* x, T* y, int N) that returns the number of overlapping elements
  // and/or a function bool rsArrayTools::areOverlapping(T* x, T* y, int N)

  y[0]   = T(1/2.) * (x[0]   + x[1]);
  y[N-1] = T(1/2.) * (x[N-2] + x[N-1]);

  //// test - leave edge values untouched:
  //y[0]   = x[0];
  //y[N-1] = x[N-1];

  for(int n = 1; n < N-1; n++)
    y[n] = T(1/3.) * (x[n-1] + x[n] + x[n+1]);

  bool preserveMean = false;  // make parameter
  typedef RAPT::rsArrayTools AR;
  if(preserveMean)
    AR::add(y, AR::mean(x, N) - AR::mean(y, N), y, N);
}
// does this function preserve the mean? -> nope! -> implement a version that subtracts the 
// difference between the means before and after
// -what, if we just leave the edge-values untouched?
// todo: implement a version that may operate in-place
// -maybe make a class with options for different handling of the ends, optional mean-preservation,
//  etc.
// -maybe the mean is not the only thing that may be useful to preserve - what about energy, for 
//  example?



// maybe implement savitzky-golay filters:
// https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
// for equidistant data, use convolution
// https://en.wikipedia.org/wiki/Local_regression


// convenience function:
template<class T>
std::vector<T> rsSmooth(const std::vector<T>& x)
{
  std::vector<T> y(x.size());
  rsSmooth(&x[0], (int) x.size(), &y[0]);
  return y;
}

// elementwise maximum of two vectors - maybe move to RAPT:
template<class T>
std::vector<T> rsMax(const std::vector<T>& x, const std::vector<T>& y)
{
  rsAssert(x.size() == y.size());
  std::vector<T> z(x.size());
  for(size_t i = 0; i < x.size(); i++)
    z[i] = rsMax(x[i], y[i]);
  return z;
}

// quick and dirty (and inefficient) implementation of an idea to find the envelope:
std::vector<double> testEnvelope1(const std::vector<double>& x)
{
  std::vector<double> t1 = x, t2 = x;
  int maxIt = 1000;
  for(int i = 1; i <= maxIt; i++)
  {
    t2 = rsSmooth(t1);

    if(i % 2 == 0)
      rsPlotVectors(x, t1);

    // two variants - figure out, which is better - or if there's any difference at all:
    //t1 = rsMax(t1, t2);
    t1 = rsMax(x, t2);
  }
  return t1;
}
// this algorithm creates an envelope that looks like a ropeway connecting the peaks - more 
// iterations tend to use less peaks and increase the tension of the rope - is should call it 
// ropeway algorithm :-)
// -maybe it makes sense to run a few ropeway iterations before searching for peaks - maybe that
//  number of required iterations can be reduced by using a longer filter kernel - maybe even IIR
//  -maybe instead of doing several iteration, we should do just one but let the user choose the 
//   length of the smoothing filter kernel - maybe a Gaussian IIR filter would be suitable
//  ->try it with the rhodes sample envelope
// -and/or maybe we should run an attack/release envelope follower with zero attack, 
//  bidirectionally - in this case - does it make a difference, which direction we run first? if so,
//  it may make sense to use an average of forward-first and backward-first application to make the 
//  algorithm invariant with respect to mirroring the data (this invariance seems desirable)


// rename to peakPickerMasks
void peakPickerShadows(const std::vector<double>& x, double shadowWidth)
{
  // Plots the intermediate signals with the shadows

  // Create time axis:
  int N = (int) x.size();
  std::vector<double> t = rsRangeLinear(0.0, double(N-1), N);

  // Create and set up peak picker:
  rsPeakPicker<double> pp;
  pp.setShadowWidths(shadowWidth, shadowWidth);

  // Create shadowed data (for plotting):
  std::vector<double> y(N), yL(N), yR(N);
  rsArrayTools::shiftToMakeMinimumZero(&x[0], N, &y[0]);
  pp.maskLeft( &t[0], &y[0], &yL[0], N);
  pp.maskRight(&t[0], &y[0], &yR[0], N);

  // Find relevant peaks:
  std::vector<int> peaksR;
  peaksR = pp.getRelevantPeaks(&t[0], &x[0], N);

  // Plot results:
  GNUPlotter plt;
  plt.setPixelSize(1200, 400);
  plt.addDataArrays(N, &t[0], &y[0]);                     // shifted data
  plt.addDataArrays(N, &t[0], &yL[0], &yR[0]);            // shadowed data
  addDataPartially(plt, t, y, peaksR);                    // relevant peaks
  plt.plot();
}

void peakPickerShadowSettings(const std::vector<double>& x, const std::vector<double>& widths)
{
  // plots envelopes resulting from various settings for the shadow width

  using VecD = std::vector<double>;

  int N = (int) x.size();
  VecD t = rsRangeLinear(0.0, double(N-1), N);

  VecD y(N);
  rsArrayTools::shiftToMakeMinimumZero(&x[0], N, &y[0]);

  rsPeakPicker<double> pp;
  std::vector<int> peaks;

  GNUPlotter plt;
  plt.setPixelSize(1200, 400);
  plt.addDataArrays(N, &t[0], &y[0]);                     // shifted data
  peaks = pp.getFinePeaks(  &t[0], &x[0], N); 
  addDataPartially(plt, t, y, peaks);
  for(size_t i = 0; i < widths.size(); i++) {
    pp.setShadowWidths(widths[i], widths[i]);
    peaks = pp.getRelevantPeaks(&t[0], &x[0], N);
    addDataPartially(plt, t, y, peaks); }
  peaks = pp.getCoarsePeaks(&t[0], &x[0], N); 
  addDataPartially(plt, t, y, peaks);
  plt.plot();
}

void showPeakPickerPlots(const std::vector<double>& x)
{
  peakPickerShadows(x, 20);
  peakPickerShadowSettings(x, { 10, 20, 30 });


  // -at n=26 (seed=65), the fine algo misses a stickout, seed 37: 149,150
  // -with seed 37, at n=194, width = 20.. a few stickouts are missed with the relevance algo
  //  -> we need to run the no-stickout algo on the fully data in any case
  //  -> this probably makes it more meaningful to always include the edge-values
  // -how would this compare to a regular bi-directional envelope follower with zero-attack?

  // todo: experiment with using the no-stickout criterion *only* starting with 0 and N-1 as the
  // initial array of peaks - this would give the coarsest possible peak-array that sastisfies the
  // no-stickout criterion
  // -maybe we could provide a higher-level "sensitivity" parameter 
  //  -if it's zero, we only get the peaks according the stickout criterion - the coarsest possible
  //   envelope
  //  -if it's 1 (or 100%), we get all peaks - the finest possible envelope
  //  -maybe this should somehow crossfade the shadowWidths from the minimum values that would 
  //   give the same result as the no-stickout criterion only and the maximum values that would 
  //   give all peaks (but how do we compute these values? this seem to be nontrivial)
  // -maybe connect the peaks with splines instead of lines
  // -maybe plot prominence related stuff - maybe plot the prominences themselves as
  //  stems

  // Notes:
  // -Can it happen that there are 3 peaks p1 > p2 > p3 at time instants t1 < t2 < t3 such that the 
  //  shadowing algorithm misses the intermediate peak p2 but p2 still sticks out from the line 
  //  that connects p1 and p3? I don't think so, because if the algo catches p1 and p3, p3 will lie
  //  below the exponetial shadow/trail emanating from p1 and that shadow will lie wholly under the
  //  line connecting p1 and p3 - so p2 must lie under the shadow (because it was missed) and 
  //  therefore even more so must lie under the line. This is because the expontial shadows are 
  //  always convex - the function segment between any two points always lies *under* the line 
  //  that directly connects the points
}

template<class T>
std::vector<T> peakSmoothabilities(
  const std::vector<T>& data, const std::vector<int>& peakIndices)
{
  std::vector<T> s(data.size());

  // should compute the number of smoothing iterations, each peak survives


  return s;
}
// maybe rename to peakScaleRanges

void peakPicker()
{
  bool peakPickerWorks = testPeakPicker();  // unit test

  typedef std::vector<double> VecD;
  typedef std::vector<int>    VecI;

  rsPeakPicker<double> picker;
  VecD x;
  VecI p;

  // test the peak-prominence algorithm with examples given here:
  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html
  // it would be totally cool, if we could just copy-and-paste the python-code form there and run
  // it in C++ - we need classes rsNumPy, rsSciPy, rsPyPlot


  rsNumPy<double> np;  // this replaces pythons "import numpy as np"
  x = np.linspace(0, 6 * np.pi, 1000);
  x = np.sin(x) + 0.6 * np.sin(2.6 * x);
  p = picker.getPeakCandidates(x);  // rename to getPeakCandidates or factor out such a function
  VecD proms = picker.peakProminences(x, p);

  // python has 8 peaks with prominences:
  // 1.24159486, 0.47840168, 0.28470524, 3.10716793, 0.284603, 0.47822491, 2.48340261, 0.47822491
  // we have 9 peaks and the last prominence is 0 - why? ...because a peak at the edge will always 
  // get zero prominence with the current implementation - it would be desirable to have a 
  // meaningful prominence for edge-peaks, too - not true anymore - that was the case with the
  // standard algo - we now have a variation

  // figure out, how we should handle the edge cases, when one or both of the loops hit the 
  // data-boundary - what would be the most desirable result then?

  //rsPlotVector(x);

  // plot the peak-marks and the prominences also...

  // try monotonic arrays, an unimodal distribution, bimodal, trimodal, try peaks near and at 
  // boundaries

  // {5,4,3,2,1} should have a peak at 0 with prominence 4, {1,2,3,4,5} should have a peak at 4
  // with prominence 4, {1,2,3,4,3,2} should have a peak at 3 with prominence 3, ....


  using AT = RAPT::rsArrayTools;

  int N = 200;

  x = rsRandomVector(N, -1.0, +1.0, 7);  // nice seeds: 4,6,7,8,37,41,65 - strange: 5,9
  AT::cumulativeSum(&x[0], &x[0], N);    // one pass seems good to create "landscapey" data
  //AT::cumulativeSum(&x[0], &x[0], N);  // two passes is too much
  //AT::add(&x[0], AT::minValue(&x[0], N), &x[0], N); // values hould be >= 0 for the shadower to work


  // todo: 
  // -test with randomized t-data, too - i.e. with non-equidistant data, because that's 
  //  actually the use-case in finding envelope peaks
  // -maybe try filters with different slopes instead of the cumulative sum - maybe slopes between
  //  0 and 9 dB/oct are reasonable (cumulative sum has 6 dB/oct) - with 12, the smooting is so 
  //  strong that we don't really see much
  // -and/or maybe use a 1-pole lowpass with some cutoff frequency - yes, cutoff and slope seem 
  //  resonable parameters for the data creation (maybe even to create mountain terrain for 
  //  graphics)?


  //x = VecD({0,0,0,8,4,4,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0});
  //x = VecD({0,0,0,10,9,9,8,8,8,7,7,7,7,6,6,6,6,6,0,0,5});  // illustrates peak wandering/drift
  //x = VecD(20); for(int i=4; i<20; i++) x[i] = 1./i;
  //x = VecD(20); for(int i=4; i<20; i++) x[i] = exp(-0.25*i);
  //x = VecD(20); x[8] = 10; x[12] = 20;  // shows split (1st iteration) and merge (2nd and 3rd)
  //x = VecD(20); x[8] = 10; x[12] = 10; x[13] = 15;
  //x = VecD(31); x[15] = 10;


  // create datasets that illustrate peak-merging and peak-splitting

  // this shows plots of the various iterations of the "ropeway" algorithm - todo: show many of them in
  // a single plot:
  //VecD yEnv1 = testEnvelope1(x);  
  // uses iterative smoothing - inefficient!

  // shows result using recursive smoothing via rsPeakTrailDragger - this should be more efficient
  // and not suffer from peak-wandering/drift
  showPeakPickerPlots(x);  // maybe this should also take an (optional) array of time-stamps
  //rsPlotVectors(x, yEnv2);



  // ...actually, using an attack-decay envelope follower, we'll see exponential decays and if it's
  // used bidirectionally (serially and/or in parallel), we may indeed get something like a cosh
  // shape - with this algo, it looks more like parabolas


  /*
  // visualize effect of iterative smoothing
  int N = (int) x.size();
  VecD y1(N), y2(N), y3(N), y4(N), y5(N), y6(N), y7(N), y8(N);
  rsSmooth(&x[0],  N, &y1[0]);
  rsSmooth(&y1[0], N, &y2[0]);
  rsSmooth(&y2[0], N, &y3[0]);
  rsSmooth(&y3[0], N, &y4[0]);
  rsSmooth(&y4[0], N, &y5[0]);
  rsSmooth(&y5[0], N, &y6[0]);
  rsSmooth(&y6[0], N, &y7[0]);
  rsSmooth(&y7[0], N, &y8[0]);
  //rsPlotVectors(x, y1, y2, y3, y4, y5);
  rsPlotVectors(x, y1, y2, y4, y8);
  //rsPlotVectors(x, y8);
  int dummy = 0;
  */
  // maybe mark the peaks and print their prominences and smoothing-robustnesses above....what 
  // could be a good name for this insensitivity-under-smoothing porperty? maybe something like 
  // robustness/resilience? keeping with analogy of landscapes, we may consider these smoothing 
  // operations as a sort of abrasion due to wind or water - we may think of the measure as a 
  // sort of abrasion-resilience...maybe durability? smoothability? whetherproofness, 
  // wheather-resistance, or just resistance...i think, perhaps durability is the best term
  // or ScaleRange?
  // maybe do not use number of iterations of MA smoothing, but decay-time for 
  // rsPeakTrailDragger - maybe compute such a value with taking the full input into account and
  // with taking only the peaks themselves into account, i.e. (0,0,4,4,4,5,4,4,4,0,0) vs
  // (0,0,0,0,0,5,0,0,0,0,0) where the latter can be implemented using the non-uniform getSample
  // function only feeding the peaks (not the zeros in between) - gives two different measures,
  // one which considers only peak-height and one which takes peak-width also into accout


  // tests the class rsPeakPicker
  // todo: 
  // -create an array of random numbers
  // -maybe smooth it (with adjustable settings)
  // -find peaks using rsPeakPicker with different settings
  // -plot the array with marks at the found peaks
  // -check, if plot looks as expected
  // -do that with various seeds - if for some settings something strange happens, note dow the
  //  settings



  // for amp-envelopes, it may make sense to use as additional criterion (in addition to be greater
  // or equal to its neighbors), that the smoothed versions of the envelope have a peak there, too
  // (maybe not exactly at the same index but very closely nearby)
  // here:
  // https://books.google.de/books?id=8h-ZDgAAQBAJ&pg=PT788&lpg=PT788&dq=finding+relevant+peaks&source=bl&ots=ZmIkfb9Zl6&sig=ACfU3U1vLoeM7mPahOTov618ijcnCDTh1w&hl=en&sa=X&ved=2ahUKEwja5Myu2bPmAhViQUEAHZjiAqQQ6AEwB3oECAkQAQ#v=onepage&q=finding%20relevant%20peaks&f=false
  // it says something about real peaks should be visible on multiple scales



  // for ideas, see:
  // https://stackoverflow.com/questions/5672095/finding-relevant-peaks-in-messy-ffts
  // https://se.mathworks.com/help/signal/ref/findpeaks.html#bufbbs1-MinPeakProminence
  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html
  // http://matlab.izmiran.ru/help/toolbox/images/morph13.html


  // thread: "I'd like to add a peak finding function to scipy.signal"
  // https://github.com/scipy/scipy/issues/8211


  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html

  // https://blog.ytotech.com/2015/11/01/findpeaks-in-python/

  // maybe we could do some post-processing on a preliminary peak-array by which we discard some 
  // of the preliminary peaks - for example based on minimum desired distance between peaks
  // or - based on peak-prominence:
  // https://en.wikipedia.org/wiki/Topographic_prominence
  // https://listsofjohn.com/PeakStats/glossary.html#Rise
  // maybe write a function getPeakProminences that takes an array of values and an array of peak
  // indices and return an array of peak prominences
  // the wikipedia articla also links to this:
  // https://en.wikipedia.org/wiki/Morse_theory

  // maybe the shape of the signal around the peak should be generally curved downward...so maybe 
  // we should look at the 2nd derivative...maybe of a smoothed signal?

  // in the thread on stackoverflow, somehwat suggest to smooth the data with this:
  // https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter

  // Ideas for peak-picking:
  // -compute peak prominences:
  //  -consider a peak in a landscape
  //  -loop over all possible directions (north, north-east, east, ...)
  //   -for chosen direction, walk until you find an elevation that is higher (or equal?) than the 
  //    elevation of the peak under consideration
  //   -along that line, note down the minimum elevation
  //  -among those noted minimum elevations, take the maximum
  //  -the difference between the peak-height and the maximum-minimum-elevation is the 
  //   peak-prominence
  //  -if the data is one-dimensional, there are only two directions to scan: left and right, in 2D
  //   there are actually uncountably infinitely many - but we don't have to deal with that case
  // -maybe define the relative peak-prominence as the peak-prominence divided by the peak-height
  //  ...maybe that's a more useful measure for peak relevance
  // -and/or try peak-prominence divided by height of global maximum
  // 
  // -another idea is to consider scale-invariance of a peak:
  //  -a relevant peak typically remains a peak even when the data is smoothed
  //  -maybe progressively smooth the data and note, how many such smoothing passes the peak 
  //   "survives" and take that as an indicator of peak relevance
  //  -during smoothing, the exact location/index may change a little - maybe allow the peak index 
  //   to be off by one with respect to the previous smoothing iteratiom in order to consider the 
  //   peak having survived the smoothing iteration
  //  -maybe use y[n] = (1/3.) * (x[n-1] + x[n] + x[n+1]) in each smoothing pass (at the endpoints,
  //   use y[0] = (1/2.) * (x[0] + x[1]), y[N-1] = (1/2.) * (x[N-1] + x[N-2])
  //   -implement that smoothing filter, operating in-place
  //  -maybe that averaging/smoothing can be combined with decimation in each smoothing stage - to
  //   save cpu-time and maybe it may even give better results?
  //  -make sure that after each smoothing pass, the mean is the same as before - check that and it 
  //   it's not the case, compute the mean before and after and add the difference to the smoothing
  //   output to reconstruct the old mean
  //  -during the algorithm, we need to track, how the peak's index wanders from one smoothing 
  //   iteration to the next (it may wander one index left or right in each iteration)...so maybe 
  //   for each peak, we need to record its path...maybe this could also be used as a relevance 
  //   measure? relevant peaks should not wander very much? ...but maybe that doesn't make sense
  //   -also, one peak may "split" into two - if we have a peak at index n in iteration i, we may
  //    find two peaks at n-1,n+1 in iteration i+1 - which one would then be considered the 
  //    offspring? maybe the higher one, but what if both have the same height? ..maybe we need to 
  //    track both offsprings?
  //   -or maybe in each iteration i we should allow an index deviation of i positions with respect
  //    to the peak's original index ...or maybe we just need sqrt(i) - try to figure out, by how 
  //    much the peak may drift after i iterations in the worst case - but what is the worst case?
  //    it should probably be a one-sided distribution - but what exact shape? i guess that depends
  //    on the filter kernel(s) used as well
  //    -i think, the scope/width in which we search for a corresponding peak should be 
  //     proportional to the width of the (convolutively) accumulated filter kernel, so sqrt
  //     seems reasonable - maybe to incrase the width by one, we should look at results after
  //     1,2,4,9,16,25,36,49...N^2 iterations
  //  -note that during the process, initially distinct peaks may merge into a single peak (due to 
  //   tolerating peak wandering) and both peaks would be considered to have "survived"
  //   -i think, that in each iteration, each peak may have up to two parents and up to two 
  //    children
  //
  // -maybe combine both ideas: track how the peak-prominences change during smoothing - a peak 
  //  whose prominence does not decrease much during smoothing is more relevant
  //
  // -as a totally different idea: in the signals considered here, peaks often occur in patterns, 
  //  for example, at regular intervals (tremolo in amp-envs) - a human observer would probably 
  //  (subconsciously) also apply some sort of pattern matching to identify the relevant peaks
  //  so maybe more sophisticated peak-picking algos should also include pattern-matching - to be
  //  used only in cases when there is a pattern...which may itself be determined by 
  //  autocorrelation? ...however, these are just some wild ideas that are overkill here - i'm
  //  just trying to contemplate, what humans probably would do subconsciously
  //
  // -additional idea: if, after finding the peaks, it happens that a dismissed peak stands out 
  //  above the linear interpolant of the kept peaks, re-include it - the envelope through the 
  //  peaks should not allow any original peak to stick out of the envelope
  //  -such a post-processing seems especially important, if some sort of minimum-peak distance
  //   is used 

  // -other algo: start with highest peak, successively add the next highest peak, iff it is 
  //  outside the min-peak-distance - wehn no other such peaks can be found, add those peaks that
  //  violate the no-stick-out condition

  // -this: http://billauer.co.il/peakdet.html searches alternatingly for minima and maxima and 
  //  accepts them only when they are some threshold ("delta") below/above the previous one
  //  ...is this similar to prominence thresholding?

  // -evaluate the implemented algorithms in terms of false positives (spurious peaks) and false 
  //  negatives (missed peaks) ...maybe use artifical spectral data of sinusoids embedded in noise
  //  and try to find the sinusoidal peak in the noise?

  // -for the amp-env problem in the debeater, i tend to think that we should consider peaks 
  //  relevant, if they survive 2 rounds of smoothing...or maybe 3? this should help against 
  //  spurios local peaks within minima....but for the little peaks on the tail, probably a 
  //  prominence based threshold is better? -> use both

  // -relationships between peak-picking and envelope-estimation:
  //  -a sufficient (but not necessary) condition for a peak to be relevant is to be also a peak
  //   of the envelope
  //  -an envelope can be estimated by connecting the found peaks

  // -i think, the numNeighbours stuf is actually the same as the minDistance parameter found in
  //  some peak-picking algos - minDistance = numNeighbours+1?

  // -for envelope estimation: how about finding also the most relevant valleys and place 
  //  datapoints the, too (the rsPeakFinder class can be used by just sign-inverting the input)

  // -maye have a function erodePlateaus that keeps only the left and right endpoint of a plateau
  //  in a peak-array - data economization...maybe economizePlateaus would be a better name

  // Strange observation: look at random walks with seeds 6 and 9 - they have a very different 
  // quality - with 9, there are far less medium-scale wiggles

}

void singleSineModelForSineSweep()
{
  // We test the algorithm to estimate the instantaneous amplitude, phase and frequency of a pure 
  // sinusoid using rsSingleSineModeler.

  // User parameters:
  int    N  =   2000;      // number of samples
  double fs =  44100;      // sample rate
  double f0 =    500;      // frequency at start
  double f1 =   5000;      // frequency at end
  double a0 =  0.25;       // amplitude at start
  double a1 =  0.50;       // amplitude at end
  double freqShape = 0.0;  // shape of freq-env, 0: exponential, 1: linear, in general: power rule
  double ampShape  = 1.0;  // dito for amp env

  // Generate arrays for (true) instantaneous frequency, omega, phase and amplitude:
  using Vec = std::vector<double>;
  using AT  = rsArrayTools;
  Vec f(N), w(N), a(N), p(N);
  AT::fillWithRange(&f[0], N, f0, f1, freqShape);
  AT::fillWithRange(&a[0], N, a0, a1,  ampShape);
  AT::scale(&f[0], &w[0], N, 2*PI/fs);
  AT::cumulativeSum(&w[0], &p[0], N);  // maybe try different integration scheme (trapezoidal, etc.)
  //rsPlotVector(f);
  //rsPlotVector(a);
  //rsPlotVector(p);

  // Create signal and quadrature component:
  Vec x(N), q(N);
  for(int n = 0; n < N; n++)
  {
    x[n] =  a[n] * sin(p[n]);
    q[n] = -a[n] * cos(p[n]);
  }
  //rsPlotVectors(x, q);

  // Create and set up the analyzer object:
  rsSingleSineModeler<double> ssm;

  // Now, we do one experiment at a time to try to retrieve the instantaneous amplitudes and 
  // phases from the signal using different algorithms. We plot the results directly and the 
  // observations are written down immediately below each experiment. Note that in many cases, its
  // actually enough to look at the amplitude estimation error. The phase estimation error will 
  // then be completely determined by that due to the identity resynthesis property of the 
  // analysis/resythesis scheme.

  // Use the offline processing function with default settings (ToDo: explicitly set up the 
  // settings before the experiment because we may change the defaults later):
  Vec aE(N), pE(N);  // E stands for "estimated"
  ssm.analyzeAmpAndPhase(&x[0], N, &aE[0], &pE[0]);
  //rsPlotVectors(x, a, aE);
  rsPlotVectors(aE-a);  // amp estimation error
  //rsPlotVectors(pE-p);  // phase estimation error
  // Observations:
  // -The amplitude errors have a transient phase of large error of around 20 samples, then it is
  //  close to zero for the lower frequencies and it becomes a bit more "noisy" for higher 
  //  frequencies, The last sample of the amp error is quite high - this looks like an edge 
  //  artifact (-> verify).
  // -The phase errors grows steadily (but in steps). I think, this is due to (un)wrapping which 
  //  is not taken into account by the analyzer, i.e. the analyzer produces a wrapped phase 
  //  whereas the original phase data is unwrapped.  
  // ToDo:
  // -Figure out the reason for the transient in the error


  // Use phaseAndAmpFormulaForward as in realtime processing (the last sample uses the backward 
  // formula):
  for(int n = 0; n < N-1; n++)
  {
    ssm.phaseAndAmpFormulaForward(x[n], x[n+1], w[n], &aE[n], &pE[n]);
    // todo: 
    // -verify, if we should use w[n] - or should it be w[n+1] or (w[n] + w[n+1])/2?
    // -try backward formula
    // -try to estimate it without the known instantaneous freq w
  }
  ssm.phaseAndAmpFormulaBackward(x[N-2], x[N-1], w[N-1], &aE[N-1], &pE[N-1]);
  //rsPlotVectors(x, a, aE);
  //rsPlotVectors(aE-a);
  rsPlotVectors(0.002*x, aE-a); // plot (attenuated) input signal along for reference
  //rsPlotVectors(pE-p);
  // Observations:
  // -This does not show the transient amplitude error but the last sample is also wrong due to
  //  edge effets (which is normal and OK).
  // -The amplitude error wiggles in a sinusoidal way with a frequency that is twice the input 
  //  sinusoid's frequency. The amplitude of the wiggle is higher (around 0.001) at the beginning 
  //  (at 500Hz) than at the end (at 5kHz - around 0.0006). This seems to have nothing to do with 
  //  the signal frequency but rather with its amplitude or the interaction of both. For a raising 
  //  freq env and raising amp env, the error decreases with a shape that looks exponentialish.
  // ToDo:
  // -Use linear shapes for both envelopes and invetsigate how that affects the shape of the error
  // -Try different combinations of rise and fall for freq and amp 

  // Estimate instantaneous amplitude using an allpass to create a quadrature component:
  rsOnePoleFilter<double, double> apf;
  apf.setSampleRate(fs);
  apf.setMode(apf.ALLPASS_BLT);
  Vec qE(N);  // estimated quadrature component
  for(int n = 0; n < N-1; n++)
  {
    apf.setCutoff(f[n]);
    qE[n] = apf.getSample(x[n]);
    aE[n] = sqrt(x[n]*x[n] + qE[n]*qE[n]);
  }
  //rsPlotVectors(x, q, qE);
  rsPlotVectors(0.002*x, aE-a);
  // Observations:
  // -After a short transient phase, the estimation of the quadrature component is really close to
  //  the desired exact quadrature signal.
  // -It seems like the estimation gets even better towards the higher frequencies at the end. Maybe 
  //  for higher freqs, the allpass state can adapt faster due to smaller feedback coeff?
  // -It seems to work well all the way up to fs/2. The same cannot be said about the two 
  //  algorithms above. So, for high, freqs, the allpass approach is really well suited. Maybe the 
  //  algos above are better suited for lower freqs where the allpass has too much inertia?
  // ToDo:
  // -Try even higher freqs.
  // -Try to time-advance/delay setting the allpass freq.
  // -Try different topologies for the allpass (DF1, DF2, TDF1, TDF2) and compare the errors.
  // -For comparing the topologies, use a more aggressive modulation of freq and/or amp, maybe even
  //  containing discontinuous jumps. This should expose the differences best.
  // -Can we somehow fix the allpass'es (supposed) inertia problems by using adifferent topology 
  //  and/or updating the allpass state on a cutoff change?
  // -Figure ou why it does not work so well for high frequencies in the context of the ResoWave 
  //  filter...maybe there's a bug in the code there? maybe a delay by one sample in qudrature 
  //  component? ...try it with a reso-split filter here (i.e. a filter that wraps the splitting
  //  off of the resonance for the ladder)


  // Now we try estimating the instantaneous phase and amplitude of the resonance of the 
  // rsResoPlitFilter. This is the important step in the rsResoWaveFilter.
  rsResoSplitFilter resoSplitter;
  // maybe move into its own experiment...maybe singleSineModelForReso. resoAnalyzer
  // ...ot maybe try such an estimation with a resonator filter first because there, we can 
  // actually compute the correct traget output (at least when we use a phasor-based 
  // implementation like rsModalFilterNonlinear) 
  // -Feed an impulse-train into a phasor-based resonator filter
  // -Compute the exact instantaneous amplitudes and phases from sine and cosine component in
  //  the filter.
  // -Try to estimate them using various etchniques, among them, the allpass approach
  // -Compare estimate to target values


  int dummy = 0;


  // ToDo: 
  // -Maybe wrap the original synthesis phase. This may also be numerically better.
  // -Instead of calling the offline analysis function analyzeAmpAndPhase, use a realtime 
  //  estimation
  // -Compare it to an approach using an allpass filter to create a quadrature component to 
  //  estimate the amplitude
  // -Try to also estimate instantaneous frequency, etc. use different algorithms
  // -Try it with some discontinuous switches in the signal's frequency and/or amplitude
  // -Use even higher frequencies, all the way up to fs/2.


  // Maybe use a continuous "exponentiality" parameter that should somehow blend/morph between
  // linear an exponential. Maybe use (in continuous terms, t = 0..1):
  //   yL(t) = y0 + t*(y1-y0)           // linear
  //   yE(t) = y0 * exp(t*log(y1/y0));  // exponential
  // combine these formuals as follows, using "a" as exponentiality parameter:
  //   y(t) = y0 * exp(a*t*log(y1/y0)) + (1-a)*t*(y1-y0)
  // ...verify if this gives a reasonbale morph between linear and exponential. At the moment, we
  // use the power law

}


std::vector<double> estimateQuadratureViaAllpass(const std::vector<double>& x, 
  const std::vector<double>& f, double fs)
{
  int N = (int)x.size();
  rsAssert((int)f.size() == N);
  std::vector<double> q(N);
  rsOnePoleFilter<double, double> apf;
  apf.setSampleRate(fs);
  apf.setMode(apf.ALLPASS_BLT);
  for(int n = 0; n < N-1; n++)
  {
    apf.setCutoff(f[n]);
    q[n] = apf.getSample(x[n]);
  }
  //rsPlotVectors(x, q);
  return q;
  // ToDo: 
  // -have a parameter for selecting the topology
}
std::vector<double> estimateInstAmpViaAllpass(const std::vector<double>& x, 
  const std::vector<double>& f, double fs)
{
  std::vector<double> q = estimateQuadratureViaAllpass(x, f, fs);
  std::vector<double> a(q.size());
  for(size_t n = 0; n < x.size(); n++)
    a[n] = sqrt(x[n]*x[n] + q[n]*q[n]);
  //rsPlotVectors(x, q, a);
  return a;
}

void singleSineModelForResoSweep()
{
  // Like the function above just that now we use a sweep of a resonator filter applied to an
  // impulse train instead of a sweep of a perfect sine oscillator. So, in this experiment, the 
  // input signal to the estimator deviates already a bit more from the assumptions that the 
  // analysis algo makes about the signal. We want to see, how much the error of the estimates 
  // increases due to that deviation. The desired target values are computed by making use of
  // a phasor based implementation of a resonator in which exact sine and cosine components are
  // readily available from which magnitude and phase can be calculated exactly.

  // User parameters:
  using Real    = double;
  using Complex = std::complex<Real>; 
  int  N   =   2000;      // number of samples
  Real fs  =  44100;      // sample rate
  Real f0  =  10000;       // frequency at start
  Real f1  =   1000;       // frequency at end
  Real fIn = fs/500;      // freq of input impulse train
  Real dec =  0.009;      // resonance decay in seconds
  Real phs =  0.0;        // resonance phase (in degrees, i think)


  // Create input impulse train:
  using Vec = std::vector<Real>;
  using AT  = rsArrayTools;
  Vec x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0.0, true);
  AT::difference(&x[0], N);
  //rsPlotVector(x);
  AT::fillWithImpulse(&x[0], N);  // test

  // Create resonator filter, create filter cutoff array, filtered signal and arrays of true 
  // instantaneous amp and phase:
  rsNonlinearModalFilter<Real, Real> mf;
  Vec f(N), a(N), p(N), y(N), q(N);
  AT::fillWithRangeExponential(&f[0], N, f0, f1);
  //rsPlotVector(f);
  for(int n = 0; n < N; n++)
  {
    mf.setModalParameters(f[n], 1.0, dec, phs, fs);
    Complex z = mf.getComplexSample(Complex(x[n], 0.0));

    y[n] =  z.imag();    // we use the sine component as output
    q[n] = -z.real();    // this is the true quadrature component
    a[n] =  abs(z);      // true instantaneous amplitude
    // todo: true inst phase

    //y[n] = mf.getSample(x[n]);
    //a[n] = mf.getInstantaneousEnvelope(x[n]);
  }
  rsPlotVectors(x, y, q, a);

  // Try to generate an estimated quadrature signal qE using y and f as inputs by allpassing y with
  // an allpass adjusted according to f:
  Vec qE = estimateQuadratureViaAllpass(y, f, fs);
  Vec aE = estimateInstAmpViaAllpass(   y, f, fs);
  rsPlotVectors(y, q, qE, a, aE);
  rsPlotVectors(0.002*y, aE-a);

  // Try to estimate the instantaneous amp via rsSingleSineModeler:
  // ...



  // Now try it with a ladder resonance:
  rsLadderFilter<double, double> ldrR, ldrN; // resonant and nonresonant ladder
  ldrR.setSampleRate(fs);
  ldrN.setSampleRate(fs);
  ldrR.setResonance(0.99);
  ldrN.setResonance(0.0);
  for(int n = 0; n < N; n++)
  {
    ldrR.setCutoff(f[n]);
    ldrN.setCutoff(f[n]);
    double yr = ldrR.getSample(x[n]);  // resonant filter output
    double yn = ldrN.getSample(x[n]);  // nonresonant filter output
    double r  = yr - yn;               // pure resonance
    y[n] = r;
  }
  qE = estimateQuadratureViaAllpass(y, f, fs);
  aE = estimateInstAmpViaAllpass(   y, f, fs);
  rsPlotVectors(y, qE, aE);




  //rosic::writeToMonoWaveFile("ResoSweep.wav", &y[0], N, fs);


  int dummy = 0;

  // Observations:
  // -When using the anti-aliased impulse-train, we get very strange results even for the exact 
  //  instantaneous amplitude - it looks wiggly. Using a single unit impulse as input does not have
  //  this problem.
  // -The estimation of instantaneous amplitude via the allpass technique seems to work well after
  //  some transient phase. ..i think, this phase might be longer for lower frequencies due to 
  //  filter inertia?
  // -The estimated allpass component seems a little bit louder than the exact one for upward 
  //  sweeps and downward sweeps alike. The error is larger at the start. Maybe this overestimation
  //  has to do with the fact, that the signal decays? ..yes - the error is a bit larger for 
  //  shorter decay times. The dependency is not very strong but it is there.
  // -The estimated amplitude for the ladder resonance looks good as well - also for high 
  //  frequencies.

  // ToDo:
  // -Maybe use a naive impulse-train as input - the anti-aliasing may make things more difficult 
  //  to analyze
  // -Write a function generateImpulseTrain(double freq, bool antiAlias) and use it here to 
  //  generate the input signal
  // -Try the other algorithms -> factor out a function that estimates instantaneous amplitude, 
  //  taking as input y and f. have different variants of that function, using different algos: 
  //  allpass, the algos from rsSingleSineModeler, etc. ...maybe include the allpass method into
  //  rsSingleSineModeler as well
  // -Try other allpass topologies
  // -Maybe try updating the allpass state somehow when changing the alppas freq to reflect the new
  //  freq. I don't know, if that makes sense or is possible -> figure out...
  // -Try the same experiments with a ladder resonance - figure out, why the resoReplaceScream 
  //  fails at high frequencies. I thought that the instantaneous amplitude and phase are 
  //  misdetected but that seems unlikely give that the detection via allpass still works nicely 
  //  for high frequencies...unless there's a bug in the resoReplacer which leads to misdetection.
  //  ...done - result looks good actually - there must be something elso to blame. Maybe it has to
  //  do with the nonlinearity? Try driving the filter harder....

  // Sidenotes:
  // -a reso-sweep from 1kHz to 4kHz with dec = 0.01 in 2000 samples sounds like a water drop
  //  ...i just foudn this by accident
}

void singleSineModel()
{
  //singleSineModelForSineSweep();
  singleSineModelForResoSweep();
  //singleSineModelForLadderSweep();
}