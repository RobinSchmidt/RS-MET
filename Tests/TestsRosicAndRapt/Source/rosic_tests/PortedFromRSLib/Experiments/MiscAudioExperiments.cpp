#include "MiscAudioExperiments.h"

void centroid()
{
  double x[10] ={ 0,0,1,2,2,1,0,0,0,0 };

  double c; 
  c = rsCentroid(x, 10);         // should be 3.5
  c = rsCentroidOfEnergy(x, 10); // should still be 3.5

  int dummy = 0;
}

// functions are now part of the library
double cubicFadeIn(double x)
{
  return x*(x*((PI/2-2)*x+(3-PI))+(PI/2));
}
double cubicFadeOut(double x)
{
  return x*x*((2-PI/2)*x+(PI/2-3))+1;
}
void cubicCrossfade()
{
  // code needs update to use GNUPlotter:
  //Plotter p;
  ////p.plotTwoFunctions(rsCubicFadeIn, rsCubicFadeOut, -0.5, 1.5);
}

int rsCeilDiv(int n, int m)  
{
  return (int) ceil(double(n) / double(m)); 
  // stupid version - maybe use something like // (n + n % m) / n
  // or 
  // k = n/m
  // if(n > k*m) k+=1;
}


void decimate()
{
  int N  = 50;   // number of samples in example signal
  int D  = 4;    // decimation factor
  int ND = N/D;
  //int ND = rsCeilDiv(N,D);
  // todo: floor(N/D) may not be the best choice - it may throw away data at the end - how about
  // ceil? ...maybe make a funtion rsCeilDiv(n, m) as opposed to the regular floor-division aka
  // integer division
  // for N=20, D=3 -the decimated signal should be at leats one sample longer, i think - maybe even
  // two - maybe we should use a value that does not discard any samples

  std::vector<double> t(N), td(N/D), x(N), xd(ND), xa(ND);

  typedef RAPT::rsArrayTools AR;

  // time axis (original and decimated):
  AR::fillWithIndex(&t[0], N);
  AR::decimate(&t[0], N, &td[0], D);   // change order: x, N, y, D

  // signal (original and decimated):
  AR::fillWithRandomValues(&x[0], N, -1.0, 1.0, 0);
  AR::decimate(&x[0], N, &xd[0], D);  // change order: x, N, y, D
  AR::decimateViaMean(&x[0], N, &xa[0], D);

  // todo: implement and try decimateViaMean

  GNUPlotter plt;
  plt.addDataArrays(N,  &t[0],  &x[0]);
  plt.addDataArrays(ND, &td[0], &xd[0]);
  plt.addDataArrays(ND, &td[0], &xa[0]);
  plt.plot();
}

template<class T>
void getPythagoreanFreqRatios(T* r, int N)
{
  r[0] = T(1);
  for(int i = 1; i < N; i++) {
    r[i] = 1.5 * r[i-1];
    while(r[i] > 2)
      r[i] /= 2;
  }
  rsHeapSort(r, N);
  // todo: implement this algo with rational numbers
}

void pythagoreanTuning()
{
  int N = 15;   // number of different notes - interesting values: 12, 53
  std::vector<double> f(N+1);
  getPythagoreanFreqRatios(&f[0], N);
  f[N] = 2;
  rsStemPlot(f);
  //rsPlotVector(f);
  // see: https://www.youtube.com/watch?v=IT9CPoe5LnM
}

void recursiveSineSweep()
{
  // We create a linearly sweeeping complex exponential from which sine and cosine sweeps can be 
  // obtained by taking real and imaginary parts. We note that a complex eponential with constant
  // frequency can be obtained by the recursion z[n] = a * z[n-1] where a is a complex number of unit
  // magnitude and the angle of a determines the rotation speed. If a = exp(j*w), the angular 
  // frequency w gives the phase increment per sample for the output sine/cosine pair. If we now let
  // a itself be time varying by a[n] = b * a[n-1] where b is another unit magnitude complex number
  // of the form b = exp(j*dw), then dw gives the frequency increment per sample and z[n] will be
  // a sine sweep with linearly increasing frequency. This can be useful in the context of the
  // Bluestein FFT and chirp-Z transform where such linear chirp/sweep signals are used to modulate 
  // the input- and output sequences.

  double fs = 10000.0;   // sample rate
  double f0 = 10.0;      // start frequency
  double df = 30.0;      // frequency increment per second

  static const int N = 10000;
  double y[N];

  double w0 = 2*PI*f0/fs;         // initial normalized radian rotation frequency
  double dw = 2*PI*df/(fs*fs);    // increment for rotation frequency
  rsComplexDbl z(1, 0);           // initial value vor z
  rsComplexDbl j(0, 1);           // imaginary unit
  rsComplexDbl a = exp(j*w0);     // initial rotation factor
  rsComplexDbl b = exp(j*dw);     // multiplier for rotation factor

  for(int n = 0; n < N; n++)
  {
    y[n] = z.real();  // as output, we take the real part
    z *= a;           // update our rotating phasor (multiply by rotation factor)
    a *= b;           // update our rotation factor
  }

  plotData(N, 0.0, 1/fs, y);

  // Observations:
  // We see a linearly sweeping sinusoid that starts at 10 Hz and ends (after 1 second) at 
  // f0 + df = 10 + 30 = 40 Hz. Periods are 0.1s at the start and 0.025s at the end.

  // todo: generalize the idea to obtain a recursive frequency modulation - the factor b should be 
  // such that b.re = cos(u), b.im = sin(u) for some u that is determined by the modulation index.
  // that way, b would always have unit magnitude and its angle would oscillate with zero mean.
}

void recursiveSineWithCubicPhase()
{
  // We try to obtain a recursive sine generator whose phase is given by a cubic polynomial. That
  // means, we want to create a signal:
  //   x(t) = cos(a0 + a1*t + a2*t^2 + a3*t^3)
  // This is useful for an efficient realtime oscillator bank to render sinusoidal models. We can
  // write this as:
  //   x(t) = Re { e^(j*(a0 + a1*t + a2*t^2 + a3*t^3)) }
  //        = Re { e^(j*a0) * e^(j*a1*t) * e^(j*a2*t^2) * e^(j*a3*t^3) }
  // where e^(j*a0) is just a constant, e^(j*a1*t) is a phasor with constant rotation frequency,
  // e^(j*a2*t^2) is a phasor with linearly increasing rotation frequency (something that is 
  // generated in recursiveSineSweep) and e^(j*a3*t^3) has quadratically increasing rotation
  // frequency. Let's call the 4 factors e0,e1,e2,e3, so we have:
  //   e0 = e^(j*a0) 
  //   e1 = e^(j*a1*t)
  //   e2 = e^(j*a2*t^2)
  //   e3 = e^(j*a3*t^3)
  // and see, how we may create each of them recursively, assume time steps in unit increments.
  //   e0[n] = e0[n-1] * p0[n] where p0[n] = 1
  //   e1[n] = e1[n-1] * p1[n] where p1[n] = e^(j*a1)
  //   e2[n] = e2[n-1] * p2[n] where p2[n] = p2[n-1] * q2[n], q2[n] ?= q2[n-1] * e^(j*a2)
  //   e3[n] = e3[n-1] * p3[n] where ...
  // e0 is constant, 
  // p1 is constant -> e1 is linear, 
  // q2 is constant -> p2 is linear -> e2 is quadratic

  //// not yet used:
  //double fs = 44100;
  //double p  = 0.25*PI;
  //double f  = 100;


  int N = 500; // number of samples

  double a0, a1, a2, a3;
  a0 =  0.5;
  a1 =  0.1;
  a2 =  0.0;
  a2 =  -0.00001;
  a3 =  0.0;

  typedef std::complex<double> Complex;
  Complex j(0, 1);           // imaginary unit
  Complex e0, e1, e2, e3;
  Complex     p1, p2, p3;
  Complex         q2, q3;
  Complex             r3;


  // initializations:
  e0 = exp(j*a0);

  p1 = exp(j*a1);
  e1 = 1;

  q2 = exp(j*a2);
  p2 = 1;
  e2 = 1;

  r3 = exp(j*a3);
  q3 = 1;
  p3 = 1;
  e3 = 1;

  Complex ep; // = e0*e1*e2*e3;  // product of the e-values

  typedef std::vector<double> Vec;
  Vec xt(N), x(N);  // target signal and actual signal
  for(int n = 0; n < N; n++)
  {
    double t = double(n);
    xt[n] = cos(a0 + a1*t + a2*t*t + a3*t*t*t);  // generate target signal

    ep  = e0*e1*e2*e3;  // product of the e-values
    x[n] = ep.real();

    // updates/recursion:
    e1 *= p1;    // update rotating phasor 1

    p2 *= q2;    // update rotation factor 2
    e2 *= p2;    // update rotating phasor 2
  }

  GNUPlotter plt;
  plt.addDataArrays(N, &xt[0]);
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &(xt-x)[0]);
  plt.plot();

  // ...doesn't work yet, when a2 != 0 - the e2 update is still wrong (and e3 is not even 
  // implemented yet)
  // it might be a numerical problem - q2 is so close to 1+0j that p2=1 all the time (for 
  // a2 = -0.000002;) - we may have to normalize the time somehow - don't assume a unit time-step 
  // but 1/fs, use user parameters: phase, freq, dFreq, ddFreq....maybe 
  // hmmm...but that probably won't help - in recursiveSineSweep, the b-coeff is also already close
  // to 1+0j - it's actually a wonder that even the linear recursive sweep works so well (hmm - how
  // well does it actually work? there's no comparison to a target signal in the test)
  // maybe try the concept itself with extended precision arithmetic, to see, if the concept would 
  // makes sense at all, when numeric issues are ignored - of course, that doesn't help in 
  // practice, but just to see, if the idea makes sense theoretically

  // hmmm...i think, for realtime additive synthesis, we may have to resort to cosine 
  // approximations (tables or polynomial) or use SIMD-iFFT, generating the output of K voices
  // simultaneously (try with rsFloat32x4 to generate 4 voices at once - can be scaled up using
  // other simd types
}


void ringModNoise()
{
  // We generate a kind of pseudo-noise by a generalized amplitude-/ringmodulation of sinusoids. Our
  // signal will be generated as:
  // y(t) = product_{k=1}^K (a * s_k + b) where s_k = sin(w_k * t + p_k)
  // for a=1,b=0 we have a pure ringmodulation, for a=b=0.5 we have a pure amplitude modulation.

  // user parameters:
  static const int K = 9;          // number of sinusoids
  static const int N = 88200;      // number of samples to generate
  double f0 = 1000.0;              // anchor frequency
  double fs = 44100.0;             // samplerate
  double a  = 1.0;                 // multiplier for sinusoid 
  double b  = 0.5;                 // constant term
  double r  = 5.0;                 // frequency ratio between f[0] and f[K] (f[K] is not actually
                                   // produced, though - we only go up to f[K-1])

  // algorithm variables:
  int k, n;                        // index for sinusoid and sample
  double f[K], w[K], p[K];         // frequencies, normalized radian frequencies, phases
  double x[N];                     // output signal

  // fill arrays with sine frequencies and phases:
  for(k = 0; k < K; k++)
  {
    f[k] = f0 * pow(r, double(k)/K);
    w[k] = 2*PI*f[k]/fs;
    p[k] = 0.0;
  }

  // generate signal:
  for(n = 0; n < N; n++)
  {
    x[n] = 1.0;
    for(k = 0; k < K; k++)
      x[n] *= a*sin(w[k]*n + p[k]) + b;

    //x[n] = (x[n]-b)/b;   // DC adjust and normalize (not sure yet about the formula)
  }

  // we have a DC offset of..

  // normalize and write to file:
  RAPT::rsArrayTools::normalize(x, N, 1.0, true);
  writeToMonoWaveFile("RingModNoise.wav",  x, N, (int) fs, 16);

  // plot signal:
  plotData(8000, 0.0, 1/fs, x);

  // Observations:
  // with a=1.0, b=0.5 all frequencies in the output signal have the same amplitude
  // with exponential frequency spacing of the sinusoids, the signal becomes more spikey as K goes 
  // up, to counteract, one can increase r
  // K=9,f0=1000,a=1,b=0.5,r=5.0: it sounds a bit like vinyl crackling noise
  // ...it sounds a lot like crackling anyway.

  int dummy = 0;
}

void slewRateLimiterLinear()
{
  double fs      = 10;   // samplerate in Hz
  double attack  = 1000; // attack time in ms
  double release = 2000; // release time in ms

  // Attack- and release times are 1 and 2 seconds respectively. With a sample-rate of 10 Hz, that
  // translates to 10 and 20 samples to ramp up and down between 0 and 1, respectively

  rsSlewRateLimiterLinearDD slewRateLimiter;

  slewRateLimiter.setSampleRate(fs);
  slewRateLimiter.setAttackTime(attack);
  slewRateLimiter.setReleaseTime(release);

  static const int N = 300;
  double t[N];  // time axis
  double x[N];  // input signal
  double y[N];  // output signal

  RAPT::rsArrayTools::fillWithIndex(t, N);
  RAPT::rsArrayTools::scale(t, N, 1.0/fs);

  RAPT::rsArrayTools::fillWithZeros(x, N);
  RAPT::rsArrayTools::fillWithZeros(y, N);

  int n;
  for(n = N/3; n < 2*N/3; n++)
    x[n] = 1.0;

  for(n = 0; n < N; n++)
    y[n] = slewRateLimiter.getSample(x[n]);

  plotData(N, t, x, y);
}

void stretchedCorrelation()
{
  // Compute a crosscorrelation value between two signals of different length, where the shorter 
  // one is first stretched (by linear interpolation) to the length of the longer one. This kind 
  // of correlation might be a good measure of waveform similarity regardless of the lengths, 
  // i.e. an array with a waveform and another array with a time compressed version of that 
  // waveform should have a "stretched" correlation value of 1 (up to error due to linear 
  // interpolation)

  static const int N1 = 120;  // length of array 1
  static const int N2 = 100;  // length of array 2
  double x1[N1], x2[N2];

  // fill x1, x2 with sinuosids such that one cycle fits exactly into the respective array:
  double w1 = 2*PI/N1;
  double w2 = 2*PI/N2;
  int n;
  for(n = 0; n < N1; n++)
    x1[n] = sin(w1*n);
  for(n = 0; n < N2; n++)
    x2[n] = sin(w2*n);

  // compute regular and stretched cross-correlation (the 2nd also with swapped arguments)
  double c, cs12, cs21;
  c    = rsCrossCorrelation(x1, N1, x2, N2);           // depends on the length ratio N1/N2
  cs12 = rsStretchedCrossCorrelation(x1, N1, x2, N2);  // should be close to 1
  cs21 = rsStretchedCrossCorrelation(x2, N2, x1, N1);  // should equal cs12

  plotData(N2, 0, 1, x1, x2);
}

void taperedFourierSeries()
{
  static const int numHarmonics = 20;
  static const int numSamples   = 5000;
  double tMin = -PI;
  double tMax = 3*PI;

  double a[numHarmonics+1];   // harmonic amplitudes without tapering
  double aL[numHarmonics+1];  // Lanczos sigma tapering factors
  double aF[numHarmonics+1];  // Fejer tapering factors
  double p[numHarmonics+1];   // phases
  double t[numSamples];
  double x[numSamples];       // signal without tapering
  double xL[numSamples];      // signal with Lanczos tapering
  double xF[numSamples];      // signal with Fejer tapering
  double f[numHarmonics+1];   // normalized frequencies for plot

  RAPT::rsArrayTools::fillWithRangeLinear(f, numHarmonics+1, 0.0, (double) numHarmonics);

  int n, k;
  a[0]  = 0.0;
  p[0]  = 0.0;
  aL[0] = 1.0;
  for(k = 1; k <= numHarmonics; k++)
  {
    a[k] = (2/PI) / k;
    aL[k] = sin(k*PI/numHarmonics) / (k*PI/numHarmonics); // Lanczos sigma factors (sinc-function?)
    p[k] = 0.0;
  }
  for(k = 0; k <= numHarmonics; k++)
    aF[k] = 1.0 - (double) k / (double) (numHarmonics+1);


  // experimental - juggle the phases:
  //p[2] = 0.25*PI;


  RAPT::rsArrayTools::fillWithRangeLinear(t, numSamples, tMin, tMax);
  double tmp;
  for(n = 0; n < numSamples; n++)
  {
    x[n]  = 0.0; 
    xL[n] = 0.0;
    xF[n] = 0.0;
    for(k = 0; k <= numHarmonics; k++)
    {
      tmp    = a[k] * sin(k*t[n] + p[k]);
      x[n]  += tmp;
      xL[n] += aL[k] * tmp;
      xF[n] += aF[k] * tmp;
    }
  }

  double peak = RAPT::rsArrayTools::maxAbs(x, numSamples);

  //plotData(numHarmonics+1, f, aL, aF);
  plotData(numSamples, t, x, xL, xF);
}

// move to filterAnalyzer
void directFormImpulseResponse(double *a, int Na, double *b, int Nb, double *h, int Nh)
{
  double x0 = 1.0;
  //rsFilter(&x0, 1, h, Nh, b, Nb, a, Na);
  RAPT::rsArrayTools::filter(&x0, 1, h, Nh, b, Nb, a, Na);
}

void transientModeling()
{
  // We try to model a transient by an impulse response of an IIR filter. Consider the simple impulse
  // response h[n] = {0.8, 0.2, 0.1}. As model filter, we use a filter of the form
  // y[n] = b0*x[n] + b2*x[n-2] - a1*y[n-1]
  // we can solve for the coefficients by means of a linear system
  // y[0] = h[0] = 0.8 = b0*0.8 + b2*0 - a1*0    -> b0 =  0.8
  // y[1] = h[1] = 0.2 = 0.8*0 + b2*0 - a1*0.8   -> a1 = -0.25
  // y[2] = h[2] = 0.1 = 0.8*0 + b2*1 + 0.25*0.2 -> b2 =  0.05

  // In general, when the impulse response to be modeled is of length N, we use feedforward 
  // coefficients 1...N/2, feedforward coeffs (N+1)/2...N-1 ->check these formulas

  double h[10] ={ 0.8, 0.2, 0.1, -0.2, -0.4, 0.2, 0.1, 0.6, -0.5, -0.1 };
  double a[5], b[5];  

  static const int Ny = 100;  // length of filter output signal
  double y[Ny];

  // model 3-term impulse response:
  a[0] =  1.0;
  a[1] = -0.25;
  b[0] =  0.8;
  b[1] =  0.0;
  b[2] =  0.05;

  directFormImpulseResponse(a, 1, b, 2, y, Ny);

  plotData(10, 0, 1, y);

  // Observations: well, yeah, the first 3 samples match, as constructed, however, i think,
  // the idea will not scale up well.
  // Maybe it's a better idea to create a pole/zero model
  // for the transient and implement that model as cascade of biquads. that way, we would have
  // control over the frequency distribution of the transient (we could scale all frequencies
  // of all poles, zeros) and the time time duration by scaling the Q-values of all poles zeros.
  // we could also broaden or narrow down the transient spectrum by contracting pole/zero 
  // frequencies from a center. So, with such a model, we could edit the transient in a macro
  // oriented way. user parameters could be: duration, frequency-scale, spectral-contraction, 
  // contraction center frequency
  // perhaps the best way to apply transformations is to convert to the s-domain first, transform
  // there and transform back to the z-domain. the filter cofficient values will depend on the
  // samplerate at which the to-be-modeled signal was sampled. it would be nice to be able to
  // store and transform them in a sample-rate independent way. using oversampling in the 
  // estimation in order to obtain more "analog" coefficients might fail because a lot of coeffs
  // might then be used to model the bandlimitation of the oversampled signal.

  int dummy = 0;
}







void windowFunctionsContinuous()
{
  // plots the windows generated by our functions which for continuous time-domain windows

  static const int N = 5000;   // number of values for plot
  double xMin   = -10.0;       // minimum value for x-axis
  double xMax   = +10.0;       // maximum value for x-axis
  double length =  16.0;       // window length

  double x[N];                 // input values
  double wHann[N];             // Hann window values
  double wHamming[N];          // Hamming window values
  double wExactBlackman[N];    // "exact" Blackman window values

  // generate the window functions:
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    wHann[n]          = rsWindowFunction::raisedCosine( x[n], length, 0.0);
    wHamming[n]       = rsWindowFunction::raisedCosine( x[n], length, 0.08);
    wExactBlackman[n] = rsWindowFunction::exactBlackman(x[n], length);
  }

  plotData(N, x, wHann, wHamming, wExactBlackman);
}

/*
// move to library (maybe)
template<class T>
void normalizeMean(T* x, int N) 
{
  T m = RAPT::rsArrayTools::mean(x, N);
  RAPT::rsArrayTools::scale(x, N, T(1)/m);
}
*/

void cosSumWindow2(double* w, int N) // 2 term
{
  double a0 = 0.5, a1 = 0.5;
  for(int n = 0; n < N; n++) {
    //double t = RAPT::rsLinToLin(double(n), -1.0, double(N), -1.0, +1.0);   // NN - good for filter design and spectrum analysis
    //double t = RAPT::rsLinToLin(double(n), 0.0, double(N-1), -1.0, +1.0);  // ZZ - good for nothing
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0);      // ZN - good for comparison with rectangular
    w[n] = a0 + a1*cos(PI*t);
  }
  RAPT::rsArrayTools::normalizeMean(w, N);
  // with the ZN mapping (first sample 0, last nonzero), the zeros are exactly at every other zero 
  // of the rectangular window - best for comparing to rectangular - but for spectral analysis or 
  // filter design NN is better (the zeros are a bit more narrowly spaced), the ZZ mapping is not 
  // useful - it spaces the zeros wider than necessarry
}
void cosSumWindow3(double* w, int N) // 3 term
{
  double a0 = 3./8, a1 = 1./2, a2 = 1./8;
  for(int n = 0; n < N; n++) {
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0); 
    w[n] = a0 + a1*cos(PI*t) + a2*cos(2*PI*t) ;
  }
  RAPT::rsArrayTools::normalizeMean(w, N);
}
void cosSumWindow4(double* w, int N) // 4 term
{
  double a0 = 5./16, a1 = 15./32, a2 = 3./16, a3 = 1./32;
  for(int n = 0; n < N; n++) {
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0); 
    w[n] = a0 + a1*cos(PI*t) + a2*cos(2*PI*t) + a3*cos(3*PI*t);
  }
  RAPT::rsArrayTools::normalizeMean(w, N);
}
void cosSumWindow5(double* w, int N) // 5 term
{
  double a0 = 35./128, a1 = 7./16, a2 = 7./32, a3 = 1./16, a4 = 1./128;
  for(int n = 0; n < N; n++) {
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0); 
    w[n] = a0 + a1*cos(PI*t) + a2*cos(2*PI*t) + a3*cos(3*PI*t) + a4*cos(4*PI*t);
  }
  RAPT::rsArrayTools::normalizeMean(w, N);
}
// these windows are probably nice for filter design because the provide a sidelobe rolloff with 
// increasing steepness
// the mainlobe width is equal to k and the sidelobe rolloff slope equal to 
// (2*k-1) * 6 dB/oct (i think)
// maybe check, if the higher order windows also can be made to sum to a constant for certain 
// overlap factors - if so, they may be useful for spectrogram analysis/resynthesis
// we require more and more even derivatives to be zero at the end of the window support (the
// odd derivatives are all zero anyway)
// todo: 
// -find coeffs for equiripple k-term windows numerically (maybe with numpy/scipy)
// -maybe we can blend between continuity and equiripple/minimax coeffs
// -maybe also allow blend between various orders (k-values) to provide a tradeoff between 
//  mainlobe width and sidelobe rejection
// they were obtained with Sage by imposing contraints on the derivative values at the window
// endpoints, for example, for the 3rd order window:
//
// var("t a0 a1 a2 a3")
// f(t)  = a0 + a1*cos(pi*t) + a2*cos(2*pi*t) + a3*cos(3*pi*t)
// f2(t) = diff(f(t), t, 2) # 2nd derivative
// f4(t) = diff(f(t), t, 4) # 4th derivative
// eq1 = f(0)  == 1         # 1st requirement
// eq2 = f(1)  == 0         # 2nd requirement
// eq3 = f2(1) == 0         # 3rd requirement
// eq4 = f4(1) == 0         # 4th requirement
// solve([eq1,eq2,eq3,eq4],[a0,a1,a2,a3])
//
// giving the result: a0 == (5/16), a1 == (15/32), a2 == (3/16), a3 == (1/32)
// don't count the constant term a0 as term - reduce all numbers in the names by 1


void windowFunctionSpectra()
{
  //int windowLength = 11;
  //int windowLength = 128;
  //int windowLength = 32;
  //int windowLength = 37;
  int windowLength = 38;

  int fftSize = 8192;
  //int fftSize = 16384;

  // create various window functions:
  typedef RAPT::rsWindowFunction WF;
  typedef WF::WindowType WT; 
  int N = windowLength;
  std::vector<double> rectangular(N), triangular(N), hanning(N), hamming(N), 
    blackman(N), blackmanHarris(N),  blackmanNutall(N), nutall(N),
    truncGauss2(N), truncGauss3(N), truncGauss4(N), truncGauss5(N), // 2,3,4,5 = 1/sigma
    salFlatTopFast3(N), salFlatTopFast4(N), salFlatTopFast5(N),
    salFlatTopMin3(N),  salFlatTopMin4(N),  salFlatTopMin5(N),
    hrsFlatTop70(N), hrsFlatTop95(N), hrsFlatTop90D(N), hrsFlatTop116D(N), hrsFlatTop144D(N),
    hrsFlatTop169D(N), hrsFlatTop196D(N), hrsFlatTop223D(N),  hrsFlatTop248D(N);

  WF::createWindow(&rectangular[0],    N, WT::rectangular, true);
  WF::createWindow(&triangular[0],     N, WT::triangularNN,  true);
  WF::createWindow(&hanning[0],        N, WT::hanningZZ,     true);
  WF::createWindow(&hamming[0],        N, WT::hamming,     true);

  WF::createWindow(&blackman[0],       N, WT::blackman,    true);
  WF::createWindow(&blackmanHarris[0], N, WT::blackmanHarris,    true);
  WF::createWindow(&blackmanNutall[0], N, WT::blackmanNutall,    true);
  WF::createWindow(&nutall[0],         N, WT::nutall,             true);

  WF::createWindow(&truncGauss2[0],     N, WT::truncatedGaussian, true, 1/2.);
  WF::createWindow(&truncGauss3[0],     N, WT::truncatedGaussian, true, 1/3.);
  WF::createWindow(&truncGauss4[0],     N, WT::truncatedGaussian, true, 1/4.);
  WF::createWindow(&truncGauss5[0],     N, WT::truncatedGaussian, true, 1/5.);

  WF::salFlatTopFast3(&salFlatTopFast3[0], N);
  WF::salFlatTopFast4(&salFlatTopFast4[0], N);
  WF::salFlatTopFast5(&salFlatTopFast5[0], N);

  WF::salFlatTopMin3(&salFlatTopMin3[0], N);
  WF::salFlatTopMin4(&salFlatTopMin4[0], N);
  WF::salFlatTopMin5(&salFlatTopMin5[0], N);
   
  WF::hrsFlatTop70(  &hrsFlatTop70[0],   N);
  WF::hrsFlatTop95(  &hrsFlatTop95[0],   N);
  WF::hrsFlatTop90D( &hrsFlatTop90D[0],  N);
  WF::hrsFlatTop116D(&hrsFlatTop116D[0], N);
  WF::hrsFlatTop144D(&hrsFlatTop144D[0], N);
  WF::hrsFlatTop169D(&hrsFlatTop169D[0], N);
  WF::hrsFlatTop196D(&hrsFlatTop196D[0], N);
  WF::hrsFlatTop223D(&hrsFlatTop223D[0], N);
  WF::hrsFlatTop248D(&hrsFlatTop248D[0], N);


  // under construction:
  std::vector<double> cosSumWnd2(N), cosSumWnd3(N), cosSumWnd4(N), cosSumWnd5(N);
  cosSumWindow2(&cosSumWnd2[0], N);
  cosSumWindow3(&cosSumWnd3[0], N);
  cosSumWindow4(&cosSumWnd4[0], N);
  cosSumWindow5(&cosSumWnd5[0], N);

  std::vector<double> chebyTweak(N), cheby20(N), cheby40(N), cheby60(N), cheby80(N), cheby100(N);
  cheby_win(&cheby20[0], N,  20);
  cheby_win(&cheby40[0], N,  40);
  cheby_win(&cheby60[0], N,  60);
  cheby_win(&cheby80[0], N,  80);
  cheby_win(&cheby100[0], N, 100);
  cheby_win(&chebyTweak[0], N, 17.5); // tweakable
  // 17.5: mainlobe-width matches rectangular window
  // 46.5: matches cosSumWnd2

  // for testing the new code to produce a chebychev window:
  std::vector<double> cheby60_2(N);
  cheby_win2(&cheby60_2[0], N, 60); 
  rsPlotVectors(cheby60, cheby60_2);  // for debug
  // ...does not work yet


  // compute chebychev window mainlobe width:
  // https://ccrma.stanford.edu/~jos/sasp/Dolph_Chebyshev_Window_Main_Lobe_Width.html
  double r  = rsDbToAmp(-17.5);   // plug attenuation here - the formula should give a width of around 2 at 17.5
  double x0 = cosh(acosh(1/r) / (N-1)); // == chebyPoly(1/r, 1/(N-1))?
  double wc = 2*acos(1/x0);
  double w0Rect = 2*PI/N; // first zero of rectangular window for reference
  //double B  = 2*wc;
  double k = N*wc/(2*PI); // see below
  // hmm...these formulas seem to compute the cutoff frequency in radians dependning on N and give
  // a too small value - we actually want a value independent of N - maybe just leave out the 
  // division by N-1? ...try it:
  //double r  = rsDbToAmp(-17.5);   // plug attenuation here
  //double x0 = cosh(acosh(1/r));  // == 1/r - this is the identity function
  //double wc = 2*acos(1/x0);
  //double B  = 2*wc;
  // no - that doesn't seem to work, for an attenuation of 17.5 dB, we should get a value of 2 but 
  // get 5.7 - more research needed - try to figure out, how the equations above were derived - 
  // maybe they can be reverse engineered to give the normalized mainlobe width in terms of DFT 
  // bins - or maybe we can from the normalized frequency wc compute the bin-index..
  // see also here - the window has impulses at its endpoints:
  // https://ccrma.stanford.edu/~jos/sasp/Example_Chebyshev_Windows_Transforms.html

  // i think, w = 2*pi*f/fs together with f = k*fs/N gives a bin index k = M*w/2pi where k is 
  // actually our desired bin-width B and w is equal to wc above...soo we get the equation
  // wc = B*2pi/N = 2*acos(1/x0) = 2*acos(1/cosh(acosh(1/r) / (N-1))) for N - is that correct?
  // maybe let the SpectrumPlotter scale the frequency axis in different ways - in particular,
  // let it go from 0...PI and see if the formula above for wc gives the right value on that axis

  // oookay - i think, for the rectangular window, the mainlobe width is defined as the first zero 
  // of the spectrum which occurs at 2*PI/N radians
  // for the chebycev window, the cutoff is measured at the point where the mainlobe crosses the
  // attenuation for the first time - that's a little bit below the first zero, but for a rough
  // computation of mainlobe width, it should be good enough - but maybe we can find a formula
  // for the zeros of the chebychev spectrum - it's defined in the freq-domain anyway - i think, 
  // we need a formula for the zeros of chebychev polynomials
  //   Tn(x) = cos(n*acos(x)) for x < 1, so we need to solve cos(n*acos(x)) = 0 for x 
  //   let u = acos(x), the solve cos(n*u) = 0 to find u = pi/2n -> x = acos(u) = acos(pi/2n)
  // ...but i actually thing, defining the width at the first crossing of the attenuation point
  // makes more sense anyway





  // maybe optionally plot the window functions themselves
  // note that gnuplot issues an error when we try to plot the window itself and immediately 
  // thereafter its spectrum, because the in-between call of the convenience function messes up
  // the datafile - we need to do one plot at a time

  typedef SpectrumPlotter<double> SP;
  typedef SP::FreqAxisUnits FU;

  SpectrumPlotter<double> plt;
  plt.setFftSize(fftSize);
  plt.setFloorLevel(-180);
  //plt.setFreqAxisUnit(FU::binIndex);
  //plt.setFreqAxisUnit(FU::normalized);
  plt.setFreqAxisUnit(FU::omega);
  //plt.setShowPhase(true);
  //plt.setZoom(); // show only low portion up to 1/zoom of the spectrum

  //rsPlotVectors(rectangular, triangular, hanning, hamming);
  //plt.plotDecibelSpectra(N, &rectangular[0], &triangular[0], &hanning[0], &hamming[0]);

  //rsPlotVectors(rectangular, blackman, blackmanHarris, blackmanNutall, nutall);
  //plt.plotDecibelSpectra(N, &rectangular[0], &blackman[0], &blackmanHarris[0], &blackmanNutall[0], &nutall[0]);

  //plt.plotDecibelSpectra(N, &rectangular[0], &truncGauss2[0], &truncGauss3[0], &truncGauss4[0], &truncGauss5[0]);

  //rsPlotVectors(rectangular, cosSumWnd2, cosSumWnd3, cosSumWnd4, cosSumWnd5); // ZN
  //plt.plotDecibelSpectra(N, &rectangular[0], &cosSumWnd2[0], &cosSumWnd3[0], &cosSumWnd4[0], &cosSumWnd5[0]);

  //rsPlotVectors(cheby20, cheby40, cheby60, cheby80, cheby100); // 1st value repeated as last (NN)
  //plt.plotDecibelSpectra(N, &cheby20[0], &cheby40[0], &cheby60[0], &cheby80[0], &cheby100[0]);
  plt.plotDecibelSpectra(N, &cheby60[0], &cheby60_2[0]);

  //rsPlotVectors(salFlatTopFast3, salFlatTopFast4, salFlatTopFast5); 
  //plt.plotDecibelSpectra(N, &salFlatTopFast3[0], &salFlatTopFast4[0], &salFlatTopFast5[0]);

  //rsPlotVectors(salFlatTopMin3, salFlatTopMin4, salFlatTopMin5); 
  //plt.plotDecibelSpectra(N, &salFlatTopMin3[0], &salFlatTopMin4[0], &salFlatTopMin5[0]);

  //rsPlotVectors(hrsFlatTop70, hrsFlatTop95, hrsFlatTop90D, hrsFlatTop116D, hrsFlatTop144D, 
  //  hrsFlatTop169D, hrsFlatTop196D, hrsFlatTop223D, hrsFlatTop248D); 
  //plt.plotDecibelSpectra(N, &hrsFlatTop70[0], &hrsFlatTop95[0], &hrsFlatTop90D[0], 
  //  &hrsFlatTop116D[0], &hrsFlatTop144D[0], &hrsFlatTop169D[0], &hrsFlatTop196D[0], 
  //  &hrsFlatTop223D[0], &hrsFlatTop248D[0]);
  // hmm, it seems like the sidelobes are always around 5-6 dB higher than the specifications says
  // not normalizing the windows doesn't change anything (seems, they already are normalized even 
  // without explicitly doing so)
  // try flatTop90D for sinusoidal analysis

  //plt.plotDecibelSpectra(N, &rectangular[0], &chebyTweak[0]);
  //plt.plotDecibelSpectra(N, &cosSumWnd2[0], &chebyTweak[0]);
};


void windowedSinc()
{
  static const int N = 5000;   // number of values
  double xMin    = -10.0;
  double xMax    = +10.0;
  double length  = 16.0;       // window length
  double stretch = 1.2;        // stretch factor of the sinc function

  double x[N];                 // input values
  double s[N];                 // normalized sinc values
  double w[N];                 // window values
  double y[N];                 // windowed sinc values
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    s[n] = rsNormalizedSinc(x[n]/stretch);
    w[n] = rsWindowFunction::cosineSquared(x[n], length);
    y[n] = rsWindowFunction::windowedSinc( x[n], length, stretch);
  }

  plotData(N, x, y, s, w);
}


void waveMorph()
{
  // We consider the unit square and prescribe function values z = f(x,y) along the boundary. We 
  // fill the interior of the square with function values that minimize the curvature of the 
  // resulting surface. ...this can be used to morph waveforms - a read phasor my pphase into the
  // x direction, the y parameter selects the waveform. ...but we could also have a y-phasor or 
  // both....

  // seems to be not very useful - it tends to a flat area in the middle


  int Nx = 41;
  int Ny = 41;



  std::vector<double> x(Nx), y(Ny);

  RAPT::rsMatrix<double> z(Nx, Ny);


  rsNoiseGenerator<double> ng;

  // init boundaries:
  int i, j;
  for(i = 0; i < Nx; i++) {
    x[i]      = double(i) / (Nx-1);
    z(i,0)    = sin(2*PI*x[i]*x[i]*x[i]);  // bottom boundary y=0
    z(i,Ny-1) = sin(2*PI*x[i]);  // top boundary y=1

    //z(i,0)    = ng.getSample();
    //z(i,Ny-1) = ng.getSample();
  }
  for(j = 0; j < Ny; j++) {
    y[j]      = double(j) / (Ny-1);
    z(0,j)    = 0;  // bottom boundary y=0
    z(Nx-1,j) = 0;  // top boundary y=1
    z(0,j)    = sin(3*PI*y[j]);  // bottom boundary y=0
    z(Nx-1,j) = sin(1*PI*y[j]);  // top boundary y=1

    //z(0,j)    = ng.getSample();
    //z(Nx-1,j) = ng.getSample();
  }


  /*
  // init interior by bilinear interpolation:
  double zl, zr, zb, zt, z1, z2;
  double wl, wr, wb, wt;  // weights
  double s;
  for(i = 1; i < Nx-1; i++)
  {
    // interpolate top and bottom row along x: 

    for(j = 1; j < Ny-1; j++) {
      //z[i][j] = 1;
      //z[i][j] = ng.getSample();

      ///z[i][j] = zb;

      //zb = (1-x[i])*z[0][0]    + x[i]*z[Nx-1][0];    // z bottom
      //zt = (1-x[i])*z[0][Ny-1] + x[i]*z[Nx-1][Ny-1]; // z top


      zl = z[0][j];
      zr = z[Nx-1][j];
      zb = z[i][0];
      zt = z[i][Ny-1];

      //z1 = (1-x[i])*zl + x[i]*zr;
      //z2 = (1-y[j])*zb + y[j]*zt;

      wl = 1-x[i];
      wr = x[i];
      wb = 1-y[j];
      wt = y[j];
      s  = wl+wr+wb+wt;
      wl /= s;
      wr /= s;
      wb /= s;
      wt /= s;


      z[i][j] = wl*zl + wr*zr + wb*zb + wt*zt; 
      // this looks wrong - maybe we should use a weighted avergae over the whole boundary
      // but the initialization it is supposed to fade away anyway - we can find some better
      // initialization later

      //z[i][j]= 0.5*(z1+z2);


      //z[i][j] = (1-x[i])*zl + x[i]*zr;

      //z[i][j] = (1-y[j])*zb + y[j]*zt;
    }
  }
  */




  for(int k = 0; k < 100; k++)
  {
    //plotMatrix(z, x, y);
    RAPT::rsMatrix<double> t = z;  // temporary
    for(i = 1; i < Nx-1; i++)
    {
      for(j = 1; j < Ny-1; j++) 
      {
        //double avg = 0.25 * (t(i-1,j-1) + t(i-1,j+1) + t(i+1,j-1) + t(i+1,j+1));

        double avg = 0.2 * (t(i-1,j-1) + t(i-1,j+1) + t(i,j) + t(i+1,j-1) + t(i+1,j+1));

        // todo: maybe try to use diagonal neighbours in the average as well (with weight 
        // 1/sqrt(2))

        double delta = z(i,j) - avg;

        z(i, j) = avg;
        // maybe we need a temporary buffer in order to not overwrite values that will still be
        // needed
      }
    }
  }



  plotMatrix(z, x, y);
  // hmm - may be not that useful - it tends to get flat in the middle

}


void waveMorph2()
{
  // At the moment, this is just a vague idea - not yet implemented

  // We morph one waveshape into another using the wave equation and treating it as a boundary 
  // value problem. This can be seen as giving a string an initial and final shape and using
  // the wave equation to figure out the in-between shapes. These in between shapes are oriented 
  // along the time axis. The whole thing can be seen as figuring out, how a string would most 
  // natually morph from one shape into another by the rules of the wave equation...right?

  // see Höhere Mathematik in Rezepten, page 946 - ah damn - no - the book specifies an initial
  // condition for shape and velocity - not an initial and final shape - nevertheless, the idea is
  // interesting - explore it further - can we do such a thing? if so, we coul perhaps select 
  // different PDEs to govern the morph? maybe this:
  // https://en.wikipedia.org/wiki/Minimal_surface
  // initial and final waveform u0(x) and u1(x) specify the boundary conditions along the two 
  // horizontal lines of the unit square. along the y-axis, we can just connect the the initial 
  // and final amplitudes of the two waveforms - or: we could use two other waveforms for that, 
  // too - that would be interesting
  // see https://pythonhosted.org/algopy/examples/minimal_surface.html

  int n = 100;          // number of spatial samples
  int m = 1000;         // number of temporal samples, m should be > n

  double h = 1. / (n+1);  // spatial stepsize
  double k = 1. / (m+1);  // temporal stepsize



  double* ut = new double[n];
  double** u;
  RAPT::rsMatrixTools::allocateMatrix(u, m, n); 
  // 1st index i time index, 2nd index space index


  // set up boundary conditions:

  PhaseModulationWaveformRenderer wr;

  wr.setCarrierRelativeFrequency(1);
  wr.setModulatorRelativeFrequency(2);
  wr.setModulationIndex(2);
  wr.renderWaveform(u[0], n);            // initial string shape
  wr.setModulatorRelativeFrequency(3);
  wr.renderWaveform(u[m-1], n);          // final string shape
  RAPT::rsArrayTools::fillWithZeros(ut, n);   // initial string velocity




  for(int i = 1; i < m-1; i++)  // loop over time
  {
    for(int j = 0; j < n; j++) // loop over space
    {

    }

  }


  //rsPlotArray(u[0], n);
  rsPlotArrays(n, u[0], u[m-1]);



 
  delete[] ut;
  RAPT::rsMatrixTools::deallocateMatrix(u, m, n);
}