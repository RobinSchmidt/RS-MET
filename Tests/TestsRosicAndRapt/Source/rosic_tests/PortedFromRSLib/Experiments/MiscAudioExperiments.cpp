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
  rsComplexDbl a = exp(j*w0);  // initial rotation factor
  rsComplexDbl b = exp(j*dw);  // multiplier for rotation factor

  for(int n = 0; n < N; n++)
  {
    y[n] = z.real();  // as output, we take the real part
    z *= a;           // update our rotating phasor (multiply by rotation factor)
    a *= b;           // update our rotation factor
  }

  plotData(N, 0.0, 1/fs, y);
  int dummy = 0;

  // Observations:
  // We see a linearly sweeping sinusoid that starts at 10 Hz and ends (after 1 second) at 
  // f0 + df = 10 + 30 = 40 Hz. Periods are 0.1s at the start and 0.025s at the end.

  // todo: generalize the idea to obtain a recursive frequency modulation - the factor b should be 
  // such that b.re = cos(u), b.im = sin(u) for some u that is determined by the modulation index.
  // that way, b would always have unit magnitude and its angle would oscillate with zero mean.
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
  RAPT::rsArray::normalize(x, N, 1.0, true);
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

  RAPT::rsArray::fillWithIndex(t, N);
  RAPT::rsArray::scale(t, N, 1.0/fs);

  RAPT::rsArray::fillWithZeros(x, N);
  RAPT::rsArray::fillWithZeros(y, N);

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

  RAPT::rsArray::fillWithRangeLinear(f, numHarmonics+1, 0.0, (double) numHarmonics);

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


  RAPT::rsArray::fillWithRangeLinear(t, numSamples, tMin, tMax);
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

  double peak = RAPT::rsArray::maxAbs(x, numSamples);

  //plotData(numHarmonics+1, f, aL, aF);
  plotData(numSamples, t, x, xL, xF);
}

// move to filterAnalyzer
void directFormImpulseResponse(double *a, int Na, double *b, int Nb, double *h, int Nh)
{
  double x0 = 1.0;
  //rsFilter(&x0, 1, h, Nh, b, Nb, a, Na);
  RAPT::rsArray::filter(&x0, 1, h, Nh, b, Nb, a, Na);
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
  RAPT::rsArray::fillWithRangeLinear(x, N, xMin, xMax);
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
  T m = RAPT::rsArray::mean(x, N);
  RAPT::rsArray::scale(x, N, T(1)/m);
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
  RAPT::rsArray::normalizeMean(w, N);
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
  RAPT::rsArray::normalizeMean(w, N);
}
void cosSumWindow4(double* w, int N) // 4 term
{
  double a0 = 5./16, a1 = 15./32, a2 = 3./16, a3 = 1./32;
  for(int n = 0; n < N; n++) {
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0); 
    w[n] = a0 + a1*cos(PI*t) + a2*cos(2*PI*t) + a3*cos(3*PI*t);
  }
  RAPT::rsArray::normalizeMean(w, N);
}
void cosSumWindow5(double* w, int N) // 5 term
{
  double a0 = 35./128, a1 = 7./16, a2 = 7./32, a3 = 1./16, a4 = 1./128;
  for(int n = 0; n < N; n++) {
    double t = RAPT::rsLinToLin(double(n), 0.0, double(N), -1.0, +1.0); 
    w[n] = a0 + a1*cos(PI*t) + a2*cos(2*PI*t) + a3*cos(3*PI*t) + a4*cos(4*PI*t);
  }
  RAPT::rsArray::normalizeMean(w, N);
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
  int windowLength = 11;
  //int windowLength = 128;

  int fftSize = 8192;
  //int fftSize = 16384;

  // create various window functions:
  typedef RAPT::rsWindowFunction WF;
  int N = windowLength;
  std::vector<double> rectangular(N), triangular(N), hanning(N), hamming(N), 
    blackman(N), blackmanHarris(N),  blackmanNutall(N), nutall(N),
    truncGauss2(N), truncGauss3(N), truncGauss4(N), truncGauss5(N); // 2,3,4,5 = 1/sigma

  WF::createWindow(&rectangular[0],    N, WF::RECTANGULAR_WINDOW, true);
  WF::createWindow(&triangular[0],     N, WF::TRIANGULAR_WINDOW,  true);
  WF::createWindow(&hanning[0],        N, WF::HANNING_WINDOW,     true);
  WF::createWindow(&hamming[0],        N, WF::HAMMING_WINDOW,     true);

  WF::createWindow(&blackman[0],       N, WF::BLACKMAN_WINDOW,    true);
  WF::createWindow(&blackmanHarris[0], N, WF::BLACKMAN_HARRIS,    true);
  WF::createWindow(&blackmanNutall[0], N, WF::BLACKMAN_NUTALL,    true);
  WF::createWindow(&nutall[0],         N, WF::NUTALL,             true);

  WF::createWindow(&truncGauss2[0],     N, WF::TRUNCATED_GAUSSIAN, true, 1/2.);
  WF::createWindow(&truncGauss3[0],     N, WF::TRUNCATED_GAUSSIAN, true, 1/3.);
  WF::createWindow(&truncGauss4[0],     N, WF::TRUNCATED_GAUSSIAN, true, 1/4.);
  WF::createWindow(&truncGauss5[0],     N, WF::TRUNCATED_GAUSSIAN, true, 1/5.);

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
  // Tn(x) = cos(n*acos(x)) for x < 1, so we need to solve cos(n*acos(x)) = 0 for x 
  // let u = acos(x), the solve cos(n*u) = 0 to find u = pi/2n -> x = acos(u) = acos(pi/2n)





  // maybe optionally plot the window functions themselves

  typedef SpectrumPlotter<double> SP;
  typedef SP::FreqAxisUnits FU;

  SpectrumPlotter<double> plt;
  plt.setFftSize(fftSize);
  //plt.setFreqAxisUnit(FU::binIndex);
  //plt.setFreqAxisUnit(FU::normalized);
  plt.setFreqAxisUnit(FU::omega);
  //plt.setZoom(); // show only low portion up to 1/zoom of the spectrum

  //plt.plotDecibelSpectra(N, &rectangular[0], &triangular[0], &hanning[0], &hamming[0]);
  //plt.plotDecibelSpectra(N, &rectangular[0], &blackman[0], &blackmanHarris[0], &blackmanNutall[0], &nutall[0]);
  //plt.plotDecibelSpectra(N, &rectangular[0], &truncGauss2[0], &truncGauss3[0], &truncGauss4[0], &truncGauss5[0]);
  //plt.plotDecibelSpectra(N, &rectangular[0], &cosSumWnd2[0], &cosSumWnd3[0], &cosSumWnd4[0], &cosSumWnd5[0]);
  //plt.plotDecibelSpectra(N, &cheby20[0], &cheby40[0], &cheby60[0], &cheby80[0], &cheby100[0]);


  plt.plotDecibelSpectra(N, &rectangular[0], &chebyTweak[0]);
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
  RAPT::rsArray::fillWithRangeLinear(x, N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    s[n] = rsNormalizedSinc(x[n]/stretch);
    w[n] = rsWindowFunction::cosineSquared(x[n], length);
    y[n] = rsWindowFunction::windowedSinc( x[n], length, stretch);
  }

  plotData(N, x, y, s, w);
}