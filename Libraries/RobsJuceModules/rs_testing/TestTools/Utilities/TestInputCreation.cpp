#include "TestInputCreation.h" 

// these will eventually go to the RAPT library:

template <class T>
inline void fillWithZeros(T *buffer, int length)
{
  for(int i = 0; i < length; i++)
    buffer[i] = T(0);
}

template <class T>
inline T rsSawWave(T x)
{
  T tmp = (T)fmod(x, 2*PI);
  if( tmp < PI )
    return T(tmp/PI);
  else
    return T((tmp/PI) - 2);
}

template <class T>
inline T rsSqrWave(T x)
{
  T tmp = (T)fmod(x, 2*PI);
  if( tmp < PI )
    return 1;
  else
    return -1;
}

template <class T>
inline T rsTriWave(T x)
{
  T tmp = (T)fmod(x, 2*PI);
  if( tmp < 0.5f*PI )
    return T(tmp/(0.5f*PI));
  else if( tmp < 1.5f*PI )
    return T(1 - ((tmp-0.5f*PI)/(0.5f*PI)));
  else
    return T(-1 + ((tmp-1.5f*PI)/(0.5f*PI)));
}

//-------------------------------------------------------------------------------------------------
/*
moved to Plotting.cpp in rs_testting
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate)
{
  for(int n = 0; n < numSamples; n++)
    timeAxis[n] = n / sampleRate;
}

void createTimeAxis(int numSamples, double *timeAxis, double sampleRate)
{
  for(int n = 0; n < numSamples; n++)
    timeAxis[n] = n / sampleRate;
}
*/

template<class T>
void createWaveform(T *x, int N, int shape, T frequency, T sampleRate, T phase, bool antiAlias)
{
  T w = (T)(2*PI*frequency/sampleRate);
  RAPT::rsArrayTools::fillWithZeros(x, N);
  switch( shape )
  {
  case 0:
  {
    if( !(frequency >= sampleRate/2 && antiAlias == true) )
    {
      for(int n = 0; n < N; n++)
        x[n] = sin(w*n + phase);
    }
  }
  break;
  case 1:
  {
    if( antiAlias == false )
    {
      for(int n = 0; n < N; n++)
        x[n] = (T) rsSawWave(T(w*n) + T(phase));
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        T a = T(-2.f / (k*PI));
        for(int n = 0; n < N; n++)
          x[n] += T(a * sin(k*(w*n+PI) + phase));
        k++;
      }
    }
  }
  break;
  case 2:
  {
    if( antiAlias == false )
    {
      for(int n = 0; n < N; n++)
        x[n] = rsSqrWave(T(w*n + phase));
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        T a = T(-4.f / (k*PI));
        for(int n = 0; n < N; n++)
          x[n] += T(a * sin(k*(w*n+PI) + phase));
        k+=2;
      }
    }
  }
  break;
  case 3:
  {
    if( antiAlias == false )
    {
      for(int n = 0; n < N; n++)
        x[n] = rsTriWave(w*n + phase);
    }
    else
    {
      int k = 1;
      T s = 1.0; // sign 
      while( k*frequency < sampleRate/2 )
      {
        T a = T(8.f / (k*k*PI*PI));
        for(int n = 0; n < N; n++)
          x[n] += s * a * sin(k*w*n + phase);
        k +=  2;
        s *= -1.0;
      }
    }
  }
  break;
  }
}
template void createWaveform(float *x, int N, int shape, float frequency, float sampleRate, 
  float phase, bool antiAlias);
template void createWaveform(double *x, int N, int shape, double frequency, double sampleRate, 
  double phase, bool antiAlias);


void createPulseWave(double *x, int N, double frequency, double dutyCycle, 
  double sampleRate, double phase, bool antiAlias)
{
  double w = 2*PI*frequency/sampleRate;
  RAPT::rsArrayTools::fillWithZeros(x, N);
  if( antiAlias == false )
  {
    for(int n=0; n<N; n++)
      x[n] = RAPT::rsPulseWave(w*n + phase, dutyCycle);
  }
  else
  {
    int k = 1;
    while( k*frequency < sampleRate/2 )
    {
      double a = 4.0 * sin(k*PI*dutyCycle) / (k*PI);
      for(int n=0; n<N; n++)
        x[n] += a * cos(k*w*n + phase);
      k++;
    }
  }
}

double sineSum(double p, double *A, double N)
{
  double y = 0.0;
  for(int h = 1; h <= N; h++)
    y += A[h-1] * sin(h*p);
  return y;
}

void createSineWave(double *x, int N, double f, double a, double fs, double p)
{
  double w = 2*PI*f/fs;
  //double p = 0;            // startphase - make this an optional parameter later
  for(int n = 0; n < N; n++)
    x[n] = a * sin(w*n+p);
}

template<class T>
void createSineWave(T *x, int N, T *f, T a, T fs)
{
  T s   = T(2*PI)/fs;     // frequency scaler
  T phi = T(0);           // instantaneous phase
  for(int n = 0; n < N; n++)
  {
    T fn = f[n];          // use temporary to allow x and f to be the same array
    x[n] = a * sin(phi);
    phi += s * fn;
  }
}
template void createSineWave(float  *x, int N, float  *f, float  a, float  fs);
template void createSineWave(double *x, int N, double *f, double a, double fs);

void createSineSweep(double* x, int N, double f1, double f2, double fs, double a)
{
  double k = 2*PI/fs;  // conversion factor from frequency to omega
  double w = k*f1;
  double p = 0;
  for(int n = 0; n < N; n++)
  {
    x[n] = a * sin(p);
    w  = k * rsLinToLin(double(n), 0.0, N-1.0, f1, f2);
    p += w;
  }
}

void createSumOfSines(double* x, int numSamples, int numSines, double fs,
  double *f, double *a, double *p)
{
  double s = 2*PI/fs;  // frequency scaler
  double phi;          // instantaneous phase;
  RAPT::rsArrayTools::fillWithZeros(x, numSamples);
  for(int k = 0; k < numSines; k++)
  {
    if(p != nullptr)
      phi = p[k];
    else
      phi = 0;
    for(int n = 0; n < numSamples; n++)
    {
      x[n] += a[k] * sin(phi);
      phi += s * f[k];
    }
  }
}

template<class T>
std::vector<T> createNoise(int numSamples, T min, T max, int seed)
{
  std::vector<T> x(numSamples);
  RAPT::rsNoiseGenerator<T> ng;
  ng.setRange(min, max);
  ng.setSeed(seed);
  ng.reset();
  for(int n = 0; n < numSamples; n++)
    x[n] = ng.getSample();
  return x;
}
template std::vector<float> createNoise(int numSamples, float min, float max, int seed);
template std::vector<double> createNoise(int numSamples, double min, double max, int seed);


template<class T>
std::vector<T> createColoredNoise(int N, T spectralSlope, int seed)
{
  // Create white noise:
  std::vector<T> y(N);
  RAPT::rsNoiseGenerator<T> ng;
  ng.setSeed(seed);
  for(int n = 0; n < N; n++)
    y[n] = ng.getSample();

  // Apply coloring:
  rosic::SlopeFilter flt;
  flt.setSlope(spectralSlope);
  for(int n = 0; n < N; n++)
    y[n] = (float)flt.getSample(y[n]);

  return y;
}
template std::vector<float> createColoredNoise(int N, float spectralSlope, int seed);
template std::vector<double> createColoredNoise(int N, double spectralSlope, int seed);

template<class T>
std::vector<T> createCrackle(int numSamples, T cutoff, int order/*, int seed*/)
{
  std::vector<T> x(numSamples);
  rsNoiseGeneratorTriModal<double> ng;
  ng.setOrder(order);
  //ng.setSeed(seed);
  ng.selectorLowpass.setSampleRate(1.0);
  ng.selectorLowpass.setCutoff(cutoff);
  for(int n = 0; n < numSamples; n++)
    x[n] = 0.5 * ng.getSample();
  return x;
}
template std::vector<double> createCrackle(int numSamples, double cutoff, int order/*, int seed*/);

template<class T>
std::vector<T> randomSampleInstants(int N, T dtMin, T dtMax, int seed)
{
  std::vector<T> t(N);
  t[0] = 0;
  RAPT::rsArrayTools::fillWithRandomValues(&t[1], N-1, dtMin, dtMax, seed);
  RAPT::rsArrayTools::cumulativeSum(&t[0], &t[0], N);
  return t;
}
template std::vector<double> randomSampleInstants(int, double, double, int);


void createRandomDataXY(double* x, double* y, int N, double dxMin, double dxMax,
  double yMin, double yMax, int seedX, int seedY)
{
  RAPT::rsNoiseGenerator<double> ngX, ngY;
  ngX.setRange(dxMin, dxMax);
  ngX.setSeed(seedX);
  ngX.reset();
  ngY.setRange(yMin, yMax);
  ngY.setSeed(seedY);
  ngY.reset();
  x[0] = ngX.getSample();
  y[0] = ngY.getSample();
  for(int i = 1; i < N; i++) {
    x[i] = ngX.getSample() + x[i-1];
    y[i] = ngY.getSample();
  }
}

std::vector<double> createSineWave(int N, double f, double fs, double a, double p)
{
  std::vector<double> x(N);
  createSineWave(&x[0], N, f, a, fs, p);
  return x;
}

void createSawWave(double *x, int N, double f, double fs, double a, int numHarmonics)
{
  if(numHarmonics == -1)
    numHarmonics = (int) ((fs/2)/f); 
  double w = 2*PI*f/fs;
  for(int k = 1; k <= numHarmonics; k++) {
    double b = -2.f * a / (k*PI);
    for(int n = 0; n < N; n++)
      x[n] += b * sin(k*(w*n+PI));
  }
}


std::vector<double> sineAndDeacyingInharmonic(int N, double f, double fs, double decay)
{
  // two partial frequencies:
  double w1 = 2*PI*f/fs;
  double w2 = w1 * (1+sqrt(5));

  // their amplitudes:
  double a1 = 1.0;
  double a2 = 1.0/3.0;
  double tau = 1/decay;  // decay time of 2nd partial (1st is steady)

  std::vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = a1 * sin(w1*n) + a2 * exp(-(n/fs)/tau) * sin(w2*n);
  return x;
}

std::vector<double> twoSinesAndDecayingDc(int N, double f, double fs, double overtoneRatio, 
  double overtoneAmplitude, double dcAmount, double dcDecay)
{
  double w = 2*PI*f/fs;
  std::vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = sin(w*n) + overtoneAmplitude * sin(overtoneRatio*w*n) + dcAmount * exp(-(n/fs)/dcDecay);
  return x;
}

std::vector<double> sawAndSquare(int N, double fs, double fSaw, double aSaw,
  double fSqr, double aSqr, bool antiAlias)
{
  std::vector<double> xSaw(N), xSqr(N), x(N);
  createWaveform(&xSaw[0], N, 1, fSaw, fs, 0., antiAlias);
  createWaveform(&xSqr[0], N, 2, fSqr, fs, 0., antiAlias);
  for(int n = 0; n < N; n++)
    x[n] = aSaw*xSaw[n] + aSqr*xSqr[n];
  return x;
}


double saw(double phi, int kMax)
{
  double sum = 0;
  for(int k = 1; k <= kMax; k++)
    sum += sin(k*phi) / k;
  return (-2/PI) * sum;
}

double sqr(double phi, int kMax)
{
  double sum = 0;
  for(int k = 1; k <= kMax; k+=2)
    sum += sin(k*phi) / k;
  return (4/PI) * sum;
}

double saw(int n, double f, double fs, int kMax)
{
  if(kMax == -1)
    kMax = (int) floor(0.5*fs/f);
  double w   = 2*PI*f/fs;
  double phi = fmod(w*n, 2*PI);
  return saw(phi, kMax);
}

double sqr(int n, double f, double fs, int kMax)
{
  if(kMax == -1)
    kMax = (int) floor(0.5*fs/f);
  double w   = 2*PI*f/fs;
  double phi = fmod(w*n, 2*PI);
  return sqr(phi, kMax);
}

void createModalPluck(double* x, int N, double key, double sampleRate)
{
  rosic::rsModalSynth ms;
  ms.setSampleRate(sampleRate);
  //ms.setAmpSlope(-3);
  ms.setAmpSlope(-6);
  ms.setDecay(500.0);           // in ms
  ms.setDecayByRatio(-100.0);   // in %
  ms.setAttack(10.0);
  ms.noteOn(key, 64.0);
  double dummy; // unused right channel output
  for(int n = 0; n < N; n++) {
    ms.getSampleFrameStereo(&x[n], &dummy);
    x[n] *= 0.5;  // maybe make the amplitude a parameter
  }
}
std::vector<double> createModalPluck(double key, double sampleRate, int length)
{
  std::vector<double> x(length);
  createModalPluck(&x[0], length, key, sampleRate);
  return x;
}

std::vector<double> createModalBellGloriosa(double sampleRate, int length)
{
  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;
  using MFB = rsModalFilterBank<double, double>;

  // Measurements were taken from
  //   https://www.youtube.com/watch?v=9ssM8cMYmbQ
  // from around 2:55, selecting the whole range where the bell rings and inspecting it with the 
  // spectrum plot in Audacity using the cursor wo figure out the frequencies. They are somewhat 
  // inexact, but we intend to sort of "quantize" them to equal temperament or to just intonation
  // anyway.
  Vec frq = {  38,    48,    81,   130,   168,   200,   243,   265,   330,   417,   489,   533,   673,   875,  1089,  1313,  1544,  1772,  2013,  2248,  2491,  2716,  2944,  3185,  3420,  3432,  3623,  3885    };
  Vec lvl = { -66.9, -64.9, -37.7, -57.9, -24.9, -20.8, -35.4, -54.3, -27.9, -36.4, -32.8, -43.6, -32.1, -35.0, -42.7, -48.2, -49.7, -62.3, -63.4, -62.5, -69.7, -63.5, -75.2, -74.7, -77.9, -86.5, -77.2, -85.0  };

  // Convert from decibels to raw amplitude;
  int M = (int) frq.size();  // M: number of modes
  Vec amp(M);
  for(int m = 0; m < M; m++)
    amp[m] = rsDbToAmp(lvl[m]);

  // Convert from frequencies in Hz to relative frequencies with respect to the fundamental. I 
  // think, 81 is the fundamental ...but not so sure - could be a subharmonic as well -> figure out!
  double fundamental = 168;
  for(int m = 0; m < M; m++)
    frq[m] /= fundamental;

  // Ad-hoc: We just make something up for the attacks and decays - later we want an algorithm to 
  // measure these:
  double decay  = 1.5;          // in seconds
  double attack = 0.1*decay;
  Vec att(M), dec(M);
  double p = 0.4;
  for(int m = 0; m < M; m++)
  {
    dec[m] = pow(m+1, -p);
    att[m] = dec[m];
  }

  // For the phases, we just use zero. Later, we want the algorithm to measure these, too:
  Vec phs(M);
  phs = rsRandomVector(M, 0, 0, 0);  // use rsZeroVector

  // Generate the sound and return it:
  int N = length;
  MFB mfb;
  mfb.setSampleRate(sampleRate);
  mfb.setReferenceFrequency(fundamental);
  mfb.setReferenceAttack(attack);
  mfb.setReferenceDecay(decay);
  mfb.setModalParameters(frq, amp, att, dec, phs);
  Vec x(N);  // the main signal
  x[0] = mfb.getSample(1.0);
  for(int n = 1; n < N; n++)
    x[n] = mfb.getSample(0.0);

  // Post-process:
  AT::normalize(&x[0], N);
  double fadeOutMs = 50;    // 50 milliseconds fade-out
  int fadeOutSamples = (int) round(sampleRate * fadeOutMs / 1000);
  RAPT::rsFadeOut(&x[0], N-fadeOutSamples, N-1);
  return x;

  // ToDo:
  // -Maybe apply a smooth fade-out over 50 ms

  // OK - it goes into the right direction but still sounds kinda wrong. The original sound a a lot 
  // brighter. I think, this is because we conflate amplitudes and decay times, attribute our 
  // measured amplitude data solely to the amplitudes and additionaly (!) apply fatser decays to 
  // higher partials. I hope that estimating the decay times will at least partially fix this 
  // problem. There's also a lot of the complexity missing - but maybe this is mostly about reverb?
  // Maybe some warble is going on, too which we don't capture. And maybe we need more of these 
  // intermediate in-between partials with lower amplitude
}

void applyVibrato(double *x, int N, double freq, double sampleRate, double depth)
{
  rosic::Vibrato vib;
  vib.setSampleRate(sampleRate);
  vib.setCycleLength(1/freq);
  vib.setDepth(depth);
  double dummy = 0;
  for(int n = 0; n < N; n++)
    vib.getSampleFrameStereo(&x[n], &dummy);
}

// maybe move to rapt
bool startsWith(const std::string& str, const std::string& pattern)
{
  int index = RAPT::rsFindFirstOccurrenceOf(
    str.c_str(), (int)str.size(), pattern.c_str(), (int) pattern.size());
  if(index == 0) 
    return true;
  return false;
}

// move to rsArrayTools:
template <class T>
int firstNonMatchElement(
  const T *buffer, int bufferLength, const T *elementsToFind, int numElements)
{
  for(int i = 0; i < bufferLength; i++)
    if(!RAPT::rsArrayTools::contains(elementsToFind, numElements, buffer[i]))
      return i;
  return -1; // all elements matched any one of the elements to find
}

double getValue(const std::string& str, const std::string& key, double defaultValue)
{
  // find the start of the substring that contains the keyed value:
  std::string ptn = key + "=";
  int index = RAPT::rsFindFirstOccurrenceOf(
    str.c_str(), (int)str.size(), ptn.c_str(), (int) ptn.size());
  if(index == -1)  
    return defaultValue;   // key was not found

  // figure out position and length of the number substring:
  index += (int) key.size() + 1;          // first character after the '='
  int length = (int) str.size() - index;  // from index to end of str
  length = firstNonMatchElement(&str[index], length, "0123456789-.e", 13); 
  if(length == -1)
    length = (int) str.size() - index;    // substring goes up to the end

  // extract number substring as c-string and convert to double:
  char* numStr = new char[length+1];
  RAPT::rsArrayTools::copy(&str[index], numStr, length);
  numStr[length] = '\0';
  double value = atof(numStr);
  delete[] numStr;
  return value;
}

std::vector<double> attackDecayEnvelope(int N, double attackSamples, double decaySamples)
{
  std::vector<double> env(N);
  rsModalFilterWithAttack<double, double> flt;
  flt.setModalParameters(0.0, 1.0, attackSamples, decaySamples, 90.0, 1.0);
  env[0] = flt.getSample(1);
  for(int n = 1; n < N; n++)
    env[n] = flt.getSample(0);
  return env;
}

std::vector<double> attackDecaySine(int N, double frequency, double amplitude, double attack, 
  double decay, double startPhase, double sampleRate)
{
  std::vector<double> y(N);
  rsModalFilterWithAttack<double, double> flt;
  flt.setModalParameters(frequency, amplitude, attack, decay, startPhase, sampleRate);
  y[0] = flt.getSample(1);
  for(int n = 1; n < N; n++)
    y[n] = flt.getSample(0);
  return y;
}

std::vector<double> getBrownZap(int numStages, double lowFreq, double highFreq, double freqShape, 
  double lowQ, double highQ, double qShape, double maxLength, int sampleRate)
{
  // Create and set up the zapper object:
  rosic::rsFlatZapper wz;
  wz.setSampleRate(sampleRate);
  wz.setNumStages(numStages);
  wz.setLowFreq(lowFreq);
  wz.setHighFreq(highFreq);
  wz.setFreqShape(freqShape);
  wz.setLowQ(lowQ);
  wz.setHighQ(highQ);
  wz.setQShape(qShape);

  // Render sample:
  int N = ceil(maxLength * sampleRate);  // Number of samples to render
  using Vec = std::vector<double>;
  Vec x(N);
  x[0] = 1;

  // Test - with noise-burst:
  //x = createNoise(N, -1.0, +1.0, 0);
  //Vec e1 = attackDecayEnvelope(N, 50,     150);  // Decay = 100-200 seems nice
  //Vec e2 = attackDecayEnvelope(N, 5000, 10000);
  //Vec e  = e1;
  //Vec e  = e1 + 0.001*e2;   // Trying to give it that snare effect - does not work well
  //rsPlotVector(e);
  //rsPlotVectors(e1, e2);
  //x = x*e;
  //rsPlotVector(x);
  // Using a noise-burst makes it sound more acoustic and natural, less electronic. But maybe this
  // doesn't count as raw-material anymore because it can be recreated from the impulse responses 
  // using convolution. Actually, the user could just play the kick sample through a convolution
  // reverb that uses a noise-burst as impulse-response

  for(int n = 0; n < N; n++)
    x[n] = wz.getSample(x[n]);
  //rsPlotVector(x);

  // Optionally plot a phase-spectrum of the allpass impulse response:
  bool plotPhase = false;  // Maybe make this a function parameter
  if(plotPhase)
  {
    SpectrumPlotter<double> sp;
    sp.setFftSize(65536);
    sp.setLogFreqAxis(true);
    sp.setSampleRate(sampleRate);
    sp.setFreqAxisUnit(SpectrumPlotter<double>::FreqAxisUnits::hertz);
    //sp.plotDecibelSpectra(N, &x[0]);  // Should be flat. Yep - it is.
    sp.plotPhaseSpectra(N, &x[0]);
    // ToDo: Maybe alternatively plot phase-delay or group-delay
  }


  // Post-process the white zap. First we turn the spectrum from white to brown by applying a first
  // order lowpass tuned somewhere below the lowest allpass tuning freq. The resulting -6 dB/oct 
  // magnitude response is ideal in the sense that the amplitude stays constant during the sweep. 
  // Then we do a bit of cleanup by applying a highpass. This removes some bumpiness that is 
  // introduced by the lowpass (it looks like we get some sort of undulating DC towards the end in
  // the lowpass - we want to remove that). Finally, the sample is shortened to remove the silence
  //  in the tail:
  rsSamplePostProcessor pp;
  pp.setSampleRate(sampleRate);
  pp.applyOnePoleLowpass( x, 0.5*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.shortenTail(x, -60.0, 0.02, 0.25/lowFreq); 
  N = x.size();
  RAPT::rsArrayTools::normalize(&x[0], N);  // implement and use sp.normalize(x);
  // Maybe a 3rd order Butterworth highpass would be better than applying a 1st order highpass 3
  // times? Try it! Maybe also try a somwhat lower cutoff for the highpass like 0.75*lowFreq.
  // Actually, it turns out that there's a better way top reduce the subsonic flutter/rumble:
  // Reduce the lowQ parameter.

  return x;
}


std::vector<double> createNamedSound(const std::string& s, double fs, int N)
{

  std::vector<double> v(N);  // vector for the signal
  double* x = &v[0];         // pointer to first sample (for convenience)

  // maybe always retrieve Freq and Amp and store in local variables for later use..


  if(startsWith(s, "SineAndDC")) {
    double freq  = getValue(s, "Freq", 200);
    double amp   = getValue(s, "Amp",  1);
    double dc    = getValue(s, "DC",   1);
    createSineWave(x, N, freq, amp, fs, 0.0);
    RAPT::rsArrayTools::add(x, dc, x, N); 
  }
  else if( startsWith(s, "Sine") )
    createSineWave(x, N, getValue(s, "Freq", 100), getValue(s, "Amp", 1), fs, 0.0);
  else if( startsWith(s, "Cosine") ) 
    createSineWave(x, N, getValue(s, "Freq", 100), getValue(s, "Amp", 1), fs, PI/2);
  else if(startsWith(s, "TwoSines")) {
    double f[2] = { getValue(s, "Freq1", 100), getValue(s, "Freq2", 200) };
    double a[2] = { getValue(s, "Amp1",  1),   getValue(s, "Amp2",  1)    };
    createSumOfSines(x, N, 2, fs, f, a);
  }
  else if(startsWith(s, "ThreeSines")) {
    double f[3] = { getValue(s, "Freq1", 100), getValue(s, "Freq2", 200), getValue(s, "Freq3", 300) };
    double a[3] = { getValue(s, "Amp1",  1),   getValue(s, "Amp2",  1),   getValue(s, "Amp3",  1)   };
    createSumOfSines(x, N, 3, fs, f, a);
  }
  else if(startsWith(s, "FourSines")) {
    double f[4] = { getValue(s, "Freq1", 100), getValue(s, "Freq2", 200),
      getValue(s, "Freq3", 300), getValue(s, "Freq4", 400) };
    double a[4] = { getValue(s, "Amp1", 1), getValue(s, "Amp2", 1), getValue(s, "Amp3", 1),
      getValue(s, "Amp4", 1) };
    createSumOfSines(x, N, 4, fs, f, a);
  }
  else if(startsWith(s, "FiveSines")) {
    double f[5] = { getValue(s, "Freq1", 100), getValue(s, "Freq2", 200), 
      getValue(s, "Freq3", 300), getValue(s, "Freq4", 400), getValue(s, "Freq5", 500) };
    double a[5] = { getValue(s, "Amp1", 1), getValue(s, "Amp2", 1), getValue(s, "Amp3", 1), 
      getValue(s, "Amp4", 1), getValue(s, "Amp5", 1) };
    createSumOfSines(x, N, 5, fs, f, a);
  }
  // maybe, it would make sense to have a general "Sines..." option that can take an arbitrary number 
  // of freqs and amps



  else if(startsWith(s, "VibratoSine")) { 
    double freq  = getValue(s, "Freq", 200);
    double rate  = getValue(s, "Rate",  10);
    double depth = getValue(s, "Depth", 10) * freq * 0.01; // 0.01 because input is in percent
    createSineWave(x, N, rate, depth, fs, 0.0);            // write sine LFO output into x
    RAPT::rsArrayTools::add(x, freq, x, N);                     // add the center freq
    createSineWave(x, N, x, 0.5, fs);                      // overwrite x by vibratoed sinewave
  }
  else if(startsWith(s, "TremoloSine")) { 
    double freq  = getValue(s, "Freq", 200);
    double rate  = getValue(s, "Rate",  10);
    double depth = getValue(s, "Depth", 10) * 0.01;        // 0.01 because input is in percent
    std::vector<double> a(N);                              // vector for instantaneous amplitude
    createSineWave(&a[0], N, rate, depth, fs, 0.0);        // write sine LFO output into a
    createSineWave(x, N, freq, 0.5, fs);                   // write non-modulated sine into x
    RAPT::rsArrayTools::add(&a[0], 1.0, &a[0], N);              // add one to mod-signal
    RAPT::rsArrayTools::multiply(x, &a[0], x, N);               // apply amp-modulation to x
  }
  else if(startsWith(s, "ModalPluck")) {
    double key = rsFreqToPitch(getValue(s, "Freq", 200));
    createModalPluck(x, N, key, fs);
  }
  else if(startsWith(s, "LowpassSaw")) {
    double freq  = getValue(s, "Freq", 200);
    double kMax  = getValue(s, "kMax", -1);
    createSawWave(x, N, freq, fs, 0.5, (int) kMax);
  }


  else rsError("Unknown sound name");

  return v;
}

//=================================================================================================

void rsSamplePostProcessor::applyOnePoleLowpass(double* x, int N, double cutoff)
{
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setSampleRate(sampleRate);
  flt.setMode(flt.LOWPASS_IIT);
  flt.setCutoff(cutoff);
  for(int n = 0; n < N; n++)
    x[n] = flt.getSample(x[n]);
}

void rsSamplePostProcessor::applyOnePoleHighpass(double* x, int N, double cutoff)
{
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setSampleRate(sampleRate);
  flt.setMode(flt.HIGHPASS_MZT);
  flt.setCutoff(cutoff);
  for(int n = 0; n < N; n++)
    x[n] = flt.getSample(x[n]);
}

void rsSamplePostProcessor::shortenTail(std::vector<double>& x, double thresholdDb, 
  double fadeOutTime, double releaseTime)
{
  RAPT::rsEnvelopeFollower<double, double> ef;
  ef.setSampleRate(sampleRate);
  ef.setAttackTime(0.0);
  ef.setReleaseTime(1000 * releaseTime);
  int N = (int) x.size();
  std::vector<double> env(N);
  for(int n = 0; n < N; n++)
    env[n] = ef.getSample(x[n]);
  //rsPlotVectors(x, env);
  // I tried to use the more advanced rsEnvelopeFollower2 but this produced total garbage results
  // for this sort of signal. The simpler works much better.

  // Find last sample that exceeds the threshold:
  double envMax   = RAPT::rsArrayTools::maxValue(&env[0], N);
  double thresh   = RAPT::rsDbToAmp(thresholdDb);
  thresh *= envMax;  // threshold should be relative
  int nCut = N-1;
  while(nCut > 0)
  {
    if(env[nCut] >= thresh)
      break;
    nCut--;
  }

  // Shorten the signal:
  N = nCut+1;
  x.resize(N);

  // Apply a smooth fade out:
  int fadeSamples = sampleRate * fadeOutTime;
  fadeSamples = rsMin(fadeSamples, N/4);
  rsFadeOut(&x[0], N-fadeSamples-1, N-1);
  //rsPlotVectors(x, env);
}
