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

template<class T>
void createWaveform(T *x, int N, int shape, T frequency, T sampleRate, T phase, bool antiAlias)
{
  T w = (T)(2*PI*frequency/sampleRate);
  RAPT::rsArray::fillWithZeros(x, N);
  switch( shape )
  {
  case 0:
  {
    if( !(frequency >= sampleRate/2 && antiAlias == true) )
    {
      for(int n=0; n<N; n++)
        x[n] = sin(w*n + phase);
    }
  }
  break;
  case 1:
  {
    if( antiAlias == false )
    {
      for(int n=0; n<N; n++)
        x[n] = (T) rsSawWave(T(w*n) + T(phase));
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        T a = T(-2.f / (k*PI));
        for(int n=0; n<N; n++)
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
      for(int n=0; n<N; n++)
        x[n] = rsSqrWave(T(w*n + phase));
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        T a = T(-4.f / (k*PI));
        for(int n=0; n<N; n++)
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
      for(int n=0; n<N; n++)
        x[n] = rsTriWave(w*n + phase);
    }
    else
    {
      int k = 1;
      T s = 1.0; // sign 
      while( k*frequency < sampleRate/2 )
      {
        T a = T(8.f / (k*k*PI*PI));
        for(int n=0; n<N; n++)
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
  RAPT::rsArray::fillWithZeros(x, N);
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

void createSineWave(double *x, int N, double f, double a, double fs)
{
  double w = 2*PI*f/fs;
  double p = 0;            // startphase - make this an optional parameter later
  for(int n = 0; n < N; n++)
    x[n] = a * sin(w*n+p);
}

void createSineWave(double *x, int N, double *f, double a, double fs)
{
  double s   = 2*PI/fs;  // frequency scaler
  double phi = 0.0;      // instantaneous phase
  for(int n = 0; n < N; n++)
  {
    x[n] = a * sin(phi);
    phi += s * f[n];
  }
}

void createSumOfSines(double* x, int numSamples, int numSines, double fs,
  double *f, double *a, double *p)
{
  double s = 2*PI/fs;  // frequency scaler
  double phi;          // instantaneous phase;
  RAPT::rsArray::fillWithZeros(x, numSamples);
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


std::vector<double> createSineWave(int N, double f, double fs)
{
  vector<double> x(N);
  createSineWave(&x[0], N, f, 1.0, fs);
  return x;
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

  vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = a1 * sin(w1*n) + a2 * exp(-(n/fs)/tau) * sin(w2*n);
  return x;
}

std::vector<double> twoSinesAndDecayingDc(int N, double f, double fs, double overtoneRatio, 
  double overtoneAmplitude, double dcAmount, double dcDecay)
{
  double w = 2*PI*f/fs;
  vector<double> x(N);
  for(int n = 0; n < N; n++)
    x[n] = sin(w*n) + overtoneAmplitude * sin(overtoneRatio*w*n) + dcAmount * exp(-(n/fs)/dcDecay);
  return x;
}

std::vector<double> sawAndSquare(int N, double fs, double fSaw, double aSaw,
  double fSqr, double aSqr, bool antiAlias)
{
  vector<double> xSaw(N), xSqr(N), x(N);
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