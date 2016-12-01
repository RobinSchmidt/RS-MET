#include "TestInputCreation.h" 

// these will eventually go to the RAPT library:

template <class T>
inline void rsFillWithZeros(T *buffer, int length)
{
  for(int i = 0; i < length; i++)
    buffer[i] = T(0);
}

inline float rsSawWave(float x)
{
  float tmp = (float)fmod(x, 2*PI_F);
  if( tmp < PI_F )
    return tmp/PI_F;
  else
    return (tmp/PI_F) - 2;
}

inline float rsSqrWave(float x)
{
  float tmp = fmod(x, 2*PI_F);
  if( tmp < PI_F )
    return 1;
  else
    return -1;
}

inline float rsTriWave(float x)
{
  float tmp = fmod(x, 2*PI_F);
  if( tmp < 0.5f*PI_F )
    return tmp/(0.5f*PI_F);
  else if( tmp < 1.5f*PI_F )
    return 1 - ((tmp-0.5f*PI_F)/(0.5f*PI_F));
  else
    return -1 + ((tmp-1.5f*PI_F)/(0.5f*PI_F));
}

//-------------------------------------------------------------------------------------------------

void createTimeAxis(int numSamples, float *timeAxis, float sampleRate)
{
  for(int n = 0; n < numSamples; n++)
    timeAxis[n] = n / sampleRate;
}

void createWaveform(float *x, int N, int shape, float frequency, float sampleRate, 
  float phase, bool antiAlias)
{
  float w = (float)(2*PI*frequency/sampleRate);
  rsFillWithZeros(x, N);
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
        x[n] = rsSawWave(w*n + phase);
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        float a = -2.f / (k*PI_F);
        for(int n=0; n<N; n++)
          x[n] += a * sin(k*(w*n+PI_F) + phase);
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
        x[n] = rsSqrWave(w*n + phase);
    }
    else
    {
      int k = 1;
      while( k*frequency < sampleRate/2 )
      {
        float a = -4.f / (k*PI_F);
        for(int n=0; n<N; n++)
          x[n] += a * sin(k*(w*n+PI_F) + phase);
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
      int   k = 1;
      float s = 1.0; // sign 
      while( k*frequency < sampleRate/2 )
      {
        float a = 8.f / (k*k*PI_F*PI_F);
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

//void createPulseWave(double *x, int N, double frequency, double dutyCycle, 
//  double sampleRate, double phase, bool antiAlias)
//{
//  double w = 2*PI*frequency/sampleRate;
//  rsFillWithZeros(x, N);
//  if( antiAlias == false )
//  {
//    for(int n=0; n<N; n++)
//      x[n] = rsPulseWave(w*n + phase, dutyCycle);
//  }
//  else
//  {
//    int k = 1;
//    while( k*frequency < sampleRate/2 )
//    {
//      double a = 4.0 * sin(k*PI*dutyCycle) / (k*PI);
//      for(int n=0; n<N; n++)
//        x[n] += a * cos(k*w*n + phase);
//      k++;
//    }
//  }
//}
//
//double sineSum(double p, double *A, double N)
//{
//  double y = 0.0;
//  for(int h = 1; h <= N; h++)
//    y += A[h-1] * sin(h*p);
//  return y;
//}
//
//void createSineWave(double *x, int N, double f, double a, double fs)
//{
//  double w = 2*PI*f/fs;
//  double p = 0;            // startphase - make this an optional parameter later
//  for(int n = 0; n < N; n++)
//    x[n] = a * sin(w*n+p);
//}
//
//void createSineWave(double *x, int N, double *f, double a, double fs)
//{
//  double s   = 2*PI/fs;  // frequency scaler
//  double phi = 0.0;      // instantaneous phase
//  for(int n = 0; n < N; n++)
//  {
//    x[n] = a * sin(phi);
//    phi += s * f[n];
//  }
//}

