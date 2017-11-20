#include "ExperimentsUtilities.h"
using namespace RSLib;

// Testsignals:

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

void createTimeAxis(int numSamples, double timeAxis[], double sampleRate)
{
  rsFillWithIndex(timeAxis, numSamples);
  rsScale(timeAxis, numSamples, 1.0/sampleRate);
}

void createWaveform(double *x, int N, int shape, double frequency, double sampleRate, 
  double phase, bool antiAlias)
{
  double w = 2*PI*frequency/sampleRate;
  rsFillWithZeros(x, N);
  switch( shape )
  {
  case SINE:
  {
    if( !(frequency >= sampleRate/2 && antiAlias == true) )
    {
      for(int n=0; n<N; n++)
        x[n] = sin(w*n + phase);
    }
  }
  break;
  case SAW:
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
        double a = -2.0 / (k*PI);
        for(int n=0; n<N; n++)
          x[n] += a * sin(k*(w*n+PI) + phase);
        k++;
      }
    }
  }
  break;
  case SQUARE:
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
        double a = -4.0 / (k*PI);
        for(int n=0; n<N; n++)
          x[n] += a * sin(k*(w*n+PI) + phase);
        k+=2;
      }
    }
  }
  break;
  case TRIANGLE:
  {
    if( antiAlias == false )
    {
      for(int n=0; n<N; n++)
        x[n] = rsTriWave(w*n + phase);
    }
    else
    {
      int    k = 1;
      double s = 1.0; // sign 
      while( k*frequency < sampleRate/2 )
      {
        double a = 8.0 / (k*k*PI*PI);
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

void createPulseWave(double *x, int N, double frequency, double dutyCycle, 
  double sampleRate, double phase, bool antiAlias)
{
  double w = 2*PI*frequency/sampleRate;
  rsFillWithZeros(x, N);
  if( antiAlias == false )
  {
    for(int n=0; n<N; n++)
      x[n] = rsPulseWave(w*n + phase, dutyCycle);
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


// Plotting:

 // new - does not yet work:
void plotData(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

void plotData(int N, double x0, double dx, double *y1, double *y2, double *y3, double *y4, 
  double *y5)
{
  double *x = new double[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + n*dx;
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
  delete[] x;
}

void plotDataLogX(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.setLogScale("x");
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

void plotVector(std::vector<double> v)
{
  int N = v.size();
  double *y = new double[N];
  for(int n = 0; n < N; n++)
    y[n] = v[n];
  GNUPlotter plt;
  plt.plotArrays(N, y);
  delete[] y;
}

// use shorter variable names here...
void plotSpectrogram(int numFrames, int numBins, double **s, double fs, int H)
{
  // fs: sample rate, H: hop size

  GNUPlotter p;

  // create time- and frequency axis and add the data:
  double tMax = H*(numFrames-1)/fs;
  double fMax =  0.5*fs*(numBins-1)/numBins;
  double *t = new double[numFrames];
  double *f = new double[numBins];
  //p.rangeLinear(t, numFrames, 0.0, H*(numFrames-1)/fs );
  //p.rangeLinear(f, numBins,   0.0, 0.5*fs*(numBins-1)/numBins); 
  p.rangeLinear(t, numFrames, 0.0, tMax);
  p.rangeLinear(f, numBins,   0.0, fMax); 


  p.setDataPrecision(4);                                  // make temp files smaller
  p.addDataMatrix(numFrames, numBins, t, f, s);

  // set up style settings:
  //p.setAxisLabels("Time in s", "Frequency in Hz", "Level in dB");
  p.setPixelSize(800, 400);
  p.addGraph("i 0 nonuniform matrix w image notitle");   

  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  p.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale


  // maybe use signal-dependent values later:
  double dbMin = -100.0;
  double dbMax =   10.0;
  p.setRange(0.0, tMax, 0.0, fMax, dbMin, dbMax);
    // doesn't seem to have any effect on the color axis (z). it seems like GNUPlot uses 
    // autoscaling for it, no matter what.

  p.plot();
  delete[] t;
  delete[] f;
}
