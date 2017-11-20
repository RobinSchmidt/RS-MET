#include "ExperimentsUtilities.h"

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
    x[n] = sin(w*n+p);
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



void plotData(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  Plotter plt;
  plt.plotData(N, x, y1, y2, y3, y4, y5);
}

void plotData(int N, double x0, double dx, double *y1, double *y2, double *y3, double *y4, 
  double *y5)
{
  double *x = new double[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + n*dx;
  Plotter plt;
  plt.plotData(N, x, y1, y2, y3, y4, y5);
  delete[] x;
}

void plotDataLogX(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  Plotter plt;
  plt.setLogScaleX(true);
  plt.plotData(N, x, y1, y2, y3, y4, y5);
}
