#ifndef RS_EXPERIMENTSUTILITIES_H
#define RS_EXPERIMENTSUTILITIES_H

#include "../Common/GNUPlotter.h"

/** Creates a sine-wave of length N (in samples) with frequency f at and amplitude a, 
at samplerate fs. */
void createSineWave(double *x, int N, double f, double a, double fs);

/** Creates a sine-wave of length N (in samples) with an instantaneous frequency given by f at each 
sample and an amplitude given by a, at samplerate fs. */
void createSineWave(double *x, int N, double *f, double a, double fs);

/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, double timeAxis[], double sampleRate);

/** Returns the value of a sum of N harmonic sines with amplitudes given by A and instantaneous 
phases given by n*p, where n is the harmonic number of the sine (starting at 1) */
double sineSum(double p, double *A, double N);

// convenience functions for interfacing with the Plotter class (they manage instantiation and 
// setup of Plotter objects for frequently used cases):
void plotData(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);

void plotData(int N, double x0, double dx, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);
  // for equidistant data, abscissa values start at x0 with increment dx

void plotDataLogX(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);


#endif
