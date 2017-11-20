#ifndef RS_EXPERIMENTSUTILITIES_H
#define RS_EXPERIMENTSUTILITIES_H

#include "../Common/GNUPlotter.h"
#include "../Common/TestUtilities.h"

// Testsignal Generation:

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

/** Synthesizes a standard waveform at the desired frequency and samplerate. The phase is 
expected in radians and there is an option to do anti-aliased synthesis (by means of adding
the sinusoidal components up to the Nyquist frequency). Shapes: 0: sine, 1: saw, 2: square, 
3: triangle. */
void createWaveform(double *x, int N, int shape, double frequency, double sampleRate, 
  double phase = 0.0, bool antiAlias = false);

void createPulseWave(double *x, int N, double frequency, double dutyCycle, double sampleRate, 
  double phase = 0.0, bool antiAlias = false);





// Plotting convenience functions:

// convenience functions for interfacing with the Plotter class (they manage instantiation and 
// setup of Plotter objects for frequently used cases):
void plotData(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);

void plotData(int N, double x0, double dx, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);
  // for equidistant data, abscissa values start at x0 with increment dx

void plotDataLogX(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL, 
  double *y4 = NULL, double *y5 = NULL);


void plotVector(std::vector<double> v);

/** Plots the magnitude spectrogram given in s against time axis t (of length numFrames) and
frequency axis f (of length numBins). */
void plotSpectrogram(int numFrames, int numBins, double **magnitudes, double sampleRate, 
  int hopSize);
  // introduce parameters to control scaling of time- and frequency axis..


#endif
