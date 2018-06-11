#ifndef RAPT_TESTINPUTCREATION_H
#define RAPT_TESTINPUTCREATION_H

#include "../RaptLibraryCode/RaptInstantiations.h"

/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate);

void createTimeAxis(int numSamples, double *timeAxis, double sampleRate);


/** Synthesizes a standard waveform at the desired frequency and samplerate. The phase is 
expected in radians and there is an option to do anti-aliased synthesis (by means of adding
the sinusoidal components up to the Nyquist frequency). Shapes: 0: sine, 1: saw, 2: square, 
3: triangle. Explicit template instantiations exist for float and double. */
template<class T>
void createWaveform(T *x, int N, int shape, T frequency, T sampleRate, T phase = 0.0, 
  bool antiAlias = false);




void createPulseWave(double *x, int N, double frequency, double dutyCycle, double sampleRate, 
  double phase = 0.0, bool antiAlias = false);

/** Creates a sine-wave of length N (in samples) with frequency f at and amplitude a, 
at samplerate fs. */
void createSineWave(double *x, int N, double f, double a, double fs);

/** Creates a sine-wave of length N (in samples) with an instantaneous frequency given by f at each 
sample and an amplitude given by a, at samplerate fs. */
void createSineWave(double *x, int N, double *f, double a, double fs);

/** Creates a sum of sine waves with given frequencies, amplitudes and phases (phases will be 
zero, if a nullptr is passed); */
void createSumOfSines(double* x, int numSamples, int numSines, double fs, double *f, double *a, 
  double *p = nullptr);

/** Returns the value of a sum of N harmonic sines with amplitudes given by A and instantaneous 
phases given by n*p, where n is the harmonic number of the sine (starting at 1) */
double sineSum(double p, double *A, double N);



/** Returns a vector of random samples uniformly distributed between min and max. */
template<class T>
std::vector<T> createNoise(int numSamples, T min, T max, int seed = 0);




#endif