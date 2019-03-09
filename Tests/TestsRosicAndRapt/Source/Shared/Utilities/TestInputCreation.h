#ifndef RAPT_TESTINPUTCREATION_H
#define RAPT_TESTINPUTCREATION_H

#include "../RaptLibraryCode/RaptInstantiations.h"

/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate);

void createTimeAxis(int numSamples, double *timeAxis, double sampleRate);

//template<class T>
inline std::vector<double> createSilence(int numSamples)
{
  std::vector<double> v;
  v.resize(numSamples);
  RAPT::rsArray::fillWithZeros(&v[0], numSamples);
  return v;
}

//template<class T>
inline std::vector<double> createImpulse(int numSamples)
{
  std::vector<double> v = createSilence(numSamples);
  v[0] = 1;
  return v;
}


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
void createSineWave(double *x, int N, double f, double a, double fs, double p = 0);

std::vector<double> createSineWave(int N, double f, double fs, double a = 1, double p = 0); // convenience function



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

void createRandomDataXY(double* x, double* y, int N, double dxMin, double dxMax, 
  double yMin, double yMax, int seedX = 0, int seedY = 1);

/** Combination of a steady and a decaying sine wave, the 2nd having a frequency-ratio of 
1+sqrt(5) = 2*goldenRatio to the first. This makes the second partial's ratio close to 3 but very 
inharmonic (the golden ratio is in some sense the most irrational number). The decay is given as 
1/tau where tau is the decay time in seconds. So, a decay parameter of 0 makes the 2nd partial 
steady, too. */
std::vector<double> sineAndDeacyingInharmonic(int N, double f, double fs, double decay = 0);

std::vector<double> twoSinesAndDecayingDc(int N, double f, double fs, double overtoneRatio, 
  double overtoneAmplitude, double dcAmount, double dcDecay);

std::vector<double> sawAndSquare(int N, double fs, double fSaw, double aSaw, 
  double fSqr, double aSqr, bool antiAlias);


double saw(double phi, int maxHarmonic);

double sqr(double phi, int maxHarmonic);

/** Creates output of an sawtooth wave at a given sample index n with given frequency and 
sample-rate. You may optionally pass a maximum harmonic number to generate, if you pass nothing,
it will be automatically selected to be the maximum harmonic possible without introducing 
aliasing. */
double saw(int n, double f, double fs, int maxHarmonic = -1);

/** @see saw - same stuff for a square wave. */
double sqr(int n, double f, double fs, int maxHarmonic = -1);


void createModalPluck(double* x, int N, double key, double sampleRate);

std::vector<double> createModalPluck(double key, double sampleRate, int length); // convenience function



#endif