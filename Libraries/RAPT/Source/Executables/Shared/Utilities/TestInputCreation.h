#ifndef RAPT_TESTINPUTCREATION_H
#define RAPT_TESTINPUTCREATION_H

#include "../../../Modules/RAPT.h"
using namespace RAPT;

/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate);

/** Synthesizes a standard waveform at the desired frequency and samplerate. The phase is 
expected in radians and there is an option to do anti-aliased synthesis (by means of adding
the sinusoidal components up to the Nyquist frequency). Shapes: 0: sine, 1: saw, 2: square, 
3: triangle. */
void createWaveform(float *x, int N, int shape, float frequency, float sampleRate, 
  float phase = 0.0, bool antiAlias = false);

// reactivate and adadpt later:
//void createPulseWave(double *x, int N, double frequency, double dutyCycle, double sampleRate, 
//  double phase = 0.0, bool antiAlias = false);
//
///** Creates a sine-wave of length N (in samples) with frequency f at and amplitude a, 
//at samplerate fs. */
//void createSineWave(double *x, int N, double f, double a, double fs);
//
///** Creates a sine-wave of length N (in samples) with an instantaneous frequency given by f at each 
//sample and an amplitude given by a, at samplerate fs. */
//void createSineWave(double *x, int N, double *f, double a, double fs);

//
///** Returns the value of a sum of N harmonic sines with amplitudes given by A and instantaneous 
//phases given by n*p, where n is the harmonic number of the sine (starting at 1) */
//double sineSum(double p, double *A, double N);




#endif