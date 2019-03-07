#pragma once


/** A class to create a sinusoidal model for quasi-harmonic signals. The assumed harmonic 
relationship between the partials allows for a different analysis algorithm that is tailored 
specifically for harmonic signals and is supposed to give better results than a more general
sinusoidal analysis algorithm that doesn't make such an assumption (if the assumption indeed holds
true for the analyzed signal, of course). The algorithm works as follows:

-pre-process audio (pitch flattening):
 -obtain cycle marks
 -obtain a time warping map that makes the cycle-marks equidistant (distance is chosen to be a 
  power of two >= the longest cycle in the input signal)
 -time-warp the signal -> all cycles have the same, FFT-friendly length
-harmonic analysis:
 -analyze each cycle of pitch flattened signal by an FFT
 -the FFT size is equal to the cycle-length -> frequencies don't need to be estimated, they are 
  known in advance
 -only amplitude and phase have to be measured
-post-process model data (account for pitch flattening):
 -move time instants of datapoints according to the inverse time-warping map
 -modify frequencies according to reciprocal of the slope of time-warping map

The so obtained model models inharmonicity, transient and noise in the input signal as fast 
variations of instantaneous frequencies and amplitudes of the partials - when the envelopes are
lowpass-filtered before resynthesis, we can resynthesize only the quasi-harmonic part of the sound.
Time domain subtraction from the original should give a residual containing inharmonic, noisy and 
transient parts. */

template<class T>
class rsHarmonicAnalyzer
{

public:



  rsHarmonicAnalyzer();



  /** \name Setup */


  /** Sets up the sample-rate. This determines the values of the frequencies that will be written
  into the model. */
  void setSampleRate(T newRate) { sampleRate = newRate; }


  /** \name Processing */

  /** Analyzes the given input sound and returns the sinusoidal model object that models the given 
  sound. */
  RAPT::rsSinusoidalModel<T> analyze(T* sampleData, int numSamples);


protected:

  /** Computes and returns an array of cycle-marks. */
  std::vector<T> findCycleMarks(T* x, int N);

  /** Used internally to fill in the data in the model at the given frame-index based on the 
  current content of our "sig" buffer member variable. The time-stamp of the frame should be passed
  by the caller - frequency, magnitude and phase data are computed from the FFT of the sig 
  buffer. */
  void fillHarmonicData(RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T timeStamp);


  //RAPT::rsPitchFlattener<T, T> flattener;
  RAPT::rsFourierTransformerRadix2<T> trafo;



  T sampleRate = 1;



  std::vector<T> sig, mag, phs;  // buffers for blocks of signal, magnitude and phase

};