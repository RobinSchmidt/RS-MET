#ifndef RS_MODALEXAMPLES_H
#define RS_MODALEXAMPLES_H

#include <iostream>
#include "SampleMapGenerator.h"

void createModalFilterExamples();
void createModalFilterBankExamples();

// Single sample creations:
void createPiano1();

// Multisample creations
void createSamplerWaveforms();
void createBassdrumPsy1Samples();

// Multisample + sfz sample-map creations:
void createBass1();
void createGong1();
void createPluck1();

// Experimental:
void testHighPluck();
// todo: make the sort sounds available (at the bottom of the .cpp file)



/** A class to generate the low-level parameters for each mode from a bunch of higher level macro
parameters. */

template<class T>
class rsModalParameterGenerator
{

public:


  rsModalParameterGenerator();

  rsModalBankParameters<T> getModalParameters();



  /** \name Frequency parameters */

  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; }

  void setFrequency(T newFrequency) { frequency = newFrequency; }

  void setInharmonicity(T newInharmonicity) { inharmonicity = newInharmonicity; }

  // void setFrequencyRatios
  // harmonic/12-tet/idealBar/etc...




  /** \name Amplitude parameters */

  void setAmplitude(T newAmplitude) { amplitude = newAmplitude; }

  void setAmpSlope1(T newSlope) { ampSlope1 = newSlope; }

  void setAmpCutoff(T newCutoff) { ampCutoff = newCutoff; }

  void setAmpSlope2(T newSlope) { ampSlope2 = newSlope; }

  void setEvenAmpScale(T newScale) { evenAmpScale = newScale; }

  void setAmpCombHarmonic(T newHarmonic) { ampCombHarmonic = newHarmonic; }

  void setAmpCombAmount(T newAmount) { ampCombAmount = newAmount; }


  /** \name Phase parameters */

  void setPhaseRandomness(T newRandomness) { phaseRandomness = newRandomness; }

  void setPhaseRandomSeed(int newSeed) { phaseRandomSeed = newSeed; }

  // void setStartPhase(T newPhase);
  // void setPhaseShiftByFreq(T newShift)
  // shift = (f-1)*shiftByFreq * 2*PI; // shift of individual harmonic
  // and/or shift = i*shiftByIndex * 2*PI;



  /** \name Decay parameters */

  void setDecay(T newDecay) { decayTime = newDecay; }

  void setDecaySlope1(T newSlope) { decaySlope1 = newSlope; }

  void setDecayCutoff(T newCutoff) { decayCutoff = newCutoff; }

  void setDecaySlope2(T newSlope) { decaySlope2 = newSlope; }

  void setEvenDecayScale(T newScale) { evenDecayScale = newScale; }

  void setDecayCombHarmonic(T newHarmonic) { decayCombHarmonic = newHarmonic; }

  void setDecayCombAmount(T newAmount) { decayCombAmount = newAmount; }


  /** \name Attack parameters */

  void setAttack(T newAttack) { attackTime = newAttack; }



  /** \name Static functions */

  /** Gives the relative mode decay time for mode with relative frequency f given a (relative)
  cutoff frequency fc and an ultimate slope of the decay-time function (with respect to f)
  determined by the exponent p
  d(f) = a / (b + (f/fc)^p) where a and b are adjusted such that d(f=1)=1 and d(f=fc)=1/2
  fc must be > 1 */
  static T modeDecayTime(T f, T fc, T p);
  // rename and use also for amplitude

  static T combAmplitude(T frequency, T notchDistance, T amount = 1, T notchOffset = 0,
    T shape = 1);


protected:

  void getFrequencies(std::vector<T>& f); 
  // maybe rename to getRelativeFrequencies or getFrequencyRatios, also assigns numPartials so it 
  // should be always called before getPhases

  void getPhases(std::vector<T>& p, const std::vector<T>& f);

  void getAmplitudes(std::vector<T>& a, const std::vector<T>& f);


  void getDecayTimes(std::vector<T>& d, const std::vector<T>& f);

  void getAttackTimes(std::vector<T>& a, const std::vector<T>& f, const std::vector<T>& d);












  // data:

  // frequency related:
  T sampleRate    = 44100;
  T frequency     = 100;
  T inharmonicity = 0.0;
  //int partialRatios  = HARMONIC;

  // amplitude related:
  T amplitude       = 1.0;
  T ampSlope1       = 1.0;
  T ampCutoff       = 2.0;    // as harmonic number, must be > 1
  T ampSlope2       = 0.0;    // additional slope above cutoff as direct exponent/power
  T evenAmpScale    = 1.0;    // amplitude scaler for even harmonics
  T ampCombHarmonic = 7.0;    // harmonic number of 1st notch
  T ampCombAmount   = 0.0;

  // phase related:
  int phaseRandomSeed  = 0;
  //int phaseRandomShape = 1;  // shape of the distribution: 1: uniform, 2: triangular, 3: parabolic
  T phaseRandomness    = 0;

  // envelope related:
  T decayTime         = 1.0;
  T decaySlope1       = 1.0;
  T decayCutoff       = 2.0;
  T decaySlope2       = 0.0;
  T evenDecayScale    = 1.0;    // decay-time scaler for even harmonics
  T decayCombHarmonic = 7.0;
  T decayCombAmount   = 0.0;

  T attackTime        = 1.0;
  T attackDecayRatioLimit = 0.99;

  int maxNumPartials = 1024;
  //int numPartials    = 0;

  std::vector<T> tmp; // for temporary values
  RAPT::rsNoiseGenerator<T> prng;
};


// move somewhere else:
void createInsertionSortSound();


#endif
