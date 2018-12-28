#ifndef RS_MODALEXAMPLES_H
#define RS_MODALEXAMPLES_H

#include <iostream>
#include "SampleMapGenerator.h"

void createModalFilterExamples();
void createModalFilterBankExamples();

// single sample creations:
void createPiano1();

// sample-map creations:
void createBass1();
void createGong1();
void createPluck1();


/** A class to generate the low-level parameters for each mode from a bunch of higher level macro
parameters. */

template<class T>
class rsModalParameterGenerator
{

public:

  rsModalBankParameters getModalParameters();
    // rsModalBankParameters should be templatized, too


  // data:

  T sampleRate    = 44100;

  T frequency     = 100;
  T amplitude     = 1.0;
  T attackTime    = 0.1;
  T decayTime     = 1.0;

  T lowpassSlope  = 0.0;    // in dB/oct
  T lowpassCutoff = 1;      // as harmonic number

  T inharmonicity = 0.0;
  T evenAmplitude = 1.0;    // amplitude scaler for even harmonics
  T evenDecay     = 1.0;    // decay-time scaler for even harmonics

  T combHarmonic  = 7.0;
  T combAmount    = 1.0;

  int maxNumPartials = 50;
};


// move somewhere else:
void createInsertionSortSound();


#endif
