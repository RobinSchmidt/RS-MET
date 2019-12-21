/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               rapt
  vendor:           RS-MET
  version:          0.0.1
  name:             Rob's Audio Processing Templates
  description:      Library of audio DSP algorithms as C++ templates
  website:          http://www.rs-met.com
  license:          Custom

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


#ifndef RAPT_H_INCLUDED
#define RAPT_H_INCLUDED


// maybe these standard library includes should go somewhere else?
//#include <ctgmath>       // template generic math?
#include <cstdio>        // to fix warning in gcc when using printf
#include <complex>       // included already by ctgmath ...but not on gcc/windows?
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
//#include <string>
#include <cstring>
#include <climits>
//#include <xmmintrin.h>
//#include <emmintrin.h>    // SSE2 - should be included by rosic
//#include <iostream>       // only for printing debug info


#include "Basics/Plotting.h"             // included by itself (not via Basics.h) to avoid wrapping into namespace
#include "Basics/Basics.h"               // type definitions, constants, functions, etc.
#include "Data/Data.h"                   // data structures like arrays, lists, etc.
#include "Math/Math.h"                   // interpolation, transforms, linear algebra, numerical analysis, ...
#include "AudioBasics/AudioBasics.h"
//#include "Music/Music.h"               // scales (pitchToFreq, ..), time signatures, sequences, patterns, notes, MIDI, etc.
#include "Filters/Filters.h"             // butterworth, ladder, biquad, elliptic, SVF, etc.
#include "Analysis/Analysis.h"           // envelope follower, smoother, pitch-detector, etc.
#include "Physics/Physics.h"             // waveguides, finite difference approximations of PDEs, etc.
//#include "Circuits/Circuits.h"         // circuit modeling
//#include "Spectral/Spectral.h"         // phase vocoder, source/filter modeling
#include "Visualization/Visualization.h" // buffers for scopes, waveform displays etc
#include "Generators/Generators.h"     // oscillator, sample player, etc. - maybe rename to sources
#include "Modulators/Modulators.h"       // ADSR, LFO, step-sequencer, breakpoint-modulator, etc.
//#include "Effects/Effects.h"           // reverb, distortion, dynamics, chorus, etc.
//#include "Framework/Framework.h"       // parameter handling, save/recall, threading, polyphony, etc.
//#include "Instruments/Instruments.h"   // full blown instruments with polyphony, state-recall, etc.

#include "Unfinished/Unfinished.h"       // code under construction - should go to the bottom!
#include "Spectral/Spectral.h"           // move up, when cycle-mark finder was moved to analysis folder

namespace RAPT
{
#include "_Deprecated/MultiArray.h"      // remove when unit tests have been updated to use new
}





// ...the ordering above should roughly reflect the dependencies (a module later in the chain may
// depend on one or more modules that come before it but not vice versa - well, we'll see if that's
// possible)

/*
// typedefs for convenience:
typedef RAPT::rsSinCosTable<double> rsSinCosTableD;

typedef RAPT::rsSmoothingFilter<double, double> rsSmoothingFilterDD;
typedef RAPT::rsLadderFilter<double, double> rsLadderDD;
typedef RAPT::rsLadderFilter<rsFloat64x2, double> rsLadderD2D;

typedef RAPT::rsRayBouncer<double> rsRayBouncerD;
typedef RAPT::rsRayBouncerDriver<double> rsRayBouncerDriverD;
*/


#endif
