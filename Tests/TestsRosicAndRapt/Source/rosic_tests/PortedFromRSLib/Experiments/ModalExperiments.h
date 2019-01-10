#ifndef RS_MODALEXPERIMENTS_H
#define RS_MODALEXPERIMENTS_H

#include "../../../Shared/Shared.h"

// maybe move to FilterExperiments

void modalFilter();
void modalFilterFreqResp(); // rename to modalWithAttack
void modalTwoModes();
void attackDecayFilter();   // remove (redundant now)
void dampedSineFilterDesign();
void biquadImpulseResponseDesign();
void modalBankTransient(); // tests the use of nonlinear feedback to produce transients

void fourExponentials();
void modalWithFancyEnv();

void modalSynthSpectra();

#endif
