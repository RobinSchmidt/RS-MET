#ifndef RS_MODALEXPERIMENTS_H
#define RS_MODALEXPERIMENTS_H

#include "../../../Shared/Shared.h"

// maybe move to FilterExperiments

void modalFilter();
void attackDecayFilter();
void dampedSineFilterDesign();
void biquadImpulseResponseDesign();

void modalBankTransient(); // tests the use of nonlinear feedback to produce transients

#endif
