#ifndef RS_PHASEVOCODEREXPERIMENTS_H
#define RS_PHASEVOCODEREXPERIMENTS_H

#include "../../../Shared/Shared.h"

void phaseRepresentation();

void grainRoundTrip();        // under construction
void plotWindows();
void spectrogramSine();

void sineParameterEstimation();


void phaseInterpolation();

void sinusoidalSynthesis1();
void sinusoidalSynthesis2();  // tests partial with negative time-stamps

void sinusoidalAnalysis1();   // a single, static sinusoid
void sinusoidalAnalysis2();   // two single, static sinusoids
void sinusoidalAnalysis3();   // a linear frequency sweep

void phaseFreqConsistency();
void harmonicDetection2Sines();
void harmonicDetection3Sines();
void harmonicDetection5Sines();
void harmonicAnalysis1();





// resynthesis:

// synthesis from generated PV data:


#endif
