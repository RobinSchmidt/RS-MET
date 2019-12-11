#ifndef RS_MODALEXPERIMENTS_H
#define RS_MODALEXPERIMENTS_H


// maybe move to FilterExperiments

void modalFilter();
void modalFilterFreqResp(); // rename to modalWithAttack
void modalTwoModes();
void attackDecayFilter();   // remove (redundant now)
void dampedSineFilterDesign();
void dampedSineFilterImpResp();
void biquadImpulseResponseDesign();
void modalBankTransient(); // tests the use of nonlinear feedback to produce transients

void fourExponentials();
void modalWithFancyEnv();

void modalSynthSpectra();

void modalDecayFit();
void modalAnalysis1();
void modalAnalysisPluck();



#endif
