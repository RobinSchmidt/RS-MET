#ifndef RS_MODALEXPERIMENTS_H
#define RS_MODALEXPERIMENTS_H


// maybe move to FilterExperiments

void twoPoleFilter();
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

// stuff for modal analysis/resynthesis:
void modalDecayFit();
void modalAnalysis1();
void modalAnalysisPluck();
void modalPartialResynthesis();



#endif
