#ifndef RS_PHASEVOCODEREXPERIMENTS_H
#define RS_PHASEVOCODEREXPERIMENTS_H


void phaseRepresentation();

void grainRoundTrip();        // under construction
void plotOverlappingWindowSum();
void spectrogramSine();
void spectrogramFilter();

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

//void harmonicBeatingDemo();   // 
void amplitudeDeBeating();
void amplitudeDeBeating2();
void harmonicDeBeating1();
void harmonicDeBeating2();     // remove beating from a partial via harmonic analysis/resynthesis framework





// resynthesis:

// synthesis from generated PV data:


#endif
