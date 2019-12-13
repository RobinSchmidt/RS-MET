#ifndef RS_ANALYSISEXPERIMENTS_H
#define RS_ANALYSISEXPERIMENTS_H

void autoCorrelation();
void autocorrelationPeakVariation();

void autoCorrelationPitchDetector();
void autoCorrelationPitchDetectorOffline();

void crossCorrelationBestMatch();

void combineFFTs(); // move to RSMath experiments


void envelopeFollower();

//void zeroCrossingPitchDetector();

void instantaneousFrequency();
void instantaneousPhase();

void maxShortTimeRMS();
void arrayRMS();

void zeroCrossingFinder();
void zeroCrossingFinder2();
void zeroCrossingFinder3();
void cycleMarkFinder();
void cycleMarkErrors();
void zeroCrossingPitchDetector();
void zeroCrossingPitchDetectorTwoTones();

void peakPicker();





#endif
