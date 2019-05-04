#ifndef RS_SAMPLINGEXPERIMENTS_H
#define RS_SAMPLINGEXPERIMENTS_H


void fadeOut();  // move to a new file SampleEditingExperiments

void resampler();
void sincResamplerAliasing();
void sincResamplerModulation();
void sincResamplerPassbandRipple();
void sincResamplerSumOfTapWeights();
void timeWarp();
void pitchDemodulation();

void phaseLockedCrossfade();
void phaseLockedCrossfade2();

void sineShift();
void sineShift2();
void pitchDetectWithSilence();

// tests with Elan's example files:
void pitchDetectA3();
void phaseLockSaxophone();
void phaseLockSaxophone2();
void autoTuneHorn();
void autoTuneHorn2();
void sylophoneCycleMarks();
void autoTuneSylophone();
void bestMatchShift();



#endif
