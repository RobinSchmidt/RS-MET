#ifndef RS_FILTEREXPERIMENTS_H
#define RS_FILTEREXPERIMENTS_H

void bandwidthScaling();
void biquadResoGain();
void butterworthEnergy();
void biDirectionalStateInit();
void biquadDesignVicanek();
void biquadTail();
void biquadModulation();
void stateVariableFilter();
void stateVariableFilterMorph();
void stateVectorFilter();
void transistorLadder();
void phonoFilterPrototypePlot();
void magnitudeMatchedOnePoleFilter();
void phonoFilterModelPlot();
void phonoFilterSimulation();
void serialParallelBlend();
void averager();
void movingAverage();
void trapezAverager();
void compareApproximationMethods();
void compareOldAndNewEngineersFilter();
void testPoleZeroMapper();
void ringingTime();
void butterworthSquaredLowHighSum();
void gaussianPrototype();
void halpernPrototype();

void maxFlatMaxSteepPrototypeM1N2();
void maxFlatMaxSteepPrototypeM2N2();
//void experimentalPrototypeM1N2();

void splitLowFreqFromDC();

void directFormFreqResp();
void ladderResonanceGain();
void ladderTransferFunction();
void ladderMultipole();
void ladderResonanceModeling();
void ladderResoShape();
void ladderThresholds();           // maybe remove - this seemed to be a dead end

void ladderFeedbackSaturation();
void ladderFeedbackSaturation2();
void ladderFeedbackSaturation3();
void ladderFeedbackSatDCGain();
void ladderFeedbackSatReso();
void ladderFeedbackSatGrowl();
void ladderFeedbackSatGrowl2();

void ladderZDF();
void ladderZDFvsUDF();

void ladderResoModulation();

void resoShapeFeedbackSat();
void resoSaturationModes();
void resoShapeGate();
void resoShapePseudoSync();
void resoSeparationNonlinear();

void resoReplace();
void resoReplacePhaseBumping();
void resoReplaceScream();

void resoWave();

void fakeResonance();
void fakeResoLowpassResponse();
void fakeResoDifferentDelays();

void quantileFilter();

template<class T, int N> void simdFilter();


#endif
