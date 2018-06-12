#ifndef RS_FILTEREXPERIMENTS_H
#define RS_FILTEREXPERIMENTS_H

#include "../../../Shared/Shared.h"

void bandwidthScaling();
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

void resoShapeFeedbackSat();
void resoSaturationModes();
void resoShapeGate();
void resoShapePseudoSync();
void resoSeparationNonlinear();

void resoReplace();
void resoReplacePhaseBumping();
void resoReplaceScream();

void fakeResonance();
void fakeResoLowpassResponse();
void fakeResoDifferentDelays();


#endif
