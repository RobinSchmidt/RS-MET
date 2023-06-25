#pragma once

#include "UnitTests/UnitTestsRapt.h"

//-------------------------------------------------------------------------------------------------
// Experiments

// Filters Experiments:

void bandpassAndNotch();
void bandSplittingTwoWay();
void bandSplittingThreeWay();
void bandSplittingThreeWay2p2z();
void bandSplittingMultiWay();
void bandSplittingTreeAlgo();
void bandSplitFreqResponses();
void complementaryFiltersIIR();
void firstOrderFilters();
void engineersFilterRingResp();
void engineersFilterFreqResps();
void ladderResonanceManipulation();
// insert movingAverage here
void nonUniformMovingAverage();
void nonUniformOnePole1();
void nonUniformOnePole2();
void nonUniformComplexOnePole();
void nonUniformAllpole();
void nonUniformBiquad();
void nonUniformBiDirectional();
void smoothingFilterOrders();
void smoothingFilterTransitionTimes();
void prototypeDesign();    // still using the old code
void poleZeroPrototype();  // new implementation
void seriesConnectionDecay();


// Generator Experiments:

void waveformFractalization();
void noise();
void noiseReverseMode();
void noiseTriModal();
void noiseWaveShaped();
void blit();
void blep();
//void blamp();
void polyBlep();
void superBlep();
void superSawDensitySweep();
void superSawStereo();
void twoPieceOsc();
void syncSweep();
void syncPhasor();
void syncPhasor2();
void syncOsc();
void dualBlepOsc();
void bouncillator();
void bouncillatorFormula();
void freqVsPhaseMod();
void rayBouncer();
void circleFractals();
void hilbertCurve();
void lindenmayer();
void snowFlake();
void triSawOsc();
void triSawOscAntiAlias();
void xoxosOsc();
void shepardTone();
void additiveEngine();
void multiplicativeSynth();
void puleWidthModulationViaTwoSaws();


// Graphics Experiments:

void colorGradientHSL();
void lineDrawing();
void lineDrawingThick();
//void lineDrawingThick2(); // obsolete
void lineJoints();
void lineTo();
void polyLineRandom();
void phaseScopeLissajous();
void splineArc();
void triangles();
void pixelCoverage();
void gradientify();
void contours();
void complexContours();
void implicitCurves();
void parametricCurve();
void spirals();   // move elsewhere Demos/Examples/Art/...
void fractal();
void differentialGeometry();
void imageScaling();

// Math Experiments:
// wait - what about the file MathExperiments.h - shouldn't they be declared there?

void ellipseLineIntersections();
void interpolatingFunction();
void linearRegression();
void multipleRegression();
void polynomialRegression();
void gaussianRegression();
void butterworthViaGaussians();

void numericOptimization();

void polynomialSinc();
void polynomialPrediction();
void probabilityLogic();
void productLogPlot();
void ratioGenerator();
void ratiosLargeLcm();
void ratiosEquidistantPowers();
void ratiosMetallic();
void sinCosTable();
void expBipolar();
void expGaussBell();
void twoParamRemap();
void numberTheoreticTrafo();


// move to something like FunExperiments (maybe copy the mandelbrot stuff over there, too)
bool groupString();
void primeAlternatingSums();
void divisibility();
void arithmeticDerivative();




// Modulator Experiments:

void attackDecayEnvelope();

//-------------------------------------------------------------------------------------------------
// Performance Tests

// Audio Performance Tests:
void filterSignConventionPerformance();
void ladderPerformance();
void stateVectorFilterPerformance();
void engineersFilterPerformance();
void turtleGraphicsPerformance();
void samplerEnginePerformance();

// Math Performance Tests:
void fftPerformance();
void matrixAdressingTest();
void simdPerformance(); 
void sinCosPerformance();

// Misc Performance Tests:
void callbackPerformance();

//-------------------------------------------------------------------------------------------------
// Unit Tests
// ...they actually have their own include file

// Data Unit Tests

