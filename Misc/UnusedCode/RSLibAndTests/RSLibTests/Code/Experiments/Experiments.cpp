#include "Misc/PlotterExperiments.h"
#include "Misc/CellularExperiments.h"

#include "RSMath/MathExperiments.h"

#include "RSAudio/AnalysisExperiments.h"
#include "RSAudio/FilterExperiments.h"
#include "RSAudio/DelayExperiments.h"
#include "RSAudio/MiscAudioExperiments.h"
#include "RSAudio/ModulatorExperiments.h"
#include "RSAudio/ResamplingExperiments.h"
#include "RSAudio/PhaseVocoderExperiments.h"
#include "RSAudio/ModalExperiments.h"
#include "RSAudio/OscillatorExperiments.h"
#include "RSAudio/PartialExtractionExperiments.h"
#include "RSAudio/PhysicsExperiments.h"
#include "RSAudio/SaturationExperiments.h"

int main(int argc, char** argv)
{
  // todo - sort these function call by category and in each category alphabetically

  //-----------------------------------------------------------------------------------------------
  // Plotter class (remove these - these are now part of the GNUPltCPP project):

  //testDataPlot2();
  //testIntVector();
  //testStyleSetup2();
  //testMultiColumn2();
  //testSurface2(); 

  // old (remove when new plotter is ready):
  //testVariableArgumentList();
  //testDataPlot1();
  //testColorSetup();
  //testStyleSetup();
  //testMultiColumn();
  //testSurface(); // does not work properly yet - we need partioned datasets



  // Cellular Automata:
  //testCellularAutomaton1();



  //-----------------------------------------------------------------------------------------------
  // RSMath:

  //binomialDistribution();
  //sineParameters();
  //bandLimitedStep();
  //cubicInterpolationNonEquidistant(); // -> turn into unit-test!
  //hyperbolicFunctions(); 
  //splineInterpolationNonEquidistant();
  //rationalInterpolation();
  //splineInterpolationAreaNormalized();
  //numericDerivative();
  //shiftPolynomial();
  //monotonicPolynomials();
  //parametricBell();
  //partialFractionExpansion();
  //partialFractionExpansion2();
  //partialFractionExpansionQuadratic();
  //dampedSineEnergy(); //->relevant for modal synthesis
  //sineIntegral();
  //logarithmQuotient();
  //stirlingNumbers();
  //bernoulliNumbers();
  //sequenceSquareRoot();
  //conicSystem();  
  //logisticMapNoise();
  //divide();
  //nestedComplexOrder2();
  //nestedComplexOrder3();
  //bigFloatErrors();


  //primeRecursion();
  //primeSieveSchmidt1();
  //primeSieveSchmidt2();
  //primeSieveAtkin();
  //primeSieve();
  //primeDistribution();
  //numberTheoreticTransform();

  //-----------------------------------------------------------------------------------------------
  // RSAudio:

  //autoCorrelation();
  //autocorrelationPeakVariation();
  //autoCorrelationPitchDetector();
  //autoCorrelationPitchDetectorOffline();
  //crossCorrelationBestMatch();
  //combineFFTs();                            // move to RSMath experiments
  cycleMarkFinder();
  //instantaneousFrequency();
  //instantaneousPhase();
  //zeroCrossingFinder();
  //zeroCrossingFinder2();
  //zeroCrossingPitchDetector();
  //zeroCrossingPitchDetectorTwoTones();

  //basicIntegerDelayLine();

  // filters:
  //bandwidthScaling();
  //stateVariableFilter();
  //stateVariableFilterMorph();
  //transistorLadder();
  //phonoFilterPrototypePlot();
  //magnitudeMatchedOnePoleFilter();
  //phonoFilterModelPlot();
  //phonoFilterSimulation();
  //serialParallelBlend();
  //averager();
  //movingAverage();
  //trapezAverager();
  //compareApproximationMethods();
  //ringingTime();
  //butterworthSquaredLowHighSum();
  //gaussianPrototype();
  //halpernPrototype();
  //maxFlatMaxSteepPrototypeM1N2();
  //maxFlatMaxSteepPrototypeM2N2();
  //splitLowFreqFromDC();

  //ladderResonanceModeling(); 
  //ladderResoShape();
  //ladderThresholds();
  //ladderFeedbackSaturation();
  //ladderFeedbackSaturation2();
  //ladderFeedbackSaturation3();
  //ladderFeedbackSatDCGain();
  //ladderFeedbackSatReso();
  //ladderFeedbackSatGrowl();  
  //ladderFeedbackSatGrowl2(); 
  //ladderZDF(); 
  //ladderZDFvsUDF();
  
  //resoSaturationModes(); 
  //resoShapeFeedbackSat();
  //resoShapeGate();
  //resoShapePseudoSync();
  //resoSeparationNonlinear();

  //resoReplace();  
  //resoReplacePhaseBumping();   
  //resoReplaceScream();

  //fakeResonance(); 
  //fakeResoLowpassResponse();
  //fakeResoDifferentDelays();



  // misc:
  //centroid();
  //cubicCrossfade();
  //slewRateLimiterLinear();
  //recursiveSineSweep();
  //ringModNoise();

  //stretchedCorrelation();
  //taperedFourierSeries();
  //transientModeling();
  //windowFunctionsContinuous();
  //windowedSinc();


  //breakpointModulator();
  //breakpointModulatorSmoothFadeOut();

  // resampling:
  //fadeOut(); // move
  //resampler();
  //sincResamplerAliasing();
  //sincResamplerModulation();
  //sincResamplerPassbandRipple();
  //sincResamplerSumOfTapWeights();
  //timeWarp();
  //pitchDemodulation();

  //phaseLockedCrossfade();
  //phaseLockedCrossfade2(); 

  //sineShift();
  //sineShift2();
  //pitchDetectWithSilence();

  // phase-vocoder:
  //phaseRepresentation();
  //grainRoundTrip();
  //plotWindows();
  //spectrogramSine();

  // tests with Elan's example files:
  //pitchDetectA3();
  //phaseLockSaxophone();
  //phaseLockSaxophone2();
  //autoTuneHorn();
  //autoTuneHorn2();
  //sylophoneCycleMarks();
  //autoTuneSylophone();
  //bestMatchShift();

  //modalFilter();
  //attackDecayFilter();
  //dampedSineFilterDesign();
  //biquadImpulseResponseDesign();

  //biDirectionalFilter();
  //sineRecreation();
  //sineWithPhaseCatchUp();
  //partialExtractionTriple();
  //partialExtractionBell();
  //partialExtractionViaBiquadTriple();
  //partialExtractionSample();

  // oscillator:
  //triSaw(); 
  //phaseShapingCurvePoly4(); 
  //phaseShapingCurvesRational();
  //phaseShaping();
  //phaseShapingSkew();

  // physical modeling
  //doublePendulum();

  // saturation:
  //powRatioParametricSigmoid();
  //cubicRatioParametricSigmoid();  // does not exist anymore?
  //quinticParametricSigmoid();  // these two...
  //septicParametricSigmoid();   // ...do not yet work as intended
  //parametricSigmoid();
  //parametricSigmoid2(); 
  //saturator();
  //sigmoidPrototypes();
  //sigmoidScaleAndShift();
  //quarticMonotonic();
  //sixticPositive();


  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";
  getchar();
  return(EXIT_SUCCESS);
}
