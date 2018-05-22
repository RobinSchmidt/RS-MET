#include "../JuceLibraryCode/JuceHeader.h"

#include "rosic_tests/rosic_CorrectnessTests.h"
using namespace rotes; // get rid of that

#include "rapt_tests/Experiments/MathExperiments.h"
#include "rapt_tests/Experiments/FilterExperiments.h"
#include "rapt_tests/Experiments/GeneratorExperiments.h"
#include "rapt_tests/Experiments/GraphicsExperiments.h"

#include "rapt_tests/PerformanceTests/MathPerformanceTests.h"
#include "rapt_tests/PerformanceTests/AudioPerformanceTests.h"
#include "rapt_tests/PerformanceTests/MiscPerformanceTests.h"

#include "rapt_tests/UnitTests/UnitTests.h"

#include "Experiments/Experiments.h"

// the new stuff proted from RSLib (todo: merge directories and files where appropriate, bring 
// everything into a consistent order):
#include "rosic_tests/PortedFromRSLib/ExamplesRSLib.h"
#include "rosic_tests/PortedFromRSLib/ExperimentsRSLib.h"
#include "rosic_tests/PortedFromRSLib/PerformanceTestsRSLib.h"
#include "rosic_tests/PortedFromRSLib/UnitTestsRSLib.h"


int main(int argc, char* argv[])
{
  //===============================================================================================
  // RAPT tests:

  //-----------------------------------------------------------------------------------------------
  // Unit tests:

  //runAllUnitTests();  // merge with eunit testsf for RSLib
  //mathUnitTests();
  //filterUnitTests();


  //-----------------------------------------------------------------------------------------------
  // Performance tests:

  //callbackPerformance();
  //matrixAdressingTest();
  //simdPerformanceFloat64x2();
  //sinCosPerformance();

  //ladderPerformance();
  //engineersFilterPerformance();
  //turtleGraphicsPerformance();

  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Prototypes (not yet in library):
  //particleBouncerExperiment();

  // Math:
  //ellipseLineIntersections();
  //expBipolar();
  //expGaussBell();
  //linearRegression();
  //productLogPlot();
  //sinCosTable();

  // Filter:
  //bandSplittingTwoWay();
  //bandSplittingMultiWay();
  //bandSplittingTreeAlgo();
  //ladderResonanceManipulation();
  //nonUniformMovingAverage();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();
  //prototypeDesign();  // old implementation - tddo: check gains of prototype filters
  //poleZeroPrototype();  // new implementation - but we don't need that

  // Physics:
  //particleForceDistanceLaw();
  //particleSystem(); 

  // Generators:
  //bouncillator();
  //bouncillatorFormula();
  //rayBouncer();
  //hilbertCurve();
  //circleFractals(); // rename to spirograph
  //lindenmayer();
  //xoxosOsc();

  // Graphics:
  //lineDrawing();
  //lineDrawingThick();
  ///lineDrawingThick2(); // obsolete
  //lineJoints();
  //lineTo();
  //polyLineRandom();
  //phaseScopeLissajous();


  //===============================================================================================
  // Tests for dragged over RSLib code:

  int dummy;
  std::string str;
  bool passed = true;

  //-----------------------------------------------------------------------------------------------
  // Unit Tests:

  //passed &= testBufferFunctions(str);
  //passed &= testCopySection(    str);
  //passed &= testMoveElements(   str);
  //passed &= testRemoveElements( str);


  //passed &= testFilterPolynomials(str);

  //passed &= testHighOrderFilter(  str);  // fails


  ////passed &= testModalFilter2(str);       // fails
  ////passed &= testModalSynth(str);         // triggers assert

  ////passed &= testNumberManipulations( str); // triggers assert (calls the two below)
  ////passed &= testDoubleIntConversions(str); // triggers same assert (called by function above)
  //passed &= testExponentExtraction(str);

  //passed &= testAutoCorrelationPitchDetector(str);

  //passed &= testSortAndSearch(str);          // calls the two below (redundant)
  //passed &= testHeapSort(str);
  //passed &= testKnuthMorrisPrattSearch(str);

  //passed &= testTypeSizes(str);

  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Math:
  //binomialDistribution();
  //sineParameters();
  //bandLimitedStep();
  //cubicInterpolationNonEquidistant();   // move to unit tests
  //hyperbolicFunctions(); 
  //splineInterpolationNonEquidistant();
  //rationalInterpolation();
  //splineInterpolationAreaNormalized();
  //numericDerivative();
  //shiftPolynomial();
  ////void stretchPolynomial();  // commented in header
  //monotonicPolynomials();
  //parametricBell();
  //partialFractionExpansion();
  //partialFractionExpansion2();
  //partialFractionExpansionQuadratic();
  //dampedSineEnergy();
  //sineIntegral();
  //logarithmQuotient();
  //stirlingNumbers();
  //bernoulliNumbers();
  //sequenceSquareRoot();
  //conicSystem();
  ////logisticMapNoise(); // takes long to compute
  //bigFloatErrors();
  //primeRecursion();
  ////primeSieveSchmidt1(); // crashes
  ////primeSieveSchmidt2(); // crashes
  //primeSieveAtkin();
  ////primeSieve();  // crashes
  //primeDistribution();
  ////numberTheoreticTransform(); // triggers assert

  // Analysis:
  //autoCorrelation();
  //autocorrelationPeakVariation();
  //autoCorrelationPitchDetector();
  //autoCorrelationPitchDetectorOffline();
  //crossCorrelationBestMatch();
  //combineFFTs(); // move to math experiments
  ////zeroCrossingPitchDetector(); // commented in header - why?
  //instantaneousFrequency(); 
  ////instantaneousPhase();  // triggers assert (there's something not yet implemented)
  //zeroCrossingFinder();
  //zeroCrossingFinder2();
  //cycleMarkFinder();
  ////zeroCrossingPitchDetector(); // triggers assert
  //zeroCrossingPitchDetectorTwoTones();

  // Delay:
  //basicIntegerDelayLine();

  // Filter:
  //bandwidthScaling();
  //stateVariableFilter();
  //stateVariableFilterMorph();
  ////transistorLadder(); // triggers assert
  //phonoFilterPrototypePlot();
  //magnitudeMatchedOnePoleFilter();
  //phonoFilterModelPlot();
  //phonoFilterSimulation();
  //serialParallelBlend();
  //averager();
  ////movingAverage();  // hangs
  ////trapezAverager();   // hangs

  //gaussianPrototype();
  //halpernPrototype();
  //compareApproximationMethods();
  //compareOldAndNewEngineersFilter();
  //testPoleZeroMapper();

  //ringingTime();
  //butterworthSquaredLowHighSum();
  //maxFlatMaxSteepPrototypeM1N2();
  //maxFlatMaxSteepPrototypeM2N2();
  ////experimentalPrototypeM1N2();  // commented in header
  //splitLowFreqFromDC();
  //ladderResonanceModeling();
  //ladderResoShape();
  //ladderThresholds();           // maybe remove - this seemed to be a dead end
  //ladderFeedbackSaturation();
  //ladderFeedbackSaturation2();
  //ladderFeedbackSaturation3();
  //ladderFeedbackSatDCGain();
  //ladderFeedbackSatReso();
  //ladderFeedbackSatGrowl();
  //ladderFeedbackSatGrowl2();
  //ladderZDF();
  //ladderZDFvsUDF();
  //resoShapeFeedbackSat();
  //resoSaturationModes();
  //resoShapeGate();
  //resoShapePseudoSync();
  //resoSeparationNonlinear();
  //resoReplace();
  //resoReplacePhaseBumping();
  //resoReplaceScream();
  //fakeResonance();
  //fakeResoLowpassResponse();
  //fakeResoDifferentDelays();



  // Misc Audio:
  //centroid();
  //cubicCrossfade();
  //recursiveSineSweep();
  //ringModNoise();
  //slewRateLimiterLinear();
  //stretchedCorrelation();
  //taperedFourierSeries();
  //transientModeling();
  //windowFunctionsContinuous();
  //windowedSinc();

  // Modal:
  //modalFilter();        // impulse response of decaying-sine filter
  //modalFilterFreqResp();  // frequency response of attack/decay-sine filter - rename
  //attackDecayFilter();  // ...hmm..almost redundant
  //dampedSineFilterDesign();
  //biquadImpulseResponseDesign();
  modalBankTransient();

  // Modulator:
  //breakpointModulator();
  //breakpointModulatorSmoothFadeOut();

  // Oscillator:
  //triSaw();
  //phaseShapingCurvePoly4();
  //phaseShapingCurvesRational();
  //phaseShaping();
  //phaseShapingSkew();

  // Partial Extraction:
  //biDirectionalFilter();    // maybe move to filter tests
  //sineRecreation();         // maybe move elsewhere
  //sineWithPhaseCatchUp();   // dito
  //partialExtractionTriple();
  //partialExtractionViaBiquadTriple();
  ////partialExtractionBell();  // crashes because sample not available
  ////partialExtractionSample();  // dito

  // Phase Vocoder:
  //phaseRepresentation();
  //grainRoundTrip();        // under construction
  //plotWindows();
  //spectrogramSine();

  // Physics:
  //doublePendulum(); // takes long

  // Resampling:
  //fadeOut();  // move to a new file SampleEditingExperiments
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
  ////// tests with Elan's example files:
  ////pitchDetectA3();
  ////phaseLockSaxophone();
  ////phaseLockSaxophone2();
  ////autoTuneHorn();
  ////autoTuneHorn2();
  ////sylophoneCycleMarks();
  ////autoTuneSylophone();
  ////bestMatchShift();

  // Saturation:
  //powRatioParametricSigmoid();
  //parametricSigmoid();
  //parametricSigmoid2();
  //quinticParametricSigmoid(); 
  //septicParametricSigmoid();
  //saturator();
  //sigmoidScaleAndShift();
  //quarticMonotonic();
  //sigmoidPrototypes();
  //sixticPositive();

  //-----------------------------------------------------------------------------------------------
  // Performance Tests:

  // Analysis:
  //testFourierTransformer(str); 
  //testAutoCorrelationPitchDetector2(str);

  // Core:
  //testFlagArray(str);

  // Math:
  //testAbsAndSign2(str);              // rename (mabye "Perf" or sth)
  //testMultinomialCoefficients2(str); // rename
  //testPrimeSieves(str);
  ////testMatrix(str);                 // triggers assert
  //testMatrixAddressing(str);

  // Misc Audio:
  //testSincInterpolator(str);

  // Modal:
  //testModalFilter2(str);
  //testModalFilterBank(str);

  //-----------------------------------------------------------------------------------------------
  // Examples:

  // Modal:
  //createInsertionSortSound();  // move somewhere else
  //createModalFilterExamples();
  //createModalFilterBankExamples(); // takes long
  // sample-map creations (they take long):
  //createBass1();
  //createGong1();
  //createPluck1();





  dummy = 0;

  //===============================================================================================
  // RoSiC tests:

  //-----------------------------------------------------------------------------------------------
  // Unit tests:

  //testAllRosicClasses();
  //testRosicAnalysis();
  //testRosicBasics();
  //testRosicFile();
  //testRosicEffects();
  //testRosicGenerators();
  //testRosicFilter();
  //testRosicNumerical();
  //testRosicMath();
  //testRosicNonRealTime();
  //testRosicOthers();

  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Analysis:

  // Basics:

  // Effects:
  //testFastGeneralizedHadamardTransform();
  //testFeedbackDelayNetwork();
  //testMultiComp();

  // File:

  // Filters:
  //testLadderFilter();
  //testModalFilter();
  //testModalFilterWithAttack();
  //testBiquadPhasePlot();
  //testFiniteImpulseResponseDesigner();
  //testConvolverPartitioned();
  //testFiniteImpulseResponseFilter();
  //testFilterAnalyzer();
  //testBiquadCascade();
  //testCrossover4Way();
  //testCrossover4Way2();
  //testSlopeFilter();
  //testPrototypeDesigner();
  //testLowpassToLowshelf();
  //testBesselPrototypeDesign();
  //testPapoulisPrototypeDesign();
  //testEngineersFilter();
  //testPoleZeroMapping();
  //highOrderFilterPolesAndZeros();

  // Genrators:
  //testOscillatorStereo();
  //testLorentzSystem();
  //testSnowflake();
  //testResetter();
  //testTurtleReverse();
  //testTurtleSource();




  //-----------------------------------------------------------------------------------------------
  // Demos:
  // ...


  //-----------------------------------------------------------------------------------------------
  // Rendering:
  // ...

  //DEBUG_HOOK;

  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";
    // If memory leaks occur even though no objects are actually created on the heap, it could mean 
    // that some class in a library module has a static data member that does a dynamic memory 
    // allocation or some global object is created somewhere that dynamically allocates memory.

  getchar();
  return(EXIT_SUCCESS);
}
