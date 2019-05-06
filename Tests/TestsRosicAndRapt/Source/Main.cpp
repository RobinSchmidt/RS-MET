//#define RS_PLOTTING  // somehow, it doesn't work to define it here - it needs to be defined
// in rapt.h - that's bad - figure out a better solution - maybe define it in the jucer-file
// when nothing else helps

#include "../JuceLibraryCode/JuceHeader.h"

// includes for unity build:
//#include "Shared/Shared.h"

#include "Experiments/Experiments.h"
#include "UnitTests/UnitTests.h"
#include "PerformanceTests/PerformanceTests.h"
#include "Misc/Misc.h"  // demos, examples, rendering, ... // todo: make unity build cpp file

// todo: move all the code into rs_testing module such that it can be compiled as a single 
// compilation unit -> faster build times for testing

#include "../../../Libraries/RobsJuceModules/romos/TestSuite/TestsMain.h"

// crash (access violation) if runAllUnitTests and envelopeDeBeating are run one after another
// it goes away when commenting out the code from
// passed &= runUnitTest(&testDifferentialEquationSystem, "rsDifferentialEquationSystem");
// to
// passed &= runUnitTest(&triangleRasterization,  "Triangle Rasterization");
// in UnitTests.cpp

// todo: maybe move some of the code into juce module rs_tests (or something) to make it more 
// convenient to used in the test-repo

int main(int argc, char* argv[])
{
  // tempoarary throw-away-code:
  //testCrossoverNewVsOld();

  //===============================================================================================
  // RAPT tests:

  //-----------------------------------------------------------------------------------------------
  // Unit tests:

  //runAllUnitTests();  // todo: merge with unit tests for RSLib

  //mathUnitTests();    // doesn't exist anymore ...it's all in runAllUnitTests now
  //filterUnitTests();  // dito (?)


  //-----------------------------------------------------------------------------------------------
  // Performance tests:

  //callbackPerformance();
  //matrixAdressingTest();
  //simdPerformance(1.0, rsFloat64x2(1.0));
  //simdPerformance(1.f, rsFloat32x4(1.f));
  //sinCosPerformance();

  //fftPerformance();
  //filterSignConventionPerformance();
  //ladderPerformance();
  //stateVectorFilterPerformance();
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
  //interpolatingFunction();
  //linearRegression();
  //metallicRatios();
  //polynomialPrediction();  // not yet implemented
  //probabilityLogic();
  //productLogPlot();
  //sinCosTable();
  //twoParamRemap();

  // Filter:
  //bandSplittingTwoWay();
  //bandSplittingMultiWay();   // turn into unit test (it currently hits an assert on fail)
  //bandSplittingTreeAlgo();
  //bandSplitFreqResponses();
  //complementaryFiltersIIR();
  //firstOrderFilters();
  //ladderResonanceManipulation();
  //nonUniformMovingAverage();
  //nonUniformOnePole();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();
  //prototypeDesign();  // old implementation - todo: check gains of prototype filters
  //poleZeroPrototype();  // new implementation - but we don't need that

  // Physics:

  //heatEquation1D();
  //particleForceDistanceLaw();
  //particleSystem();

  // Generators:
  //blit();
  //blep();
  //polyBlep();
  superBlep();
  //syncSweep();
  //syncPhasor();
  //syncPhasor2();
  //syncOsc();
  //dualBlepOsc();
  //bouncillator();
  //bouncillatorFormula();
  //freqVsPhaseMod();
  //rayBouncer();
  //hilbertCurve();
  //circleFractals(); // rename to spirograph
  //lindenmayer();
  //triSawOsc();
  //triSawOscAntiAlias();
  //xoxosOsc();

  // Graphics:
  //lineDrawing();
  //lineDrawingThick();
  ///lineDrawingThick2(); // obsolete
  //lineJoints();
  //lineTo();
  //polyLineRandom();
  //phaseScopeLissajous();
  //splineArc();
  //triangles();
  //pixelCoverage();


  // just for fun:
  //groupString();
  //primeAlternatingSums();
  //divisibility();
  //arithmeticDerivative();

  // third party code:
  //sampleTailExtenderTest();


  //===============================================================================================
  // Tests for dragged over RSLib code:

  //int dummy;
  //std::string str;    // causes a memleak
  bool passed = true;

  //-----------------------------------------------------------------------------------------------
  // Unit Tests:

  //passed &= testBufferFunctions(str);
  //passed &= testCopySection(    str);
  //passed &= testMoveElements(   str);
  //passed &= testRemoveElements( str);


  //passed &= testFilterPolynomials(str);

  //passed &= testHighOrderFilter(  str);  // fails


  //passed &= testModalFilter2(str);
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
  //bandMatrix();                       // under construction
  //pentaDiagnonalMatrix();
  //pentaDiagnonalMatrix2();
  //minSqrdDifsForFixSums();
//  minSqrdCurvForFixSums();
  //binomialDistribution();
  //sineParameters();
  //bandLimitedStep();
  //naturalCubicSpline();
  //naturalCubicSpline2();
  //cubicInterpolationNonEquidistant();   // move to unit tests
  //hyperbolicFunctions();
  //splineInterpolationNonEquidistant();
  //rationalInterpolation();
  //splineInterpolationAreaNormalized();
  //numericDiffAndInt();
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
//  envelopeFollower();
  ////zeroCrossingPitchDetector(); // commented in header - why?
  //instantaneousFrequency();
  ////instantaneousPhase();  // triggers assert (there's something not yet implemented)
  //maxShortTimeRMS();
  arrayRMS();
  //zeroCrossingFinder();
  //zeroCrossingFinder2();
  //zeroCrossingFinder3();
  //cycleMarkFinder();
  //cycleMarkErrors();
  ////zeroCrossingPitchDetector(); // triggers assert
  //zeroCrossingPitchDetectorTwoTones();


  // Delay:
  //basicIntegerDelayLine();

  // Filter:
  //bandwidthScaling();
  //butterworthEnergy();
  //stateVariableFilter();
  //stateVariableFilterMorph();
  //stateVectorFilter();   // just a stub, at the moment
  //biquadModulation();   // compares modulation properties of various biquad structures
  ////transistorLadder(); // triggers assert
  //phonoFilterPrototypePlot();
  //magnitudeMatchedOnePoleFilter();
  //phonoFilterModelPlot();
  //phonoFilterSimulation();
  //serialParallelBlend();
  //averager();
  //movingAverage();
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
  //ladderZDFvsUDF();   // compares modulation properties
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

  // Modal Filters/Synthesis:
  //modalFilter();        // impulse response of decaying-sine filter
  //modalFilterFreqResp();  // frequency response of attack/decay-sine filter - rename
  //modalTwoModes();
  //attackDecayFilter();  // ...hmm..almost redundant
  //dampedSineFilterDesign();
  //biquadImpulseResponseDesign();
  //modalBankTransient();
  //fourExponentials();  // weighted sum of 4 exponential envelopes - for shaping mode envelope
  //modalWithFancyEnv();
  //modalSynthSpectra();

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
//  windowFunctionSpectra();
  //windowedSinc();

  // Modulator:
  //breakpointModulator();
  //breakpointModulatorSmoothFadeOut();
  //triSawModulator();

  // Oscillator:
  //triSaw();
  //phaseShapingCurvePoly4();
  //phaseShapingCurvesRational();
  //phaseShaping();
  //phaseShapingSkew();

  // Partial Extraction:
  //biDirectionalFilter();    // maybe move to filter tests
  //envelopeDeBeating();
  //sineRecreation();         // maybe move elsewhere
  //sineWithPhaseCatchUp();   // dito
  //partialExtractionTriple();
  //partialExtractionViaBiquadTriple();
  ////partialExtractionBell();  // crashes because sample not available
  ////partialExtractionSample();  // dito

  // Phase Vocoder:
  //phaseRepresentation();
  //grainRoundTrip();        // under construction
  //spectrogramSine();
  //sineParameterEstimation();
  //plotWindows();   // rename to plotOverlappingWindowSum
  //phaseInterpolation();
  //sinusoidalSynthesis1();
  //sinusoidalSynthesis2();
  //sinusoidalAnalysis1();
  //sinusoidalAnalysis2();
  //sinusoidalAnalysis3();
  //phaseFreqConsistency();
  //harmonicDetection2Sines();
  //harmonicDetection3Sines();
  //harmonicDetection5Sines();
  //harmonicAnalysis1();


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

  ////// tests with Elan's example files (they don't work unless the files are available):
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
  //testModalFilter2(str);    // this is actually a unit test
  //testModalSynth(str);      // this too ...move elsewhere
  //testModalFilter3(str);
  //testModalFilterBank(str);

  //-----------------------------------------------------------------------------------------------
  // Examples:

  // Modal:
  //createInsertionSortSound();  // move somewhere else
  //createModalFilterExamples();
  //createModalFilterBankExamples(); // takes long
  //createPiano1();

  // sample-map creations (they take long):
  //createBass1();
  //createGong1();
  //createPluck1();

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

  //===============================================================================================
  // Modular:

  //runModularUnitTests();
  //runModularPerformanceTests(true);  // produces a memleak
  //testModularCodeGenerator();
  //runModularInteractiveTests();  // triggers assert due to plotting code
  romos::moduleFactory.clearRegisteredTypes(); // avoids memleak in unit tests

  // important atomic modules for performance tests:
  // Biquad: pure code, atomic module, wired model
  // Phasor

  //runModularTests(); // we need to make a .h file with the declarations


  //DEBUG_HOOK;
  //int* test = new int;  // uncomment, to see, if it fires correctly
  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected (pre exit of main()) !!! \n";
    //std::cout << "\n\n!!! Memory leaks detected !!! \n";
  // If memory leaks occur even though no objects are actually created on the heap, it could mean
  // that some class in a library module has a static data member that does a dynamic memory
  // allocation or some global object is created somewhere that dynamically allocates memory.

  // todo: fix the memory leak - i guess it's in rosic - test it by building the project with
  // rapt only - doesn't work - try to figure out, where heap memory is allocated..

  // set a debug-breakpoint _malloc_dbg in debug_heap.cpp and run the program - it gets called a
  // bunch of times from the startup code - skipping through these, later there will be a call
  // from the offending memory allocating code from our codebase
  // ..it seems to come from romos - compiling with rapt and rosic only doesn't produce memleaks

  // ok - its in:
  // void BlitIntegratorInitialStates::createStateValueTables()
  // deleteStateValueTables(); is called after the memleaks were detected hmmm..
  // solution: don't use global objects that freely lie around, instead use a
  // GlobalData class which is a singleton and encapsulates all sorts of data that should be
  // globally accessible
  // write a unit test for blit-saw (maybe there is one already) and change BlitSaw module to use
  // that - then, before checking for memory leak, call GlobalData::cleanUp (there may be a
  // GlobalData::initialize method as well that computes all the data/tables once and for all
  // or: allocate memory and compute the data-tables only when they are needed the first time
  // (lazy initialization)...can also be used for blep-tables, etc. - anything that needs globally
  // constant tables
  // maybe class ProcessingStatus can be extended for that - it would fit well there, too


  getchar();
  return(EXIT_SUCCESS);
}

// ToDo:
// -fix access violation in rsPrimeFactors - done
// -use more efficient implementation for rsPowInt
// -check that fabs or rsAbs is used everywhere where floating point numbers can occurr
//  -maybe use rsAbs preferably because it may also be used for modular integers