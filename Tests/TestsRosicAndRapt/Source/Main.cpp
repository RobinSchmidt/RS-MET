//#define RS_PLOTTING  // somehow, it doesn't work to define it here - it needs to be defined
// in rapt.h - that's bad - figure out a better solution - maybe define it in the jucer-file
// when nothing else helps


#include "../JuceLibraryCode/JuceHeader.h"

// includes for unity build:
//#include "Shared/Shared.h"


#include "rapt_tests/RaptTests.h"


// get rid of these includes - the best would be, to move all that stuff into the rs_testing juce
// module:
#include "Experiments/Experiments.h"


#include "rosic_tests/UnitTestsRosic.h"
using namespace rotes;  // Get rid of this !


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

int main(int argc, char* argv[])
{
  // Here, a lot of experimentation and test functions are called. It's for research and 
  // development. Most of the time, most function calls are commented out - the idea is to 
  // uncomment one at a time while working on a particular problem. We may also call the driver
  // routines for the unit tests here.


  // temporary throw-away-code:
  //testCrossoverNewVsOld();

  // todo: 
  // -with the increasing number of tests and experiments, it gets increasingly hard to find a 
  //  particular one - we need some better organization and order
  // -get rid of the distinction between testing classes from rapt and rosic - that makes it 
  //  confusing and hard to find a particular test -> merge the tests
  // -maybe make even just a single include file for all rapt tests
  // -maybe split this into Demos and Research - the demos shall remain in the main RS-MET 
  //  codebase, the research stuff should go into the research repo, eventually, 
  //  research-experiments may be propagated up into main repo
  // -maybe "override" malloc for debug builds to initialize the allocated memory with garbage as
  //  explained here https://www.youtube.com/watch?v=RoVD6zlftF0 this will help to detect bugs 
  //  related to uninitialized memory




  //===============================================================================================
  // Unit Tests:

  bool ok = true;
  //ok &= runUnitTestsRapt();
  //ok &= runUnitTestsRosic();
  //ok = ok;  // dummy instruction for setting a debug breakpoint here, if needed

  // The allpass unit test currently fails because I changed the implementation of 
  // rsSparsePolynomial to accomodate for coming new infrastruture to handle inexact floating 
  // point comparisons. The unit tests now implicitly do some exact float comparisions where 
  // inexact ones are needed and therefore fails. I also added a comment about this in this
  // commit:  84a9b2baf591da5c41d5d2dd31b310d6f1ead34f
  // What needs to be done to make it pass again is to use rsSparsePolynomial<double, double>
  // instead of rsSparsePolynomial<double> (which defaults the 2nd template parameter (for TTol, 
  // the tolerance) to rsEmptyType which results in the exact floating point comparisons. We need
  // to instantiate with TTol = double as well and then set up a suitable value for the tolerance.
  // I think, the last version that passes the test is the commit 
  // 5d5e4a5ee809c57ddc0be24d6f30a3845ab295ac with comment "Give rsSparsePolynomial a second 
  // template parameter TTol and a member of..." from 2025/03/28 or some commit very near that
  // one.
  // Done?

  // See:
  // Using Floating-point in C++: What Works, What Breaks, and Why - Egor Suvorov - CppCon 2025
  // https://www.youtube.com/watch?v=m83TjrB6wYw  at 23:50
  // for a possible way to fix the currently failing float-string-float roundtrip unit test. Use
  // the constant std::numeric_limits<float>::max_digits10  ...and likewise for double. This 
  // constant defines how many decimal digits we need for exact roundtrips. Dont's confuse it with
  // digits10 which (I think) is the precision limit for exact string-float-string roundtrips.
  // Also interesting: at around 53 min, he mentions that std::complex can be significantly slower
  // than a handwritten class


  // ToDo: let the functions take an integer argument that specifies the "level" of exhaustiveness
  // of testing. 0: should be able to do all tests in 5 seconds, 1: 20 seconds, 2: 80 seconds etc.
  // ...we may run very exhaustive tests that may take hours - but we don't want to run them as 
  // often as the quick tests

  //mathUnitTests();    // doesn't exist anymore ...it's all in runAllUnitTests now
  //filterUnitTests();  // dito (?)



  //===============================================================================================
  // Tests for the test infrastructure:

  //printTests();   // Prints numbers, matrices, etc.




  //===============================================================================================
  // RAPT Experiments

  //-----------------------------------------------------------------------------------------------
  // Performance tests:

  //callbackPerformance();
  //matrixAdressingTest();
  //simdPerformance();
  //sinCosPerformance();
  //fftPerformance();

  //filterSignConventionPerformance();
  //ladderPerformance();
  //stateVectorFilterPerformance();
  //engineersFilterPerformance();
  //turtleGraphicsPerformance();
  //samplerEnginePerformance();

  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Prototypes (not yet in library):
  //particleBouncerExperiment();

  // Math:

  // Linear Algebra:
  //determinant();
  //characteristicPolynomial();
  //testSubSpaces();               // todo: move to unit tests
  //testSigularValueDecomp();      // dito
  //linearIndependence();
  //eigenstuff();
  //iterativeLinearSolvers();
  //linearSolverPrecision();

  // Geometry:
  //ellipseLineIntersections();

  // Function Evaluation:
  //expBipolar();
  //expGaussBell();
  //fmodTest();                      // Tests for wrap-around functions
  //gaussBellProduct();
  //polynomialSinc();                // approximation of sinc by polynomial constructed from roots
  //productLogPlot();
  //sinCosTable();                   // new polynomial approximation is called from unit tests
  //softBitCrush();
  //twoParamRemap();
  //unitIntervalMap();
  //powerIterator();
  //gaussianIterator();
  //expPolyIterator();
  //reciprocalIterator();  // rename to multiStepSolverIVP (IVP: initial value problem)
  // It implements prototypes of Adams-Bashforth, Adams-Moulton, BDF methods and more using the 
  // ODE for 1/x as example problem (I think)

  // Interpolation:
  //linearFractionalInterpolation();  // aka LinFrac
  //monotonicInterpolation();
  //interpolatingFunction();

  // Curve Fitting:
  //linearRegression();
  //multipleRegression();
  //polynomialRegression();
  //gaussianRegression();
  //butterworthViaGaussians();

  // Numeric Calculus
  //numericOptimization();
  //numericMinimization1D();
  //numericRootFinding1D();
  //iteratedNumDiff();

  //polynomialPrediction();     // not yet implemented
  //probabilityLogic();
  //ratioGenerator();
  //ratiosLargeLcm();
  //ratiosEquidistantPowers();
  //ratiosMetallic();
  //numberTheoreticTrafo();   // move to unit tests!
  //numberTheoreticTrafoModuli();

  // Misc math:
  //mathErrorsTest();                // Tests for behaviors in case of math errors, stub



  // Filter:
  //bandpassAndNotch();          // trying to build a bank of complementary bandpass/notch pairs
  //bandSplittingTwoWay();
  //bandSplittingThreeWay();
  //bandSplittingThreeWay2p2z();
  //bandSplittingMultiWay();   // turn into unit test (it currently hits an assert on fail)
  //bandSplittingTreeAlgo();
  //bandSplitFreqResponses();
  //biDirectionalStateInit();
  //biquadDesignVicanek();         // maybe rename to biquadMatchVicanek
  //biquadStability();
  //biquadTail();
  //biquadModulation();            // compares modulation properties of various biquad structures
  //brickwallDeRinging();         // stub
  //complementaryFiltersIIR();
  //ringingResponses();
  //engineersFilterFreqResps();
  //engineersFilterFreqRespsMeasured();  // Measured freq responses - the computated ones may lie!
  //engineersFilterMethodsComparison();
  //firstOrderFilters();
  //ladderResonanceManipulation();
  //nonUniformMovingAverage();
  //nonUniformOnePole1();
  //nonUniformOnePole2();
  //nonUniformComplexOnePole();
  //nonUniformAllpole(); // rename - it's not restricted to allpoles anymore
  //nonUniformBiquad();
  //nonUniformBiDirectional();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();
  //prototypeDesign();  // old implementation - todo: check gains of prototype filters
  //poleZeroPrototype();  // new implementation - but we don't need that
  //seriesConnectionDecay();
  //quantileFilter();
  //hilbertFilter();
  //simdFilter<float, 4>();  // doesn't work
  //subBandFilter();

  // Physics:
  //doublePendulum(); // takes long
  //heatEquation1D();
  //waveEquation1D();
  //rectangularMembrane();
  //rectangularRoom();
  //particleForceDistanceLaw();
  //particleSystem();
  //quantumSpinMeasurement();  // move to unit tests
  //quantumSpinEntanglement();
  //quantumGates();            // move to unit tests
  //quantumComputer();
  //quantumSpinEvolution();
  //quantum3StateSystem();
  //quantumParticle();
  //tennisRacket();
  //tennisRacket2();
  //tennisRacket3();
  //tennisRacketFreq();

  // Generators:
  //waveformFractalization();
  //noise();
  //noiseReverseMode();
  //noiseTriModal();
  //noiseWaveShaped();
  //blit();
  //blep();
  //polyBlep();
  //superBlep();
  //superSawDensitySweep();
  //superSawStereo();
  //twoPieceOsc();
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
  //snowFlake();
  //triSawOsc();
  //triSawOscAntiAlias();
  //xoxosOsc();
  //shepardTone();
  //additiveEngine();
  //multiplicativeSynth();
  //pulseWidthModulationViaTwoSaws();  // just a stub at the moment
  //flatZapper();  // Various experiments
  //freqSweeper();
  //sineSweepBassdrum();

  // Modulators:
  //attackDecayEnvelope();

  // Graphics:
  //colorGradientHSL();
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

  // Plotting:
  //gradientify();
  //contours();
  //complexContours();
  //implicitCurves();
  //parametricCurve();
  //spirals();                 // move to Rendering.cpp in research codebase
  //fractal();
  //differentialGeometry();

  // Image processing:
  //imageScaling();

  // Just for fun (TODO: move to research repo):
  //groupString();
  //primeAlternatingSums();
  //divisibility();
  //arithmeticDerivative();






  // third party code:
  //sampleTailExtenderTest();


  // Math:
  //bandMatrix();                           // under construction
  //pentaDiagnonalMatrix();
  //pentaDiagnonalMatrix2();
  //minSqrdDifsForFixSums();
  //minSqrdCurvForFixSums();
  //binomialDistribution();
  //sineParameters();
  //bandLimitedStep();

  //chebychevInterpolant();
  //cubicSplines();
  //cubicInterpolationNonEquidistant();     // move to unit tests
  //hyperbolicFunctions();
  //splineInterpolationNonEquidistant();
  //rationalInterpolation();
  //splineInterpolationAreaNormalized();

  //numericDifferentiation();                // num. dif. on a function object
  //numericIntegration();                   // a.k.a. numeric "quadrature"
  //nonUniformArrayDiffAndInt();            // differentiation and integration of sampled data
  //testNonUniformInvertibleDiff();         // algo that produces reciprocals when swapping x and y
  //uniformArrayDiffAndInt();               // under construction

  //derivativeFormulas();                   // still empty
  //vertexMeshGradient();                   // implementations are in MeshExperiments.cpp
  //vertexMeshHessian();
  //vertexMeshLaplacian();

  //convolvePolynomials();                  // obsolete
  //convolvePiecewise();
  //shiftPolynomial();
  ////void stretchPolynomial();             // commented in header
  //monotonicPolynomials();
  //mixedPolynomialRoots();

  //parametricBell();
  //partialFractionExpansion();
  //partialFractionExpansion2();
  //partialFractionExpansion3();
  //partialFractionExpansionQuadratic();
  //dampedSine();
  //sineIntegral();
  //logarithmQuotient();
  //stirlingNumbers();
  //bernoulliNumbers();
  //bernoulliPolynomials();
  //sequenceSquareRoot();
  //conicSystem();
  ////logisticMapNoise();                   // takes long to compute
  //bigFloatErrors();
  //primeRecursion();
  ////primeSieveSchmidt1(); // crashes
  ////primeSieveSchmidt2(); // crashes
  //primeSieveAtkin();
  ////primeSieve();  // crashes
  //primeDistribution();
  //numberTheoreticTransform();           // triggers assert
  //variousFunctions();
  //functionOperators();

  // Analysis:
  //autoCorrelation();
  //autocorrelationPeakVariation();
  //autoCorrelationPitchDetector();
  //autoCorrelationPitchDetectorOffline();
  //crossCorrelationBestMatch();
  //combineFFTs(); // move to math experiments
  //envelopeFollower();
  ////zeroCrossingPitchDetector(); // commented in header - why?
  //instantaneousFrequency();
  ////instantaneousPhase();  // triggers assert (there's something not yet implemented)
  //maxShortTimeRMS();
  //arrayRMS();
  //peakFinder();                    // find peaks with subsample precision
  peakSmoother();
  //zeroCrossingFinder();
  //zeroCrossingFinder2();
  //zeroCrossingFinder3();
  //cycleMarkFinder();
  //cycleMarkErrors();
  //zeroCrossingPitchDetector(); // triggers assert
  //zeroCrossingPitchDetectorTwoTones();
  //ropewayAlgo();
  //peakPicker();
  //singleSineModel();
  //singleSineCycleWobbles();    // maybe rename - sineFromWub, sineFromWavelet
  //multiSineCycleWobbles();
  //multiHalfCycleWobbles();
  //sineFromDecayingSines();

  // Delay, allpass, reverb stuff:
  //delayLines();                    // Also includes allpass, universal comb, comb vs modal bank stuff
  //twoPoleAllpassDelays();
  //dampedCombAllpasses();
  //feedbackDelayNetworks();
  //allpassFDN();                  // Under construction

  // Filter:
  //bandwidthScaling();
  //biquadResoGainToQ();           // investigate relation beween filter Q and resonance gain
  //butterworthEnergy();
  //onePoleFilterSimper();           // Stub - maybe get rid of it!
  //sallenKeyFilterSimper();
  //stateVariableFilters();          // Tests with various SVF implementations
  //stateVectorFilter();           // Stub


  //transistorLadder();    // triggers assert
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
  //compareOldAndNewEngineersFilter(); // is this obsolete now? i think, there is no "old" version anymore
  //testPoleZeroMapper();

  //ringingTime();
  //butterworthSquaredLowHighSum();
  //maxFlatMaxSteepPrototypeM1N2();
  //maxFlatMaxSteepPrototypeM2N2();
  ////experimentalPrototypeM1N2();  // commented in header
  //splitLowFreqFromDC();
  //directFormFreqResp();
  //ladderResonanceGain();
  //ladderTransferFunction();
  //ladderMultipole();
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
  //ladderZDFvsUDF();       // compares cutoff modulation properties
  //ladderResoModulation();
  //resoShapeFeedbackSat();
  //resoSaturationModes();
  //resoShapeGate();
  //resoShapePseudoSync();
  //resoSeparationNonlinear();
  //resoReplace();
  //resoReplacePhaseBumping();
  //resoReplaceScream();
  //resoWave();
  //fakeResonance();
  //fakeResoLowpassResponse();
  //fakeResoDifferentDelays();
  //samplerFilters();

  // Modal Filters/Synthesis:
  //twoPoleFilter();
  //modalFilter();        // impulse response of decaying-sine filter
  //modalFilterFreqResp();  // frequency response of attack/decay-sine filter - rename
  //attackDecayFilter();  // ...hmm..almost redundant
  //modalTwoModes();
  //dampedSineFilterDesign();
  //dampedSineFilterImpResp();
  //biquadImpulseResponseDesign();
  //modalBankTransient();
  //fourExponentials();  // weighted sum of 4 exponential envelopes - for shaping mode envelope
  //modalWithFancyEnv();
  //modalSynthSpectra();
  //modalDecayFit();
  //modalAnalysis1();
  //modalAnalysisPluck();
  //modalPartialResynthesis();
  //modalAnalysisGloriosa();
//  modalReverb();               // Under construction

  // Misc Audio:
  //centroid();
  //cubicCrossfade();
  //decimate();
  //pythagoreanTuning();    // move to some music/music-theory section
  //recursiveSineSweep();
  //recursiveSineWithCubicPhase();
  //ringModNoise();
  //slewRateLimiterLinear();
  //slewRateLimiterPolynomial();   // experimental, does not yet work as intended
  //stretchedCorrelation();
  //taperedFourierSeries();
  //transientModeling();
  //windowFunctionsContinuous();
  //windowFunctionSpectra(); // todo: try bump-function and piecewise window using integrated bump tapers
  //windowedSinc();
  //waveMorph();  // under construction
  //squareToSaw();    // shapes square-wave into saw-wave



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
  //phaseShapingLinFrac();           // stub
  //zeroDelayFeedbackPhaseMod();

  // Partial Extraction:
  //biDirectionalFilter();    // maybe move to filter tests
  //beatingSines();
  //envelopeDeBeating();
  //sineRecreation();               // maybe move elsewhere
  //sineRecreationBandpassNoise();
  //sineWithPhaseCatchUp();       // dito
  //partialExtractionTriple();
  //partialExtractionViaBiquadTriple();
  ////partialExtractionBell();  // crashes because sample not available
  ////partialExtractionSample();  // dito

  // Phase Vocoder (split into Spectrogram and SineModel):
  //phaseRepresentation();
  //grainRoundTrip();        // under construction
  //spectrogramSine();
  //spectrogramFilter();
  //sineParameterEstimation();
  //plotOverlappingWindowSum();
  //phaseInterpolation();
  //sinusoidalSynthesisDC();
  //sinusoidalSynthesis1();
  //sinusoidalSynthesis2();
  //sinusoidalAnalysis1();
  //sinusoidalAnalysis2();  // something fails terribly here
  //sinusoidalAnalysis3();
  //phaseFreqConsistency();

  //harmonicDetection2Sines();
  //harmonicDetection3Sines();
  //harmonicDetection5Sines();
  //harmonicAnalysis1();

  //amplitudeDeBeating();
  //amplitudeDeBeating2();
  //harmonicDeBeating1();
  //harmonicDeBeating2();


  // Resampling:
  //fadeOut();  // move to a new file SampleEditingExperiments
  //resampler();
  //resamplerDelay();
  //sincResamplerAliasing();
  //sincResamplerModulation();
  //sincResamplerPassbandRipple();
  //sincResamplerSumOfTapWeights();
  //timeWarp();
  //pitchFlattening();
  //phaseLockedCrossfade();
  //phaseLockedCrossfade2();
  //pitchDetectWithSilence();

  // Matching:
  //sineShift();
  //sineShift2();
  //amplitudeMatch();
  //amplitudeMatch2();


  ////// tests with Elan's example files (they don't work unless the files are available):
  ////pitchDetectA3();
  ////phaseLockSaxophone();
  ////phaseLockSaxophone2();
  ////autoTuneHorn();
  ////autoTuneHorn2();
  ////sylophoneCycleMarks();
  ////autoTuneSylophone();
  ////bestMatchShift();
  // move them into the test repo and add the relevant sample files there (if i still can find
  // them, that is)

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
  //sigmoidConvergenceRates();
  //sixticPositive();
  hilbertDistortion();
  //adHocTapeEmuIdeas();
  //tapeEmulation();                  // Jatin Chowdhury's tape hysteresis algorithm


  // Distortion:




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
  //testModalFilter3(str);
  //testModalFilterBank(str);

  //-----------------------------------------------------------------------------------------------
  // Examples:


  //===============================================================================================
  // RoSiC tests:


  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Analysis:
  //testOscilloscopeBuffer();

  // Basics:
  //testWindowFunctions();
  //testInterpolation();

  // Reverb:
  // allpassDisperser();                // Rename to allpassDiffusor
  // allpassDelay();
  // allpassDelayChain();
  // allpassDelayChainVsNest();
  // These 4 are actually experiments but they ended up in file (an namespace) for unit tests. 
  // Move them to a better place. Maybe to DelayExperiments.h/cpp

  // feedbackDelayNetwork();            // writes wave file
  algoVerb();                        // writes wave file



  // Spectral effects:
  //spectralFilter();                  // Maybe move into a file for spectral processors
  //formantShifter();
  //spectralShifter();


  // File:

  // Math:
  //testLinLogEquationSolver();

  // Filters:
  //testLadderFilter();
  //testModalFilter();
  //testModalFilterWithAttack();
  //testBiquadPhasePlot();
  //testFilterAnalyzer();
  //testBiquadCascade();
  //testSlopeFilter();
  //testPrototypeDesigner();
  //testLowpassToLowshelf();
  //testBesselPrototypeDesign();
  //testPapoulisPrototypeDesign();
  //testEngineersFilter();
  //testPoleZeroMapping();
  //highOrderFilterPolesAndZeros();
  //testCrossover4Way();
  //testCrossover4Way2();
  //testFiniteImpulseResponseDesigner();
  //testConvolverPartitioned();         // is unit test
  //testFiniteImpulseResponseFilter();  // is unit test

  // Genrators:
  //testOscillatorStereo();
  //testLorentzSystem();              // it's spelled Lorenz - without the t!
  //testCombustionEngine();         // stub
  //testSnowflake();
  //testResetter();
  //testTurtleReverse();
  //testTurtleSource();
  //testSamplerEngine();

  // Modulators:
  //testConsecutiveExponentialDecay();

  // Others:
  //testSlewRateLimiterLinear();






  //-----------------------------------------------------------------------------------------------
  // Demos:
  // ...


  //-----------------------------------------------------------------------------------------------
  // Rendering:

  // Modal:
  //createInsertionSortSound();  // move somewhere else
  createModalFilterExamples();
  createModalFilterBankExamples(); // takes long
  //createPiano1();

  // The new renering scripts for creating sample content for the sfz engine:
  //createMiscSamples();
  //createAllpassDrums();
  //createSamplerWaveforms();

  // Older sample-map creations based on modal synthesis (they take long):
  createBass1();
  //createGong1();
  //createBell1();
  //createPluck1();
  //testHighPluck();
  // ToDo: create from the same sample-sets also soundfonts with 1,2,3,4,6 samples per octave via
  // key-crossfading (the default is 12 per octave, i.e. 1 sample per key) - compare them to find 
  // the best trade-off between size and quality (probably 3 or 4?)


  //===============================================================================================
  // Modular:

  //runModularUnitTests();             // MUST run before performance tests (or access violation)
  //runModularPerformanceTests(true);  // produces a memleak unless we call clearRegisteredTypes() 
  //testModularCodeGenerator();
  //runModularInteractiveTests();  // triggers assert due to plotting code
  romos::moduleFactory.clearRegisteredTypes(); // avoids memleak in unit tests

  // important atomic modules for performance tests:
  // Biquad: pure code, atomic module, wired model
  // Phasor

  //runModularTests(); // we need to make a .h file with the declarations


  //DEBUG_HOOK;
  //int* leakTest = new int;  // Uncomment to see, if the memleak test fires correctly to rule out
                              // false negatives.
  if(detectMemoryLeaks())
  {
    std::cout << "\n\n!!! Memory leaks detected (pre exit of main()) !!! \n";
    //getchar();
  }
  // If memory leaks occur even though no objects are actually created on the heap, it could mean
  // that some class in a library module has a static data member that does a dynamic memory
  // allocation or some global object is created somewhere that dynamically allocates memory.

  // ToDo:
  //
  // -Maybe #define _CRTDBG_MAP_ALLOC in debug on in in Viusla Studio. This should give more
  //  detailed leak information. See:
  //  https://github.com/surge-synthesizer/surge/issues/7478
  //  https://forum.juce.com/t/conflict-with--crtdbg-map-alloc-visual-studio/6040/10
  //  https://forum.juce.com/t/need-memory-leak-detection-with-filename/10399
  //  https://forum.juce.com/t/memory-leak-tips/1109/9
  //  https://forum.juce.com/t/leak-detector/6516/6
  //  https://forum.juce.com/t/memory-leak-tool/13330/13
  //
  // -Check out these tools and maybe start using them:
  //  -Microsoft Application Verifier (AppVerifier)
  //   https://learn.microsoft.com/en-us/windows-hardware/drivers/devtest/application-verifier
  //  -Valgrind:
  //   https://valgrind.org/docs/manual/quick-start.html
  //  -Visual Leak Detector:
  //   https://marketplace.visualstudio.com/items?itemName=ArkadyShapkin.VisualLeakDetectorforVisualC
  //  -Address Sanitizer:
  //   https://learn.microsoft.com/de-de/cpp/sanitizers/asan?view=msvc-170

  // [Done?]:
  // ToDo: fix the memory leak - i guess it's in rosic - test it by building the project 
  // with rapt only - doesn't work - try to figure out, where heap memory is allocated.
  //
  // set a debug-breakpoint _malloc_dbg in debug_heap.cpp and run the program - it gets called a
  // bunch of times from the startup code - skipping through these, later there will be a call
  // from the offending memory allocating code from our codebase
  // ..it seems to come from romos - compiling with rapt and rosic only doesn't produce memleaks
  //
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
  // Having the getchar here is more convenient for running the unit tests. Having it wrapped in
  // the "if(detectMemoryLeaks())" conditional makes more sense for experiments.


  return(EXIT_SUCCESS);
}

// ToDo:
// -use more efficient implementation for rsPowInt
// -check that fabs or rsAbs is used everywhere where floating point numbers can occurr
//  -maybe use rsAbs preferably because it may also be used for modular integers