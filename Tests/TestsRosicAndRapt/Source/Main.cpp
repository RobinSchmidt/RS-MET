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

  //runAllUnitTests();
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
  //prototypeDesign();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();

  // Physics:
  //particleForceDistanceLaw();
  //particleSystem(); 

  // Generators:
  //bouncillator();
  //bouncillatorFormula();
  //rayBouncer();
  //hilbertCurve();
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

  passed &= testBufferFunctions(str);
  passed &= testCopySection(    str);
  passed &= testMoveElements(   str);
  passed &= testRemoveElements( str);

  passed &= testFilterPolynomials(str);
  //passed &= testHighOrderFilter(  str);  // fails

  //passed &= testModalFilter2(str);       // fails
  //passed &= testModalSynth(str);         // triggers assert

  //passed &= testNumberManipulations( str); // triggers assert (calls the two below)
  //passed &= testDoubleIntConversions(str); // triggers same assert (called by function above)
  passed &= testExponentExtraction(str);

  passed &= testAutoCorrelationPitchDetector(str);

  passed &= testSortAndSearch(str);          // calls the two below (redundant)
  passed &= testHeapSort(str);
  passed &= testKnuthMorrisPrattSearch(str);

  passed &= testTypeSizes(str);

  //-----------------------------------------------------------------------------------------------
  // Experiments:

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




  //-----------------------------------------------------------------------------------------------
  // Performance Tests:

  //-----------------------------------------------------------------------------------------------
  // Examples:


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
