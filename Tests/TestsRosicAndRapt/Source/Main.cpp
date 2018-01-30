#include "../JuceLibraryCode/JuceHeader.h"

#include "rosic_tests/rosic_CorrectnessTests.h"
using namespace rotes;

#include "rapt_tests/Experiments/MathExperiments.h"
#include "rapt_tests/Experiments/FilterExperiments.h"
#include "rapt_tests/Experiments/GeneratorExperiments.h"
#include "rapt_tests/Experiments/GraphicsExperiments.h"

#include "rapt_tests/PerformanceTests/MathPerformanceTests.h"
#include "rapt_tests/PerformanceTests/AudioPerformanceTests.h"

#include "rapt_tests/UnitTests/UnitTests.h"


#include "Experiments/Experiments.h"


int main(int argc, char* argv[])
{
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


  //===============================================================================================
  // RAPT tests:

  //-----------------------------------------------------------------------------------------------
  // Unit tests:

  //runAllUnitTests();
  //mathUnitTests();
  //filterUnitTests();


  //-----------------------------------------------------------------------------------------------
  // Performance tests:

  //matrixAdressingTest();
  //simdPerformanceFloat64x2();
  //sinCosPerformance();

  //ladderPerformance();
  //engineersFilterPerformance();


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
  xoxosOsc();

  // Graphics:
  //lineDrawing();
  //lineDrawingThick();
  ///lineDrawingThick2(); // obsolete
  //lineJoints();
  //lineTo();
  //polyLineRandom();
  //phaseScopeLissajous();

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
