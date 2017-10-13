#include "../JuceLibraryCode/JuceHeader.h"

#include "rosic_tests/rosic_CorrectnessTests.h"
using namespace rotes;

#include "rapt_tests/Experiments/MathExperiments.h"
#include "rapt_tests/Experiments/FilterExperiments.h"
#include "rapt_tests/Experiments/GeneratorExperiments.h"
#include "rapt_tests/Experiments/GraphicsExperiments.h"

#include "rapt_tests/PerformanceTests/MathPerformanceTests.h"



#include "Experiments/Experiments.h"


//#include "Shared\Shared.h" // temporary - to figure out where the emory leak comes from

int main(int argc, char* argv[])
{
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
  // Performance tests:

  //matrixAdressingTest();
  //sinCosPerformance();


  //-----------------------------------------------------------------------------------------------
  // Experiments:

  // Prototypes (not yet in library):
  //particleBouncerExperiment();

  // Math:
  //ellipseLineIntersections();
  //expGaussBell();
  //linearRegression();
  //productLogPlot();
  //sinCosTable();

  // Filter:
  //ladderResonanceManipulation();
  //nonUniformMovingAverage();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();

  // Physics:
  //particleForceDistanceLaw();
  //particleSystem(); 

  // Generators:
  //rayBouncer();
  rayBouncer1D();

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

  // somehow, we always get memory leaks, even if we do nothing - figure out why - does some header
  // create a global object which leaks?
  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";

  getchar();
  return(EXIT_SUCCESS);
}
