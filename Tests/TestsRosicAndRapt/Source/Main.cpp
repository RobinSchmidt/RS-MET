#include "../JuceLibraryCode/JuceHeader.h"

#include "rosic_tests/rosic_CorrectnessTests.h"
using namespace rotes;

#include "rapt_tests/Experiments/MathExperiments.h"
#include "rapt_tests/Experiments/FilterExperiments.h"
#include "rapt_tests/Experiments/GeneratorExperiments.h"
#include "rapt_tests/Experiments/GraphicsExperiments.h"

#include "rapt_tests/PerformanceTests/MathPerformanceTests.h"

#include "Experiments/Experiments.h"

/*
// temporary - to figure out where the emory leak comes from (so we may can comment out all 
// headers above):
#include <crtdbg.h>
#include <iostream>
inline bool detectMemoryLeaks()
{
#ifdef _MSC_VER
  return (_CrtDumpMemoryLeaks() == 1);
#else
  return false;
#endif
}
*/

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
  // Performance tests:

  //matrixAdressingTest();
  //sinCosPerformance();


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
  //ladderResonanceManipulation();
  //nonUniformMovingAverage();
  //smoothingFilterOrders();
  //smoothingFilterTransitionTimes();

  // Physics:
  //particleForceDistanceLaw();
  //particleSystem(); 

  // Generators:
  //bouncillator();
  //bouncillatorFormula();
  //rayBouncer();
  //xoxosOsc();

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
