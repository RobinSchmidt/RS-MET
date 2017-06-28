#include "../JuceLibraryCode/JuceHeader.h"

#include "rosic_tests/rosic_CorrectnessTests.h"
using namespace rotes;

int main (int argc, char* argv[])
{

  // Unit tets:

  //testAllRosicClasses();
  //testRosicAnalysis();
  //testRosicBasics();
  //testRosicFile();
  //testRosicEffects();
  //testRosicGenerators();
  testRosicFilter();
  //testRosicNumerical();
  //testRosicMath();
  //testRosicNonRealTime();
  //testRosicOthers();

  // Performance tests:
  // ...

  // Experiments:
  //particleBouncerExperiment();
  // ...

  // Demos:
  // ...

  // Rendering:
  // ...




  //DEBUG_HOOK;

  return 0; // if we end here without interruption by a debug-break, all tests have passed
}
