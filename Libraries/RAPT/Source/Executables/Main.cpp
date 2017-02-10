#include "Demos/Demos.h"
#include "Experiments/Experiments.h"
#include "PerformanceTests/PerformanceTests.h"
#include "UnitTests/UnitTests.h"

int main(int argc, char** argv)
{
  // Select, if you want to run demos, experiments, unit tests or performance tests:

  //runDemos();
  //runExperiments();
  runPerformanceTests();
  //runUnitTests();



  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";

  getchar();
  return(EXIT_SUCCESS);
}