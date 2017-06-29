// This file and all the included ones below (and their cpp files) are obsolete - they are still
// kept here as reference as long as the new files are not complete

#include "Demos/Demos.h"
#include "Experiments/Experiments.h"
#include "PerformanceTests/PerformanceTests.h"
#include "UnitTests/UnitTests.h"

int main(int argc, char** argv)
{
  // Select, if you want to run demos, experiments, unit tests or performance tests:

  //runDemos();
  runExperiments();
  //runPerformanceTests();
  //runUnitTests();



  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";

  getchar();
  return(EXIT_SUCCESS);
}