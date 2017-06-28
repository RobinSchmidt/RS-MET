#include "PerformanceTestTools.h"

void printPerformanceTestResult(const char *name, double numCycles)
{
  std::cout << name;
  std::cout << ": ";
  std::cout << numCycles;
  std::cout << "\n";

  // maybe the formatting can be refined...
}