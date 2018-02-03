#include "PerformanceTestTools.h"

void printPerformanceTestResult(const char *name, double numCycles)
{
  std::cout << name;
  std::cout << ": ";
  std::cout << numCycles;
  std::cout << "\n";

  // maybe the formatting can be refined...
}

void print(const char *name, int x)
{
  std::cout << name;
  std::cout << ": ";
  std::cout << x;
  std::cout << "\n";
}
// code-duplication - make template print, instantiate for double and int...maybe printLine


template<class T>
void printMemoryOccupation(const char *name, T& object)
{
  std::cout << name;
  std::cout << ": ";
  std::cout << sizeof(object);
  std::cout << "\n";
}

