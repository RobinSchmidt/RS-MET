#include "PerformanceTestTools.h"

void printPerformanceTestResult(const char *name, double numCycles)
{
  printPerformanceTestResult(std::string(name), numCycles);
  /*
  std::cout << name;
  std::cout << ": ";
  std::cout << numCycles;
  std::cout << "\n";
  */
}

void printPerformanceTestResult(const std::string& name, double numCycles)
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

//=================================================================================================

void PerformanceAnalyzer::init()
{
  clearTests();
  clearInputSizes();
  numRuns = 30;
  smallOutlier = -std::numeric_limits<double>::infinity();
  largeOutlier =  std::numeric_limits<double>::infinity();
  rawData.clear();
  means.clear();
  variances.clear();
}

void PerformanceAnalyzer::runTests()
{
  initResultArray();
  for(size_t i = 0; i < tests.size(); i++) {
    for(size_t j = 0; j < inputSizes.size(); j++) {
      for(size_t k = 0; k < numRuns; k++) {
        cpuCounter.init();
        (*tests[i])(inputSizes[j]);  // run test i with input size (indexed by) j
        double cycles = cpuCounter.getNumCyclesSinceInit();
        rawData[i][j][k] = cycles; }}}
  int dummy = 0;
}

std::string PerformanceAnalyzer::getReport()
{
  std::string report;

  //...

  return report;
}

void PerformanceAnalyzer::initResultArray()
{
  rawData.clear();
  rawData.resize(tests.size());
  for(size_t i = 0; i < tests.size(); i++) {
    rawData[i].resize(inputSizes.size());
    for(size_t j = 0; j < inputSizes.size(); j++)
      rawData[i][j].resize(numRuns);
  }
}