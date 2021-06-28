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
      for(int k = 0; k < numRuns; k++) {
        cpuCounter.init();
        (*tests[i])((int)inputSizes[j]);  // run test i with input size (indexed by) j
        double cycles = (double)cpuCounter.getNumCyclesSinceInit();
        rawData[i][j][k] = cycles; }}}

  // Factor out into processResults:
  removeOutliers();
  computeMeansAndVariances();
  // Maybe we should not remove them but rather flag the outliers? Then we could not take them 
  // into account in the statistical analysis but still plot them in a scatter plot, perhaps with
  // a different color or different symbol like: normal data point: filled circle, outlier: cross
}

std::string PerformanceAnalyzer::getReport()
{
  rsError("Not yet implemented");

  std::string report;

  // formatting variables:
  //int lineWidth     = 80; // make user parameter
  //int columnWidth   = 6;
  //int maxNameLength = 6; // preliminary

  //size_t i, j;

  report += "Mean Values\n";

  //std::vector<std::vector<double>> means = getMeans();

  //std::vector<std::vector<double>> getMeans();
  //...

  return report;
}

void PerformanceAnalyzer::plotResults()
{
  rsError("Not yet implemented");
}
// may this be better handled by another, higher level class?


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

void PerformanceAnalyzer::removeOutliers()
{
  // for each test and input size, we must sort the array of results, find the median (center 
  // value) - maybe they should be stored in an array, too - but maybe not - they are trivial to 
  // recompute, once the rawData is sorted

  // i think, we should use inf for invalid datapoints, sort the raw data, make a copy, and
  // of each array, just cut off the inf-values at the end via resize in the copy - that is the
  // cleaned up data

  // or maybe use interquartile range:
  // https://en.wikipedia.org/wiki/Quartile
  // https://en.wikipedia.org/wiki/Interquartile_range

  // or use values between user adjustable quantiles, for example above the 10% quantile and
  // below the 90% quantile

}

void PerformanceAnalyzer::computeMeansAndVariances()
{


  // maybe refactor into computeMeans and computeVariances
}

/*

-Maybe specifying the size of the input is not enough in all cases. Maybe sometimes we want to test
 different datasets of the same size? Maybe we should have a more general API to pass datasets 
 instead of sizes? How is the data to be generated anyway? ...maybe the caller should pass 
 functions that generate the data? Or maybe we should just store a bunch of void-pointers that
 point to the datasets and it's up to the user-provided function to cast and interpret them?

*/