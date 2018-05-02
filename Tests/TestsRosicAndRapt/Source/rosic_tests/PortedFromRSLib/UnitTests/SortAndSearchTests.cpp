#include "SortAndSearchTests.h"

bool testSortAndSearch(std::string &reportString)
{
  std::string testName = "rsSortAndSearch";
  bool testResult = true;

  testResult &= testHeapSort(reportString);
  testResult &= testKnuthMorrisPrattSearch(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testHeapSort(std::string &reportString)
{
  std::string testName = "rsHeapSort";
  bool testResult = true;

  int numTests = 20;
  static const int length = 100;
  int testArray[length];
  for(int i=0; i<numTests; i++)
  {
    rsFillWithRandomValues(testArray, length, -1000, +1000, 1);
    rsHeapSort(testArray, length);
    testResult &= rsIsSortedAscending(testArray, length);

    // check odd lengths by just sorting the subarray up to length-1:
    rsFillWithRandomValues(testArray, length, -1000, +1000, 1);
    rsHeapSort(testArray, length-1);
    testResult &= rsIsSortedAscending(testArray, length-1);
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testKnuthMorrisPrattSearch(std::string &reportString)
{
  std::string testName = "rsKnuthMorrisPrattSearch";
  bool testResult = true;

  char text[31]   = "abababbbaababbaabababaabababaa";  // 30 + '\0'
  char pattern[5] = "abab";                            //  4 + '\0'

  std::vector<int> matches = rsFindAllOccurencesOf(text, 30, pattern, 4);

  testResult &= matches.size() == 7;
  testResult &= matches[0] == 0;
  testResult &= matches[1] == 2;
  testResult &= matches[2] == 9;
  testResult &= matches[3] == 15;
  testResult &= matches[4] == 17;
  testResult &= matches[5] == 22;
  testResult &= matches[6] == 24;

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


