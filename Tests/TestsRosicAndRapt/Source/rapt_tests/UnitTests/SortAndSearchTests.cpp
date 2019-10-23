#include "SortAndSearchTests.h"

bool testHeapSort()
{
  bool testResult = true;

  int numTests = 20;
  static const int length = 100;
  int testArray[length];
  for(int i=0; i<numTests; i++)
  {
    RAPT::rsArray::fillWithRandomValues(testArray, length, -100, +100, 1);
    rsHeapSort(testArray, length);
    testResult &= rsArray::isSortedAscending(testArray, length);

    // check odd lengths by just sorting the subarray up to length-1:
    RAPT::rsArray::fillWithRandomValues(testArray, length, -100, +100, 1);
    rsHeapSort(testArray, length-1);
    testResult &= rsArray::isSortedAscending(testArray, length-1);
  }

  return testResult;
}

bool testKnuthMorrisPrattSearch()
{
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

  return testResult;
}

bool testBinarySearch()
{
  bool testResult = true;

  static const int N = 10;  // length of example array

  int a[N] = {0,1,2,3,4,5,6,7,8,9};

  using AR = RAPT::rsArray;

  int i;

  // do this in a loop i = 0..N-1:
  i = AR::binarySearch(a, N, 7); 
  testResult &= i == 7;



  return testResult;
}


bool testSortAndSearch()
{
  bool testResult = true;

  testResult &= testHeapSort();
  testResult &= testKnuthMorrisPrattSearch();
  testResult &= testBinarySearch();

  return testResult;
}
