#include "SortAndSearchTests.h"

bool testHeapSort()
{
  bool testResult = true;

  int numTests = 20;
  static const int length = 100;
  int testArray[length];
  for(int i=0; i<numTests; i++)
  {
    RAPT::rsArrayTools::fillWithRandomValues(testArray, length, -100, +100, 1);
    rsHeapSort(testArray, length);
    testResult &= rsArrayTools::isSortedAscending(testArray, length);

    // check odd lengths by just sorting the subarray up to length-1:
    RAPT::rsArrayTools::fillWithRandomValues(testArray, length, -100, +100, 1);
    rsHeapSort(testArray, length-1);
    testResult &= rsArrayTools::isSortedAscending(testArray, length-1);
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


template<class F1, class F2>
bool testSplitIndex(int* array, int length, F1 indexToValue, F2 valueToIndex)
{
  bool testResult = true;
  for(int subLength = 0; subLength <= length; subLength++) 
  {
    for(int index = 0; index < subLength; index++) 
    {
      int searchedValue = indexToValue(index);
      int foundIndex = RAPT::rsArrayTools::findSplitIndex(array, subLength, searchedValue);
      testResult &= foundIndex == valueToIndex(searchedValue);
    }
  }
  return testResult;
}
bool testSplitIndex()
{
  bool testResult = true;

  static const int length = 10;  // length of example array

  int a[length] = {0,1,2,3,4,5,6,7,8,9};

  using AR = RAPT::rsArrayTools;

  auto identity = [](int i){ return i; };  // converts indices to values via the identity function...
  //rsAssert(isInverseFunction(identity, identity, 0, length-1, 1));
  AR::fill(a, length, identity);
  testResult &= testSplitIndex(a, length, identity, identity);

  auto timesTwo = [](int i){ return 2*i; };
  auto divByTwo = [](int i){ return i/2; };
  //rsAssert(isInverseFunction(timesTwo, divByTwo, 0, length-1, 1));
  AR::fill(a, length, timesTwo);
  testResult &= testSplitIndex(a, length, timesTwo, divByTwo);

  auto times3plus1 = [](int i){ return 3*i+1; }; 
  auto minus1divBy3 = [](int i){ return (i-1)/3; };
  //rsAssert(isInverseFunction(times3plus1, minus1divBy3, 0, length-1, 1));
  AR::fill(a, length, times3plus1);
  testResult &= testSplitIndex(a, length, times3plus1, minus1divBy3);

  int foundIndex;
  std::vector<int> b = {1,2,2,2,4,4,4,4,5,6};
  foundIndex  = AR::findSplitIndex(&b[0], (int)b.size(), 4);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndex(&b[0], (int)b.size(), 3);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndex(&b[0], (int)b.size(), 2);
  testResult &= foundIndex == 1;

  std::vector<double> c = {1.,2.,2.,2.,4.,4.,4.,4.,5.,6.};
  foundIndex  = AR::findSplitIndex(&c[0], (int)c.size(), 3.9);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndex(&c[0], (int)c.size(), 4.0);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndex(&c[0], (int)c.size(), 4.1);
  testResult &= foundIndex == 8;
  foundIndex  = AR::findSplitIndex(&c[0], (int)c.size(), 2.1);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndexClosest(&c[0], (int)c.size(), 3.9);
  testResult &= foundIndex == 4;
  foundIndex  = AR::findSplitIndexClosest(&c[0], (int)c.size(), 2.1);
  testResult &= foundIndex == 3;

  // the rsAssert(isInverseFunction...) are commented out because they were only needed for 
  // verifying *once*, that the index-to-value and value-to-index functions actually are inverses 
  // of each other as intended. they are not part of the unit test for the search-algo - they are
  // used to test, that we actually generate proper test-inputs

  return testResult;
}

bool testSortAndSearch()
{
  bool testResult = true;

  testResult &= testHeapSort();
  testResult &= testKnuthMorrisPrattSearch();
  testResult &= testSplitIndex();

  return testResult;
}
