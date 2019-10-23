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

  static const int length = 10;  // length of example array

  int a[length] = {0,1,2,3,4,5,6,7,8,9};

  using AR = RAPT::rsArray;



  for(int subLength = 0; subLength <= length; subLength++) {  // what about sublength == 0?
    for(int searchedValue = 0; searchedValue < subLength; searchedValue++) {
      int foundIndex = AR::binarySearch(a, subLength, searchedValue);
      testResult &= foundIndex == searchedValue; // the indices equal the values in array a
    }
  }
  // what if searchedValue >= any_of(a)


  AR::fillWithRangeLinear(a, length, 0, 18); // 0,2,4,6,...16,18
  for(int subLength = 0; subLength <= length; subLength++) {
    for(int searchedValue = 0; searchedValue < subLength; searchedValue++) {
      int foundIndex = AR::binarySearch(a, subLength, searchedValue);
      testResult &= foundIndex == searchedValue/2;
    }
  }

  // todo: test with value scaler*index + offset for different scalers and offsets
  // maybe factor out the loop into a function that takes an array, length, searchedValue and 
  // (lambda)function that computes the target index from the value (in the first case, that would
  // be the iedentity, in the 2nd, x/2 etc.



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
