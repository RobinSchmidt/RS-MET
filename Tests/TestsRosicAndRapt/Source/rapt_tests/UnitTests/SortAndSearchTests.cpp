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



template<class F1, class F2>
bool testBinSearch(int* array, int length, F1 indexToValue, F2 valueToIndex)
{
  bool testResult = true;
  for(int subLength = 0; subLength <= length; subLength++) 
  {
    for(int index = 0; index < subLength; index++) 
    {
      int searchedValue = indexToValue(index);
      int foundIndex = RAPT::rsArray::binarySearch(array, subLength, searchedValue);
      testResult &= foundIndex == valueToIndex(searchedValue);
    }
  }
  return testResult;
}

template<class T, class F>
void fill(T* a, int N, F indexToValue)
{
  for(int i = 0; i < N; i++)
    a[i] = indexToValue(i);
}
// move to rsArray
// ...has the std library something similar for containers? if so, conform to its syntax

template<class T, class F1, class F2>
T applyComposedFunction(T x, F1 innerFunction, F2 outerFunction)
{
  return outerFunction(innerFunction(x));
}

template<class T, class F>
bool mapsToItself(T x, F f)
{
  return x == f(x);
}

/** Checks, if the 2nd function is the inverse function of the first for the given input argument 
x. */
template<class T, class F1, class F2>
bool mapsBack(T x, F1 forwardFunction, F2 maybeInverseFunction)
{
  return x == applyComposedFunction(x, forwardFunction, maybeInverseFunction);
  //return x == maybeInverseFunction(forwardFunction(x));
}
// maybe rename to isFunctionLocallyInverse

template<class T, class F1, class F2>
bool isInverseFunction(F1 forwardFunc, F2 maybeInverseFunc, T minValue, T maxValue, T increment)
{
  T value = minValue;
  while(value < maxValue) {
    if( !mapsBack(value, forwardFunc, maybeInverseFunc) )
      return false;
    value += increment;
  }
  return true;
}


bool testBinarySearch()
{
  bool testResult = true;

  static const int length = 10;  // length of example array

  int a[length] = {0,1,2,3,4,5,6,7,8,9};

  using AR = RAPT::rsArray;

  /*
  // these two loops are now redundant with the calls to testBinSearch below:
  for(int subLength = 0; subLength <= length; subLength++) {
    for(int searchedValue = 0; searchedValue < subLength; searchedValue++) {
      int foundIndex = AR::binarySearch(a, subLength, searchedValue);
      testResult &= foundIndex == searchedValue; // the indices equal the values in array a
    }
  }
  AR::fillWithRangeLinear(a, length, 0, 18); // 0,2,4,6,...16,18
  for(int subLength = 0; subLength <= length; subLength++) {
    for(int searchedValue = 0; searchedValue < subLength; searchedValue++) {
      int foundIndex = AR::binarySearch(a, subLength, searchedValue);
      testResult &= foundIndex == searchedValue/2;
    }
  }
  */

  // todo: test with value scaler*index + offset for different scalers and offsets
  // maybe factor out the loop into a function that takes an array, length, searchedValue and 
  // (lambda)function that computes the target index from the value (in the first case, that would
  // be the iedentity, in the 2nd, x/2 etc.

  auto identity = [](int i){ return i; };  // converts indices to values via the identity function...
  rsAssert(isInverseFunction(identity, identity, 0, length-1, 1));
  fill(a, length, identity);
  testResult &= testBinSearch(a, length, identity, identity);

  auto timesTwo = [](int i){ return 2*i; };
  auto divByTwo = [](int i){ return i/2; };
  rsAssert(isInverseFunction(timesTwo, divByTwo, 0, length-1, 1));
  fill(a, length, timesTwo);
  testResult &= testBinSearch(a, length, timesTwo, divByTwo);

  // the rsAssert(isInverseFunction...) may be commented out - they are onyl needed for verifying 
  // *once*, that the index-to-value and value-to-index functions actually are inverses of each 
  // other as intended. they are not part of the unit test for the search-algo -c they are used to 
  // test, that we actually generate proper test-inputs



  // what, if we exchange the roles?
  //testResult &= testBinSearch(a, length, divByTwo, timesTwo);
  // ...it fails, because in the integers, multiplying by 2 does not undo division by 2, 
  // for example 3 != (3/2)*2 = 2
  //rsAssert(isInverseFunction(divByTwo, timesTwo, 0, length-1, 1));


  auto f1 = [](int i){ return 3*i+1; };    // give better names ...times3plus1
  auto f2 = [](int i){ return (i-1)/3; };  // ...minus1divBy3
  //bool isInverse = isInverseFunction(f1, f2, 0, length-1, 1);
  rsAssert(isInverseFunction(f1, f2, 0, length-1, 1));
  fill(a, length, f1);
  testResult &= testBinSearch(a, length, f1, f2);

  // find invertible functions in the integers that do interesting things

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
