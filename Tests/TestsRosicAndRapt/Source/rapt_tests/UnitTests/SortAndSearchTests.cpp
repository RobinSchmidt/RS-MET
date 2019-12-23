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
bool testFindSplitIndex(int* array, int length, F1 indexToValue, F2 valueToIndex)
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
bool testFindSplitIndex()
{
  bool r = true;   // test result

  static const int length = 10;  // length of example array

  int a[length] = {0,1,2,3,4,5,6,7,8,9};

  using AR = RAPT::rsArrayTools;

  auto identity = [](int i){ return i; };  // converts indices to values via the identity function...
  //rsAssert(isInverseFunction(identity, identity, 0, length-1, 1));
  AR::fill(a, length, identity);
  r &= testFindSplitIndex(a, length, identity, identity);

  auto timesTwo = [](int i){ return 2*i; };
  auto divByTwo = [](int i){ return i/2; };
  //rsAssert(isInverseFunction(timesTwo, divByTwo, 0, length-1, 1));
  AR::fill(a, length, timesTwo);
  r &= testFindSplitIndex(a, length, timesTwo, divByTwo);

  auto times3plus1 = [](int i){ return 3*i+1; }; 
  auto minus1divBy3 = [](int i){ return (i-1)/3; };
  //rsAssert(isInverseFunction(times3plus1, minus1divBy3, 0, length-1, 1));
  AR::fill(a, length, times3plus1);
  r &= testFindSplitIndex(a, length, times3plus1, minus1divBy3);

  int i;
  std::vector<int> b = {1,2,2,2,4,4,4,4,5,6};
  i = AR::findSplitIndex(&b[0], (int)b.size(), 0); r &= i == 0;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 1); r &= i == 0;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 2); r &= i == 1;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 3); r &= i == 4;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 4); r &= i == 4;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 5); r &= i == 8;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 6); r &= i == 9;
  i = AR::findSplitIndex(&b[0], (int)b.size(), 7); r &= i == 10;

  // test with floating point numbers:
  std::vector<double> c = {1.,2.,2.,2.,4.,4.,4.,4.,5.,6.};
  i = AR::findSplitIndex(&c[0], (int)c.size(), 3.9); r &= i == 4;
  i = AR::findSplitIndex(&c[0], (int)c.size(), 4.0); r &= i == 4;
  i = AR::findSplitIndex(&c[0], (int)c.size(), 4.1); r &= i == 8;
  i = AR::findSplitIndex(&c[0], (int)c.size(), 2.1); r &= i == 4;
  i = AR::findSplitIndexClosest(&c[0], (int)c.size(), 3.9); r &= i == 4;
  i = AR::findSplitIndexClosest(&c[0], (int)c.size(), 2.1); r &= i == 3;
  i = AR::findSplitIndexClosest(&c[0], (int)c.size(), 6.0); r &= i == 9;
  i = AR::findSplitIndexClosest(&c[0], (int)c.size(), 7.0); r &= i == 9; //



  // the rsAssert(isInverseFunction...) are commented out because they were only needed for 
  // verifying *once*, that the index-to-value and value-to-index functions actually are inverses 
  // of each other as intended. they are not part of the unit test for the search-algo - they are
  // used to test, that we actually generate proper test-inputs

  return r;
}

template<class T>
T findFloatPosition(const T* a, const int N, const T v)
{
  rsAssert(N >= 2, "array must have at least 2 elements");


  int i = RAPT::rsArrayTools::findSplitIndex(a, N, v);

  // If all values in a are less than v, extrapolate rightward:
  if(i == N)
  {
    return rsInterpolateLinear(a[i-2], a[i-1], T(i-2), T(i-1), v); // optimize!
    //...extrapolate and return early
  }

  // If all values in a are greater than v, extrapolate leftward:
  if(i == 0 && a[0] > v)
  {
    return rsInterpolateLinear(a[0], a[1], T(0), T(1), v); // optimize!
    //...extrapolate and return early
  }

  // If there is a run of same values v following the found index, we need to take the middle of
  // that run - find upper index of such a run:
  if(a[i] == v)
  {
    int iu = i;
    while(iu < N-1 && a[iu+1] == v)
      iu++;
    return T(0.5) * T(i+iu);
  }

  // If a[i] > v, we need to interpolate the position:
  if(a[i] > v)
  {
    // we can assume here that i > 0 because i==0 && a[i] > v was handled above

  }


  // ...

  // If 0 is returned


  T iFlt = T(i);  // preliminary




  return iFlt;
}

bool testFindFloatPosition()
{
  bool r = true;   // test result

  double a[9] = {1,2,3,4,5,6,7,8,9};  // array
  double p;                           // position of value

  p = findFloatPosition(a, 9, -1.0); r &= p == -2.0;
  p = findFloatPosition(a, 9,  0.0); r &= p == -1.0;
  p = findFloatPosition(a, 9,  1.0); r &= p ==  0.0;
  p = findFloatPosition(a, 9,  2.0); r &= p ==  1.0;
  p = findFloatPosition(a, 9,  5.0); r &= p ==  4.0;
  p = findFloatPosition(a, 9,  8.0); r &= p ==  7.0;
  p = findFloatPosition(a, 9,  9.0); r &= p ==  8.0;
  p = findFloatPosition(a, 9, 10.0); r &= p ==  9.0;
  p = findFloatPosition(a, 9, 11.0); r &= p == 10.0;

  // a = {1,2,3,4,5,5,7,8,9}
  a[5] = 5.0; p = findFloatPosition(a, 9, 5.0);  r &= p == 4.5;

  // a = {1,2,3,4,5,5,5,8,9}
  a[6] = 5.0; p = findFloatPosition(a, 9, 5.0);  r &= p == 5.0;





  return r;
}


bool testSortAndSearch()
{
  bool testResult = true;

  testResult &= testHeapSort();
  testResult &= testKnuthMorrisPrattSearch();
  testResult &= testFindSplitIndex();
  testResult &= testFindFloatPosition();

  return testResult;
}
