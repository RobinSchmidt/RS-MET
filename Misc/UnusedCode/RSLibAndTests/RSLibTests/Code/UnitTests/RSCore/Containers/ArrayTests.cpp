#include "ArrayTests.h"

bool testArray(std::string &reportString)
{
  std::string testName = "rsArray";
  bool testResult = true;

  testResult &= testArrayAppend(reportString);
  testResult &= testArrayInsert(reportString);
  testResult &= testArrayRemove(reportString);
  testResult &= testArrayGrowAndShrink(reportString);
  testResult &= testArrayMisc(reportString);

  testResult &= testFlagArray(reportString);

  testResult &= testInfiniteDataStream(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testArrayAppend(std::string &reportString)
{
  std::string testName = "rsArrayAppend";
  bool testResult = true;

  rsArray<int> array1;
  array1.appendElement(1); // [1]
  array1.appendElement(2); // [1 2]
  array1.appendElement(3); // [1 2 3]
  array1.appendElement(7); // [1 2 3 7]
  array1.appendElement(3); // [1 2 3 7 3]
  array1.appendElement(5); // [1 2 3 7 3 5]
  array1.appendElement(2); // [1 2 3 7 3 5 2]
  array1.appendElement(6); // [1 2 3 7 3 5 2 6]
  array1.appendElement(9); // [1 2 3 7 3 5 2 6 9]
  array1.appendElement(4); // [1 2 3 7 3 5 2 6 9 4], 10 elements

  testResult &= ( array1.getNumElements() == 10 );
  testResult &= ( array1.findElement(7)   == 3  );
  testResult &= ( array1.findElement(3)   == 2  );

  array1.appendIfNotAlreadyThere(7);
  testResult &= ( array1.getNumElements() == 10 );

  rsArray<int> array2;
  array2.appendElement(1); // [1]
  array2.appendElement(2); // [1 2]
  array2.appendElement(3); // [1 2 3]
  array2.appendElement(7); // [1 2 3 7]
  array2.appendElement(3); // [1 2 3 7 3]
  array2.appendElement(5); // [1 2 3 7 3 5]

  testResult &= ( array1 != array2 );

  array2.appendElement(2); // [1 2 3 7 3 5 2]
  array2.appendElement(6); // [1 2 3 7 3 5 2 6]
  array2.appendElement(9); // [1 2 3 7 3 5 2 6 9]
  array2.appendElement(4); // [1 2 3 7 3 5 2 6 9 4], 10 elements

  testResult &= ( array1 == array2 );

  testResult &= ( array1.hasElement(7) );
  array1.removeElementByValue(7);
  testResult &= ( !array1.hasElement(7) );

  rsArray<int> array3;
  array3.appendArray(array2);
  array3.appendArray(array1);
  testResult &= ( array3.getNumElements() == array2.getNumElements() + array1.getNumElements() );


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testArrayInsert(std::string &reportString)
{
  return true;
}

bool testArrayRemove(std::string &reportString)
{
  std::string testName = "rsArrayRemove";
  bool testResult = true;

  int cArrayInt[10] = {1,3,3,4,5,5,7,8,9,4};
  rsArray<int> completeArray(cArrayInt, 10);

  // create array of elements to remove or keep:
  int cSomeElements[2] = {7,3};
  int cSomeRemoved[7]  = {1,4,5,5,8,9,4};
  int cSomeKept[3]     = {3,3,7};
  rsArray<int> someElements(cSomeElements, 2);

  // remove and check:
  rsArray<int> someRemovedTarget(cSomeRemoved, 7);
  rsArray<int> someRemovedActual = completeArray;
  someRemovedActual.removeMatchingElements(someElements);
  testResult &= ( someRemovedActual == someRemovedTarget );

  // keep and check:
  rsArray<int> someKeptTarget(cSomeKept, 3);
  rsArray<int> someKeptActual = completeArray;
  someKeptActual.keepOnlyMatchingElements(someElements);
  testResult &= ( someKeptActual == someKeptTarget );

  // remove range and check:
  int cRangeRemoved[6] = {1,3,7,8,9,4};
  rsArray<int> rangeRemovedTarget(cRangeRemoved, 6);
  rsArray<int> rangeRemovedActual = completeArray;
  rangeRemovedActual.removeRange(2, 5);
  testResult &= ( rangeRemovedActual == rangeRemovedTarget );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testArrayGrowAndShrink(std::string &reportString)
{
  std::string testName = "rsArrayGrowAndShrink";
  bool testResult = true;

  rsArray<int> testArray;
  testResult &= ( testArray.getNumAllocatedElements() == 0 );
  testArray.appendElement(1);
  testResult &= ( testArray.getNumAllocatedElements() == 1 );
  testArray.appendElement(2);
  testResult &= ( testArray.getNumAllocatedElements() == 2 );
  testArray.appendElement(3);
  testResult &= ( testArray.getNumAllocatedElements() == 4 );
  testArray.appendElement(4);
  testResult &= ( testArray.getNumAllocatedElements() == 4 );
  testArray.appendElement(5);
  testResult &= ( testArray.getNumAllocatedElements() == 8 );
  testArray.removeElementByValue(5);
  testResult &= ( testArray.getNumAllocatedElements() == 4 );
  testArray.removeElementByValue(4);
  testResult &= ( testArray.getNumAllocatedElements() == 4 );
  testArray.removeElementByValue(3);
  testResult &= ( testArray.getNumAllocatedElements() == 2 );
  testArray.removeElementByValue(2);
  testResult &= ( testArray.getNumAllocatedElements() == 1 );
  testArray.removeElementByValue(1);
  testResult &= ( testArray.getNumAllocatedElements() == 1 );  // it does not shrink to entirely zero
  testArray.appendElement(1);
  testResult &= ( testArray.getNumAllocatedElements() == 1 );
  testArray.appendElement(2);
  testResult &= ( testArray.getNumAllocatedElements() == 2 );
  testArray.clear();
  testResult &= ( testArray.getNumAllocatedElements() == 0 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testArrayMisc(std::string &reportString)
{
  std::string testName = "rsArrayMisc";
  bool testResult = true;

  int cArray[10] = {0,1,2,3,4,5,6,7,8,9};
  rsArray<int> testArray(cArray, 10);

  rsArray<int> subArray = testArray.getSubArray(3, 6);
  int cTargetSubArray[4] = {3,4,5,6};
  rsArray<int> targetSubArray(cTargetSubArray, 4);
  testResult &= ( subArray == targetSubArray );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFlagArray(std::string &reportString)
{
  std::string testName = "rsFlagArray";
  bool testResult = true;

  rsFlagArray a(1000);
  testResult &= a.areAllFlagsFalse();
  a.setAllTrue();  
  testResult &= a.areAllFlagsTrue();
  testResult &= a.getNumTrueFlags() == 1000;
  a.setFlagFalse(987);
  a.setFlagFalse(463);
  a.setFlagFalse(123);
  testResult &= a.getNumTrueFlags() == 997;
  a.setFlagTrue(123);
  testResult &= a.getNumTrueFlags() == 998;
  a.setFlagTrue(987);
  testResult &= a.getNumTrueFlags() == 999;
  testResult &= a.isFlagTrue(123);
  testResult &= a.isFlagFalse(463);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testInfiniteDataStream(std::string &reportString)
{
  std::string testName = "rsInfiniteDataStream";
  bool testResult = true;

                  // 0  1  2  3  4  5  6  7 - indices
  int testData[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  static const int bN = 20; // length of buffer b
  int b[20];                // output buffer

  rsInfiniteDataStream<int> stream(testData, 8);
  //stream.setInputDataAddress(testData, 8);

  // check single value access via get-function:
  testResult &= stream.getValue(-1) == 0;
  testResult &= stream.getValue( 0) == 1;
  testResult &= stream.getValue( 7) == 8;
  testResult &= stream.getValue( 8) == 0;

  // check single value read access via index operator:
  testResult &= stream[-1] == 0;
  testResult &= stream[0]  == 1;
  testResult &= stream[7]  == 8;
  testResult &= stream[8]  == 0;

  // check single value write access via index operator:
  stream[-1] = 10; testResult &= stream[-1] ==  0; // out-of-bounds write should have no effect
  stream[8]  = 10; testResult &= stream[8]  ==  0; // same here
  stream[0]  = 10; testResult &= stream[0]  == 10; // within-bounds write should have effect
  stream[0]  = 1;                                  // restore original value
  
  // check very small output buffer-sizes (1, 0):
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, -2, 1);
  testResult &= b[0] == 0;
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);
  stream.getBuffer(b,  5, 1);
  testResult &= b[0] == 6;
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);
  stream.getBuffer(b,  9, 1);
  testResult &= b[0] == 0;
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, -2, 0);
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);
  stream.getBuffer(b,  5, 0);
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);
  stream.getBuffer(b,  9, 0);
  testResult &= rsIsFilledWithValue(&b[1], bN-1, 1000);  
  
  // check pre-padding:
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, -2, 5);
  testResult &= b[0] == 0;
  testResult &= b[1] == 0;
  testResult &= b[2] == 1;
  testResult &= b[3] == 2;
  testResult &= b[4] == 3;
  testResult &= rsIsFilledWithValue(&b[5], bN-5, 1000);

  // check post-padding:
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, 5, 5);
  testResult &= b[0] == 6;
  testResult &= b[1] == 7;
  testResult &= b[2] == 8;
  testResult &= b[3] == 0;
  testResult &= b[4] == 0;
  testResult &= rsIsFilledWithValue(&b[5], bN-5, 1000);

  // check pre- and post-padding:
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, -2, 15);
  testResult &= b[0]  == 0;
  testResult &= b[1]  == 0;
  testResult &= b[2]  == 1;
  testResult &= b[3]  == 2;
  testResult &= b[4]  == 3;
  testResult &= b[5]  == 4;
  testResult &= b[6]  == 5;
  testResult &= b[7]  == 6;
  testResult &= b[8]  == 7;
  testResult &= b[9]  == 8;
  testResult &= b[10] == 0;
  testResult &= b[11] == 0;
  testResult &= b[12] == 0;
  testResult &= b[13] == 0;
  testResult &= b[14] == 0;
  testResult &= rsIsFilledWithValue(&b[15], bN-15, 1000);

   // check no padding:
  rsFillWithValue(b, 20, 1000);
  stream.getBuffer(b, 2, 3);
  testResult &= b[0] == 3;
  testResult &= b[1] == 4;
  testResult &= b[2] == 5;
  testResult &= rsIsFilledWithValue(&b[3], bN-3, 1000);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

