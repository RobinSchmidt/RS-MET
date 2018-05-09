#include "BufferFunctionTests.h"

bool testBufferFunctions(std::string &reportString)
{
  std::string testName = "rsFunctionTemplates";
  bool testResult = true;

  testResult &= testCopySection(reportString);
  //testResult &= testRemoveElements(reportString);
  //testResult &= testMoveElements(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


/*
// can be deleted (is already copied to RSLib):
template<class T>
void rsCopySection(T *source, int sourceLength, T *destination, int copyStart, int copyLength)
{
  int cl, pl1, pl2;  // actual copy-, pre-padding-, post-padding-lengths
  if( copyStart >= 0 )
  {
    // copying:
    cl = rsMin(copyLength, sourceLength-copyStart);
    rsCopyBuffer(&source[copyStart], destination, cl);

    // post-padding:
    pl2 = copyLength-cl;
    rsFillWithZeros(&destination[cl], pl2);
  }
  else
  {
    // pre-padding:
    pl1 = rsMin(-copyStart, copyLength);
    rsFillWithZeros(destination, pl1);

    // copying:
    cl = rsMin(copyLength-pl1, sourceLength);
    rsCopyBuffer(source, &destination[pl1], cl);

    // post-padding:
    pl2 = copyLength-cl-pl1;
    rsFillWithZeros(&destination[pl1+cl], pl2);
  }
}
*/

bool testCopySection(std::string &reportString)
{
  std::string testName = "rsCopySection";
  bool testResult = true;

  static const int Na = 10;            // length of input buffer
  int a[Na] = {1,2,3,4,5,6,7,8,9,10};  // input buffer

  static const int Nb = 20;            // (maximum) length of output buffer
  int b[Nb];
  int n;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, Na, b, 2, 3);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5;
  for(n = 3; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, 2, 3);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5;
  for(n = 3; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, 2, 4);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, -2, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, -2, 7);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2 && b[4] == 3 && b[5] == 4 && 
                b[6] == 5;
  for(n = 7; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, -2, 8);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2 && b[4] == 3 && b[5] == 4 && 
                b[6] == 5 && b[7] == 0;
  for(n = 8; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, -4, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, -3, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 1;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, 4, 4);
  testResult &= b[0] == 5 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  RAPT::rsArray::fillWithValue(b, Nb, -1);
  RAPT::rsArray::copySection(a, 5, b, 5, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

// these two functions are apparently not yet complete:
bool testMoveElements(  std::string &reportString)
{
  std::string testName = "rsMoveElements";
  bool testResult = true;

  static const int length = 10;
  int b[length] = {1,6,7,1,4,5,7,6,1,2};

  RAPT::rsArray::rightShift(b, length, 2);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
bool testRemoveElements(std::string &reportString)
{
  std::string testName = "rsRemoveElements";
  bool testResult = true;

  /*
  static const int length = 10;
  int testBuffer1[length] = {1,6,7,1,4,5,7,6,1,2};
  int testBuffer2[length] = {0,0,0,0,0,0,0,0,0,0};
  int matchBuffer[3]      = {1,6,7};
  int numMatches = rsCopyIfMatching(testBuffer1, testBuffer2, length, matchBuffer, 3);
  int numNonMatches = rsCopyIfNotMatching(testBuffer1, testBuffer2, length, matchBuffer, 3);

  numNonMatches = rsCopyIfNotMatching(testBuffer1, testBuffer1, length, matchBuffer, 3);
   */

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
