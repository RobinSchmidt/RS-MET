#include "BufferFunctionTests.h"



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

bool testCopySection()
{
  bool testResult = true;

  static const int Na = 10;            // length of input buffer
  int a[Na] = {1,2,3,4,5,6,7,8,9,10};  // input buffer

  static const int Nb = 20;            // (maximum) length of output buffer
  int b[Nb];
  int n;

  using AR = RAPT::rsArrayTools;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, Na, b, 2, 3);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5;
  for(n = 3; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, 2, 3);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5;
  for(n = 3; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, 2, 4);
  testResult &= b[0] == 3 && b[1] == 4 && b[2] == 5 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, -2, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, -2, 7);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2 && b[4] == 3 && b[5] == 4 && 
                b[6] == 5;
  for(n = 7; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, -2, 8);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 1 && b[3] == 2 && b[4] == 3 && b[5] == 4 && 
                b[6] == 5 && b[7] == 0;
  for(n = 8; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, -4, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, -3, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 1;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, 4, 4);
  testResult &= b[0] == 5 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  AR::fillWithValue(b, Nb, -1);
  AR::copySection(a, 5, b, 5, 4);
  testResult &= b[0] == 0 && b[1] == 0 && b[2] == 0 && b[3] == 0;
  for(n = 4; n < Nb; n++)
    testResult &= b[n] == -1;

  return testResult;
}


bool testReverse()
{
  bool  r   = true;                 // test result
  using Vec = std::vector<int>;
  using Arr = RAPT::rsArrayTools;
  auto  rev = [&](Vec& v) { Arr::reverse(&v[0], (int)rsSize(v)); }; // reversal function

  Vec a;
  a = {1};         rev(a); r &= a == Vec({1});
  a = {1,2};       rev(a); r &= a == Vec({2,1});
  a = {1,2,3};     rev(a); r &= a == Vec({3,2,1});
  a = {1,2,3,4};   rev(a); r &= a == Vec({4,3,2,1});
  a = {1,2,3,4,5}; rev(a); r &= a == Vec({5,4,3,2,1});

  // tests with array-sizes <= 0 - they should not touch the data at all:
  Arr::reverse(&a[0],  0); r &= a == Vec({5,4,3,2,1});
  Arr::reverse(&a[0], -1); r &= a == Vec({5,4,3,2,1});
  Arr::reverse(&a[0], -2); r &= a == Vec({5,4,3,2,1});
  Arr::reverse(&a[0], -3); r &= a == Vec({5,4,3,2,1});
  Arr::reverse(&a[0], -4); r &= a == Vec({5,4,3,2,1});

  return r;
}


// these two functions are apparently not yet complete:
bool testMoveElements()
{
  bool testResult = true;

  static const int length = 10;
  int b[length] = {1,6,7,1,4,5,7,6,1,2};

  RAPT::rsArrayTools::rightShift(b, length, 2);

  return testResult;
}
bool testRemoveElements()
{
  bool testResult = true;

  typedef RAPT::rsArrayTools AR;
  static const int length = 10;
  int testBuffer1[length] = {1,6,7,1,4,5,7,6,1,2};
  int testBuffer2[length] = {0,0,0,0,0,0,0,0,0,0};
  int matchBuffer[3]      = {1,6,7};
  int numMatches    = AR::copyIfMatching(   testBuffer1, testBuffer2, length, matchBuffer, 3);
  int numNonMatches = AR::copyIfNotMatching(testBuffer1, testBuffer2, length, matchBuffer, 3);
  numNonMatches     = AR::copyIfNotMatching(testBuffer1, testBuffer1, length, matchBuffer, 3);

  // todo: verify that the function did the right thing - after copyIfMatching, the target buffer
  // should be: 1,6,7,1,7,6,1 and after copyIfNotMatching: 4,5,2

  return testResult;
}





bool testBufferFunctions()
{
  bool testResult = true;

  testResult &= testCopySection();
  testResult &= testReverse();
  testResult &= testRemoveElements();
  testResult &= testMoveElements();

  return testResult;
}
