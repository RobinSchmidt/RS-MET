#include "ComplexTests.h"

// \todo rename to BigNumberTests

bool testArbitraryPrecision(std::string &reportString)
{
  std::string testName = "rsArbitraryPrecision";
  bool testResult = true;

  testResult &= testDigitArithmetic(reportString);
  testResult &= testBaseChange(reportString);
  testResult &= testBigInt(reportString);
  testResult &= testBigFloat(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


bool testDigitDivision2()
{
  static const int N = 2; // number of digits
  rsUint32 a[N];           // dividend
  rsUint32 b[N];           // divisor
  rsUint32 q[N];           // a / b (quotient)
  rsUint32 r[N];           // a % b (remainder);

  bool result = true;

  // 87 / 5 = 17, r = 2
  rsBigNumber::digitsFromInt(87, a, N);
  rsBigNumber::digitsFromInt( 5, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==1 && q[1]==7;
  result &= r[0]==0 && r[1]==2;

  // 99 / 1 = 99, r = 0
  rsBigNumber::digitsFromInt(99, a, N);
  rsBigNumber::digitsFromInt( 1, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==9 && q[1]==9;
  result &= r[0]==0 && r[1]==0;

  // 99 / 5 = 19, r = 4
  rsBigNumber::digitsFromInt(99, a, N);
  rsBigNumber::digitsFromInt( 5, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==1 && q[1]==9;
  result &= r[0]==0 && r[1]==4;

  // 99 / 23 = 4, r = 7
  rsBigNumber::digitsFromInt(99, a, N);
  rsBigNumber::digitsFromInt(23, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==0 && q[1]==4;
  result &= r[0]==0 && r[1]==7;

  return result;
}
bool testDigitDivision5()
{
  static const int N = 5; // number of digits
  rsUint32 a[N];           // dividend
  rsUint32 b[N];           // divisor
  rsUint32 q[N];           // a / b (quotient)
  rsUint32 r[N];           // a % b (remainder);

  bool result = true;

  // 04537 / 00003 = 1512, r = 1
  rsBigNumber::digitsFromInt( 4537, a, N);
  rsBigNumber::digitsFromInt(    3, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==0 && q[1]==1 && q[2]==5 && q[3]==1 && q[4]==2;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==1;

  // 99999 / 1 = 99999, r = 0
  rsBigNumber::digitsFromInt(99999, a, N);
  rsBigNumber::digitsFromInt(    1, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==9 && q[1]==9 && q[2]==9 && q[3]==9 && q[4]==9;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==0;

  // 99999 / 2 = 49999, r = 1
  rsBigNumber::digitsFromInt(99999, a, N);
  rsBigNumber::digitsFromInt(    2, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==4 && q[1]==9 && q[2]==9 && q[3]==9 && q[4]==9;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==1;

  // 99999 / 9 = 11111, r = 0
  rsBigNumber::digitsFromInt(99999, a, N);
  rsBigNumber::digitsFromInt(    9, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==1 && q[1]==1 && q[2]==1 && q[3]==1 && q[4]==1;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==0;

  // 07000 / 00029 = 241, r = 11
  rsBigNumber::digitsFromInt( 7000, a, N);
  rsBigNumber::digitsFromInt(   29, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==0 && q[1]==0 && q[2]==2 && q[3]==4 && q[4]==1;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==1 && r[4]==1;

  // 03000 / 00099 = 30, r = 30
  rsBigNumber::digitsFromInt( 3000, a, N);
  rsBigNumber::digitsFromInt(   99, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==0 && q[1]==0 && q[2]==0 && q[3]==3 && q[4]==0;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==3 && r[4]==0;

  // 03142 / 00017 = 184, r = 14
  rsBigNumber::digitsFromInt( 3142, a, N);
  rsBigNumber::digitsFromInt(   17, b, N);
  rsBigNumber::divide(a, b, N, q, r);
  result &= q[0]==0 && q[1]==0 && q[2]==1 && q[3]==8 && q[4]==4;
  result &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==1 && r[4]==4;

  return result;
}
bool exhaustiveIntegerDivisionTest(int numDigits)
{
  // divides all numbers with the given (maximum) number of digits) by all other such numbers and
  // checks the results

  static const int maxNumDigits = 4;
  rsAssert( numDigits >= 1 && numDigits <= maxNumDigits );

  int upperLimit = 9;
  for(int k = 1; k < numDigits; k++)
    upperLimit = 10 * upperLimit + 9;

  rsUint32 q[maxNumDigits];  // dividend and quotient (we divide in-place)
  rsUint32 r[maxNumDigits];  // divisor and remainder
  rsUint64 qi, ri, qit, rit; // quotient and remainder and target avlue as int

  bool result = true;

  for(int i = 0; i <= upperLimit; i++)
  {
    for(int j = 1; j <= upperLimit; j++) // start at 1 to avoid division by zero
    {
      rsBigNumber::digitsFromInt(i, q, numDigits);
      rsBigNumber::digitsFromInt(j, r, numDigits);
      rsBigNumber::divide(       q, r, numDigits, q, r);
      qi  = rsBigNumber::digitsToInt(q, numDigits);
      ri  = rsBigNumber::digitsToInt(r, numDigits);
      qit = i / j;
      rit = i % j;
      result &= qi == qit;
      result &= ri == rit;
      rsAssert(result == true);
    }
  }

  return result;
}

double roundLikeInDigitArithmetic(double x, rsUint64 base = 10)
{
  double xi  = floor(x);
  double xf  = x - xi;
  double tol = 1.0e-13;
  if( xf < 0.5 - tol )
    return xi;
  else if( xf > 0.5 + tol )
    return xi+1;
  else
  {
    if( rsBigNumber::shouldRoundUpAtMidway((rsUint32) xi, base) )
      return xi+1;
    else
      return xi;
  }
}
void getDesiredRoundedSequence(double x, rsUint32 *d, rsUint32 numDigits)
{
  double tmpF = 0.0;
  double tmpQ = x;
  for(rsUint32 i = 0; i < numDigits; i++)
  {
    if( i == numDigits-1 )
    {
      tmpF = roundLikeInDigitArithmetic(tmpQ);
      //rsAssert( tmpF < 10.0 );
      // todo: if tmpF == 10.0, we need to perfom a carry operation on the desired digits
    }
    else
      tmpF = floor(tmpQ);
    tmpQ = 10*(tmpQ - tmpF);
    d[i] = (rsUint32) tmpF;
  }
}
bool exhaustiveFractionDivisionTest(int numDigits)
{
  static const int maxNumDigits = 4;
  rsAssert( numDigits >= 1 && numDigits <= maxNumDigits );
  int upperLimit = 9;
  for(int k = 1; k < numDigits; k++)
    upperLimit = 10 * upperLimit + 9;


  int a, b;                    // operands as integers
  rsUint32 q[2*maxNumDigits];  // dividend and quotient (we divide in-place)
  rsUint32 r[2*maxNumDigits];  // divisor and remainder
  rsUint32 d[2*maxNumDigits];  // desired rounded digit sequence
  double qd;                   // quotient (times 10^numDigits) as double
  //double factor = pow(10, numDigits);

  //rsUint64 qi, ri, qit, rit;   // quotient and remainder and target avlue as int

  bool result = true;

  for(a = 0; a <= upperLimit; a++)
  {
    for(b = 1; b <= upperLimit; b++) // start at 1 to avoid division by zero
    {
      rsBigNumber::digitsFromInt(a, q, numDigits);
      rsBigNumber::digitsFromInt(b, r, numDigits);
      rsBigNumber::divideFractions(q, r, q, numDigits);
      qd = ((double) a / (double) b); //  /  pow(10.0, numDigits-1);

      //if( a == 2 && b == 21 )
      //  int dummy = 0;

      getDesiredRoundedSequence(qd / pow(10.0, numDigits-1), d, 2*numDigits);

      if( !rsContains(d, 2*numDigits, (rsUint32)10) )
      {
        // the condition above is necesarry because getDesiredRoundedSequence sometimes may
        // put a 10 into the array, where there should be a zero and carrying into the preceding
        // digits - the function does not implment the carry, so we must leave some gaps in the
        // tests - todo: fix this
        for(int i = 0; i < 2*numDigits; i++)
          result &= q[i] == d[i];
      }
      rsAssert( result == true );
    }
  }

  return result;
}

bool testDigitArithmetic(std::string &reportString)
{
  std::string testName = "rsDigitArithmetic";
  bool testResult = true;

  static const int N = 4;       // number of digits
  //int base      = 10;           // base for the number system
  rsUint32 a[N] = {9, 6, 4, 8}; // left operand
  rsUint32 b[N] = {3, 1, 5, 7}; // right operand
  rsUint32 s[N+1];              // a + b (sum)
  rsUint32 d[N];                // a - b (difference)
  rsUint32 p[2*N];              // a * b (product)
  rsUint32 q[20];               // a / b (quotient)
  rsUint32 r[N];                // a % b (remainder);

  // test conversion to/from int:
  testResult &= rsBigNumber::digitsToInt(a, N) == 9648;
  rsBigNumber::digitsFromInt(3472, r, N);
  testResult &= r[0]==3 && r[1]==4 && r[2]==7 && r[3]==2;

  // 9648 + 3157 = 12 805
  s[0] = rsBigNumber::add(a, b, N, &s[1]);
  testResult &= s[0]==1 && s[1]==2 && s[2]==8 && s[3]==0 && s[4]==5;

  // 9648 * 3157 = 30 458 736
  //rsDigitMul(a, b, N, p);
  rsBigNumber::multiply(a, N, b, N, p);
  testResult &= p[0]==3 && p[1]==0 && p[2]==4 && p[3]==5 && p[4]==8 && p[5]==7 && p[6]==3 &&
                p[7]==6;

  // 964 * 31 = 29 884
  rsBigNumber::multiply(a, 3, b, 2, p);
  testResult &= p[0]==2 && p[1]==9 && p[2]==8 && p[3]==8 && p[4]==4;

  // 9648 - 3157 = 6491 (one borrow):
  testResult &= rsBigNumber::subtract(a, b, N, d) == 0;
  testResult &= d[0]==6 && d[1]==4 && d[2]==9 && d[3]==1;

  // 9035 - 7268 = 1767 (one borrow from zero)
  a[0]=9; a[1]=0; a[2]=3; a[3]=5;
  b[0]=7; b[1]=2; b[2]=6; b[3]=8;
  testResult &= rsBigNumber::subtract(a, b, N, d) == 0;
  testResult &= d[0]==1 && d[1]==7 && d[2]==6 && d[3]==7;

  // 9035 - 7 = 9028
  testResult &= rsBigNumber::subtract(a, 7, N, d) == 0;
  testResult &= d[0]==9 && d[1]==0 && d[2]==2 && d[3]==8;

  // 22 - 46 is negative - the result is invalid
  a[0]=2; a[1]=2;
  b[0]=4; b[1]=6;
  testResult &= rsBigNumber::subtract(a, b, 2, d) == 1;
  testResult &= d[0]==7 && d[1]==6;

  // 2 / 1 = 2, remainder = 0
  a[0] = 2;
  b[0] = 1;
  rsBigNumber::divide(a, b, 1, q, r, 10);
  testResult &= q[0] == 2;
  testResult &= r[0] == 0;

  // 2 / 2 = 1, remainder = 0
  a[0] = 2;
  b[0] = 2;
  rsBigNumber::divide(a, b, 1, q, r, 10);
  testResult &= q[0] == 1;
  testResult &= r[0] == 0;

  // 8 / 4 = 2, remainder = 0
  a[0] = 8;
  b[0] = 4;
  rsBigNumber::divide(a, b, 1, q, r, 10);
  testResult &= q[0] == 2;
  testResult &= r[0] == 0;

  // 20 / 10 = 2, remainder = 0, in base 21
  a[0] = 20;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 21);
  testResult &= q[0] == 2;
  testResult &= r[0] == 0;

  // 3142 / 617 = 5, remainder = 57
  a[0]=3; a[1]=1; a[2]=4; a[3]=2;
  b[0]=0; b[1]=6; b[2]=1; b[3]=7;
  rsBigNumber::divide(a, b, N, q, r);
  testResult &= q[0]==0 && q[1]==0 && q[2]==0 && q[3]==5;
  testResult &= r[0]==0 && r[1]==0 && r[2]==5 && r[3]==7;

  // test this case for in-place division such that the pointers a and b are used for operands and
  // results (in both possible orders):
  a[0]=3; a[1]=1; a[2]=4; a[3]=2;
  b[0]=0; b[1]=6; b[2]=1; b[3]=7;
  rsBigNumber::divide(a, b, N, a, b);
  testResult &= a[0]==0 && a[1]==0 && a[2]==0 && a[3]==5;
  testResult &= b[0]==0 && b[1]==0 && b[2]==5 && b[3]==7;
  a[0]=3; a[1]=1; a[2]=4; a[3]=2;
  b[0]=0; b[1]=6; b[2]=1; b[3]=7;
  rsBigNumber::divide(a, b, N, b, a);
  testResult &= b[0]==0 && b[1]==0 && b[2]==0 && b[3]==5;
  testResult &= a[0]==0 && a[1]==0 && a[2]==5 && a[3]==7;

  // 200 / 101 = 1, remainder = 99
  a[0]=0; a[1]=2; a[2]=0; a[3]=0;
  b[0]=0; b[1]=1; b[2]=0; b[3]=1;
  rsBigNumber::divide(a, b, N, q, r);
  testResult &= q[0]==0 && q[1]==0 && q[2]==0 && q[3]==1;
  testResult &= r[0]==0 && r[1]==0 && r[2]==9 && r[3]==9;

  // 1 / 100 = 0, remainder = 1
  a[0]=0; a[1]=0; a[2]=0; a[3]=1;
  b[0]=0; b[1]=1; b[2]=0; b[3]=0;
  rsBigNumber::divide(a, b, N, q, r);
  testResult &= q[0]==0 && q[1]==0 && q[2]==0 && q[3]==0;
  testResult &= r[0]==0 && r[1]==0 && r[2]==0 && r[3]==1;

  // 3142 / 17 = 184, remainder = 14
  a[0]=3; a[1]=1; a[2]=4; a[3]=2;
  b[0]=0; b[1]=0; b[2]=1; b[3]=7;
  rsBigNumber::divide(a, b, N, q, r);
  testResult &= q[0]==0 && q[1]==1 && q[2]==8 && q[3]==4;
  testResult &= r[0]==0 && r[1]==0 && r[2]==1 && r[3]==4;

  // 8596 / 1347 = 6, remainder=514
  a[0]=8; a[1]=5; a[2]=9; a[3]=6;
  b[0]=1; b[1]=3; b[2]=4; b[3]=7;
  rsBigNumber::divide(a, b, N, q, r);
  testResult &= q[0]==0 && q[1]==0 && q[2]==0 && q[3]==6;
  testResult &= r[0]==0 && r[1]==5 && r[2]==1 && r[3]==4;

  // test digit division in base 2^8:
  a[0] = 20;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 256);
  testResult &= q[0] == 2;
  testResult &= r[0] == 0;

  // test digit division in base 2^32:
  a[0] = 21;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 4294967296);
  testResult &= q[0] == 2;
  testResult &= r[0] == 1;
  a[0] = 10;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 4294967296);
  testResult &= q[0] == 1;
  testResult &= r[0] == 0;
  a[0] = 20;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 4294967296);
  testResult &= q[0] == 2;
  testResult &= r[0] == 0;
  a[0] = 321573290;
  b[0] = 10;
  rsBigNumber::divide(a, b, 1, q, r, 4294967296);
  testResult &= q[0] == 32157329;
  testResult &= r[0] == 0;

  // test a couple of more cases for division:
  testResult &= testDigitDivision2();
  testResult &= testDigitDivision5();

  // test division for fractions:

  // 0.87 / 0.31 = 02.806451612903226
  a[0] = 8; a[1] = 7;
  b[0] = 3; b[1] = 1;
  rsBigNumber::divideFractions(a, 2, b, 2, q, 4);
  testResult &= q[0]==0 && q[1]==2 && q[2]==8 && q[3]==1;
  rsBigNumber::divideFractions(a, 2, b, 2, q, 15);
  testResult &= q[ 0]==0 && q[ 1]==2 && q[ 2]==8 && q[ 3]==0 && q[ 4]==6 && q[ 5]==4 && q[ 6]==5
             && q[ 7]==1 && q[ 8]==6 && q[ 9]==1 && q[10]==2 && q[11]==9 && q[12]==0 && q[13]==3
             && q[14]==2;

  // 0.871 / 0.31 = 02.809677419354839
  a[2] = 1;
  rsBigNumber::divideFractions(a, 3, b, 2, q, 5);
  testResult &= q[0]==0 && q[1]==2 && q[2]==8 && q[3]==1 && q[4]==0;
  rsBigNumber::divideFractions(a, 3, b, 2, q, 15);
  testResult &= q[ 0]==0 && q[ 1]==2 && q[ 2]==8 && q[ 3]==0 && q[ 4]==9 && q[ 5]==6 && q[ 6]==7
             && q[ 7]==7 && q[ 8]==4 && q[ 9]==1 && q[10]==9 && q[11]==3 && q[12]==5 && q[13]==4
             && q[14]==8;

  // 0.87 / 0.311 = 002.797427652733119
  b[2] = 1;
  rsBigNumber::divideFractions(a, 2, b, 3, q, 5);
  testResult &= q[0]==0 && q[1]==0 && q[2]==2 && q[3]==8 && q[4]==0;
  rsBigNumber::divideFractions(a, 2, b, 3, q, 15);
  testResult &= q[ 0]==0 && q[ 1]==0 && q[ 2]==2 && q[ 3]==7 && q[ 4]==9 && q[ 5]==7 && q[ 6]==4
             && q[ 7]==2 && q[ 8]==7 && q[ 9]==6 && q[10]==5 && q[11]==2 && q[12]==7 && q[13]==3
             && q[14]==3;

  testResult &= exhaustiveIntegerDivisionTest(1);
  testResult &= exhaustiveIntegerDivisionTest(2);
  //testResult &= exhaustiveIntegerDivisionTest(3); // takes several seconds
  //testResult &= exhaustiveIntegerDivisionTest(4); // takes several minutes
    // uncomment the latter 2 tests only, when the division algo was modified

  testResult &= exhaustiveFractionDivisionTest(1);
  testResult &= exhaustiveFractionDivisionTest(2);
   //testResult &= exhaustiveFractionDivisionTest(3);


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testBaseChange(std::string &reportString)
{
  std::string testName = "BaseChange";
  bool t = true;  // testResult

  rsUint32 x[10];
  rsArray<rsUint32> y;
  int i; // loop counter

  // 27_10 = 11011_2 (27 in base 10 is 110011 in base 2):
  x[0] = 2; x[1] = 7;
  y = rsBigNumber::changeBaseForInteger(x, 2, 10, 2);
  t &= y.getNumElements() == 5;
  t &= y[0] == 1; t &= y[1] == 1; t &= y[2] == 0; t &= y[3] == 1; t &= y[4] == 1;

  // 28_10 = 11100_2:
  x[0] = 2; x[1] = 8;
  y = rsBigNumber::changeBaseForInteger(x, 2, 10, 2);
  t &= y.getNumElements() == 5;
  t &= y[0] == 1; t &= y[1] == 1; t &= y[2] == 1; t &= y[3] == 0; t &= y[4] == 0;

  // 97_10 = 141_8 = 61_16 = 1100001_2:
  x[0] = 9; x[1] = 7;
  y = rsBigNumber::changeBaseForInteger(x, 2, 10, 8);
  t &= y.getNumElements() == 3;
  t &= y[0] == 1; t &= y[1] == 4; t &= y[2] == 1;
  y = rsBigNumber::changeBaseForInteger(x, 2, 10, 16);
  t &= y.getNumElements() == 2;
  t &= y[0] == 6; t &= y[1] == 1;
  y = rsBigNumber::changeBaseForInteger(x, 2, 10, 2);
  t &= y.getNumElements() == 7;
  t&=y[0]==1;t&=y[1]==1;t&=y[2]==0;t&=y[3]==0;t&=y[4]==0;t&=y[5]==0;t&=y[6]==1;
  // use t &= y[0]==1 && y[1]==1 && ...

  // 9357806241_10 = 2(767871649)_4294967296
  x[0]=9;x[1]=3;x[2]=5;x[3]=7;x[4]=8;x[5]=0;x[6]=6;x[7]=2;x[8]=4;x[9]=1;
  y = rsBigNumber::changeBaseForInteger(x, 10, 10, 4294967296);
  t &= y.getNumElements() == 2;
  t &= y[0] == 2; t &= y[1] == 767871649;

  // let's denote numbers in bases >= 36 like: 2'56'2'456'23'45'6'234 as strings

  static const int nz = 20;  // number of digits for z
  rsUint32 z[nz];            // == x in another base

  // test comparison to 1/2 in bases 10, 8 and 5:
  x[0] = 6; x[1] = 0; x[2] = 0; x[3] = 0;  // 0.6000_10
  t &= rsBigNumber::compareToOneHalf(x, 4, 10) == +1;
  x[0] = 4; x[1] = 9; x[2] = 9; x[3] = 9; // 0.4999_10
  t &= rsBigNumber::compareToOneHalf(x, 4, 10) == -1;
  x[0] = 5; x[1] = 0; x[2] = 0; x[3] = 1; // 0.5001_10
  t &= rsBigNumber::compareToOneHalf(x, 4, 10) == +1;
  x[0] = 5; x[1] = 0; x[2] = 0; x[3] = 0; // 0.5000_10
  t &= rsBigNumber::compareToOneHalf(x, 4, 10) ==  0;

  x[0] = 5; x[1] = 0; x[2] = 0; x[3] = 0;  // 0.5000_8
  t &= rsBigNumber::compareToOneHalf(x, 4, 8) == +1;
  x[0] = 3; x[1] = 7; x[2] = 7; x[3] = 7; // 0.3777_8 = 3/8+7/8^2+7/8^3+7/8^4 = 0.499755859375...
  t &= rsBigNumber::compareToOneHalf(x, 4, 8) == -1;
  x[0] = 4; x[1] = 0; x[2] = 0; x[3] = 1; // 0.4001_8
  t &= rsBigNumber::compareToOneHalf(x, 4, 8) == +1;
  x[0] = 4; x[1] = 0; x[2] = 0; x[3] = 0; // 0.4000_8 == 1/2
  t &= rsBigNumber::compareToOneHalf(x, 4, 8) ==  0;

  x[0] = 3; x[1] = 0; x[2] = 0; x[3] = 0;  // 0.3000_5 = 3/5 = 0.6
  t &= rsBigNumber::compareToOneHalf(x, 4, 5) == +1;
  x[0] = 2; x[1] = 4; x[2] = 4; x[3] = 4;  // 0.2444_5 = 2/5 + 4/5^2 + 4/5^3 + 4/5^4 = 0.5984
  t &= rsBigNumber::compareToOneHalf(x, 4, 5) == +1;
  x[0] = 2; x[1] = 2; x[2] = 2; x[3] = 3;  // 0.2223_5 = 2/5 + 2/5^2 + 2/5^3 + 3/5^4 = 0.5008
  t &= rsBigNumber::compareToOneHalf(x, 4, 5) == +1;
  x[0] = 2; x[1] = 2; x[2] = 2; x[3] = 2;  // 0.2223_5 = 2/5 + 2/5^2 + 2/5^3 + 2/5^4 = 0.4992
  t &= rsBigNumber::compareToOneHalf(x, 4, 5) == -1;
  x[0] = 2; x[1] = 2; x[2] = 1; x[3] = 4;  // 0.2214_5 = 2/5 + 2/5^2 + 1/5^3 + 4/5^4 = 0.4944
  t &= rsBigNumber::compareToOneHalf(x, 4, 5) == -1;


  // test rounding:
  x[0] = 6; x[1] = 8; x[2] = 7; x[3] = 5; x[4] = 3;  // 0.68753_10 -> 0.68800_10
  t &= rsBigNumber::round(x, 5, 3, 10) == false;
  t &= x[0] == 6; t &= x[1] == 8; t &= x[2] == 8; t &= x[3] == 0; ; t &= x[4] == 0;
  x[0] = 6; x[1] = 8; x[2] = 7; x[3] = 4; x[4] = 7;  // 0.68747_10 -> 0.68700_10
  t &= rsBigNumber::round(x, 5, 3, 10) == false;
  t &= x[0] == 6; t &= x[1] == 8; t &= x[2] == 7; t &= x[3] == 0; ; t &= x[4] == 0;
  x[0] = 6; x[1] = 8; x[2] = 7; x[3] = 5; x[4] = 0;  // 0.68750_10 -> 0.68800_10 (up to even)
  t &= rsBigNumber::round(x, 5, 3, 10) == false;
  t &= x[0] == 6; t &= x[1] == 8; t &= x[2] == 8; t &= x[3] == 0; ; t &= x[4] == 0;
  x[0] = 6; x[1] = 8; x[2] = 6; x[3] = 5; x[4] = 0;  // 0.68650_10 -> 0.68700_10 (down to even)
  t &= rsBigNumber::round(x, 5, 3, 10) == false;
  t &= x[0] == 6; t &= x[1] == 8; t &= x[2] == 6; t &= x[3] == 0; ; t &= x[4] == 0;
  x[0] = 9; x[1] = 9; x[2] = 9; x[3] = 5; x[4] = 0;   // 0.99950_10 -> 1.00000_10 (up to even)
  t &= rsBigNumber::round(x, 5, 3, 10) == true; // overflow
  t &= x[0] == 0; t &= x[1] == 0; t &= x[2] == 0; t &= x[3] == 0; ; t &= x[4] == 0;

  // test base-change for fractions:
  // 0.6875_10 = 0.54_8 = 0.B_16 = 0.1011_2:
  x[0] = 6; x[1] = 8; x[2] = 7; x[3] = 5;
  rsBigNumber::changeBaseForFraction(x, 4, 10, z, nz, 8);
  t &= z[0] == 5; t &= z[1] == 4; for(i = 2; i < nz; i++) t&= z[i] == 0;
  rsBigNumber::changeBaseForFraction(x, 4, 10, z, nz, 16);
  t &= z[0] == 11; for(i = 1; i < nz; i++) t&= z[i] == 0;
  rsBigNumber::changeBaseForFraction(x, 4, 10, z, nz, 2);
  t&=z[0]==1; t&=z[1]==0; t&=z[2]==1; t&=z[3]==1; for(i=4; i<nz; i++) t&=z[i]==0;

  // use a z-array that's exactly long enough:
  rsBigNumber::changeBaseForFraction(x, 4, 10, z, 4, 2);
  t&=z[0]==1; t&=z[1]==0; t&=z[2]==1; t&=z[3]==1;

  // use a too short z-array (z should be rounded up):
  rsBigNumber::changeBaseForFraction(x, 4, 10, z, 3, 2);
  t&=z[0]==1; t&=z[1]==1; t&=z[2]==0;

  // 0.1_10 = 0.000110011001100..._2, non-terminating (periodic), but expect the last 2 digits to
  // break the pattern because of the rounding rule:
  x[0] = 1;
  rsBigNumber::changeBaseForFraction(x, 1, 10, z, nz, 2);
  t&=z[0]==0;t&=z[18]==1;t&=z[19]==0;
  for(i=1;i<=13;i+=4){t&=z[i]==0;t&=z[i+1]==0;t&=z[i+2]==1;t&=z[i+3]==1;}

  // 0.200_6 = 0.3333333..._10
  x[0] = 2; x[1] = 0; x[2] = 0;
  rsBigNumber::changeBaseForFraction(x, 3, 6, z, nz, 10);
  for(i=0; i<nz; i++) t&=z[i]==3;

  // 0.400_6 = 0.6666666..._10
  x[0] = 4; x[1] = 0; x[2] = 0;
  rsBigNumber::changeBaseForFraction(x, 3, 6, z, nz, 10);
  t&=z[nz-1]==7; for(i=0; i<nz-1; i++) t&=z[i]==6;

  appendTestResultToReport(reportString, testName, t);
  return t;
}

bool testBigInt(std::string &reportString)
{
  std::string testName = "rsBigInt";
  bool testResult = true;

  rsUint32 buf[4]; //

  //rsBigInt i1(0, 1), i2(0, 2), i3(0, 3), i4(0, 4); // 4 different sizes

  rsUint32 testDigits[8] = {0,0,0,2,1,4,7,3};
  rsBigInt i1, i2, i3, i4;

  rsInt64 base    = 4294967296;
  rsInt64 largest = base-1;

  // test constructor that creates a number form a digit array:
  i1 = rsBigInt(&testDigits[3], 4, false);
  testResult &= i1.getNumDigits() == 4;
  i1.copyDigitsToBuffer(buf);
  testResult &= buf[0] == 2;
  testResult &= buf[1] == 1;
  testResult &= buf[2] == 4;
  testResult &= buf[3] == 7;
  i1 = rsBigInt(&testDigits[0], 7, false);
  testResult &= i1.getNumDigits() == 4;
  i1.copyDigitsToBuffer(buf);
  testResult &= buf[0] == 2;
  testResult &= buf[1] == 1;
  testResult &= buf[2] == 4;
  testResult &= buf[3] == 7;
  i1 = rsBigInt(&testDigits[0], 3, false);
  testResult &= i1.getNumDigits() == 1;
  i1.copyDigitsToBuffer(buf);
  testResult &= buf[0] == 0;

  // test largest 1-digit numbers:
  i1.setValue(largest);
  i1.copyDigitsToBuffer(buf);
  testResult &= i1.getNumDigits() == 1;
  testResult &= buf[0] == (rsUint32) largest;
  testResult &= i1.isNegative() == false;
  i1.setValue(-largest);
  i1.copyDigitsToBuffer(buf);
  testResult &= i1.getNumDigits() == 1;
  testResult &= buf[0] == (rsUint32) largest;
  testResult &= i1.isNegative() == true;

  // test smallest 2 digit numbers:
  i1.setValue(largest+1);
  i1.copyDigitsToBuffer(buf);
  testResult &= i1.getNumDigits() == 2;
  testResult &= buf[0] == 1;
  testResult &= buf[1] == 0;
  testResult &= i1.isNegative()  == false;
  i1.setValue(-(largest+1));
  i1.copyDigitsToBuffer(buf);
  testResult &= i1.getNumDigits() == 2;
  testResult &= buf[0] == 1;
  testResult &= buf[1] == 0;
  testResult &= i1.isNegative()  == true;

  // test some other 2 digit number:
  i1.setValue(largest+5);
  i1.copyDigitsToBuffer(buf);
  testResult &= i1.getNumDigits() == 2;
  testResult &= buf[0] == 1;
  testResult &= buf[1] == 4;
  testResult &= i1.isNegative()  == false;

  // test main comparison operators ">" and "==" (all others are based on these):
  testResult &= (rsBigInt( 5) > rsBigInt( 5)) == false;
  testResult &= (rsBigInt( 5) > rsBigInt( 3)) == true;
  testResult &= (rsBigInt( 3) > rsBigInt( 5)) == false;
  testResult &= (rsBigInt(-5) > rsBigInt(-5)) == false;
  testResult &= (rsBigInt(-5) > rsBigInt(-3)) == false;
  testResult &= (rsBigInt(-3) > rsBigInt(-5)) == true;
  testResult &= (rsBigInt(-5) > rsBigInt( 3)) == false;
  testResult &= (rsBigInt( 3) > rsBigInt(-5)) == true;
  testResult &= ( rsBigInt(largest+1) >  rsBigInt(largest  )) == true;
  testResult &= ( rsBigInt(largest  ) >  rsBigInt(largest+1)) == false;
  testResult &= (-rsBigInt(largest+1) > -rsBigInt(largest  )) == false;
  testResult &= (-rsBigInt(largest  ) > -rsBigInt(largest+1)) == true;

  testResult &= (rsBigInt( 5) == rsBigInt( 5)) == true;
  testResult &= (rsBigInt( 5) == rsBigInt(-5)) == false;
  testResult &= (rsBigInt(-5) == rsBigInt( 5)) == false;
  testResult &= (rsBigInt( 5) == rsBigInt( 3)) == false;

  // test addition and "==" operator, we assume that creating a BigInt from an int64 works
  // correctly - the test above should fail, if this is not the case:
  i1 = rsBigInt(largest) + rsBigInt(0);
  testResult &= i1 == rsBigInt(largest);
  i1 = rsBigInt(largest) + rsBigInt(1);
  testResult &= i1 == rsBigInt(largest+1);
  i1 = rsBigInt(largest) + rsBigInt(5);
  testResult &= i1 == rsBigInt(largest+5);
  i1 = rsBigInt(largest) + rsBigInt(-1);
  testResult &= i1 == rsBigInt(largest-1);
  i1 = rsBigInt(largest+1) + rsBigInt(-1);
  testResult &= i1 == rsBigInt(largest);
  i1 = rsBigInt(largest) + rsBigInt(-5);
  testResult &= i1 == rsBigInt(largest-5);
  i1 =  rsBigInt(-5) + rsBigInt(largest);
  testResult &= i1 == rsBigInt(largest-5);

  i1 = rsBigInt( 5) + rsBigInt( 3); testResult &= i1 == rsBigInt( 8); //  5 +  3 =  8
  i1 = rsBigInt( 5) + rsBigInt(-3); testResult &= i1 == rsBigInt( 2); //  5 + -3 =  2
  i1 = rsBigInt(-5) + rsBigInt( 3); testResult &= i1 == rsBigInt(-2); // -5 +  3 = -2
  i1 = rsBigInt(-5) + rsBigInt(-3); testResult &= i1 == rsBigInt(-8); // -5 + -3 = -8
  i1 = rsBigInt( 3) + rsBigInt( 5); testResult &= i1 == rsBigInt( 8); //  3 +  5 =  8
  i1 = rsBigInt( 3) + rsBigInt(-5); testResult &= i1 == rsBigInt(-2); //  3 + -5 = -2
  i1 = rsBigInt(-3) + rsBigInt( 5); testResult &= i1 == rsBigInt( 2); // -3 +  5 =  2
  i1 = rsBigInt(-3) + rsBigInt(-5); testResult &= i1 == rsBigInt(-8); // -3 + -5 = -8
  i1 = rsBigInt(-5) + rsBigInt( 5); testResult &= i1 == rsBigInt( 0); // -5 +  5 =  0
  i1 = rsBigInt( 5) + rsBigInt(-5); testResult &= i1 == rsBigInt( 0); //  5 + -5 =  0

  // test subtraction:
  i1 = rsBigInt( 5) - rsBigInt( 3); testResult &= i1 == rsBigInt( 2); //  5 -  3 =  2
  i1 = rsBigInt( 5) - rsBigInt(-3); testResult &= i1 == rsBigInt( 8); //  5 - -3 =  8
  i1 = rsBigInt(-5) - rsBigInt( 3); testResult &= i1 == rsBigInt(-8); // -5 -  3 = -8
  i1 = rsBigInt(-5) - rsBigInt(-3); testResult &= i1 == rsBigInt(-2); // -5 - -3 = -2
  i1 = rsBigInt( 3) - rsBigInt( 5); testResult &= i1 == rsBigInt(-2); //  3 -  5 = -2
  i1 = rsBigInt( 3) - rsBigInt(-5); testResult &= i1 == rsBigInt( 8); //  3 - -5 =  8
  i1 = rsBigInt(-3) - rsBigInt( 5); testResult &= i1 == rsBigInt(-8); // -3 -  5 = -8
  i1 = rsBigInt(-3) - rsBigInt(-5); testResult &= i1 == rsBigInt( 2); // -3 - -5 =  2
  i1 = rsBigInt(-3) - rsBigInt( 3); testResult &= i1 == rsBigInt(-6); // -3 -  3 = -6
  i1 = rsBigInt( 3) - rsBigInt(-3); testResult &= i1 == rsBigInt( 6); //  3 - -3 =  6

  // test multiplication:
  i1 = rsBigInt(largest) * rsBigInt(largest);
  testResult &= i1.getNumDigits() == 2;
  i1.copyDigitsToBuffer(buf);
  testResult &= buf[0] == largest-1;
  testResult &= buf[1] == 1;
  i1 = rsBigInt(largest+5) * rsBigInt(10);
  testResult &= i1 == rsBigInt(((rsInt64)largest+5)*10);

  // test conversion from string:
  i1 = rsBigInt::fromString("1234");
  testResult &= i1 == rsBigInt(1234);
  i1 = rsBigInt::fromString("-1234");
  testResult &= i1 == rsBigInt(-1234);
  i1 = rsBigInt::fromString("001234");
  testResult &= i1 == rsBigInt(1234);
  i1 = rsBigInt::fromString("-001234");
  testResult &= i1 == rsBigInt(-1234);
  i1 = rsBigInt::fromString("4294967295");
  testResult &= i1 == rsBigInt(4294967295);
  i1 = rsBigInt::fromString("4294967296");
  testResult &= i1 == rsBigInt(4294967296);

  // test division and modulo, examples were created with maxima using, for example:
  // divide(7437112567129538, 167481):
  i1 = rsBigInt::fromString("7437112567129538");
  i2 = rsBigInt::fromString("167481");
  i1.divMod(i2, i3, i4);
  testResult &= i3 == rsBigInt::fromString("44405709108");
  testResult &= i4 == rsBigInt::fromString("12590");
  i1 = rsBigInt::fromString("3215732907234665911873258216741872582173456128732762176216234");
  i2 = rsBigInt::fromString("23183651827465253431");
  i1.divMod(i2, i3, i4);
  testResult &= i3 == rsBigInt::fromString("138706918615170248522905946156808970120913");
  testResult &= i4 == rsBigInt::fromString("6408014379918113731");
  i3 = i1 / i2;
  i4 = i1 % i2;
  testResult &= i3 == rsBigInt::fromString("138706918615170248522905946156808970120913");
  testResult &= i4 == rsBigInt::fromString("6408014379918113731");
  i1 = rsBigInt(35) / rsBigInt(10);
  testResult &= i1 == rsBigInt::fromString("3");   //  35 /  10 =  3
  i1 = rsBigInt(35) % rsBigInt(10);
  testResult &= i1 == rsBigInt::fromString("5");   //  35 %  10 =  5
  i1 = rsBigInt(35) / rsBigInt(-10);
  testResult &= i1 == rsBigInt::fromString("-3");  //  35 / -10 = -3
  i1 = rsBigInt(35) % rsBigInt(-10);
  testResult &= i1 == rsBigInt::fromString("5");   //  35 % -10 =  5
  i1 = rsBigInt(-35) / rsBigInt(10);
  testResult &= i1 == rsBigInt::fromString("-3");  // -35 /  10 = -3
  i1 = rsBigInt(-35) % rsBigInt(10);
  testResult &= i1 == rsBigInt::fromString("-5");  // -35 %  10 = -5
  i1 = rsBigInt(-35) / rsBigInt(-10);
  testResult &= i1 == rsBigInt::fromString("3");   // -35 / -10 =  3
  i1 = rsBigInt(-35) % rsBigInt(-10);
  testResult &= i1 == rsBigInt::fromString("-5");  // -35 % -10 = -5

  // test conversion to a string:
  rsString s1("3215732907234665911873258216741872582173456128732762176216234");
  rsString s2;
  i1 = rsBigInt::fromString(s1);
  s2 = i1.toString();
  testResult &= s2 == s1;
  s1 = rsString("-3215732907234665911873258216741872582173456128732762176216234");
  i1 = rsBigInt::fromString(s1);
  s2 = i1.toString();
  testResult &= s2 == s1;

  // create a rsBigInt from a hexadecimal representation, where the digits above 9 are defined as
  // A=10, B=11, C=12, D=13, E=14, F=15 (lowercase should also work)
  // B7F3A02C_16 = 11*16^7 + 7*16^6 + 15*16^5 + 3*16^4 + 10*16^3 + 0*16^2 + 2*16^1 + 12*16^0
  //             = 3086196780
  testResult &= rsBigInt::fromString("B7F3A02C", 16) == 3086196780;
  testResult &= rsBigInt::fromString("b7f3a02c", 16) == 3086196780;

  s1 = rsBigInt(3086196780).toString(16);
  testResult &= rsBigInt(3086196780).toString(16) == rsString("B7F3A02C");

  // ...from other bases:
  testResult &= rsBigInt::fromString("1100001", 2)   == 97;
  testResult &= rsConvertIntegerIntoBase("97", 10, 2) == rsString("1100001");

  // test string conversion to other bases:
  i1 = rsBigInt::fromString("97");
  //....

  // test math functions:
  testResult &= rsPow(rsBigInt(7), 13) == rsBigInt(96889010407); // 7^13 = 96 889 010 407
  i1 = rsPow(rsBigInt(4294967296, 10), 3);
  testResult &= i1 == rsBigInt::fromString("79228162514264337593543950336", 10, 10);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testBigFloatConversion(std::string &reportString)
{
  std::string testName = "rsBigFloatStringConversion";
  bool testResult = true;

  //rsInt64  base = 4294967296;

  double xd, yd;     // values before and after conversion cycle
  //double ulp;      // unit in the last place
  double ulp2;       // half of unit in the last place
  double err;        // absolute error in conversion cycle
  double tol;        // error tolerance
  //double exc;      // excess error
  rsBigFloat x, y;
  rsUint32 f[100];   // buffer for the fraction
  rsString s;
  rsUint64 B;        // base
  rsUint32 N;        // number of digits
  rsUint32 P;        // precision (in bits);


  // base conversion:

  f[0] = 0; f[1] = 1; f[2] = 2; // x = 0.012_3 = 0.185185185..._10 -> 0.1851852
  x.setData(f, 3, 3, 0, false, true);

  y = x.toBase(10, 1); y.copyDigitsToBuffer(f);
  testResult &= f[0]==2;
  y = x.toBase(10, 2); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==9;
  y = x.toBase(10, 3); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==8 && f[2]==5;
  y = x.toBase(10, 4); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==8 && f[2]==5 && f[3]==2;
  y = x.toBase(10, 5); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==8 && f[2]==5 && f[3]==1 && f[4]==9;
  y = x.toBase(10, 6); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==8 && f[2]==5 && f[3]==1 && f[4]==8 && f[5]==5;
  y = x.toBase(10, 7); y.copyDigitsToBuffer(f);
  testResult &= f[0]==1 && f[1]==8 && f[2]==5 && f[3]==1 && f[4]==8 && f[5]==5 && f[6]==2;

  // x = 0.121e-4_3 = 3^(-4) * (1/3^1 + 2/3^2 + 1/3^3) = 16_10/2187_10
  //   = 0.0073159579332419..._10 = 0.73159579332419e-2..._10:
  f[0] = 1; f[1] = 2; f[2] = 1; x.setData(f, 3, 3, -4, false, true);
  y = x.toBase(10, 1); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7;
  y = x.toBase(10, 2); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3;
  y = x.toBase(10, 3); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3 && f[2]==2;
  y = x.toBase(10, 4); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3 && f[2]==1 && f[3]==6;
  y = x.toBase(10, 5); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3 && f[2]==1 && f[3]==6 && f[4]==0;
  y = x.toBase(10, 6); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3 && f[2]==1 && f[3]==5 && f[4]==9 && f[5]==6;
  y = x.toBase(10, 7); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == -2;
  testResult &= f[0]==7 && f[1]==3 && f[2]==1 && f[3]==5 && f[4]==9 && f[5]==5 && f[6]==8;

  // x = 0.121e0_3 = 3^0 * (1/3^1 + 2/3^2 + 1/3^3) = 16_10/27_10 = 0.592592592..._10:
  f[0] = 1; f[1] = 2; f[2] = 1; x.setData(f, 3, 3, 0, false, true);
  y = x.toBase(10, 1); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==6;
  y = x.toBase(10, 2); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9;
  y = x.toBase(10, 3); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9 && f[2]==3;
  y = x.toBase(10, 4); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9 && f[2]==2 && f[3]==6;
  y = x.toBase(10, 5); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9 && f[2]==2 && f[3]==5 && f[4]==9;
  y = x.toBase(10, 6); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9 && f[2]==2 && f[3]==5 && f[4]==9 && f[5]==3;
  y = x.toBase(10, 7); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 0;
  testResult &= f[0]==5 && f[1]==9 && f[2]==2 && f[3]==5 && f[4]==9 && f[5]==2 && f[6]==6;

  // x = 0.121e3_3 = 3^3 * (1/3^1 + 2/3^2 + 1/3^3) = 16_10 = 0.16e2_10:
  f[0] = 1; f[1] = 2; f[2] = 1; x.setData(f, 3, 3, 3, false, true);
  y = x.toBase(10, 1); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==2;
  y = x.toBase(10, 2); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6;
  y = x.toBase(10, 3); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6 && f[2]==0;
  y = x.toBase(10, 4); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6 && f[2]==0 && f[3]==0;
  y = x.toBase(10, 5); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6 && f[2]==0 && f[3]==0 && f[4]==0;
  y = x.toBase(10, 6); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6 && f[2]==0 && f[3]==0 && f[4]==0 && f[5]==0;
  y = x.toBase(10, 7); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 2;
  testResult &= f[0]==1 && f[1]==6 && f[2]==0 && f[3]==0 && f[4]==0 && f[5]==0 && f[6]==0;

  // x = 0.121e9_3 = 3^9 * (1/3^1 + 2/3^2 + 1/3^3) = 11664_10 = 0.11664e5_10:
  f[0] = 1; f[1] = 2; f[2] = 1; x.setData(f, 3, 3, 9, false, true);
  y = x.toBase(10, 1); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1;
  y = x.toBase(10, 2); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==2;
  y = x.toBase(10, 3); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==1 && f[2]==7;
  y = x.toBase(10, 4); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==1 && f[2]==6 && f[3]==6;
  y = x.toBase(10, 5); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==1 && f[2]==6 && f[3]==6 && f[4]==4;
  y = x.toBase(10, 6); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==1 && f[2]==6 && f[3]==6 && f[4]==4 && f[5]==0;
  y = x.toBase(10, 7); y.copyDigitsToBuffer(f); testResult &= y.getExponent() == 5;
  testResult &= f[0]==1 && f[1]==1 && f[2]==6 && f[3]==6 && f[4]==4 && f[5]==0 && f[6]==0;

  f[0] = 0; f[1] = 0; f[2] = 2; f[3] = 0; f[4] = 2; f[5] = 1;
  x.setData(f, 6, 3, 0, false, true);
  x = x.toBase(10, 7);  // 0.002021_3 = 0.083676268861454..._10 -> 0.08367627
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 8 &&f[1]== 3 &&f[2]== 6 &&f[3]== 7 &&f[4]== 6 &&f[5]== 2 &&f[6]== 7;
  testResult &= x.getExponent() == -1;

  f[0] = 0; f[1] = 0; f[2] = 1;  // 0.001_2 = 0.125_10
  x.setData(f, 3, 2, 0, false, true);
  x = x.toBase(10, 3);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 5;
  x = x.toBase(65536, 3);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 8192 &&f[1]== 0 &&f[2]== 0;
  x = x.toBase(1048576, 3);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 131072 &&f[1]== 0 &&f[2]== 0;
  x = x.toBase(1073741824, 3); // 2^30
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 134217728 &&f[1]== 0 &&f[2]== 0;  // 2^27
  x = x.toBase(4294967296, 3); // 2^32
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 536870912 &&f[1]== 0 &&f[2]== 0;  // 2^29

  x = x.toBase(10, 3);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 5;

  f[0] = 2; f[1] = 3;  // 0.23_10
  x.setData(f, 2, 10, 0, false, true);
  x = x.toBase(4294967296, 8);
  x = x.toBase(3489660929, 8);
  x = x.toBase(10, 8);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 2 &&f[1]== 3 &&f[2]== 0 &&f[3]== 0 &&f[4]== 0 &&f[5]== 0 &&f[6]== 0;
  testResult &= x.getExponent() == 0;

  // conversion from/to built-in double precision numbers:

  x = 0.25;
  x.fromDouble(0.25, 4, 2);
  x = x.toBase(10, 8);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 2 &&f[1]== 5 &&f[2]== 0 &&f[3]== 0 &&f[4]== 0 &&f[5]== 0 &&f[6]== 0;

  x = 0.23;
  x = x.toBase(10, 8);
  x.copyDigitsToBuffer(f);
  testResult &= f[0]== 2 &&f[1]== 3 &&f[2]== 0 &&f[3]== 0 &&f[4]== 0 &&f[5]== 0 &&f[6]== 0;
  s = x.toString(5, 10);
  testResult &= s == rsString("0.23000e0");

  x  = 0.125;  // 0.10100000e-2_2
  y = x.toBase(2, 3);
  y.copyDigitsToBuffer(f);
  testResult &= f[0]== 1 &&f[1]== 0 && f[2]== 0;
  y = x.toBase(10, 5);
  testResult &= y == "0.125";
  xd = x.toDouble();
  testResult &= xd == 0.125;
  x  = 1.0;
  xd = x.toDouble();
  testResult &= xd == 1.0;
  x  = 2.0;
  xd = x.toDouble();
  testResult &= xd == 2.0;
  x  = 0.123;
  xd = x.toDouble();
  testResult &= xd == 0.123;

  x.fromDouble(2.0, 100, 10);
  xd = x.toDouble();
  testResult &= xd == 2.0;

  x.fromDouble(4.0, 100, 10);
  xd = x.toDouble();
  testResult &= xd == 4.0;

  x.fromDouble(4.0, 3, 10);
  xd = x.toDouble();
  testResult &= xd == 4.0;

  x.fromDouble(7.99, 3, 10);
  xd = x.toDouble();
  testResult &= xd == 7.99;

  x.fromDouble(8.0, 3, 10);
  xd = x.toDouble();
  testResult &= xd == 8.0;

  x = "9";
  //y = x.toBase(2, 5);
  //y.copyDigitsToBuffer(f);
  xd = x.toDouble();
  testResult &= xd == 9.0;

  xd  = rsBigFloat::precisionInBits(53, 2);  // test - should be 53


  xd  = 2.5;
  P   = 3;    // number of bits required to represent 2.5_10 = 0.101e2

  tol = RS_EPS(double);
  tol = xd * RS_EPS(double);

  B = 2; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  x.copyDigitsToBuffer(f);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2;

  B = 3; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  //x.copyDigitsToBuffer(f); // = 0.22e1_3 = 8/3 = 2.6666..._10
  //y = x.toBase(2, 53);
  //y.copyDigitsToBuffer(f);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 5; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  x.copyDigitsToBuffer(f); // = 0.23e1_5 = 13/5 = 2.6_10
  //y = x.toBase(2, 53);
  //y.copyDigitsToBuffer(f);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 7; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 16; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  //x.copyDigitsToBuffer(f);
  //y = x.toBase(2, 53);
  //y.copyDigitsToBuffer(f);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 32; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  //x.copyDigitsToBuffer(f);
  //y = x.toBase(2, 53);
  //y.copyDigitsToBuffer(f);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 100; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 1000000000; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 2147483648; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 4294967295; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 4294967296; N = rsBigFloat::numRequiredDigits(P, B); x.fromDouble(xd, N, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 4294967296; x.fromDouble(xd, 2, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  B = 4000000000; x.fromDouble(xd, 2, B);
  ulp2 = 0.5*x.getUnitInTheLastPlace(); yd = x.toDouble(); err = fabs(yd-xd);
  testResult &= err <= ulp2+tol;

  x.fromDouble(2.5, 10, 1000000000);
  xd = x.toDouble();
  testResult &= xd == 2.5;

  x.fromDouble(2.5, 10, 4294967296);
  xd = x.toDouble();
  testResult &= xd == 2.5;

  x.fromDouble(2.5, 10, 4000000000);
  xd = x.toDouble();
  testResult &= xd == 2.5;

  // initialization from a decimal string:

  x = "0";       testResult &= x.isZero();
  x = "00";      testResult &= x.isZero();
  x = "0.0";     testResult &= x.isZero();
  x = "0.00";    testResult &= x.isZero();
  x = "00.000";  testResult &= x.isZero();

  x = "-0043.7856e-03";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 6;
  testResult &= x.getExponent()  == -1;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 7 &&f[3]== 8 &&f[4]== 5 &&f[5]== 6;

  x = "43.7856e-3";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 6;
  testResult &= x.getExponent()  == -1;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 7 &&f[3]== 8 &&f[4]== 5 &&f[5]== 6;

  x = "43.7856e3";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 6;
  testResult &= x.getExponent()  == 5;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 7 &&f[3]== 8 &&f[4]== 5 &&f[5]== 6;

  x = "43.7856e+3";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 6;
  testResult &= x.getExponent()  == 5;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 7 &&f[3]== 8 &&f[4]== 5 &&f[5]== 6;

  x = "-00430.7856e-03";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 7;
  testResult &= x.getExponent()  == 0;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 0 &&f[3]== 7 &&f[4]== 8 &&f[5]== 5 &&f[6]== 6;

  x = "-00.004307856e-03";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 7;
  testResult &= x.getExponent()  == -5;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 0 &&f[3]== 7 &&f[4]== 8 &&f[5]== 5 &&f[6]== 6;

  x = "-.004307856e-03";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 7;
  testResult &= x.getExponent()  == -5;
  testResult &= f[0]== 4 &&f[1]== 3 &&f[2]== 0 &&f[3]== 7 &&f[4]== 8 &&f[5]== 5 &&f[6]== 6;

  x = "1272";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 4;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "-1272";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 4;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "001272";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 4;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "-001272";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 4;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "-001272e-6";  //e=-2
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == true;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == -2;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "67.125";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 5;
  testResult &= x.getExponent()  == 2;
  testResult &= f[0]== 6 &&f[1]== 7 &&f[2]== 1 &&f[3]== 2 &&f[4]== 5;

  x = "0.563";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 3;
  testResult &= x.getExponent()  == 0;
  testResult &= f[0]== 5 &&f[1]== 6 &&f[2]== 3;

  x = "+000.0025";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 2;
  testResult &= x.getExponent()  == -2;
  testResult &= f[0]== 2 &&f[1]== 5;

  x = "127200";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 6;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "1272.00";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 4;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;

  x = "12720.00";
  x.copyDigitsToBuffer(f);
  testResult &= x.isNegative()   == false;
  testResult &= x.getNumDigits() == 4;
  testResult &= x.getExponent()  == 5;
  testResult &= f[0]== 1 &&f[1]== 2 &&f[2]== 7 &&f[3]== 2;


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testBigFloatArithmetic(std::string &reportString)
{
  std::string testName = "BigFloatArithmetic";
  bool testResult = true;

  rsBigFloat x, x1, x2, x3, x4;
  rsUint32 f[10];  // buffer for the fraction

  x1 = "1.34567";
  x2 = "0.0243256";
  x3 = x2 + x1;      // 1.3699956
  testResult &= x3 == "1.37000";

  x1 = "0.317";
  x2 = "0.574";

  x3 = x1 + x2;
  testResult &= x3 == "0.891";
  x3 = x2 + x1;
  testResult &= x3 == "0.891";
  x4 = x3 - x2;
  testResult &= x4 == x1;
  x4 = -x2 + x3;
  testResult &= x4 == x1;

  x4 = x3 + x1;
  testResult &= x4 == "1.21"; // 1.208
  x4 = x1 + x3;
  testResult &= x4 == "1.21";

  x2 = "0.000501";
  x3 = x1 + x2;
  testResult &= x3 == "0.318";  // 0.317501
  x3 = x2 + x1;
  testResult &= x3 == "0.318";

  x2 = "0.000500";
  x3 = x1 + x2;
  testResult &= x3 == "0.318";  // 0.317500 -> round to even
  x3 = x2 + x1;
  testResult &= x3 == "0.318";

  x2 = "0.000499";
  x3 = x1 + x2;
  testResult &= x3 == "0.317";  // 0.317499
  x3 = x2 + x1;
  testResult &= x3 == "0.317";

  x1 = "0.316";
  x2 = "0.000500";
  x3 = x1 + x2;
  testResult &= x3 == "0.316";  // 0.316500 -> round to even
  x3 = x2 + x1;
  x3.copyDigitsToBuffer(f);
  testResult &= x3 == "0.316";

  x1 = "0.999";
  x3 = x1 + x2;  // 1.0
  testResult &= x3 == "1.0";

  x3 = x1 + x1 = "2.0";         // 1.998
  testResult &= x3 == "2.0";

  x1 = "1.3";
  x2 = "1.52";
  x3 = x2 + x1;
  testResult &= x3 == "2.82";

  // force denormalized calculation
  x1 = "0.537201";
  x1.setExponent(rsBigFloat::minExponent);
  x2 = "0.536137";
  x2.setExponent(rsBigFloat::minExponent);
  x3 = x1 - x2;
  x3.copyDigitsToBuffer(f);
  x4 = x3 + x2;
  testResult &=  x4 == x1;
  //testResult &= x3 == "0.001064";

  // test multiply-function with 0.522 * 0.32 = 0.16182 and different precisions for the result:
  x1 = "0.522";
  x2 = "0.31";
  x3 = rsBigFloat((rsUint32) 3, (rsUint64) 10);
  rsBigFloat::multiply(x1, x2, x3);
  testResult &= x3 == "0.162";
  rsBigFloat::multiply(x2, x1, x3);
  testResult &= x3 == "0.162";
  x3 = rsBigFloat((rsUint32) 4, (rsUint64) 10);
  rsBigFloat::multiply(x1, x2, x3);
  testResult &= x3 == "0.1618";
  rsBigFloat::multiply(x2, x1, x3);
  testResult &= x3 == "0.1618";
  x3 = rsBigFloat((rsUint32) 5, (rsUint64) 10);
  rsBigFloat::multiply(x1, x2, x3);
  testResult &= x3 == "0.16182";
  rsBigFloat::multiply(x2, x1, x3);
  testResult &= x3 == "0.16182";
  x3 = rsBigFloat((rsUint32) 6, (rsUint64) 10);
  rsBigFloat::multiply(x1, x2, x3);
  testResult &= x3 == "0.16182";
  rsBigFloat::multiply(x2, x1, x3);
  testResult &= x3 == "0.16182";

  // test the multiplication-operator:
  x1 = "0.523";
  x2 = "0.147";
  x3 = x1 * x2;
  testResult &= x3 == "0.0769";  // = 0.076881

  // test the divide-function:
  x1 = "0.523";
  x2 = "0.31";
  x3 = rsBigFloat((rsUint32) 3, (rsUint64) 10);
  rsBigFloat::divide(x1, x2, x3);
  testResult &= x3 == "1.69";          // 1.68709677419... -> 1.69
  rsBigFloat::divide(x2, x1, x3);
  testResult &= x3 == "0.593";         // 0.59273422562... -> 0.593
  x3 = rsBigFloat((rsUint32) 4, (rsUint64) 10);
  rsBigFloat::divide(x1, x2, x3);
  testResult &= x3 == "1.687";
  rsBigFloat::divide(x2, x1, x3);
  testResult &= x3 == "0.5927";
  x3 = rsBigFloat((rsUint32) 5, (rsUint64) 10);
  rsBigFloat::divide(x1, x2, x3);
  testResult &= x3 == "1.6871";
  rsBigFloat::divide(x2, x1, x3);
  testResult &= x3 == "0.59273";
  x3 = rsBigFloat((rsUint32) 6, (rsUint64) 10);
  rsBigFloat::divide(x1, x2, x3);
  testResult &= x3 == "1.68710";
  rsBigFloat::divide(x2, x1, x3);
  testResult &= x3 == "0.592734";
  x3 = rsBigFloat((rsUint32) 7, (rsUint64) 10);
  rsBigFloat::divide(x1, x2, x3);
  testResult &= x3 == "1.687097";
  rsBigFloat::divide(x2, x1, x3);
  testResult &= x3 == "0.5927342";

  // test the division-operator:
  x1 = "0.523";
  x2 = "0.31";
  x4 = x2 / x1;
  testResult &= x4 == "0.593";    // 0.59273422562... -> 0.593
  x3 = x1 / x2;
  testResult &= x3 == "1.69";     // 1.68709677419... -> 1.69
  x4 = x3 / x2;
  testResult &= x4 == "5.45";    // 5.451612903225806... -> 5.45

  // try overflow: 0.999 + 0.0999 = 1.0989 -> 1.01

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testBigFloatFunctions(std::string &reportString)
{
  std::string testName = "rsBigFloatFunctions";
  bool testResult = true;

  rsUint64 B = 10;  // base
  rsUint32 N = 50;  // number of digits

  rsString s;
  //double xd;      // desired value

  rsBigFloat x, y, z;

  // test square-root:
  x.fromDouble(2.0, N, B);
  y = rsSqrt(x);
  s = y.toString(N, (rsUint32)B);
  testResult &= s == rsString("0.14142135623730950488016887242096980785696718753769e1");
    // result: "0.14142135623730950488016887242096980785696718753770e1"
    // correct: 0.1414213562373095048801688724209698078569671875376948073176679737990733e1

  //x.fromDouble(2.0, p, B);
  //y.fromDouble(4.0, p, B);
  //z = rsAGM(x, y);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testBigFloat(std::string &reportString)
{
  std::string testName = "rsBigFloat";
  bool testResult = true;

  //double agm = rsAGM(4.0, 2.0);

  testResult &= testBigFloatArithmetic(reportString);
  testResult &= testBigFloatConversion(reportString);
  //testResult &= testBigFloatFunctions( reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

