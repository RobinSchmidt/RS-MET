
bool testIntAbs(std::string &reportString)
{
  std::string testName = "IntAbs";
  bool r = true; // result

  rsInt64 i64;
  i64 = rsAbs( 4294967294LL); r &= i64 == 4294967294;
  i64 = rsAbs( 4294967295LL); r &= i64 == 4294967295;
  i64 = rsAbs( 4294967296LL); r &= i64 == 4294967296;
  i64 = rsAbs( 2147483648LL); r &= i64 == 2147483648;
  i64 = rsAbs(-4294967294LL); r &= i64 == 4294967294;
  i64 = rsAbs(-4294967295LL); r &= i64 == 4294967295;
  i64 = rsAbs(-4294967296LL); r &= i64 == 4294967296;
  i64 = rsAbs(-2147483648LL); r &= i64 == 2147483648;

  //appendTestResultToReport(reportString, testName, r);
  return r;
}


int rsSumOfProducts1(int n, int m)
{
  int sum = 0;
  for(int k = 0; k < n; k++)
  {
    int prod = 1;
    for(int i = 1; i <= m; i++)
      prod *= k+i;
    sum += prod;
  }
  return sum;
}
// naive version

int rsSumOfProducts2(int n, int m)
{
  int sum = 0;
  for(int k = 0; k < n; k++)
    sum += rsProduct(k+1, k+m);
  return sum;
}
// using rsProduct, needed an index-shift in the product

int rsSumOfProducts(int n, int m)
{
  return rsProduct(n, n+m) / (m+1);
}
// without the outer loop, see: https://www.youtube.com/watch?v=-2PA7SbWoJ0&t=27m
// move to IntegerFunctions.h, but we use different meaning for n - here, n is the number of terms
// whereas in the video, it's the upper summation index for a sum starting at 0, so there are n+1
// terms in total in the sums in the video
// how about a formula for sum_i prod_j (i+j) for i=iMin..iMax and j=jMin..jMax...the iMin..iMax
// can be realized by subtracting two evaluations of our existing function:
//   rsSumOfProducts(iMax, m) - rsSumOfProducts(iMin, m)
// but what about the inner loop?

// is there also a formula for products of sums?


int rsSumOfProducts1(int iMin, int iMax, int jMin, int jMax)
{
  int sum = 0;
  for(int i = iMin; i <= iMax; i++)
  {
    int prod = 1;
    for(int j = jMin; j <= jMax; j++)
      prod *= i+j;
    sum += prod;
  }
  return sum;
}
// naive version - try to find a formula that avoids the outer loop

bool testSumAndProduct()
{
  bool r = true;

  int y;
  y = rsSum(    3, 7); r &= y == 25;   // 3+4+5+6+7 = 25
  y = rsProduct(3, 7); r &= y == 2520; // 3*4*5*6*7 = 2520

  // what if we use negative numbers? maybe it doesn't work - so maybe we should use unsigned
  // integers?


  // todo: implement this formula for a rsSumOfProducts function:
  // https://www.youtube.com/watch?v=-2PA7SbWoJ0&t=27m
  //   sum_{k=0}^n prod_{i=1}^m (k+i)
  // = sum_{k=0}^n prod_{i=k+1}^{k+m} i
  // = (1/(m+1)) prod_{i=1}^{m+1} (n+i)
  // = (1/(m+1)) prod_{i=n+1}^{n+m+1} i

  // should take n,m as parameters - can this be generalized to use different start-values for k
  // and i?

  // n = 4, m = 3: 3 terms and each term has 3 factors:
  // 1*2*3 + 2*3*4 + 3*4*5 + 4*5*6 = 210
  y = rsSumOfProducts1(4, 3); r &= y == 210;
  y = rsSumOfProducts2(4, 3); r &= y == 210;
  y = rsSumOfProducts( 4, 3); r &= y == 210;

  // n = 3, m = 3: 3 terms and each term has 3 factors:
  // 1*2*3 + 2*3*4 + 3*4*5 = 90
  y = rsSumOfProducts1(3, 3); r &= y == 90;
  y = rsSumOfProducts2(3, 3); r &= y == 90;
  y = rsSumOfProducts( 3, 3); r &= y == 90;

  // 5 terms with 4 factors each:
  // 1*2*3*4 + 2*3*4*5 + 3*4*5*6 + 4*5*6*7 + 5*6*7*8 = 3024
  y = rsSumOfProducts1(5, 4); r &= y == 3024;
  y = rsSumOfProducts2(5, 4); r &= y == 3024;
  y = rsSumOfProducts( 5, 4); r &= y == 3024;


  y = rsSumOfProducts1(2, 5, 1, 3);
  // (2+1)*(2+2)*(2+3) + (3+1)*(3+2)*(3+3) + (4+1)*(4+2)*(4+3) + (5+1)*(5+2)*(5+3) =
  // 3*4*5 + 4*5*6 + 5*6*7 + 6*7*8 = 726

  y = rsSumOfProducts1(2, 5, 4, 6);
  // (2+4)*(2+5)*(2+6) + (3+4)*(3+5)*(3+6) + (4+4)*(4+5)*(4+6) + (5+4)*(5+5)*(5+6)
  // 6*7*8 + 7*8*9 + 8*9*10 + 9*10*11 = 2550


  return r;
}

bool testBinomialCoefficients(std::string &reportString)
{
  std::string testName = "rsBinomialCoefficients";
  bool testResult = true;

  using uint = unsigned int;
  static const uint nMax = 20;

  double B1[nMax+1]; rsArrayTools::fillWithNaN(B1, nMax+1);
  double B2[nMax+1]; rsArrayTools::fillWithNaN(B2, nMax+1);
  for(uint n = 0; n <= nMax; n++)
  {
    rsNextPascalTriangleLine(B1, B1, n);
    rsPascalTriangleLine(B2, n);
    for(uint k = 0; k <= n; k++)
    {
      uint c1 = rsBinomialCoefficient(      n, k);
      uint c2 = rsBinomialCoefficientUpTo20(n, k);
      testResult &= (c1 == c2);
      double c3 = B1[k];
      double c4 = B2[k];
      testResult &= (double)c1 == c3;
      testResult &= (double)c1 == c4;
    }
  }

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testMultinomialCoefficients(std::string &reportString)
{
  std::string testName = "rsMultinomialCoefficients";
  bool testResult = true;

  static const rsUint32 mMax = 5;  // maximum number of indices for the k-values
  static const rsUint32 nMax = 12; // maximum value for the sum of the k-values to be tested
  rsUint32 k[mMax];

  // example from http://en.wikipedia.org/wiki/Multinomial_theorem - number of distinct
  // permutations of the letters of the word "MISSISSIPPI":
  rsUint32 m, coeff;
  m    =  4;  // number of distinct letters (== number of indices for k)
  k[0] =  1;  // number of Ms
  k[1] =  4;  // number of Is
  k[2] =  4;  // number of Ss
  k[3] =  2;  // number of Ps
  coeff = rsMultinomialCoefficientUpTo12(k, m);
  testResult &= (coeff == 34650);
  coeff = rsMultinomialCoefficient(k, m);
  testResult &= (coeff == 34650);

  // these nested loops compute all multinomial coefficients for m=2, m=3, m=4, m=5 by means of two
  // different algorithms (naive and optimized) and compare the results:
  for(rsUint32 k1 = 0; k1 <= nMax; k1++)
  {
    k[0] = k1;
    for(rsUint32 k2 = 0; k2 <= nMax; k2++)
    {
      k[1] = k2;
      if( rsArrayTools::sum(k, 2) <= nMax )
        testResult &= rsMultinomialCoefficient(k, (rsUint32)2) == rsMultinomialCoefficientUpTo12(k, (rsUint32)2);
      for(rsUint32 k3 = 0; k3 <= nMax; k3++)
      {
        k[2] = k3;
        if( rsArrayTools::sum(k, 3) <= nMax )
          testResult &= rsMultinomialCoefficient(k, (rsUint32)3) == rsMultinomialCoefficientUpTo12(k, (rsUint32)3);
        for(rsUint32 k4 = 0; k4 <= nMax; k4++)
        {
          k[3] = k4;
          if( rsArrayTools::sum(k, 4) <= nMax )
            testResult &= rsMultinomialCoefficient(k, (rsUint32)4) == rsMultinomialCoefficientUpTo12(k, (rsUint32)4);
          for(rsUint32 k5 = 0; k5 <= nMax; k5++)
          {
            k[4] = k5;
            if( rsArrayTools::sum(k, 5) <= nMax )
              testResult &= rsMultinomialCoefficient(k, (rsUint32)5) == rsMultinomialCoefficientUpTo12(k, (rsUint32)5);
          }
        }
      }
    }
  }

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testMultinomialFormula(std::string &reportString)
{
  std::string testName = "rsMultinomialFormula";
  bool testResult = true;

  // actually, this is just an experiment, if we can verify the multinomial theorem by means of
  // an example
  static const rsUint32 m = 4;
  rsInt32 x[m] = {5, -1, 3, -2};
  //rsUint32 n   = 6;
  rsInt32 n   = 6;
  int t = rsPowInt(rsArrayTools::sum(x, m), n); // target value: (5 - 1 + 3 - 2)^6
  int sum = 0;                             // accumulator for the terms

  // accumulate the 4-fold sum:
  for(int k1 = 0; k1 <= n; k1++)
  {
    for(int k2 = 0; k2 <= n; k2++)
    {
      for(int k3 = 0; k3 <= n; k3++)
      {
        for(int k4 = 0; k4 <= n; k4++)
        {
          if( k1 + k2 + k3 + k4 == n )
          {
            int k[m]; k[0] = k1; k[1] = k2; k[2] = k3; k[3] = k4;
            sum += rsMultinomialCoefficient(k, 4) * rsPowInt(x[0], k[0]) * rsPowInt(x[1], k[1])
                                                  * rsPowInt(x[2], k[2]) * rsPowInt(x[3], k[3]);
          }
        }
      }
    }
  }
  testResult &= sum == t;

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testLeviCivita(std::string &reportString)
{
  std::string testName = "rsLeviCivita";
  bool testResult = true;

  // test Levi-Civita symbol in 2 dimensions (exhaustively):
  int i2[2];
  i2[0] = 1; i2[1] = 1; testResult &= rsLeviCivita(i2, 2) ==  0;
  i2[0] = 1; i2[1] = 2; testResult &= rsLeviCivita(i2, 2) == +1;
  i2[0] = 2; i2[1] = 1; testResult &= rsLeviCivita(i2, 2) == -1;
  i2[0] = 2; i2[1] = 2; testResult &= rsLeviCivita(i2, 2) ==  0;

  // test Levi-Civita symbol in 3 dimensions (exhaustively):
  int i3[3];
  i3[0] = 1; i3[1] = 1; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 1; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 1; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 2; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 2; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 2; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) == +1;
  i3[0] = 1; i3[1] = 3; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 1; i3[1] = 3; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) == -1;
  i3[0] = 1; i3[1] = 3; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 1; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 1; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 1; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) == -1;
  i3[0] = 2; i3[1] = 2; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 2; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 2; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 3; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) == +1;
  i3[0] = 2; i3[1] = 3; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 2; i3[1] = 3; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 1; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 1; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) == +1;
  i3[0] = 3; i3[1] = 1; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 2; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) == -1;
  i3[0] = 3; i3[1] = 2; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 2; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 3; i3[2] = 1; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 3; i3[2] = 2; testResult &= rsLeviCivita(i3, 3) ==  0;
  i3[0] = 3; i3[1] = 3; i3[2] = 3; testResult &= rsLeviCivita(i3, 3) ==  0;

  // test Levi-Civita symbol in 20 dimensions probabilistically by swapping two randomly selected
  // indices at a time and checking if the result of the functions changes the sign with each of
  // these swaps:
  int i20[20];
  rsArrayTools::fillWithRangeLinear(i20, 20, 1, 20);
  int target = 1;
  rsRandomUniform(1.0, 20.0, 1);
  for(int i = 1; i <= 100; i++)
  {
    testResult &= rsLeviCivita(i20, 20) == target;
    int k1 = (int) rsRound(rsRandomUniform(0.5001, 20.4999));
    int k2 = (int) rsRound(rsRandomUniform(0.5001, 20.4999));
    if( k1 != k2 )
    {
      rsSwap(i20[k1-1], i20[k2-1]);
      target *= -1;
    }
  }

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

int rsSignedMod(int a, int b)
{
  return a % b; // this perhaps only works with msc? -> test! nope - works in gcc, too
}
// maybe move to rapt ..compare to srem (used in the elliptic function code, i think)

// gcd, lcm, egcd, ...
bool testGcdAndCo()
{
  bool r = true;

  int y;

  // test behavior of signed modulo operation:
  y = rsSignedMod(+23, +10); r &= y == +3;
  y = rsSignedMod(+23, -10); r &= y == +3;
  y = rsSignedMod(-23, +10); r &= y == -3;
  y = rsSignedMod(-23, -10); r &= y == -3;


  y = rsGcd(12, 15); r &= y == 3;
  y = rsLcm(36, 48); r &= y == 144;

  return r;
}


bool testIntegerFunctions()
{
  std::string testName = "rsIntegerFunctions";
  std::string dummy;
  bool testResult = true;

  testResult &= testIntAbs(dummy);
  testResult &= testSumAndProduct();
  testResult &= testBinomialCoefficients(dummy);
  testResult &= testMultinomialCoefficients(dummy);
  testResult &= testMultinomialFormula(dummy);
  testResult &= testLeviCivita(dummy);
  testResult &= testGcdAndCo();

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
