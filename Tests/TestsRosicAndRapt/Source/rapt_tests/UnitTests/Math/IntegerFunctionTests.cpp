
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
  for(int k = 0; k <= n; k++)
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
  for(int k = 0; k <= n; k++)
    sum += rsProduct(k+1, k+m);
  return sum;
}
// using rsProduct, needed an index-shift in the product

int rsSumOfProducts3(int n, int m)
{
  return rsProduct(n+1, n+m+1) / (m+1);
}
// without the outer loop, see: https://www.youtube.com/watch?v=-2PA7SbWoJ0&t=27m
// move to IntegerFunctions.h


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

  // 1*2*3 + 2*3*4 + 3*4*5 + 4*5*6 = 210
  y = rsSumOfProducts1(2, 4);  // 504
  // seems like the 1st parameter is the number of terms minus 1 and the 2nd is the number of 
  // factors - more natural would be, if the first is the number of terms itself

  // n = 2, m = 3:
  // 1*2*3 + 2*3*4 + 3*4*5 = 90
  y = rsSumOfProducts1(2, 3);
  y = rsSumOfProducts2(2, 3);
  y = rsSumOfProducts3(2, 3);


  y = rsSumOfProducts1(3, 5); // 1*2*3*4 + 2*3*4*5 + 3*4*5*6 + 4*5*6*7 + 5*6*7*8 = 3024
  y = rsSumOfProducts2(3, 5); // 24      + 120     + 360     + 840     + 1680
  y = rsSumOfProducts3(3, 5);
  // hmm - they all return the same values, as they should - but what actually do they compute? the
  // 3,5 case gives a different number (10080) than i expect (3024), but in all 3 cases the same 
  // number -> figure out, what's going on...

  // are the parameters reversed?
  y = rsSumOfProducts1(5, 3); // nope!

  return r;
}

bool testBinomialCoefficients(std::string &reportString)
{
  std::string testName = "rsBinomialCoefficients";
  bool testResult = true;

  unsigned int nMax = 20;
  unsigned int c1, c2;
  for(unsigned int n = 0; n <= nMax; n++)
  {
    for(unsigned int k = 0; k <= n; k++)
    {
      c1 = rsBinomialCoefficient(      n, k);
      c2 = rsBinomialCoefficientUpTo20(n, k);
      testResult &= (c1 == c2);
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

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}