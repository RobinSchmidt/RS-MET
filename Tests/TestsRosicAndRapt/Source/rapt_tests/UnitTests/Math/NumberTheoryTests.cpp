#include "NumberTheoryTests.h"

bool testNumberTheory()
{
  std::string testName = "NumberTheory";
  std::string dummy;
  bool testResult = true;

  testResult &= testPrimeTableGeneration(dummy);
  testResult &= testPrimeFactorization(  dummy);
  testResult &= testNumberTheoryMisc(    dummy);

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testPrimeTableGeneration(std::string &reportString)
{
  std::string testName = "PrimeTableGeneration";
  bool testResult = true;






  // create table of primes up to 1000 and check against this target-table:
  static const rsUint32 np = 168;
  rsUint32 tp[np] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
                     103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
                     199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
                     313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
                     433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
                     563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
                     673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
                     811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
                     941,947,953,967,971,977,983,991,997};

  rsUint32 i;
  std::vector<rsUint32> pa;
  rsFindPrimesUpTo(pa, (rsUint32)1000);
  testResult &= pa.size() == (rsInt32)np;
  for(i = 0; i < np; i++)
    testResult &= pa[i] == tp[i];

/*  // commented because doesn't compile with gcc on windows - will cause this test to fail -
    // uncomment and fix underlying compilation problem

  rsUint32 p[np];
  rsFillPrimeTable(p, np, (rsUint32)29);
  //rsFillPrimeTable(p, np, 32);
  for(i = 0; i < np; i++)
    testResult &= p[i] == tp[i];

  // assume rsFindPrimesUpTo is correct - check alternative version against it:
  rsUint32 maxPrime = 1000000;
  rsFindPrimesUpTo(pa, maxPrime);
  rsUint32 numPrimes = pa.getNumElements();

  rsUint32 *p2 = new rsUint32[numPrimes];
  rsFillPrimeTable(p2, numPrimes, 1024);
  for(i = 0; i < numPrimes; i++)
  {
    testResult &= p2[i] == pa[i];
    rsAssert(testResult == true);
  }
  delete[] p2;
*/

  // suppose, we want to find all primes between 20000 and 21000
  // the table with number, multiples of which should be crossed out should have its largest
  // element >= 21000-20000 = 1000
  // 131 * 157 = 20567
  // or 41*43 = 1763 -> find all primes between 1700 and 1800

  // just for fun:
  //rsFindPrimesUpTo(pa, (rsUint32)1000000);
  //rsFindPrimesUpTo(pa, (rsUint32)10000000);
  //rsFindPrimesUpTo(pa, (rsUint32)4294967295);
  //for(int i = 0; i < pa.getNumElements(); i++)
  //  printf("%d %s", pa[i], " ");

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

// reconstructs a numkber from its prime factorization (which can be obtained by rsPrimeFactors) 
// ...maybe move to RAPT library:
rsUint32 rsPrimeProduct(std::vector<rsUint32>& factors, std::vector<rsUint32>& exponents)
{
  rsAssert(factors.size() == exponents.size());
  rsUint32 r = 1;
  for(size_t i = 0; i < factors.size(); i++)
    r *= rsPowInt(factors[i], exponents[i]); // todo: switch algo for rsPowInt...see implementation comment
  return r;
}

bool testPrimeFactorization(std::string &reportString)
{
  std::string testName = "PrimeFactorization";
  bool testResult = true;


  // factor all numbers from 0 to N and reconstruct them (todo: include negative numbers, too):
  std::vector<rsUint32> p, f, e;
  rsUint32 N = 100;
  for(rsUint32 i = 0; i <= N; i++) {
    rsPrimeFactors(i, f, e);
    rsUint32 j = rsPrimeProduct(f, e); // crashes for i == 5
    testResult &= j == i;
  }

  /*
  // code below is factored out because it takes an unreasonably long time. i think, it doesn't 
  // hang but legitimately takes a long time to factor these larger primes since i fixed the bug
  // with taking the integer square-root twice in rsPrimeFactors ...figure out...

  // factor some larger example numbers:
  rsUint32 x = 507996720; // = 2^4 * 3^2 * 5^1 * 7^3 * 11^2 * 17^1
  rsPrimeFactors(x, f, e);
  testResult &= f.size() == 6;
  testResult &= f[0] ==  2 && e[0] == 4;
  testResult &= f[1] ==  3 && e[1] == 2;
  testResult &= f[2] ==  5 && e[2] == 1;
  testResult &= f[3] ==  7 && e[3] == 3;
  testResult &= f[4] == 11 && e[4] == 2;
  testResult &= f[5] == 17 && e[5] == 1;

  x = 431985125; // = 5^3 * 11^2 * 13^4
  rsPrimeFactors(x, f, e);
  testResult &= f.size() == 3;
  testResult &= f[0] ==  5 && e[0] == 3;
  testResult &= f[1] == 11 && e[1] == 2;
  testResult &= f[2] == 13 && e[2] == 4;
  */

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


// template instantiation-wrappers to make the template-functions usable via function-pointers:
rsInt32 modularInverseInst(rsInt32 x, rsInt32 m)
{
  return rsModularInverse(x, m);
}
rsInt32 primeModularInverseInst(rsInt32 x, rsInt32 m)
{
  return rsPrimeModularInverse(x, m);
}
rsInt32 primeModularInverseInst2(rsInt32 x, rsInt32 m)
{
  return rsPrimeModularInverse2(x, m);
}
// test-function for different implementations of the modular inverse:
bool testModularInverse(rsInt32 (*pModularInverse)(rsInt32, rsInt32), bool testNonPrimeModuli)
{
  // values coprime to modulus m should have an inverse, otherwise the function should return 0

  bool testResult = true;

  rsInt32 m, x, y;
  m =  7;
  x =  1; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  x =  2; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  x =  3; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  x =  4; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  x =  5; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  x =  6; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  m =  3;
  x = 35; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;

  if( testNonPrimeModuli == true )
  {
    m = 12;
    x =  0; y = pModularInverse(x, m); testResult &= y == 0;
    x =  1; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  2; y = pModularInverse(x, m); testResult &= y == 0;
    x =  3; y = pModularInverse(x, m); testResult &= y == 0;
    x =  4; y = pModularInverse(x, m); testResult &= y == 0;
    x =  5; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  6; y = pModularInverse(x, m); testResult &= y == 0;
    x =  7; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  8; y = pModularInverse(x, m); testResult &= y == 0;
    x =  9; y = pModularInverse(x, m); testResult &= y == 0;
    x = 10; y = pModularInverse(x, m); testResult &= y == 0;
    x = 11; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 12; y = pModularInverse(x, m); testResult &= y == 0;
    x = 13; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 14; y = pModularInverse(x, m); testResult &= y == 0;
    x = 15; y = pModularInverse(x, m); testResult &= y == 0;
    x = 16; y = pModularInverse(x, m); testResult &= y == 0;
    x = 17; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 18; y = pModularInverse(x, m); testResult &= y == 0;
    x = 19; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 20; y = pModularInverse(x, m); testResult &= y == 0;
    x = 21; y = pModularInverse(x, m); testResult &= y == 0;
    x = 22; y = pModularInverse(x, m); testResult &= y == 0;
    x = 23; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;

    // todo: test x <= 0 - in this case, the inversion function should add m to x unti it's
    // >= 0

    m = 15;
    x =  1; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  2; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  3; y = pModularInverse(x, m); testResult &= y == 0;
    x =  4; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  5; y = pModularInverse(x, m); testResult &= y == 0;
    x =  6; y = pModularInverse(x, m); testResult &= y == 0;
    x =  7; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  8; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x =  9; y = pModularInverse(x, m); testResult &= y == 0;
    x = 10; y = pModularInverse(x, m); testResult &= y == 0;
    x = 11; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 12; y = pModularInverse(x, m); testResult &= y == 0;
    x = 13; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
    x = 14; y = pModularInverse(x, m); testResult &= (x * y) % m == 1;
  }
  return testResult;
}

bool testNumberTheoryMisc(std::string &reportString)
{
  std::string testName = "NumberTheoryMisc";
  bool testResult = true;

  rsInt32 a, b, g, x, y;

  // extended Euclid algorithm:
  x = 42; y = 30; rsEGCD(x, y, a, b, g);        // -2*42 + 3*30 = -84 + 90 = 6 = gcd(42, 30)
  testResult &= a == -2 && b == 3 && g == 6;
  x = 30; y = 42; rsEGCD(x, y, a, b, g);        // 3*30 - 2*42 = 90 - 84 = 6 = gcd(30, 42)
  testResult &= a == 3 && b == -2 && g == 6;
  x = 1; y = 0; rsEGCD(x, y, a, b, g);          // 1*1 + 0*0 = 1 + 0 = 1 = gcd(1, 0)
  testResult &= a == 1 && b == 0 && g == 1;
  x = 3; y = 0; rsEGCD(x, y, a, b, g);          // 1*3 + 0*0 = 3 + 0 = 3 = gcd(3, 0)
  testResult &= a == 1 && b == 0 && g == 3;

  // modular power:
  testResult &= rsModularPow(3, 10, 0) == 1;
  testResult &= rsModularPow(3, 10, 1) == 3;
  testResult &= rsModularPow(3, 10, 2) == 9;
  testResult &= rsModularPow(3, 10, 3) == 7;
  testResult &= rsModularPow(3, 10, 4) == 1;
  testResult &= rsModularPow(3, 10, 5) == 3;
  testResult &= rsModularPow(3, 10, 6) == 9;
  testResult &= rsModularPow(3, 10, 7) == 7;

  // modular inversion:
  testResult &= testModularInverse(&modularInverseInst,       true);
  testResult &= testModularInverse(&primeModularInverseInst,  false);
  testResult &= testModularInverse(&primeModularInverseInst2, false);

  // chinese remainder theorem:
  rsInt32 remainders[3] = {2, 2, 6};
  rsInt32 moduli[3]     = {3, 5, 7};
  rsInt32 R = rsChineseRemainderTheorem(remainders, moduli, 3);
  testResult &= R == 62;
  moduli[0]     = 25; moduli[1]     = 27;  moduli[2]     = 28;
  remainders[0] =  9; remainders[1] = 21;  remainders[2] =  2;
  R = rsChineseRemainderTheorem(remainders, moduli, 3);
  testResult &= R == 534;

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
