#include "MathTests.h"

#undef min

inline double absFast(double x)
{
  static const unsigned long long mask = 0x7FFFFFFFFFFFFFFFULL; // binary: 011111...
  unsigned long long tmp = *(reinterpret_cast<unsigned long long*>(&x)) & mask;
  return *(reinterpret_cast<double*>(&tmp));
}
// in the performance test, it turned out that standard-fabs is actually faster.

void testAbsAndSign2(std::string &reportString)
{
  // Compares the performance of the standard library functions fabs(double x) and
  // rsAbsFast(double x)

  static const int N = 10000;   // number of tests
  double x[N];
  double y[N];
  RAPT::rsArrayTools::fillWithRandomValues(x, N, -N, N, 0);
  int n;

  ::PerformanceCounterTSC counter;
  counter.init();
  for(n = 0; n < N; n++)
    y[n] = fabs(x[n]);
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("standard fabs", cycles / N);  // 0.5 cycles (?)

  counter.init();
  for(n = 0; n < N; n++)
    y[n] = absFast(x[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("absFast", cycles / N);      // 1.5 cycles - actually slower

  counter.init();
  for(n = 0; n < N; n++)
    y[n] = copysign(x[n], 1.0);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("copysign", cycles / N);       // 35 cycles

  counter.init();
  for(n = 0; n < N; n++)
    y[n] = rsSign(x[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsSign", cycles / N);           // 0.43 cycles
}

void testMultinomialCoefficients2(std::string &reportString)
{
  static const rsUint32 mMax = 5;  // maximum number of indices for the k-values
  static const rsUint32 nMax = 12; // maximum value for the sum of the k-values to be tested
  rsUint32 k[mMax];
  rsUint32 result;

  ::PerformanceCounterTSC counter;
  counter.init();

  // we use the same nested loop computation as in the unit test - refer to the comments there for
  // more info:
  for(rsUint32 k1 = 0; k1 <= nMax; k1++)
  {
    k[0] = k1;
    for(rsUint32 k2 = 0; k2 <= nMax; k2++)
    {
      k[1] = k2;
      if( RAPT::rsArrayTools::sum(k, 2) <= nMax )
        result = rsMultinomialCoefficient(k, (rsUint32)2);
      for(rsUint32 k3 = 0; k3 <= nMax; k3++)
      {
        k[2] = k3;
        if( RAPT::rsArrayTools::sum(k, 3) <= nMax )
          result = rsMultinomialCoefficient(k, (rsUint32)3);
        for(rsUint32 k4 = 0; k4 <= nMax; k4++)
        {
          k[3] = k4;
          if( RAPT::rsArrayTools::sum(k, 4) <= nMax )
            result = rsMultinomialCoefficient(k, (rsUint32)4);
          for(rsUint32 k5 = 0; k5 <= nMax; k5++)
          {
            k[4] = k5;
            if( RAPT::rsArrayTools::sum(k, 5) <= nMax )
              result = rsMultinomialCoefficient(k, (rsUint32)5);
          }
        }
      }
    }
  }
  double cycles = (double) counter.getNumCyclesSinceInit();

  printPerformanceTestResult("MultinomialCoefficients", cycles);
}

// from:
// https://code.google.com/p/primesieve/wiki/HowTo
const int L1D_CACHE_SIZE = 32768;
//typedef rsUint64 int64_t; // inserted by Robin Schmidt
//typedef rsInt64 int64_t; // inserted by Robin Schmidt
/// Generate primes using the segmented sieve of Eratosthenes.
/// This algorithm uses O(n log log n) operations and O(sqrt(n)) space.
/// @param limit         Sieve primes <= limit.
/// @param segment_size  Size of the sieve array in bytes.
///
void segmented_sieve(rsInt64 limit, int segment_size = L1D_CACHE_SIZE)
{
  int sqrt = (int) std::sqrt((double) limit);
  rsInt64 count = (limit < 2) ? 0 : 1;
  rsInt64 s = 2;
  rsInt64 n = 3;

  // vector used for sieving
  std::vector<char> segment(segment_size);

  // generate small primes <= sqrt
  std::vector<char> is_prime(sqrt + 1, 1);
  for (int i = 2; i * i <= sqrt; i++)
    if (is_prime[i])
      for (int j = i * i; j <= sqrt; j += i)
        is_prime[j] = 0;

  std::vector<int> primes;
  std::vector<int> next;

  for (rsInt64 low = 0; low <= limit; low += segment_size)
  {
    std::fill(segment.begin(), segment.end(), 1);

    // current segment = interval [low, high]
    rsInt64 high = std::min(low + segment_size - 1, limit);

    // store small primes needed to cross off multiples
    for (; s * s <= high; s++)
      if (is_prime[s])
      {
        primes.push_back((int) s);
          next.push_back((int)(s * s - low));
      }

    // segmented sieve of Eratosthenes
    for (std::size_t i = 1; i < primes.size(); i++)
    {
      int j = next[i];
      for (int k = primes[i] * 2; j < segment_size; j += k)
        segment[j] = 0;
      next[i] = j - segment_size;
    }
    for (; n <= high; n += 2)
      if (segment[n - low])
      {
        int p = (int)n; // added by robin schmidt - this is how the primes are retrieved
        count++;
      }
  }
  std::cout << count << " primes found." << std::endl;
}
// this prime sieve outperforms mine by a factor 2 - but i'm confused how to use it to actually
// pass the generated primes as return value.

void testPrimeSieves(std::string &reportString)
{
  //static const rsUint64 maxPrime = 100000000;
  static const rsUint64 maxPrime = 10000000;
  //static const rsUint64 maxPrime = 1000;
    // we should measure for different values of maxPrime - the algorithms may scale differently

  ::PerformanceCounterTSC counter;
  double cyclesPerNumber;

  //rsArrayTools<rsUint64> pa;
  std::vector<rsUint64> pa;
  counter.init();
  rsFindPrimesUpTo(pa, maxPrime);
  rsUint32 numPrimes = (int)pa.size();
  cyclesPerNumber = (double) counter.getNumCyclesSinceInit() / maxPrime;
  printPerformanceTestResult("Prime Sieve (Joerg Arndt)", cyclesPerNumber);

  rsUint32 *p = new rsUint32[numPrimes];
  //rsUint64 *p = new rsUint64[numPrimes];
  counter.init();
  //rsFillPrimeTable(p, numPrimes, 1024);
  //rsFillPrimeTable(p, numPrimes, 2048);
  //rsFillPrimeTable(p, numPrimes, 4096);
  //rsFillPrimeTable(p, numPrimes, 8192);
  //rsFillPrimeTable(p, numPrimes, 16384);
  rsFillPrimeTable(p, numPrimes, 32768);  // seems to be the sweet spot
  //rsFillPrimeTable(p, numPrimes, 65536);
  cyclesPerNumber = (double) counter.getNumCyclesSinceInit() / maxPrime;
  printPerformanceTestResult("Prime Sieve (Robin Schmidt)", cyclesPerNumber);
  delete[] p;

  counter.init();
  segmented_sieve(maxPrime, 16384);
  cyclesPerNumber = (double) counter.getNumCyclesSinceInit() / maxPrime;
  printPerformanceTestResult("Prime Sieve (Kim Walisch)", cyclesPerNumber);


  // seems mine is faster for smaller value of maxPrime and Arndt's is faster for greater values
  // with a crossover point around 20000000 for a buffer-size of 1024, using larger buffer-sizes
  // the crossover can be pushed further up
}

std::string matrixString(const rsMatrixDbl& A)
{
  /*
  rsString s;
  s += "Matrix(";
  s += rsString(A.getNumRows());
  s += "x";
  s += rsString(A.getNumColumns());
  s += "):";
  return s.getAsStdString(); // rename into toStdString
  */
  std::string s;
  s += "Matrix(";
  s += to_string(A.getNumRows());
  s += "x";
  s += to_string(A.getNumColumns());
  s += "):";
  return s;
}
void runMatrixTest(std::string &reportString, int numRows, int numColumns, int numRuns = 20)
{
  rsMatrixDbl A(numRows, numColumns);
  A.randomizeElements(-1.0, 1.0);

  ::PerformanceCounterTSC counter;
  counter.init();
  for(int i = 1; i <= numRuns; i++)
    A = 2.0 * ((A*trans(A))*A + A);  // A := 2 * ((A*A^T)*A + A)
  double cycles = (double) counter.getNumCyclesSinceInit();

  double cyclesPerElement = cycles / (numRows*numColumns*numRuns);

  rsAssert(false); // printing below doesn't compile
  //printPerformanceTestResult(matrixString(A), cyclesPerElement);
}
void testMatrix(std::string &reportString)
{
                                       // without lazy-copy | with lazy-copy
                                       // Release | Debug   | Release | Debug
  runMatrixTest(reportString,  2,  3); // 4500    | 14000   | 5400    | 13000
  runMatrixTest(reportString, 20, 30); // 1200    |  2000   |  900    |  3000
  runMatrixTest(reportString, 60, 90); //         |         | 2000    |  7000
}

void testMatrixAddressing(std::string &reportString)
{
  // Compares two matrix addressing schemes: using a flat array with pointer arithmetic vs. using
  // an pointer-to-pointer array. The test is to copy matrix values from one matrix into another.

  int N = 500;    // number of rows
  int M = 700;    // number of columns
  int i, j;      // row and column indices

  // allocate flat matrices and row-pointers, fill a-matrix with random values:
  double *af = new double[N*M];   // a1 matrix as flat array
  double *bf = new double[N*M];
  double **a = new double*[N];    // a1 row pointers
  double **b = new double*[N];
  for(i = 0; i < N; i++){
    a[i] = &af[i*M];
    b[i] = &bf[i*M];
    RAPT::rsArrayTools::fillWithRandomValues(af, N*M, 1.0, 2.0, 0);
  }

  // measure copying a into b via pointer-to-pointer access:
  ::PerformanceCounterTSC counter;
  counter.init();
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++)
      b[i][j] = a[i][j];
  }
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("pointer-to-pointer Matrix", cycles / (N*M));

  // measure copying a into b via pointer arithmetic:
  counter.init();
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++)
      bf[M*i+j] = af[M*i+j];
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("pointer-arithmetic Matrix", cycles / (N*M));

  // Results - cycles per element-access (in release build):
  // matrix size:            5x7 20x70  70x20  30x50  50x70  300x500
  // pointer-to-pointer:     20  11     10     9.1    12     22
  // pointer-arithmetic:     12   2      2     2.3    3.4    14

  // Conclusion:
  // Accessing matrix elements via pointer arithmetic is faster than using a pointer array
  // (i.e. pointer-to-pointer). Add to that, that the pointer array needs additional storage
  // space - the relative increase in storage-size is given by
  // sizeof(size_t) / (N*sizeof(ElementType)) - the conclusion is that using pointer arithmetic
  // is better - both, in terms of speed and storage space. The advantages seem greatest for
  // medium sized matrices (a couple of thousands of elements), but even for smaller and larger
  // matrices, pointer arithmetic beats pointer arrays.

  delete[] af, bf, a, b;
}

