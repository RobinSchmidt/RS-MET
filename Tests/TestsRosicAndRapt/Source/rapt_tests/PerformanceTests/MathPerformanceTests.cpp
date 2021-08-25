using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

void fftPerformance()
{
  // mostly to test the PerformanceAnalyzer class


  // create dummy data to perform the tests on:
  //std::vector<int> sizes = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
  //std::vector<int> sizes = { 1, 2, 4, 8, 16, 32, 64, 128};
  std::vector<int> sizes = { 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
  size_t numSizes = sizes.size();
  int maxSize = sizes[numSizes-1];
  //std::vector<double> noise = createNoise(maxSize, -1.0, 1.0);
  std::vector<std::complex<double>> buf(maxSize);

  // define functors to be passed to performance analyzer:
  std::function<void(int)> dft = [&](int N){ rsDFT(&buf[0], N); };
  std::function<void(int)> fft = [&](int N){ rsFFT(&buf[0], N); };

  // set up the performance analyzer:
  PerformanceAnalyzer pa;
  pa.addTest(&dft, "DFT");
  pa.addTest(&fft, "FFT");
  pa.setTestInputSizes(sizes);
  // maybe instead of an array of input sizes, it should take an array of pointers to datasets

  // run tests and print report:
  pa.runTests();
  std::string report = pa.getReport();
  std::cout << report;
  pa.plotResults();  // todo...should plot means, variances, min/max, median and the raw data
                     // as scatter-plot

  // ToDo:
  // -Test various FFT routines, including FFTW and commercial ones. Switch the availability of 
  //  these in the library by #defines, a la RS_USE_GPL, RS_USE_MKL, RS_USE_IPP, etc., so client
  //  code can decide at compile time the (combination of) licensing and thereby make certain parts
  //  of the library un/available
  // -Maybe plot means with error bars for variances, maybe also plot mins/maxes/medians/modes. 
  //  That should actually be done in some sort of performance test framework.  Maybe a class 
  //  rsPerformance tester that takes a std::function which it is supposed to measure. Then, it 
  //  takes a bunch of of tests and computes statistics (mean, median, min, max, variance, etc.).
  //  Maybe, before taking mean and variance, outliers should be are removed, where the ourlierness
  //  is defined by being some factor below and above the median. And/or do scatter plots showing
  //  all the data raw.
}


void matrixAdressingTest()
{
  // Compares two matrix addressing schemes: using a flat array with pointer arithmetic vs. using
  // an pointer-to-pointer array. The test is to copy matrix values from one matrix into another.

  //typedef float Data;
  typedef double Data;

  size_t N = 70;     // number of rows
  size_t M = 20;     // number of columns
  size_t i, j;       // row and column indices

  // allocate flat matrices and row-pointers, fill a-matrix with random values:
  Data *af = new Data[N*M];   // a1 matrix as flat array
  Data *bf = new Data[N*M];
  Data **a = new Data*[N];    // a1 row pointers
  Data **b = new Data*[N];
  for(i = 0; i < N; i++){
    a[i] = &af[i*M];
    b[i] = &bf[i*M];
  }
  rsArrayTools::fillWithRandomValues(af, int(N*M), -1.0, +1.0, 0);

  // measure copying a into b via pointer-to-pointer access:
  PerformanceCounterTSC counter;
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

  // Hmm...having done these measurements the other day again, i got different results - now the
  // pointer-to-pointer version performing better. What now? ...also, the numbers are ridiculously
  // high anyway :-O ...whatever, i should use pointer arithmetic and flat storage for
  // compatibility with lapack routines anyway

  //delete[] af, bf, a, b; // nope - gives memory leak - we need to delete them all separately
  delete[] af;
  delete[] bf;
  delete[] a;
  delete[] b;
}

// this doesn't work because the functions are all inlined, so we can't get a function pointer
//inline void testPerformance(rsFloat64x2 (*func) (rsFloat64x2 x), const char* name, int N, rsFloat64x2 x)
//{
//  ProcessorCycleCounter counter;
//  double cycles;
//  double k = 1.0/(2*N);
//  counter.init();
//  for(int n = 0; n < N; n++)  x = func(x);
//  cycles = (double)counter.getNumCyclesSinceInit();
//  dontOptimize(&x);
//  printPerformanceTestResult(name, k*cycles);
//}

template<class TScalar, class TVector>
void simdPerformance(TScalar scl, TVector vec, const char* dataTypeName)
{
  static const int N = 1000;  // number of vector operations
  // maybe we should somehow make sure that the compiler doesn't know this, so it can't unroll
  // the loops?

  // ToDo:
  // -Split the file into:
  //  -arithmetic operators
  //  -comparison operators
  //  -math functions
  // -Create a data structure to record the measurements of several runs, so we may later do some
  //  statistical analysis - in particular: remove outliers, compute mean, median, etc.
  // -maybe create seperate project for the development of the simd class to reduce recompilation 
  //  times
  // -maybe write the data into a file and/or draw some plots

  TScalar zeroS(0);
  TScalar oneS(1);
  TScalar accuS(0);
  TVector zeroV(0);
  TVector oneV(1);
  TVector accuV(0);

  ::ProcessorCycleCounter counter;
  //PerformanceCounter counter;
  double cycles;
  double k = 1.0/N;
  int n;

  using STR = std::string;
  std::cout <<  STR("SIMD performance test for ") + dataTypeName + STR("\n");

  // Addition:
  // scalar = scalar + scalar:
  counter.init(); for(n = 0; n < N; n++) accuS = accuS + oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(accuS); printPerformanceTestResult("s1 = s1 + s2", k*cycles);

  // vector = scalar + vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneS + accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = s1 + v1", k*cycles);

  // vector = vector + scalar:
  counter.init(); for(n = 0; n < N; n++) accuV = accuV + oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v1 + s2", k*cycles);

  // vector = vector + vector:
  counter.init(); for(n = 0; n < N; n++) accuV = accuV + oneV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v1 + v2", k*cycles);


  // Subtraction:
  // scalar = scalar - scalar:
  counter.init(); for(n = 0; n < N; n++) accuS = accuS - oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(accuS); printPerformanceTestResult("s1 = s1 - s2", k*cycles);

  // vector = scalar - vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneS - accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = s1 - v1", k*cycles);

  // vector = vector - scalar:
  counter.init(); for(n = 0; n < N; n++) accuV = accuV - oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v1 - s2", k*cycles);

  // vector = vector - vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneV - accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v2 - v1", k*cycles);


  // Multiplication:
  // scalar = scalar * scalar:
  counter.init(); for(n = 0; n < N; n++) accuS = accuS * oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(accuS); printPerformanceTestResult("s1 = s1 * s2", k*cycles);

  // vector = scalar * vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneS * accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = s1 * v1", k*cycles);

  // vector = vector * scalar:
  counter.init(); for(n = 0; n < N; n++) accuV = accuV * oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v1 * s2", k*cycles);

  // vector = vector * vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneV * accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v2 * v1", k*cycles);


  // Division:
  // scalar = scalar / scalar:
  counter.init(); for(n = 0; n < N; n++) accuS = accuS / oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(accuS); printPerformanceTestResult("s1 = s1 / s2", k*cycles);

  // vector = scalar / vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneS / accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = s1 / v1", k*cycles);

  // vector = vector / scalar:
  counter.init(); for(n = 0; n < N; n++) accuV = accuV / oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v1 / s2", k*cycles);

  // vector = vector / vector:
  counter.init(); for(n = 0; n < N; n++) accuV = oneV / accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = v2 / v1", k*cycles);


  // Unary minus:
  // vector = -scalar
  counter.init(); for(n = 0; n < N; n++)  accuV = -accuS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = -s1    ", k*cycles);

  // vector = -vector
  counter.init(); for(n = 0; n < N; n++)  accuV = -accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV); printPerformanceTestResult("v1 = -v1    ", k*cycles);





  TVector x = 10;

  // Unary math functions:

  // Clip:
  //counter.init(); for(n = 0; n < N; n++) x = rsClip(x, TScalar(-1), TScalar(1));
  counter.init(); for(n = 0; n < N; n++) x = rsClip(x, TVector(-1), TVector(1));
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("clip", k*cycles);

  // Abs:
  counter.init(); for(n = 0; n < N; n++) x = rsAbs(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("abs ", k*cycles);

  // Sign:
  counter.init(); for(n = 0; n < N; n++) x = rsSign(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sign", k*cycles);

  // Sqrt:
  counter.init(); for(n = 0; n < N; n++) x = rsSqrt(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sqrt", k*cycles);

  // Exp:
  counter.init(); for(n = 0; n < N; n++) x = rsExp(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("exp ", k*cycles);

  // Log:
  counter.init(); for(n = 0; n < N; n++) x = rsLog(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("log ", k*cycles);

  // Sin:
  counter.init(); for(n = 0; n < N; n++) x = rsSin(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sin ", k*cycles);

  // Cos:
  counter.init(); for(n = 0; n < N; n++) x = rsCos(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("cos ", k*cycles);

  // Tan:
  counter.init(); for(n = 0; n < N; n++) x = rsTan(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("tan ", k*cycles);

  std::cout << "\n";


  //testPerformance(&rsAbs, "rsAbs ", N, -10.0);
  // maybe factor into function: testFunctionApplication (should perhaps be inline, so as to not
  // disturb the results by function call overhead)

  // Results:
  // approximate relative costs of operations: we set the cost of 1 addition = 1

  // rsFloat64x2:
  // add: 1, sub: 1.5, mul: 1, div: 9.5, unary minus: 1.5
  // interesting: mul costs the same as add but sub is more expensive
  // is this true also for scalar double?
  // rsClip: 1, rsSign: 3.5, rsAbs: 4
  // rsSqrt: 12.5, rsExp: 10, rsLog: 7, rsSin: 15, rsCos: 15, rsTan: 15

  // rsFloat32x4:
}

// Test instantiations for old versions:
template void simdPerformance(double, rsFloat64x2, const char*);
template void simdPerformance(float, rsFloat32x4, const char*);

// Test instantiations for new versions:
template void simdPerformance(float, rsSimdVector<float, 4>, const char*);
// does not yet compile, due to elementary math functions not yet defined


void simdPerformance()
{
  // Tests for old versions:
  //simdPerformance(1.0, rsFloat64x2(1.0), "rsFloat64x2");
  simdPerformance(1.f, rsFloat32x4(1.f), "rsFloat32x4");

  // Tests for new versions:
  //simdPerformance(1.f, rsSimdVector<float,  1>(1.f), "rsSimdVector<float, 1>");
  //simdPerformance(1.f, rsSimdVector<float,  2>(1.f), "rsSimdVector<float, 2>");
  simdPerformance(1.f, rsSimdVector<float,  4>(1.f), "rsSimdVector<float, 4>");
  simdPerformance(1.f, rsSimdVector<float,  8>(1.f), "rsSimdVector<float, 8>");
  //simdPerformance(1.f, rsSimdVector<float, 16>(1.f), "rsSimdVector<float, 16>");

  //simdPerformance(1.0, rsSimdVector<double,  1>(1.0), "rsSimdVector<double, 1>");
  //simdPerformance(1.0, rsSimdVector<double,  2>(1.0), "rsSimdVector<double, 2>");
  //simdPerformance(1.0, rsSimdVector<double,  4>(1.0), "rsSimdVector<double, 4>");
  //simdPerformance(1.0, rsSimdVector<double,  8>(1.0), "rsSimdVector<double, 8>");
  //simdPerformance(1.0, rsSimdVector<double, 16>(1.0), "rsSimdVector<double, 16>");

  // Expectations:
  // -Operations on rsSimdVector<float, 1> should use the same number of CPU cycles as the
  //  scalar operations. Code should boild down to a zero-cost abstraction of wrapping multiple 
  //  scalar operations into one
  // -operations on rsSimdVector<float, 4> should use the same amount of cpu as the old rsFloat32x4
  //  implementation.
  // -operations on rsSimdVector<float, 8> should use roughly twice as much cpu as on 
  //  rsSimdVector<float, 4> when no AVX is used. It should map to a zero-cost abstraction of using
  //  two __m128 variables
  //
  // ToDo:
  // -create a test porject that somehow lets us automatically do unit tests in different 
  //  configurations, i.e. with different macros defined...maybe run the same tests multiple times, 
  //  progressively defining more and more macros, like SSE, SSE2, AVX, AVX2, etc.
  // -maybe keep the old implementations around as prototype and reference, even when the new ones
  //  are ready for production
  // -compare the perfomance of the abstractions with those of the native types, like 
  //  rsSimdVector<float, 4> with using __m128 directly, etc.

  // Observations:
  // -for rsSimdVector<float, 8>, scalar + vector is a lot more expensive than vector + scalar
  //  -> implement both using the faster algo - it's commutative anyway
}


template<class T>
inline T rsPitchToFreqViaPow(T pitch)
{
  //return T(8.1757989156437073336828122976033 * exp(0.057762265046662109118102676788181 * pitch));
  return T(440)*( pow(T(2), (pitch-T(69))/T(12)) ); // naive, slower but numerically more precise
}

void sinCosPerformance()
{
  // Maybe rename to mathFuncsPerformance and add tests for other elementary functions here, such
  // as exp, log, tan, tanh, ... maybe also some special functions like tgamma

  //static const int N = 10000;  // number of values
  //static const int N = 128;  // number of values
  static const int N = 1024;  // number of values
  float xMin = float(-2*PI);
  float xMax = float(+2*PI);

  float  x[N],  ySin[N],  yCos[N];
  double xD[N], ySinD[N], yCosD[N];
  rsArrayTools::fillWithRandomValues(x, N, xMin, xMax, 0);
  //PerformanceCounterQPC counter;   // gives implausible result
  //PerformanceCounterPMC counter;   // trigger exception "No symbol file loaded"
  PerformanceCounterTSC counter;     // ...the good old Agner Fog implementation actually works
  int n;
  double cycles;

  int tableSize = 1024;
  rsSinCosTable<float>  table(tableSize);
  rsSinCosTable<double> tableD(tableSize);

  // measure cost of sin/cos standard library functions:

  // sin(float), cos(float):
  counter.init();
  for(n = 0; n < N; n++)
  {
    ySin[n]  = sin(x[n]);
    yCos[n]  = cos(x[n]);
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Standard library, float, sin,cos", cycles / N);

  // sin(doubl), cos(double):
  counter.init();
  for(n = 0; n < N; n++)
  {
    ySinD[n] = sin(x[n]);
    yCosD[n] = cos(x[n]);
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Standard library, double, sin,cos", cycles / N);


  // measure cost of rsSinCosTable<float> using rounding:
  counter.init();
  for(n = 0; n < N; n++)
    table.getValuesRounded(x[n], &ySin[n], &yCos[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, rounded, float", cycles / N);

  // measure cost of rsSinCosTable<float> using linear interpolation:
  counter.init();
  for(n = 0; n < N; n++)
    table.getValuesLinear(x[n], &ySin[n], &yCos[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, linear, float", cycles / N);

  // measure cost of rsSinCosTable<float> using cubic interpolation:
  counter.init();
  for(n = 0; n < N; n++)
    table.getValuesCubic(x[n], &ySin[n], &yCos[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, cubic, float", cycles / N);

  // measure cost of rsSinCosTable<double> using linear interpolation:

  rsArrayTools::fillWithRandomValues(xD, N, xMin, xMax, 0);
  counter.init();
  for(n = 0; n < N; n++)
    tableD.getValuesLinear(xD[n], &ySinD[n], &yCosD[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, linear, double", cycles / N);

  // measure cost of rsSinCos1:
  counter.init();
  for(n = 0; n < N; n++)
    rsSinCos1(x[n], &ySinD[n], &yCosD[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsSinCos1, double", cycles / N);

  counter.init();
  for(n = 0; n < N; n++)
    rsSinCos2(x[n], &ySinD[n], &yCosD[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsSinCos2, double", cycles / N);

  counter.init();
  for(n = 0; n < N; n++)
    rsSinCosApprox4(x[n], &ySinD[n], &yCosD[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsSinCosApprox4, double", cycles / N);


  // rsPitchToFreq:
  counter.init();
  for(n = 0; n < N; n++) 
    ySinD[n] = rsPitchToFreq(x[n]);  // rename ySinD to something more generic like y1D
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsPitchToFreq, double", cycles / N);

  counter.init();
  for(n = 0; n < N; n++) 
    ySinD[n] = rsPitchToFreqViaPow(x[n]);  // rename ySinD to something more generic like y1D
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("rsPitchToFreqViaPow, double", cycles / N);


  dontOptimize(&ySin);
  dontOptimize(&yCos);
  dontOptimize(&ySinD);
  dontOptimize(&yCosD);



  int dummy = 0;



  // Conclusions:
  // -rsSinCos1 is slower than using std::sin/cos, and the table gives results that can't be real
  //  (i get values around 0.3 cycles per sin/cos pair). Something must be wrong with the test 
  //  code.

  // ToDo:
  // -try to remove the code duplication
}
