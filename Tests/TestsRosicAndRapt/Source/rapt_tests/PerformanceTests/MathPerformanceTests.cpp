using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

void fftPerformance()
{
  // mostly to test the PerformanceAnalyzer class


  // create dummy data to perform the tests on:
  //std::vector<int> sizes = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
  std::vector<int> sizes = { 1, 2, 4, 8, 16, 32, 64, 128};
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

  // run tests and print report:
  pa.runTests();
  std::string report = pa.getReport();
  std::cout << report;

  // maybe plot means with error bars for variances, maybe also plot mins/maxes/medians/modes
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
void simdPerformance(TScalar scl, TVector vec)
{
  static const int N = 5000;  // number of vector operations

  TScalar zeroS = 0;
  TScalar oneS  = 1;
  TScalar accuS = 0;
  TVector zeroV = 0;
  TVector oneV  = 1;
  TVector accuV = 0;

  ::ProcessorCycleCounter counter;
  //PerformanceCounter counter;
  double cycles;
  double k = 1.0/(2*N);
  int n;

  // Print the number of cycles per scalar addition - in the case of vector types, we expect to see
  // the number to be a factor 2 smaller than in the case of scalar types (because we get two
  // scalar additions for each vector addition):

  // scalar = scalar + scalar:
  counter.init();
  for(n = 0; n < 2*N; n++)
    accuS = accuS + oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(accuS);
  printPerformanceTestResult("scl1 = scl1 + scl2", k*cycles);

  // vector = vector + vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = accuV + oneV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = vec1 + vec2", k*cycles);

  // vector = vector + scalar:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = accuV + oneS;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = vec1 + scl2", k*cycles);

  // vector = scalar + vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = oneS + accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = scl1 + vec1", k*cycles);

  // vector = scalar - vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = oneS - accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = scl1 - vec1", k*cycles);

  // vector = vector - vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = oneV - accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = vec2 - vec1", k*cycles);

  // vector = vector * vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = oneV * accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = vec2 * vec1", k*cycles);

  // vector = vector / vector:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = oneV / accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = vec2 / vec1", k*cycles);

  // unary minus:
  counter.init();
  for(n = 0; n < N; n++)
    accuV = -accuV;
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&accuV);
  printPerformanceTestResult("vec1 = -vec1      ", k*cycles);


  TVector x = 10;

  // clip:
  counter.init(); for(n = 0; n < N; n++) x = rsClip(x, -1, 1);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("clip", k*cycles);

  // abs:
  counter.init(); for(n = 0; n < N; n++) x = rsAbs(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("abs ", k*cycles);

  // sign:
  counter.init(); for(n = 0; n < N; n++) x = rsSign(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sign", k*cycles);

  // sqrt:
  counter.init(); for(n = 0; n < N; n++) x = rsSqrt(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sqrt", k*cycles);

  // exp:
  counter.init(); for(n = 0; n < N; n++) x = rsExp(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("exp ", k*cycles);

  // log:
  counter.init(); for(n = 0; n < N; n++) x = rsLog(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("log ", k*cycles);

  // sin:
  counter.init(); for(n = 0; n < N; n++) x = rsSin(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("sin ", k*cycles);

  // cos:
  counter.init(); for(n = 0; n < N; n++) x = rsCos(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("cos ", k*cycles);

  // tan:
  counter.init(); for(n = 0; n < N; n++) x = rsTan(x);
  cycles = (double)counter.getNumCyclesSinceInit();
  dontOptimize(&x); printPerformanceTestResult("tan ", k*cycles);



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
template void simdPerformance(double, rsFloat64x2);
template void simdPerformance(float, rsFloat32x4);


void rsSinCos1(double x, double* s, double* c)
{
  // it seems like this code is slower than std::cos/sin
  // taken from: http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/
  // low precision version:

  // always wrap input angle to -PI..PI:
  if (x < -3.14159265)
    x += 6.28318531;
  else
    if (x >  3.14159265)
      x -= 6.28318531;

  // compute sine:
  if (x < 0)
    *s = 1.27323954 * x + 0.405284735 * x * x;
  else
    *s = 1.27323954 * x - 0.405284735 * x * x;

  // compute cosine: sin(x + PI/2) = cos(x)
  x += 1.57079632;
  if (x >  3.14159265)
    x -= 6.28318531;
  if(x < 0)
    *c = 1.27323954 * x + 0.405284735 * x * x;
  else
    *c = 1.27323954 * x - 0.405284735 * x * x;
}

void sinCosPerformance()
{
  //static const int N = 10000;  // number of values
  static const int N = 128;  // number of values
  float xMin = 0.f;
  float xMax = float(2*PI);

  float x[N], ySin[N], yCos[N];
  rsArrayTools::fillWithRandomValues(x, N, xMin, xMax, 0);
  PerformanceCounterQPC counter;
  int n;

  rsSinCosTable<float>  table(1024);
  rsSinCosTable<double> tableD(1024);

  // measure cost of sin/cos standard library functions:
  counter.init();
  for(n = 0; n < N; n++)
  {
    ySin[n] = sin(x[n]);
    yCos[n] = cos(x[n]);
  }
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Standard library, float", cycles / N);

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
  double xD[N], ySinD[N], yCosD[N];
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

  // Conclusion:
  // rsSinCos1 is slower than using std::sin/cos, and the table gives results that can't be real
  // (i get values around 0.3 cycles per sin/cos pair). Something must be wrong with the test code.
}
