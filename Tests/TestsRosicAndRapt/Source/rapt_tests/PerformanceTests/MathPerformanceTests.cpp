#include "MathPerformanceTests.h"

void matrixAdressingTest()
{
  // Compares two matrix addressing schemes: using a flat array with pointer arithmetic vs. using
  // an pointer-to-pointer array. The test is to copy matrix values from one matrix into another.

  typedef float Data;
  //typedef double Data;

  int N = 500;     // number of rows
  int M = 200;     // number of columns
  int i, j;       // row and column indices

  // allocate flat matrices and row-pointers, fill a-matrix with random values:
  Data *af = new Data[N*M];   // a1 matrix as flat array
  Data *bf = new Data[N*M];
  Data **a = new Data*[N];    // a1 row pointers
  Data **b = new Data*[N];
  for(i = 0; i < N; i++){
    a[i] = &af[i*M];
    b[i] = &bf[i*M];
  }
  ArrayTools::rsFillWithRandomValues(af, N*M, -1.0, +1.0, 0);

  // measure copying a into b via pointer-to-pointer access:
  ProcessorCycleCounter counter;
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

  //delete[] af, bf, a, b; // nope - gives memory leak - we need to delete them all separately
  delete[] af;
  delete[] bf;
  delete[] a;
  delete[] b;
}

void sinCosPerformance()
{
  static const int N = 5000;  // number of values
  float xMin = 0.f;
  float xMax = float(2*PI);

  float x[N], ySin[N], yCos[N];
  ArrayTools::rsFillWithRandomValues(x, N, xMin, xMax, 0);
  ProcessorCycleCounter counter;
  int n;

  // measure cost of sin/cos standard library functions:
  counter.init();
  for(n = 0; n < N; n++)
  {
    ySin[n] = sin(x[n]);
    yCos[n] = cos(x[n]);
  }
  double cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Standard library, float", cycles / N);

  // measure cost of rsSinCosTable<float> using linear interpolation:
  rsSinCosTableF table(1024);
  counter.init();
  for(n = 0; n < N; n++)
    table.getValuesLinear(x[n], &ySin[n], &yCos[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, linear, float", cycles / N);

  // measure cost of rsSinCosTable<float> using rounding:
  counter.init();
  for(n = 0; n < N; n++)
    table.getValuesRounded(x[n], &ySin[n], &yCos[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, rounded, float", cycles / N);




  double xD[N], ySinD[N], yCosD[N];
  ArrayTools::rsFillWithRandomValues(xD, N, xMin, xMax, 0);

  // measure cost of rsSinCosTable<double> using linear interpolation:
  rsSinCosTableD tableD(1024);
  counter.init();
  for(n = 0; n < N; n++)
    tableD.getValuesLinear(xD[n], &ySinD[n], &yCosD[n]);
  cycles = (double) counter.getNumCyclesSinceInit();
  printPerformanceTestResult("Table, linear, double", cycles / N);
}