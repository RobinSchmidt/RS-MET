
using namespace RAPT;  // get rid

bool testBandDiagonalSolver()
{
  //std::string testName = "BandDiagonalSolver";
  bool testResult = true;

  static const int N  = 9;
  static const int kl = 3;
  static const int ku = 2;
  //static const int M  = kl+ku+1;  // # rows in band storare format

  // We solve the following 9x9 system A*x = b for x:
  // 11 12 13 00 00 00 00 00 00     1     74
  // 21 22 23 24 00 00 00 00 00     2    230
  // 31 32 33 34 35 00 00 00 00     3    505
  // 41 42 43 44 45 46 00 00 00     4    931
  // 00 52 53 54 55 56 57 00 00  *  5 = 1489
  // 00 00 63 64 65 66 67 68 00     6   2179
  // 00 00 00 74 75 76 77 78 79     7   3001
  // 00 00 00 00 85 86 87 88 89     8   3055
  // 00 00 00 00 00 96 97 98 99     9   2930

  // create right-hand-side vector b and target values for solution vector x
  double b[N] = { 74,230,505,931,1489,2179,3001,3055,2930 };  // the right hand side in A*x = b
  double x[N] = { 1,2,3,4,5,6,7,8,9 };                        // target values for x

  // create solver object and set up coefficient matrix:
  rsBandDiagonalSolver<double> solver(N, kl, ku);
  for(int k = -kl; k <= ku; k++) {
    for(int n = 0; n < N - abs(k); n++) {
      double val = (n+1)*10 + n+1;
      if(k < 0)
        val += -k*10;
      else
        val += k;
      solver.setDiagonalElement(k, n, val); // sets up the coeff matrix
    }
  }

  // solve the system with the 3 different solver routines:
  typedef rsBandDiagonalSolver<double>::Algorithm ALGO;
  double xGbsvxx[N], xGbsvx[N], xGbsv[N];                // results of the 3 algorithms
  solver.setAlgorithm(ALGO::gbsvxx);
  solver.solve(xGbsvxx, b, 1);
  solver.setAlgorithm(ALGO::gbsvx);
  solver.solve(xGbsvx, b, 1);
  solver.setAlgorithm(ALGO::gbsv);
  solver.solve(xGbsv, b, 1);

  // compute maximum errors in the 3 solutions and check if it is below some thresholds:
  double errGbsv   = rsArrayTools::maxDeviation(x, xGbsv, N);   // rename to maxDistance
  double errGbsvx  = rsArrayTools::maxDeviation(x, xGbsvx, N);
  double errGbsvxx = rsArrayTools::maxDeviation(x, xGbsvxx, N);
  testResult &= errGbsv    < 2e-9;
  testResult &= errGbsvx   < 3e-10;
  testResult &= errGbsvxx == 0.0;

  // This test fails with gcc. There, all 3 errors have the same value of around 5e-10. The
  // resulting vectors are all equal, as if the same algorithm was used in all 3 computations.
  // However, we indeed do call gbsv, gbsvx, gbsvxx respectively. It doesn't seem to have anything
  // to do with compiling for 64 vs 32 bit - the test passes with msc when compiling for 32 bit.

  return testResult;
}

bool testMatrix2x2()
{
  std::string testName = "Matrix2x2";
  bool testResult = true;
  typedef rsLinearAlgebra LA;


  // Test compuatation of eigenvalues and eigenvectors:

  double ev, ex, ey, d;
  bool nrm = false;   // normalize eigenvectors to unit length

  double tol = 1.e-16;

  ev = LA::eigenvalue2x2_1(-4.0, 6.0, -3.0, 5.0);                // -1
  LA::eigenvector2x2_1(    -4.0, 6.0, -3.0, 5.0, &ex, &ey, nrm); // k*(2,1) - is (1,0.5), so k=0.5
  testResult &= rsIsCloseTo(ev,   -1.0, tol);
  testResult &= rsIsCloseTo(ex/ey, 2.0, tol);

  ev = LA::eigenvalue2x2_2(-4.0, 6.0, -3.0, 5.0);                // 2
  LA::eigenvector2x2_2(    -4.0, 6.0, -3.0, 5.0, &ex, &ey, nrm); // k*(1,1), is (1,1), so k=1
  testResult &= rsIsCloseTo(ev,    2.0, tol);
  testResult &= rsIsCloseTo(ex/ey, 1.0, tol);

  ev = LA::eigenvalue2x2_1(4.0, 1.0, 6.0, 3.0);
  LA::eigenvector2x2_1(    4.0, 1.0, 6.0, 3.0, &ex, &ey, nrm);
  testResult &= rsIsCloseTo(ev,     1.0,  tol);
  testResult &= rsIsCloseTo(ex/ey, -1/3., tol);

  ev = LA::eigenvalue2x2_2(4.0, 1.0, 6.0, 3.0);
  LA::eigenvector2x2_2(    4.0, 1.0, 6.0, 3.0, &ex, &ey, nrm);
  testResult &= rsIsCloseTo(ev,    6., tol);
  testResult &= rsIsCloseTo(ex/ey, 1/2., tol);

  // try special cases with b = 0:

  // b = 0, a < d:

  ev = LA::eigenvalue2x2_1(-4.0, 0.0, -3.0, 5.0);                // -4
  LA::eigenvector2x2_1(    -4.0, 0.0, -3.0, 5.0, &ex, &ey, nrm); // (1,1/3)
  testResult &= rsIsCloseTo(ev,   -4.0, tol);
  testResult &= rsIsCloseTo(ex/ey,  3., tol);

  ev = LA::eigenvalue2x2_2(-4.0, 0.0, -3.0, 5.0);                // 5
  LA::eigenvector2x2_2(    -4.0, 0.0, -3.0, 5.0, &ex, &ey, nrm); // (0,1)
  testResult &= rsIsCloseTo(ev, 5.0, tol)
             && rsIsCloseTo(ex, 0.0, tol)
             && rsIsCloseTo(ey, 1.0, tol); // maybe when we normalize, check only for ey > 0


  // b = 0, a > d:

  ev = LA::eigenvalue2x2_1(-4.0, 0.0, -3.0, -5.0);                // -5
  LA::eigenvector2x2_1(    -4.0, 0.0, -3.0, -5.0, &ex, &ey, nrm); // (0,1)
  testResult &= rsIsCloseTo(ev, -5.0, tol)
             && rsIsCloseTo(ex,  0.0, tol)
             && rsIsCloseTo(ey,  1.0, tol);

  ev = LA::eigenvalue2x2_2(-4.0, 0.0, -3.0, -5.0);                // -4
  LA::eigenvector2x2_2(    -4.0, 0.0, -3.0, -5.0, &ex, &ey, nrm); // (1,-3)
  testResult &= rsIsCloseTo(ev,    -4.0,  tol);
  testResult &= rsIsCloseTo(ex/ey, -1/3., tol);


  // b = 0, a = d:

  ev = LA::eigenvalue2x2_1(-4.0, 0.0, -3.0, -4.0);                // -4
  LA::eigenvector2x2_1(    -4.0, 0.0, -3.0, -4.0, &ex, &ey, nrm); // (0,1)
  testResult &= rsIsCloseTo(ev, -4.0, tol)
             && rsIsCloseTo(ex,  0.0, tol)
             && rsIsCloseTo(ey,  1.0, tol);

  ev = LA::eigenvalue2x2_2(-4.0, 0.0, -3.0, -4.0);                // -4
  LA::eigenvector2x2_2(    -4.0, 0.0, -3.0, -4.0, &ex, &ey, nrm); // (0,1)
  testResult &= rsIsCloseTo(ev, -4.0, tol)
             && rsIsCloseTo(ex,  0.0, tol)
             && rsIsCloseTo(ey,  1.0, tol);

  // todo: try the special cases with complex numbers - we may also need tolerances in the
  // comparisons
  // the current implementaion will fall into the wrong branch when real parts of a and d
  // are equal

  typedef rsMatrix2x2<double> Mat;
  rsMatrix2x2<double> A(1,2,3,4), B(5,6,7,8), C, D;
  C = A+B; testResult &= C == Mat(  6,  8, 10, 12);
  C = A-B; testResult &= C == Mat( -4, -4, -4, -4);
  C = A*B; testResult &= C == Mat( 19, 22, 43, 50);
  C = B*A; testResult &= C == Mat( 23, 34, 31, 46);
  C = Mat::commutator(A,B); testResult &= C == Mat( -4, -12, 12, 4); // A*B - B*A
  d = A.getDeterminant(); testResult &= d == -2;
  d = B.getDeterminant(); testResult &= d == -2;
  C = A.getInverse(); testResult &= C == Mat(-2, 1, 3./2, -1./2);
  C = B.getInverse(); testResult &= C == Mat(-4, 3, 7./2, -5./2);
  C = A/B; testResult &= C == Mat( 3, -2, 2, -1);
  C = B/A; testResult &= C == Mat(-1,  2, -2, 3);
  C = A*A*A*A*A;
  D = A.getPower(5);
  testResult &= C == D; // maybe also do a comparison with tolerance

  // test getPower with a matrix with 2 equal eigenvalues:
  // https://web.calpoly.edu/~brichert/teaching/oldclass/f2002217/solutions/solutions9.pdf
  // we must have a = d and b*c = 0
  A = Mat({5,0,3,5});
  C = A*A*A*A*A;
  D = A.getPower(5);
  testResult &= C == D; // maybe also do a comparison with tolerance


  /*
  // now, let's try some complex matrices
  typedef std::complex<double> Cmplx;
  typedef rsMatrix2x2<Cmplx>   MatC;
  typedef rsVector2D<Cmplx>    VecC;
  std::complex<double> one(1,0), zero(0,0), i(0,1);
  MatC Z(one, zero, zero, -one); // the Pauli z matrix
  Cmplx c1, c2;  // eigenvalues
  VecC  C1, C2;  // eigenvectors
  c1 = Z.eigenvalue1();    // -1
  c2 = Z.eigenvalue2();    // (1,0)
  C1 = Z.eigenvector1();   // +1
  C2 = Z.eigenvector2();   // (0,1)
  // eigenvectors are swapped

  // try it with the real version of pauli Z

  typedef rsVector2D<double> Vec;
  A = Mat(1,0,0,1);
  double r1, r2;
  Vec R1, R2;
  r1 = A.eigenvalue1();    // -1
  r2 = A.eigenvalue2();    // (1,0)
  R1 = A.eigenvector1();   // +1
  R2 = A.eigenvector2();   // (0,1)

  // what? actually sage says that (0,1) belongs to +1 and (1,0) belongs to -1 indeed:
  //A = matrix([[1, 0], [0, -1]])
  //A.eigenvectors_right()
  */

  return testResult;
}

bool testLinearSystem2x2()
{
  std::string testName = "LinearSystem2x2";
  bool testResult = true;
  typedef rsLinearAlgebra LA;

  double x[2];
  double y[2]    = {17, 39};
  double A[2][2] = {{1, 2},
                    {3, 4}};

  LA::rsSolveLinearSystem2x2(A, x, y);
  testResult &= (x[0] == 5.0);
  testResult &= (x[1] == 6.0);

  return testResult;
}

bool testLinearSystem3x3()
{
  std::string testName = "LinearSystem3x3";
  bool testResult = true;

  double x[3];
  double y[3]    = {-4, 15, 11};
  double A[3][3] = {{ 1, 2, -3},
                    {-5, 7,  2},
                    {-2, 2,  3}};

  rsLinearAlgebra::rsSolveLinearSystem3x3(A, x, y);

  testResult &= (x[0] == 1.0);
  testResult &= (x[1] == 2.0);
  testResult &= (x[2] == 3.0);

  return testResult;
}

bool testLinearSystemViaGauss()
{
  std::string testName = "LinearSystemViaGauss";
  bool testResult = true;

  typedef std::complex<double> Complex;

  Complex x[3];
  Complex y[3]    = { Complex(-2,-6), Complex(-19,+16), Complex(+44,-37)};
  Complex A[3][3] = {{Complex(+1,+2), Complex( +2,-1),  Complex( -3, +1)},
                     {Complex(+4,-2), Complex( -6,+1),  Complex( +5, -2)},
                     {Complex(-9,+1), Complex( +7,-2),  Complex( +8, +2)}};
  Complex *pA[3];
  for(int i = 0; i < 3; i++)
    pA[i] = &A[i][0]; // we need a pointer-to-pointer for the solver

  rsLinearAlgebra::rsSolveLinearSystem(pA, x, y, 3);

  Complex xt[3] = {Complex(-2,+1), Complex(+3,-2), Complex(+1,-1)};
  double tol = 1.e-14;
  testResult &= abs(x[0]-xt[0]) < tol;
  testResult &= abs(x[1]-xt[1]) < tol;
  testResult &= abs(x[2]-xt[2]) < tol;

  return testResult;
}

bool testGaussJordanInversion()
{
  std::string testName = "GaussJordanInversion";
  bool testResult = true;

  double A[3][3] = {{2,  1, 4},
                    {3, 10, 3},
                    {1,  5, 1}};
  double *pA[3];
  for(int i = 0; i < 3; i++)
    pA[i] = &A[i][0];

  rsLinearAlgebra::rsInvertMatrix(pA, 3);

  // create the correct inverse matrix:
  double Ai[3][3] = {{-0.5, +1.9, -3.7},
                     {+0.0, -0.2, +0.6},
                     {+0.5, -0.9, +1.7}};

  // compare - note that we must use the pointer pA instead of A itself because the inversion
  // routine may exchange rows by just swapping the pointers to the rows:
  double tol = 1.e-14;
  testResult &= abs(pA[0][0]-Ai[0][0]) < tol;
  testResult &= abs(pA[0][1]-Ai[0][1]) < tol;
  testResult &= abs(pA[0][2]-Ai[0][2]) < tol;
  testResult &= abs(pA[1][0]-Ai[1][0]) < tol;
  testResult &= abs(pA[1][1]-Ai[1][1]) < tol;
  testResult &= abs(pA[1][2]-Ai[1][2]) < tol;
  testResult &= abs(pA[2][0]-Ai[2][0]) < tol;
  testResult &= abs(pA[2][1]-Ai[2][1]) < tol;
  testResult &= abs(pA[2][2]-Ai[2][2]) < tol;

  return testResult;
}

bool testTridiagonalSystem()
{
  std::string testName = "TridiagonalSystem";
  bool testResult = true;

  double xt[5];
  double x[5];
  double y[5]    =  {-4, +3, +1, -5, -2};
  double A[5][5] = {{+2, -1,  0,  0,  0},
                    {+3, +5, +3,  0,  0},
                    { 0, +2, -3, -5,  0},
                    { 0,  0, +5, -2, +6},
                    { 0,  0,  0, -3, -7}};
  double *pA[5];
  for(int i = 0; i < 5; i++)
    pA[i] = &A[i][0];

  // solve via standard procedure to obtain the target solution:
  rsLinearAlgebra::rsSolveLinearSystem(pA, xt, y, 5);

  // solve via special tridiagonal procedure:
  double lower[4] = {+3, +2, +5, -3};     // lower diagonal of A
  double main[5]  = {+2, +5, -3, -2, -7}; // main  diagonal of A
  double upper[4] = {-1, +3, -5, +6};     // upper diagonal of A
  rsLinearAlgebra::rsSolveTridiagonalSystem(lower, main, upper, y, x, 5);

  // compare results:
  double tol = 1.e-14;
  for(int i = 0; i < 5; i++)
    testResult &= fabs(x[i]-xt[i]) < tol;

  return testResult;
}


bool testSquareMatrixTranspose()
{
  std::string testName = "SquareMatrixTranspose";
  bool testResult = true;

  // create test-matrix:
  double A[5][5] = {{ 1,  2,  3,  4,  5},
                    { 6,  7,  8,  9, 10},
                    {11, 12, 13, 14, 15},
                    {16, 17, 18, 19, 20},
                    {21, 22, 23, 24, 25}};

  // create transposed test-matrix:
  double T[5][5] = {{ 1,  6, 11, 16, 21},
                    { 2,  7, 12, 17, 22},
                    { 3,  8, 13, 18, 23},
                    { 4,  9, 14, 19, 24},
                    { 5, 10, 15, 20, 25}};

  // test (in-place) transposition of matrix A:
  double *pA[5];
  for(int i = 0; i < 5; i++)
    pA[i] = &A[i][0];


  rsArrayTools::transposeSquareArray(pA, 5);

  for(int i = 0; i < 5; i++)
  {
    for(int j = 0; j < 5; j++)
      testResult &= A[i][j] == T[i][j];
  }

  return testResult;
}



// todo: Gram-Schmidt orthogonalization of a set of basis vectors


bool testMatrixVectorMultiply()
{
  std::string testName = "MatrixVectorMultiply";
  bool testResult = true;

  // create test-matrix and vector, allocate result vector:
  double A[3][3] = {{1, 2, 3},
                    {4, 5, 6},
                    {7, 8, 9}};
  double *pA[3];
  for(int i = 0; i < 3; i++)
    pA[i] = &A[i][0];
  double x[3] =  {1, 2, 3};
  double y[3];

  // 3x3-matrix times 3-vector:
  // |1 2 3| * |1| = |14|
  // |4 5 6|   |2|   |32|
  // |7 8 9|   |3|   |50|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::matrixVectorMultiply(pA, x, y, 3, 3);
  testResult &= y[0] == 14;
  testResult &= y[1] == 32;
  testResult &= y[2] == 50;

  // transposed 3x3-matrix times 3-vector:
  // |1 2 3|^T * |1| = |1 4 7| * |1| = |30|
  // |4 5 6|     |2|   |2 5 8|   |2|   |36|
  // |7 8 9|     |3|   |3 6 9|   |3|   |42|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 3, 3);
  testResult &= y[0] == 30;
  testResult &= y[1] == 36;
  testResult &= y[2] == 42;

  // 3x2-matrix times 2-vector:
  // |1 2| * |1| = | 5|
  // |4 5|   |2|   |14|
  // |7 8|         |23|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::matrixVectorMultiply(pA, x, y, 3, 2);
  testResult &= y[0] == 5;
  testResult &= y[1] == 14;
  testResult &= y[2] == 23;

  // transposed 2x3-matrix times 2-vector:
  // |1 2 3|^T * |1| = |1 4| * |1| = | 9|
  // |4 5 6|     |2|   |2 5|   |2|   |12|
  //                   |3 6|         |15|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 2, 3);
  testResult &= y[0] == 9;
  testResult &= y[1] == 12;
  testResult &= y[2] == 15;

  // 2x3-matrix times 3-vector:
  // |1 2 3| * |1| = |14|
  // |4 5 6|   |2|   |32|
  //           |3|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::matrixVectorMultiply(pA, x, y, 2, 3);
  testResult &= y[0] == 14;
  testResult &= y[1] == 32;
  testResult &= y[2] == -1;

  // transposed 3x2-matrix times 3-vector:
  // |1 2|^T * |1| = |1 4 7| * |1| = |30|
  // |4 5|     |2|   |2 5 8|   |2|   |36|
  // |7 8|     |3|             |3|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 3, 2);
  testResult &= y[0] == 30;
  testResult &= y[1] == 36;
  testResult &= y[2] == -1;

  return testResult;
}

bool testMatrixMultiply3x3()
{
  // This function tests only some special cases where P either equals N or M using matrices where
  // no index-range exceeds 3. It is called from the more general test function
  // testMatrixMultiply3x3 internally.

  bool testResult = true;

  double A[3][3] = {{1, 2, 3},
                    {4, 5, 6},
                    {7, 8, 9}};
  double B[3][3] = {{9, 8, 7},
                    {6, 5, 4},
                    {3, 2, 1}};
  double C[3][3];
  double *pA[3], *pB[3], *pC[3];
  for(int i = 0; i < 3; i++)
  {
    pA[i] = &A[i][0];
    pB[i] = &B[i][0];
    pC[i] = &C[i][0];
  }

  // 2x3-matrix times 3x2 matrix (unused fields of result matrix should retain init-value -1):
  // |1 2 3| * |9 8| = |30 24|
  // |4 5 6|   |6 5|   |84 69|
  //           |3 2|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiply(pA, pB, pC, 2, 3, 2);
  testResult &= C[0][0] == 30 && C[0][1] == 24 && C[0][2] == -1;
  testResult &= C[1][0] == 84 && C[1][1] == 69 && C[1][2] == -1;
  testResult &= C[2][0] == -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x2-matrix times 2x3 matrix:
  // |1 2| * |9 8 7| = | 21 18 15|
  // |4 5|   |6 5 4|   | 66 57 48|
  // |7 8|             |111 96 81|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiply(pA, pB, pC, 3, 2, 3);
  testResult &= C[0][0] ==  21 && C[0][1] == 18 && C[0][2] == 15;
  testResult &= C[1][0] ==  66 && C[1][1] == 57 && C[1][2] == 48;
  testResult &= C[2][0] == 111 && C[2][1] == 96 && C[2][2] == 81;

  // transposed 3x2-matrix times 3x2 matrix:
  // |1 2|^T * |9 8| = |1 4 7| * |9 8| = |54 42|
  // |4 5|     |6 5|   |2 5 8|   |6 5|   |72 57|
  // |7 8|     |3 2|             |3 2|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 3, 2, 2);
  testResult &= C[0][0] == 54 && C[0][1] == 42 && C[0][2] == -1;
  testResult &= C[1][0] == 72 && C[1][1] == 57 && C[1][2] == -1;
  testResult &= C[2][0] == -1 && C[2][1] == -1 && C[2][2] == -1;

  // transposed 2x3-matrix times 2x3 matrix:
  // |1 2 3|^T * |9 8 7| = |1 4| * |9 8 7| = |33 28 23|
  // |4 5 6|     |6 5 4|   |2 5|   |6 5 4|   |48 41 34|
  //                       |3 6|             |63 54 45|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 2, 3, 3);
  testResult &= C[0][0] == 33 && C[0][1] == 28 && C[0][2] == 23;
  testResult &= C[1][0] == 48 && C[1][1] == 41 && C[1][2] == 34;
  testResult &= C[2][0] == 63 && C[2][1] == 54 && C[2][2] == 45;

  // 2x3-matrix times transposed 2x3 matrix:
  // |1 2 3| * |9 8 7|^T = |1 2 3| * |9 6| = | 46 28|
  // |4 5 6|   |6 5 4|     |4 5 6|   |8 5|   |118 73|
  //                                 |7 4|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 2, 3, 2);
  testResult &= C[0][0] ==  46 && C[0][1] == 28 && C[0][2] == -1;
  testResult &= C[1][0] == 118 && C[1][1] == 73 && C[1][2] == -1;
  testResult &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x2-matrix times transposed 3x2 matrix:
  // |1 2| * |9 8|^T = |1 2| * |9 6 3| = | 25 16  7|
  // |4 5|   |6 5|     |4 5|   |8 5 2|   | 76 49 22|
  // |7 8|   |3 2|     |7 8|             |127 82 37|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 3, 2, 3);
  testResult &= C[0][0] ==  25 && C[0][1] == 16 && C[0][2] ==  7;
  testResult &= C[1][0] ==  76 && C[1][1] == 49 && C[1][2] == 22;
  testResult &= C[2][0] == 127 && C[2][1] == 82 && C[2][2] == 37;

  // transposed 2x3-matrix times transposed 3x2 matrix:
  // |1 2 3|^T * |9 8|^T = |1 4| * |9 6 3| = |41 26 11|
  // |4 5 6|     |6 5|     |2 5|   |8 5 2|   |58 37 16|
  //             |3 2|     |3 6|             |75 48 21|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 2, 3, 3);
  testResult &= C[0][0] == 41 && C[0][1] == 26 && C[0][2] == 11;
  testResult &= C[1][0] == 58 && C[1][1] == 37 && C[1][2] == 16;
  testResult &= C[2][0] == 75 && C[2][1] == 48 && C[2][2] == 21;

  // transposed 3x2-matrix times transposed 2x3 matrix:
  // |1 2|^T * |9 8 7|^T = |1 4 7| * |9 6| = | 90 54|
  // |4 5|     |6 5 4|     |2 5 8|   |8 5|   |114 69|
  // |7 8|                           |7 4|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 3, 2, 2);
  testResult &= C[0][0] ==  90 && C[0][1] == 54 && C[0][2] == -1;
  testResult &= C[1][0] == 114 && C[1][1] == 69 && C[1][2] == -1;
  testResult &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x3-matrix times 3x3-matrix in place multiplication:
  // |1 2 3| * |9 8 7| = | 30  24 18|
  // |4 5 6|   |6 5 4|   | 84  69 54|
  // |7 8 9|   |3 2 1|   |138 114 90|
  rsMatrixTools::copyMatrix(pA, pC, 3, 3);
  rsMatrixTools::matrixInPlaceMultiply(pC, pB, 3, 3);
  testResult &= C[0][0] ==  30 && C[0][1] ==  24 && C[0][2] == 18;
  testResult &= C[1][0] ==  84 && C[1][1] ==  69 && C[1][2] == 54;
  testResult &= C[2][0] == 138 && C[2][1] == 114 && C[2][2] == 90;

  // 2x3-matrix times 3x3-matrix in place multiplication:
  // |1 2 3| * |9 8 7| = | 30  24 18|
  // |4 5 6|   |6 5 4|   | 84  69 54|
  //           |3 2 1|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::copyMatrix(pA, pC, 2, 3);
  rsMatrixTools::matrixInPlaceMultiply(pC, pB, 2, 3);
  testResult &= C[0][0] ==  30 && C[0][1] ==  24 && C[0][2] == 18;
  testResult &= C[1][0] ==  84 && C[1][1] ==  69 && C[1][2] == 54;
  testResult &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  return testResult;
}

bool testMatrixMultiply()
{
  std::string testName = "MatrixMultiply";
  bool testResult = testMatrixMultiply3x3();

  double A[4][4] = {{  2,   3,   5,   7},
                    { 11,  13,  17,  19},
                    { 23,  29,  31,  37},
                    { 41,  43,  47,  53}};
  double B[4][4] = {{ 59,  61,  67,  71},
                    { 73,  79,  83,  89},
                    { 97, 101, 103, 107},
                    {109, 113, 127, 131}};
  double C[4][4];
  double *pA[4], *pB[4], *pC[4];
  for(int i = 0; i < 4; i++)
  {
    pA[i] = &A[i][0];
    pB[i] = &B[i][0];
    pC[i] = &C[i][0];
  }

  // 2x3 matrix times 3x4 matrix:
  // | 2  3  5| * |59  61  67  71| = | 822  864  898  944|
  // |11 13 17|   |73  79  83  89|   |3247 3415 3567 3757|
  //              |97 101 103 107|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiply(pA, pB, pC, 2, 3, 4);
  testResult &= C[0][0] ==  822 && C[0][1] ==  864 && C[0][2] ==  898 && C[0][3] ==  944;
  testResult &= C[1][0] == 3247 && C[1][1] == 3415 && C[1][2] == 3567 && C[1][3] == 3757;
  testResult &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  testResult &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // 2x3 matrix times transposed 3x4 matrix:
  // | 2  3  5| * | 59  61  67|^T = | 636  798 1012 1192|
  // |11 13 17|   | 73  79  83|     |2581 3241 4131 4827|
  //              | 97 101 103|
  //              |109 113 127|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 2, 3, 4);
  testResult &= C[0][0] ==  636 && C[0][1] ==  798 && C[0][2] == 1012 && C[0][3] == 1192;
  testResult &= C[1][0] == 2581 && C[1][1] == 3241 && C[1][2] == 4131 && C[1][3] == 4827;
  testResult &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  testResult &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // transposed 3x2 matrix times 3x4 matrix:
  // | 2  3|^T * |59  61  67  71| = |3152 3314 3416 3582|
  // |11 13|     |73  79  83  89|   |3939 4139 4267 4473|
  // |23 29|     |97 101 103 107|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 3, 2, 4);
  testResult &= C[0][0] == 3152 && C[0][1] == 3314 && C[0][2] == 3416 && C[0][3] == 3582;
  testResult &= C[1][0] == 3939 && C[1][1] == 4139 && C[1][2] == 4267 && C[1][3] == 4473;
  testResult &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  testResult &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // transposed 3x2 matrix times transposed 4x3 matrix:
  // | 2  3|^T * | 59  61  67|^T = |2330 2924 3647 4382|
  // |11 13|     | 73  79  83|     |2913 3653 4691 5479|
  // |23 29|     | 97 101 103|
  //             |109 113 127|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 3, 2, 4);
  testResult &= C[0][0] == 2330 && C[0][1] == 2924 && C[0][2] == 3674 && C[0][3] == 4382;
  testResult &= C[1][0] == 2913 && C[1][1] == 3653 && C[1][2] == 4591 && C[1][3] == 5479;
  testResult &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  testResult &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  return testResult;
}

bool testChangeOfBasis()
{
  std::string testName = "ChangeOfBasis";
  bool testResult = true;

  // basis A and B (the columns represent the basis vectors)
  double A[3][3] = {{-2,  3,  4},
                    {-3,  2, -5},
                    { 2, -3,  2}};

  double B[3][3] = {{ 2,  1, -2},
                    { 1, -3,  2},
                    {-4, -2,  3}};
  double AT[3][3], BT[3][3]; // transposed bases (rows are basis vectors)
  double C[3][3],  CT[3][3]; // change-of-base matrices (column-, and row-case)

  // row-pointers:
  double *pA[3], *pB[3], *pC[3], *pAT[3], *pBT[3], *pCT[3];
  for(int i = 0; i < 3; i++)
  {
    pA[i]  = &A[i][0];
    pB[i]  = &B[i][0];
    pC[i]  = &C[i][0];
    pAT[i] = &AT[i][0];
    pBT[i] = &BT[i][0];
    pCT[i] = &CT[i][0];
  }

  // obtain the row-wise bases:
  rsArrayTools::transposeSquareArray(pA, pAT, 3);
  rsArrayTools::transposeSquareArray(pB, pBT, 3);

  // coordinates of vector v in basis A:
  double va[3] = {1, -2, 3};

  // compute coordinates of v in basis B:
  double vb[3];
  rsLinearAlgebra::rsChangeOfBasisColumnWise(pA, pB, va, vb, 3);

  // compute canonical coordinates of v from its coordinates in va and vb and compare:
  double vea[3], veb[3];
  rsMatrixTools::matrixVectorMultiply(pA, va, vea, 3, 3);
  rsMatrixTools::matrixVectorMultiply(pB, vb, veb, 3, 3);
  testResult &= rsArrayTools::equal(vea, veb, 3);

  // compute change-of-base matrix to go from basis A to B:
  rsLinearAlgebra::rsChangeOfBasisMatrixColumnWise(pA, pB, pC, 3);

  // compute coordinates of v in B again, this time using the change-of-base matrix, compare to
  // previously computed coordinates:
  double vb2[3];
  rsMatrixTools::matrixVectorMultiply(pC, va, vb2, 3, 3);
  testResult &= rsArrayTools::almostEqual(vb, vb2, 3, 1.e-14);

  // tests for row-based representations of bases A and B:
  rsLinearAlgebra::rsChangeOfBasisRowWise(pAT, pBT, va, vb, 3);
  testResult &= rsArrayTools::almostEqual(vb, vb2, 3, 1.e-14);

  rsLinearAlgebra::rsChangeOfBasisMatrixRowWise(pAT, pBT, pCT, 3);
  rsMatrixTools::matrixVectorMultiply(pCT, va, vb2, 3, 3);
  testResult &= rsArrayTools::almostEqual(vb, vb2, 3, 1.e-14);

  testResult &= rsMatrixTools::areMatricesApproximatelyEqual(pC, pCT, 3, 3, 1.e-14);
   // the change-of-base matrix doesn't care about, how the base-vectors are represented in some
   // matrix-representation, so it should be the same in both cases

  return testResult;
}

// i think, i found a rule: iff C = A^-1 * B, then C^-1 = B^-1 * A - is this a known theorem?
// it can easily be found from A*x = B*y (where x,y and either A or B are assumed to be known) and
// the solving for x and y respectively and defining C as conversion matrix from x to y and C-^1
// the inverse conversion, can also be verified by E = C * C^-1 = A^-1 * B * B^-1 * A

template<class T>
void solveDiagonalSystem(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B)
{
  int M = X.getNumColumns();  // number of required solution vectors
  int N = A.getNumRows();     // number of elements in each solution vector
  for(int k = 0; k < M; k++)
    for(int i = N-1; i >= 0; i--)
      X(i, k) = B(i, k) / A(i, i);
}
// adapted from solveUpperDiagonalSystem

bool testNullSpace()
{
  bool r = true;

  using Matrix = rsMatrix<double>;
  using LA     = rsLinearAlgebraNew;
  double tol   = 1.e-14;

  // maybe it makes sense to work with rational numbers rather than floating point

  Matrix I, A, z;
  int rank = 0;

  // Karpf. pg 142:
  z = Matrix(3, 1, { 0,0,0 });
  A = Matrix(3, 3, {1,2,3, 4,5,6, 7,8,9});
  //A = Matrix(3, 3, {7,8,9, 1,2,3, 4,5,6});    // reorder rows to avoid swapping
  rank = LA::makeTriangular(A, z);
  //LA::makeTriangularNoPivot(A, B);
  r &= rank == 2;

  // Karpf. pg 141:
  I = Matrix(3, 3, {1,0,0, 0,1,0, 0,0,1});
  A = Matrix(3, 3, {1,1,1, 1,2,4, 2,3,5});
  rank = LA::makeTriangular(A, z);
  r &= rank == 2;

  // see Karpf pg 77,80,129,142,411

  return r;
}

// tests the new implementation
bool testLinearSystemViaGauss2()
{
  bool r = true;

  using Vector = std::vector<double>;
  using Matrix = rsMatrix<double>;
  using LA     = rsLinearAlgebraNew;
  double tol   = 1.e-14;

  //rsMatrix<double> A(3, 3, { 1,2,3, 4,5,6, 7,8,9 }); // this matrix is singular
  Matrix A(3, 3, { 2,1,4, 3,10,3, 1,5,1 });
  Matrix tmpA = A;                               // because algo destroys the original A
  Vector x({1,2,3});
  Vector b  = A * x;                            //       A * x = b
  Vector x2 = LA::solve(tmpA, b);               // solve A * x = b for x
  r &= RAPT::rsAreVectorsEqual(x, x2, tol);
  tmpA = A, b = A*x;                             // restore destroyed tmp and b

  // try it with 3x2 solution matrix - figure out if the X,B rhs matrices/vectors may in general be
  // the same (it works when inverting - but this is the special case where B is the identity - in
  // general, it probably won't work)

  // check matrix inversion:
  Matrix At(3, 3, {-0.5,1.9,-3.7, 0.0,-0.2,0.6, 0.5,-0.9,1.7 }); // target matrix
  Matrix Ai = LA::inverse(A);
  r &= Ai.equals(At, tol);

  // solve some systems with 3x3 matrix and 3x2 solution vector:
  Matrix X(3, 2, {1,4, 2,5, 3,6});
  Matrix B = A * X;
  r &= B == Matrix(3, 2, {16,37, 32,80, 14,35});
  tmpA = A;
  Matrix tmpB = B;
  Matrix X2(3,2);
  LA::solve(tmpA, X2, tmpB);  // destroys tmpA and tempB
  r &= X2.equals(X, tol);

  // Sage:
  // A = matrix([[2,1,4],[3,10,3],[1,5,1]])
  // X = matrix([[1,4],[2,5],[3,6]])
  // A, X, A*X
  //
  // gives:
  // [ 2  1  4]  [1 4]  [16 37]
  // [ 3 10  3]  [2 5]  [32 80]
  // [ 1  5  1], [3 6], [14 35]

  // test behavior of makeTriangular in case of singular matrices


  // figure out, when it can be used in place, i.e. X == B

  return r;
}

bool testPowerIteration()
{
  bool ok = true;

  // ...

  return ok;
}


bool testLinearAlgebra()
{
  std::string reportString = "LinearAlgebra"; // dummy-string - delete later
  bool ok = true;

  // Basics:
  ok &= testSquareMatrixTranspose();
  ok &= testMatrixVectorMultiply();
  ok &= testMatrixMultiply();

  // Direct solvers:
  //ok &= testBandDiagonalSolver();  // fails with gcc
  ok &= testMatrix2x2();
  ok &= testLinearSystem2x2();
  ok &= testLinearSystem3x3();
  ok &= testLinearSystemViaGauss();
  ok &= testGaussJordanInversion();
  ok &= testTridiagonalSystem();
  ok &= testChangeOfBasis();
  ok &= testNullSpace();
  ok &= testLinearSystemViaGauss2();

  // Iterative solvers:


  return ok;
}
