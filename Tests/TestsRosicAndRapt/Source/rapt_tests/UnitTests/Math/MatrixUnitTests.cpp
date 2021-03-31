//=================================================================================================
// Tests for old matrix implementation:

bool testMatrixScalarOperations()
{
  bool testResult = true;

  // A = |15 20 26|
  //     |12 18 24|
  rsMatrixDbl A(2, 3);
  A.set(0,0, 15); A.set(0,1, 20); A.set(0,2, 26);
  A.set(1,0, 12); A.set(1,1, 18); A.set(1,2, 24);

  // A = A + 2 = |17 22 28|
  //             |14 20 26|
  A = A + 2.0;
  testResult &= A(0,0) == 17 && A(0,1) == 22 && A(0,2) == 28;
  testResult &= A(1,0) == 14 && A(1,1) == 20 && A(1,2) == 26;

  // A = A - 2 = |15 20 26|
  //             |12 18 24|
  A = A - 2.0;
  testResult &= A(0,0) == 15 && A(0,1) == 20 && A(0,2) == 26;
  testResult &= A(1,0) == 12 && A(1,1) == 18 && A(1,2) == 24;

  // A = 2 + A = |17 22 28|
  //             |14 20 26|
  A = 2.0 + A;
  testResult &= A(0,0) == 17 && A(0,1) == 22 && A(0,2) == 28;
  testResult &= A(1,0) == 14 && A(1,1) == 20 && A(1,2) == 26;

  // A = 30 - A = |13  8  2|
  //              |16 10  4|
  A = 30.0 - A;
  testResult &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  testResult &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A = 2.0 * A = |26 16  4|
  //               |32 20  8|
  A = 2.0 * A;
  testResult &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  testResult &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A = A / 2.0 = |13  8  2|
  //               |16 10  4|
  A = A / 2.0;
  testResult &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  testResult &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A = A * 2.0 = |26 16  4|
  //               |32 20  8|
  A = A * 2.0;
  testResult &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  testResult &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A /= 2.0 = |13  8  2|
  //            |16 10  4|
  A /= 2.0;
  testResult &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  testResult &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A *= 2.0 = |26 16  4|
  //            |32 20  8|
  A *= 2.0;
  testResult &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  testResult &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A += 2.0 = |28 18  6|
  //            |34 22 10|
  A += 2.0;
  testResult &= A(0,0) == 28 && A(0,1) == 18 && A(0,2) ==  6;
  testResult &= A(1,0) == 34 && A(1,1) == 22 && A(1,2) == 10;

  // A -= 2.0 = |26 16  4|
  //            |32 20  8|
  A -= 2.0;
  testResult &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  testResult &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  return testResult;
}
bool testMatrixAddSubTransEqual()
{
  bool testResult = true;

  // intitalize some test matrices:
  // A = |2  3  5|, B = |17 19 23|
  //     |7 11 13|      |29 31 37|
  rsMatrixDbl A(2, 3);
  A.set(0,0, 2); A.set(0,1,  3); A.set(0,2,  5);
  A.set(1,0, 7); A.set(1,1, 11); A.set(1,2, 13);

  rsMatrixDbl B(2, 3);
  B.set(0,0, 17); B.set(0,1, 19); B.set(0,2, 23); 
  B.set(1,0, 29); B.set(1,1, 31); B.set(1,2, 37);

  rsMatrixDbl C;

  // check, if C is properly initialized:
  testResult &= C.getNumRows()    == 0;
  testResult &= C.getNumColumns() == 0;
  //testResult &= C.mFlat == nullptr;
  //testResult &= C.m     == nullptr;




  // test unary minus:
  // C = -A = |-2  -3  -5|
  //          |-7 -11 -13|
  C = -A;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) == -2 && C(0,1) ==  -3 && C(0,2) ==  -5;
  testResult &= C(1,0) == -7 && C(1,1) == -11 && C(1,2) == -13;

  // test transposition:
  // C = A^T = |2  7|
  //           |3 11|
  //           |5 13|
  C = trans(A);
  testResult &= C.getNumRows()    == 3;
  testResult &= C.getNumColumns() == 2;
  testResult &= C(0,0) ==  2 && C(0,1) ==  7;
  testResult &= C(1,0) ==  3 && C(1,1) == 11;
  testResult &= C(2,0) ==  5 && C(2,1) == 13;

  // test in-place transpostion and equality operator:
  // C = C^T == A
  C.transpose();
  testResult &= C == A;

  // test in/equality operators:
  C.transpose();
  testResult &= (A ==  B) == false;
  testResult &= (A !=  B) == true;
  testResult &= (A ==  C) == false;
  testResult &= (A !=  C) == true;

  // C = A + B = |19 22 28|
  //             |36 42 50|
  C = A + B;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) == 19 && C(0,1) == 22 && C(0,2) == 28;
  testResult &= C(1,0) == 36 && C(1,1) == 42 && C(1,2) == 50;

  // C += A = |21 25 33| 
  //          |43 53 63|
  C += A;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) == 21 && C(0,1) == 25 && C(0,2) == 33;
  testResult &= C(1,0) == 43 && C(1,1) == 53 && C(1,2) == 63;

  // C -= A = |19 22 28|
  //          |36 42 50|
  C -= A;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) == 19 && C(0,1) == 22 && C(0,2) == 28;
  testResult &= C(1,0) == 36 && C(1,1) == 42 && C(1,2) == 50;

  // C = B - A = |15 16 18|
  //             |22 20 24|
  C = B - A;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) == 15 && C(0,1) == 16 && C(0,2) == 18;
  testResult &= C(1,0) == 22 && C(1,1) == 20 && C(1,2) == 24;

  return testResult;
}
bool testMatrixMul()
{
  bool testResult = true;

  // intitalize some test matrices:
  // A = |2  3  5|, B = |17 19 23 29|
  //     |7 11 13|      |31 37 41 43|
  //                    |47 53 59 61|
  rsMatrixDbl A(2, 3);
  A.set(0,0, 2); A.set(0,1,  3); A.set(0,2,  5);
  A.set(1,0, 7); A.set(1,1, 11); A.set(1,2, 13);

  rsMatrixDbl B(3, 4);
  B.set(0,0, 17); B.set(0,1, 19); B.set(0,2, 23); B.set(0,3, 29);
  B.set(1,0, 31); B.set(1,1, 37); B.set(1,2, 41); B.set(1,3, 43);
  B.set(2,0, 47); B.set(2,1, 53); B.set(2,2, 59); B.set(2,3, 61);

  rsMatrixDbl C;

  // C = A * B = | 362  414  464  492|
  //             |1071 1229 1379 1469|
  C = A * B;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 4;
  testResult &= C(0,0) ==  362 && C(0,1) ==  414 && C(0,2) ==  464 && C(0,3) ==  492;
  testResult &= C(1,0) == 1071 && C(1,1) == 1229 && C(1,2) == 1379 && C(1,3) == 1469;

  // same with the accumulative multiplier:
  C  = A;
  C *= B;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 4;
  testResult &= C(0,0) ==  362 && C(0,1) ==  414 && C(0,2) ==  464 && C(0,3) ==  492;
  testResult &= C(1,0) == 1071 && C(1,1) == 1229 && C(1,2) == 1379 && C(1,3) == 1469;

  // B2 = B * B^T = |2020 3420  4932|
  //                |3420 5860  8460|
  //                |4932 8460 12220|
  rsMatrixDbl B2 = B * trans(B);
  testResult &= B2.getNumRows()    == 3;
  testResult &= B2.getNumColumns() == 3;
  testResult &= B2(0,0) == 2020 && B2(0,1) == 3420 && B2(0,2) ==  4932;
  testResult &= B2(1,0) == 3420 && B2(1,1) == 5860 && B2(1,2) ==  8460;
  testResult &= B2(2,0) == 4932 && B2(2,1) == 8460 && B2(2,2) == 12220;

  // C = A * B2 = | 38960  66720  96344|
  //              |115876 198380 286444|
  C = A * B2;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) ==  38960 && C(0,1) ==  66720 && C(0,2) ==  96344;
  testResult &= C(1,0) == 115876 && C(1,1) == 198380 && C(1,2) == 286444;
  C  = A;
  C *= B2;
  testResult &= C.getNumRows()    == 2;
  testResult &= C.getNumColumns() == 3;
  testResult &= C(0,0) ==  38960 && C(0,1) ==  66720 && C(0,2) ==  96344;
  testResult &= C(1,0) == 115876 && C(1,1) == 198380 && C(1,2) == 286444;

  return testResult;
}
/*
template<class T>
T square(T x)
{
  return x*x;
}
*/
bool testMatrixArithmeticOld()
{
  bool ok = true;

  ok &= testMatrixScalarOperations();
  ok &= testMatrixAddSubTransEqual();
  ok &= testMatrixMul();

  // intitalize some test matrices:
  // A = |2  3  5|
  //     |7 11 13|
  // 
  rsMatrixDbl A(2, 3);
  A.set(0,0, 2); A.set(0,1,  3); A.set(0,2,  5);
  A.set(1,0, 7); A.set(1,1, 11); A.set(1,2, 13);

  A.applyFunction(&square<double>);
  ok &= A(0,0) ==  4 && A(0,1) ==   9 && A(0,2) ==  25;
  ok &= A(1,0) == 49 && A(1,1) == 121 && A(1,2) == 169;   

  return ok;
}

//=================================================================================================
// Tests for array based matrices:

bool testSquareMatrixTranspose()
{
  bool ok = true;

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

  for(int i = 0; i < 5; i++) {
    for(int j = 0; j < 5; j++)
      ok &= A[i][j] == T[i][j]; }

  return ok;
}

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









//=================================================================================================
// Tests for new matrix implementation:

// compares data in the given matrix view to that in the given vector and return true, if they are
// equal
template<class T>
bool hasData(const rsMatrixView<T>& m, const std::vector<T>& v)
{
  if(m.getSize() != (int) v.size())
    return false;
  return RAPT::rsArrayTools::equal(m.getDataPointerConst(), &v[0], m.getSize());
}

bool testMatrixView()
{
  std::string testName = "MatrixView";
  bool r = true;  // test result

  typedef rsMatrixView<double> MatrixView;
  typedef std::vector<double> Vec;

  double A6[6] = { 1,2,3,4,5,6 }; // array of 6 elements

  // 1 2 3
  // 4 5 6
  MatrixView m(2, 3, A6);
  r &= m(0,0) == 1; r &= m(0,1) == 2; r &= m(0,2) == 3;
  r &= m(1,0) == 4; r &= m(1,1) == 5; r &= m(1,2) == 6;

  // 1 2 
  // 3 4 
  // 5 6
  m.setShape(3, 2);
  r &= m(0,0) == 1; r &= m(0,1) == 2; 
  r &= m(1,0) == 3; r &= m(1,1) == 4;
  r &= m(2,0) == 5; r &= m(2,1) == 6; 

  // test elementary row- and column operations:

  // Create 4x3 matrix:
  // 11 12 13
  // 21 22 23
  // 31 32 33
  // 41 42 43
  double A12[12] = { 11,12,13, 21,22,23, 31,32,33, 41,42,43 };
  MatrixView m43(4, 3, A12);
  r &= hasData(m43, Vec({ 11,12,13, 21,22,23, 31,32,33, 41,42,43 }));

  m43.swapRows(1, 2);
  r &= hasData(m43, Vec({ 11,12,13, 31,32,33, 21,22,23, 41,42,43 }));
  m43.addWeightedRowToOther(2, 1, 2.0);
  r &= hasData(m43, Vec({ 11,12,13, 73,76,79, 21,22,23, 41,42,43 }));
  m43.scaleRow(2, 2.0);
  r &= hasData(m43, Vec({ 11,12,13, 73,76,79, 42,44,46, 41,42,43 }));

  m43.swapColumns(1,2);
  r &= hasData(m43, Vec({ 11,13,12, 73,79,76, 42,46,44, 41,43,42 }));
  m43.addWeightedColumnToOther(0, 1, -1.0);
  r &= hasData(m43, Vec({ 11,2,12, 73,6,76, 42,4,44, 41,2,42 }));
  m43.scaleColumn(1, 0.5);
  r &= hasData(m43, Vec({ 11,1,12, 73,3,76, 42,2,44, 41,1,42 }));






  // Test row- and column major storage and submatrix addressing - we use this 7x9 matrix:
  //
  //    0  1  2  3  4  5  6  7  8  
  //  0 11 12 13 14 15 16 17 18 19
  //  1 21 22 23 24 25 26 27 28 29
  //  2 31 32 33 34 35 36 37 38 39
  //  3 41 42 43 44 45 46 47 48 49
  //  4 51 52 53 54 55 56 57 58 59
  //  5 61 62 63 64 65 66 67 68 69
  //  6 71 72 73 74 75 76 77 78 79
  //
  // which is stored in memory like this:
  //
  // 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29... 
  // 11 12 13 14 15 16 17 18 19 21 22 23 24 25 26 27 28 29 31 32 33 34 35 36 37 38 39 41 42 43...
  // 11 21 31 41 51 61 71 12 22 32 42 52 62 72 13 23 33 43 53 63 73 14 24 34 44 54 64 74 15 25...
  //
  // 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62
  // 44 45 46 47 48 49 51 52 53 54 55 56 57 58 59 61 62 63 64 65 66 67 68 69 71 72 73 74 75 76 77 78 79 
  // 35 45 55 65 75 16 26 36 46 56 66 76 17 27 37 47 57 67 77 18 28 38 48 58 68 78 19 29 39 49 59 69 79
  //
  // the top-row is the flat array index, the 2nd the stored numbers in row-major, the 3rd row the 
  // numbers in column-major format. 

  // This code is inactive because it uses another version of rsMatrix - but i have reverted to the
  // version that doesn't support switching between row-major and column major storage - see below

  /*
  double A63[63];                    // storage space for a 7x9 matrix
  MatrixView m79r(7, 9, A63);
  r &=  m79r.isStorageRowMajor();    // by default, we should get row-major storage
  r &= !m79r.isStorageColumnMajor();
  r &= m79r.isStorageContiguous();
  int i, j;
  for(i = 0; i < 7; i++)
    for(j = 0; j < 9; j++)
      m79r(i, j) = 10*(i+1) + j+1;   

  // verify some random samples from the flat array:
  r &= A63[16] == 28 && A63[28] == 42 && A63[51] == 67;

  // test submatrixView - for row major-case


  MatrixView m79c(7, 9, A63, false); // false should indicate "no row-major storage"
  r &= !m79c.isStorageRowMajor();
  r &=  m79c.isStorageColumnMajor();
  r &= m79r.isStorageContiguous();

  // fill the matrix again - this time thd numbers got to different places
  for(i = 0; i < 7; i++)
    for(j = 0; j < 9; j++)
      m79r(i, j) = 10*(i+1) + j+1;
  //r &= A63[16] == 33 && A63[28] == 15 && A63[51] == 38;
  */

  /* On addressing submatrices:
  
  Let's and try to address the 3x5 submatrix that starts at index-pair 2,3, i.e. at the number 34 
  and ends at index-pair 4,7, i.e at number 58 in a row-major storage format:

  row-major means: 
  rowStride = numCols = 9
  colStride = 1
  start = 2*rowStride + 3*colStride = 2*9 + 3*1 = 21

  and at flat index 21, we find the 34 - as it should be - this is A(0,0) in our 
  submatrix. The number 46 is found at flat-index 32 and this is element A(1,2) in our 
  submatrix.

  let's take the formula:
  index = start + 1*rowStride + 2*colStride = 21 + 1*9 + 2*1 = 32
  ..yup - the formula works!

  now let's do it in column-major addressing:
  rowStride = 1
  colStride = numRows = 7
  start = 2*rowStride + 3*colStride = 2*1 + 3*7 = 23

  and sure enough, at flat index 23, we find the number 34 in our column-major format. Addressing the
  46 as A(1,2) again gives the falt index:

  index = start + 1*rowStride + 2*colStride = 23 + 1*1 + 2*7 = 38

  ...aaaand yes! at flat index 38, we find our 46. That means, if we address matrix entries by the 
  formula: start + rowIndex*rowStride + colIndex*colStride, we area able to provide row -and 
  column-major storage and convenient access to submatrices without having to copy any data! :-)
  So, using the address computation i*rowStride + j*colStride allows for row- and column-major
  storage as well as addressing submatrices (we would have to use a different dataPointer - one
  that points to the start of the submatrix. row- and column-stride will then be different from
  numCols or numRows. But the disadvangtae is that for submatrices, the storage is not contiguous
  anymore, so we cant use functions from rsArrayTools for things like addition, finding the 
  maximum, etc. - so i don't know, if it's worth it - that's why i have reverted the code to before
  these things above were implemented. The modified version is in a commit from 
  17th Jan 2020 - just in case, i decide to continue this route  */

  return r;
}

bool testMatrixOperators()
{
  bool r = true;  // test result
  using Mat = rsMatrix<double>;
  using Vec = std::vector<double>;

  Mat A(3, 3, { 2,1,4, 3,10,3, 1,5,1 });
  Mat X(3, 2, {1,4, 2,5, 3,6});
  Mat B = A * X;
  r &= B == Mat(3, 2, {16,37, 32,80, 14,35});

  // test matrix-vector multiplication:
  Vec x, y; x = Vec({1,2,3});
  A = Mat(2, 3, {1,2,3, 4,5,6});  y = A*x; r &= y == Vec({14,32});
  A = Mat(3, 2, {1,2, 3,4, 5,6}); y = x*A; r &= y == Vec({22,28});
  
  return r;
}

bool testMatrixAlloc() // rename to testMatrixAllocationAndArithmetic
{
  bool testResult = true;
  using Matrix = rsMatrix<double>;
  using Vector = std::vector<double>;
  int& allocs  = Matrix::numHeapAllocations;  // to count allocations
  allocs = 0;

  // A = |1 2 3|
  //     |4 5 6|
  Matrix A(2, 3, {1.,2.,3., 4.,5.,6.}); // calls rsMatrix(int, int, std::vector<T>&&)
  testResult &= A(0,0) == 1 &&  A(0,1) == 2 && A(0,2) == 3;
  testResult &= A(1,0) == 4 &&  A(1,1) == 5 && A(1,2) == 6;
  testResult &= allocs == 1;

  // B = |1 2| 
  //     |3 4|
  //     |5 6|
  Matrix B(3, 2, {1.,2.,3., 4.,5.,6.}); // calls rsMatrix(int, int, std::vector<T>&&)
  testResult &= B(0,0) == 1 &&  B(0,1) == 2;
  testResult &= B(1,0) == 3 &&  B(1,1) == 4;
  testResult &= B(2,0) == 5 &&  B(2,1) == 6;
  testResult &= allocs == 2;

  // matrix-vector multiplication:
  Vector x({1,2,3});
  Vector y = A * x;
  testResult &= y == Vector({14,32});

  // matrix multiplication:
  Matrix C = A*B; 
  // calls: 
  //   operator*(const rsMatrix<T>&)
  //     rsMatrix(int numRows, int numColumns)
  //     rsMatrix(rsMatrix&& B)
  testResult &= allocs == 3;
  testResult &= C(0,0) == 22 &&  C(0,1) == 28;
  testResult &= C(1,0) == 49 &&  C(1,1) == 64;

  C = A;  // calls copy assigment operator
  testResult &= allocs == 4;
  testResult &= C == A;

  A = A;  // self copy assignment - should not re-allocate
  testResult &= allocs == 4;

  Matrix D = A+A;
  testResult &= allocs == 5;
  testResult &= D == Matrix(2, 3, {2.,4.,6., 8.,10.,12.});
  testResult &= allocs == 6;
  testResult &= D == A+A;
  testResult &= allocs == 7;

  C = B*A;  // calls move assignment operator
  testResult &= allocs == 8;   // fails here
  testResult &= C(0,0) ==  9 &&  C(0,1) == 12 && C(0,2) == 15;
  testResult &= C(1,0) == 19 &&  C(1,1) == 26 && C(1,2) == 33;
  testResult &= C(2,0) == 29 &&  C(2,1) == 40 && C(2,2) == 51;


  Matrix E(A+D);  // calls move constructor
  testResult &= allocs == 9;
  testResult &= E == A+D;  // temporary A+D needs allocation
  testResult &= allocs == 10; 

  Matrix F(E);   // calls copy constructor
  testResult &= allocs == 11;
  testResult &= F == E;   // no allocation here
  testResult &= allocs == 11;

  // hmm - these here call move/copy constructors and not move/copy assigment operators - why
  // ..ah - it's because it's not a re-assignment
  Matrix G = A+D; // calls move constructor
  testResult &= allocs == 12;
  testResult &= G == A+D;  // temporary A+D needs allocation
  testResult &= allocs == 13;

  Matrix H = G;   // calls copy constructor
  testResult &= allocs == 14;
  testResult &= H == G;
  testResult &= allocs == 14;


  Matrix I, J;  // no allocations yet
  testResult &= allocs == 14;
  I = A+D;      // move assignment
  testResult &= allocs == 15;
  J = I;        // copy assignment
  testResult &= allocs == 16;


  // multiplication with a scalar:
  I = 2.0 * A;
  testResult &= allocs == 17;
  testResult &=  I == Matrix(2, 3, {2.,4.,6., 8.,10.,12.});
  testResult &= allocs == 18;

  J = A * 2.0;
  testResult &= allocs == 19;
  testResult &= I == J;
  testResult &= allocs == 19;

  // in-place addition and subtraction via +=, -=:
  J += A;
  testResult &= allocs == 19;
  testResult &= J == 3.0 * A;
  testResult &= allocs == 20;
  J -= A;
  testResult &= allocs == 20;
  testResult &= J == 2.0 * A;
  testResult &= allocs == 21;

  // the *= operator can't operate in-place:
  J = A;    // should not re-allocate because J already has the right shape
  testResult &= allocs == 21;
  J *= B;   // allocates
  testResult &= allocs == 22;
  testResult &= J == A*B;
  testResult &= allocs == 23;

  // in-place scaling by a scalar:
  J = A;    
  testResult &= allocs == 24;
  J *= 2.0;
  testResult &= allocs == 24;
  testResult &= J == 2.0 * A;
  testResult &= allocs == 25;
  J /= 2.0;
  testResult &= allocs == 25;
  testResult &= J == A;
  testResult &= allocs == 25;


  C = Matrix(2, 4, {1,2,3,4, 5,6,7,8});     // 2x4
  testResult &= allocs == 26;
  D = Matrix(4, 2, {1,2, 3,4, 5,6, 7,8});   // 4x2
  testResult &= allocs == 27;

  E = B*C;                                  // 3x2 * 2x4 = 3x4
  testResult &= allocs == 28;
  testResult &= E.getNumRows()    == 3;
  testResult &= E.getNumColumns() == 4;
  testResult &= E(0,0) == 11 && E(0,1) == 14 && E(0,2) == 17 && E(0,3) == 20;
  testResult &= E(1,0) == 23 && E(1,1) == 30 && E(1,2) == 37 && E(1,3) == 44;
  testResult &= E(2,0) == 35 && E(2,1) == 46 && E(2,2) == 57 && E(2,3) == 68;

  E = D*A;
  testResult &= allocs == 29;               // 4x2 * 2x3 = 4x3
  testResult &= E.getNumRows()    == 4;
  testResult &= E.getNumColumns() == 3;
  testResult &= E(0,0) ==  9 && E(0,1) == 12 && E(0,2) == 15;
  testResult &= E(1,0) == 19 && E(1,1) == 26 && E(1,2) == 33;
  testResult &= E(2,0) == 29 && E(2,1) == 40 && E(2,2) == 51;
  testResult &= E(3,0) == 39 && E(3,1) == 54 && E(3,2) == 69;

  // move constructor:
  Matrix K(std::move(E));
  testResult &= allocs == 29;
  testResult &= K.getNumRows() == 4 && K.getNumColumns() == 3;
  testResult &= E.getNumRows() == 0 && E.getNumColumns() == 0;
  testResult &= K.getDataVectorConst().size() == 12;
  testResult &= E.getDataVectorConst().size() == 0;
  testResult &= K.getDataPointerConst() != nullptr;
  testResult &= E.getDataPointerConst() == nullptr;

  // move assignment operator:
  E = std::move(K);
  testResult &= allocs == 29;
  testResult &= E.getNumRows() == 4 && E.getNumColumns() == 3;
  testResult &= K.getNumRows() == 0 && K.getNumColumns() == 0;
  testResult &= E.getDataVectorConst().size() == 12;
  testResult &= K.getDataVectorConst().size() == 0;
  testResult &= E.getDataPointerConst() != nullptr;
  testResult &= K.getDataPointerConst() == nullptr;

  // transposition:
  E = A;
  testResult &= allocs == 30;
  testResult &= E.getNumRows() == 2 && E.getNumColumns() == 3;
  E.transpose();
  testResult &= allocs == 31;
  testResult &= E == Matrix(3, 2, {1,4, 2,5, 3,6});
  testResult &= allocs == 32;

  Matrix S4(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});
  testResult &= allocs == 33;
  S4.transpose();
  testResult &= allocs == 33;
  testResult &= S4 == Matrix(4, 4, {1,5,9,13, 2,6,10,14, 3,7,11,15, 4,8,12,16});
  testResult &= allocs == 34;

  Matrix V3(3, 1, {1,2,3});     // a 3D column-vector...
  testResult &= allocs == 35;
  V3.transpose();               // ...becomes a 3D row-vector
  testResult &= allocs == 35;
  testResult &= V3 == Matrix(1, 3, {1,2,3});
  testResult &= allocs == 36;

  // factory functions:
  E = Matrix::identity(4);
  testResult &= allocs == 37;
  testResult &= E == Matrix(4, 4, {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1});
  testResult &= allocs == 38;

  E = Matrix::zero(3, 4);
  testResult &= allocs == 39;
  testResult &= E == Matrix(3, 4, {0,0,0,0, 0,0,0,0, 0,0,0,0});
  testResult &= allocs == 40;

  // swapping two matrices:
  RAPT::rsSwapNaive(A, B);
  testResult &= allocs == 43;   // naive swapping causes 3 heap allocations...
  RAPT::rsSwapViaMove(A, B);
  testResult &= allocs == 43;   // ...whereas swapping via std::move causes none
  RAPT::rsSwap(A, B);
  testResult &= allocs == 43;   // rsSwap calls rsSwapViaMove






  // create diagonal matrix



  // todo:
  // -maybe add/subtract scalars (in-place and out-of-place)....but maybe not: mathematically, the
  //  space of MxN matrices is a vector space and addition of a scalar and a vector is not a thing
  //  in vector spaces




  return testResult;
}

bool testKroneckerProduct()
{
  bool res = true;
  typedef rsMatrix<double> Matrix;

  int& allocs  = Matrix::numHeapAllocations;  // to count allocations
  allocs = 0;

  Matrix A(2, 3, {1,2,3,4,5,6});
  res &= allocs == 1;
  Matrix B(4, 5, {1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20});
  res &= allocs == 2;

  Matrix C = Matrix::getKroneckerProduct(A, B);
  res &= allocs == 3;
  res &= C.getNumRows() == 8 && C.getNumColumns() == 15;

  Matrix T(8, 15,   // target Kronecker product
    { 1, 2, 3, 4, 5, 2, 4, 6, 8, 10, 3,  6,  9, 12, 15,
      6, 7, 8, 9,10,12,14,16,18, 20,18, 21, 24, 27, 30,
     11,12,13,14,15,22,24,26,28, 30,33, 36, 39, 42, 45,
     16,17,18,19,20,32,34,36,38, 40,48, 51, 54, 57, 60,
      4, 8,12,16,20, 5,10,15,20, 25, 6, 12, 18, 24, 30,
     24,28,32,36,40,30,35,40,45, 50,36, 42, 48, 54, 60,
     44,48,52,56,60,55,60,65,70, 75,66, 72, 78, 84, 90,
     64,68,72,76,80,80,85,90,95,100,96,102,108,114,120});
  res &= allocs == 4;
  res &= C == T;

  // the following sage code produces the result:
  // A = matrix([[1,2,3],[4,5,6]])
  // B = matrix([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]])
  // C = A.tensor_product(B)
  // A, B, C

  return res;
}

bool testSparseMatrix()
{
  bool res = true;

  using Vec = std::vector<float>;

  // We create a matrix of the form:
  // 
  //      |1 1 0 0 0 0 0 0|
  //  A = |0 0 1 1 0 0 0 0|
  //      |0 0 0 0 1 1 0 0|
  //      |0 0 0 0 0 0 1 1|
  //
  // so, it takes an 8-vector as input and produces a 4-vector as output in which each element is
  // the sum of two consecutive elements from the input.

  rsSparseMatrix<float> A(4, 8); // still the zero matrix
  Vec x({1,2,3,4,5,6,7,8});      // input vector
  Vec y({1,2,3,4});              // output vector
  A.product(&x[0], &y[0]);  res &= y == Vec({0,0,0,0}); 

  // Now build up the actual matrix in a kind of "random" manner. The nonzero (i.e. one) elements 
  // are at positions: (0,0),(0,1),(1,2),(1,3),(2,4),(2,5),(3,6),(3,7):
  A.set(1, 2, 1.f);
  A.set(2, 5, 1.f);
  A.set(0, 0, 1.f);
  A.set(0, 1, 1.f);
  A.set(3, 6, 1.f);
  A.set(1, 3, 1.f);
  A.set(2, 4, 1.f);
  A.set(3, 7, 1.f);

  // Test the multiplication with the new matrix again:
  A.product(&x[0], &y[0]);  res &= y == Vec({3,7,11,15});

  res &= A(0, 0) == 1.f;
  res &= A(0, 1) == 1.f;
  res &= A(0, 2) == 0.f;
  res &= A(0, 3) == 0.f;
  // ...

  // todo: test replacing elements, also with zero (in which case they should get removed)

  return res;
}

// move to LinearAlgebraUnitTests:
bool testSparseMatrixSolvers()
{
  bool res = true;

  using Vec = std::vector<float>;
  using Mat = rsSparseMatrix<float>;

  //      |7 1 2|
  //  A = |3 8 4|  ...must be strictly diagoally dominant for Gauss-Seidel to converge
  //      |5 2 9|

  int N = 3;
  Mat A(N, N);
  A.set(0, 0, 7.f);  A.set(0, 1, 1.f);  A.set(0, 2, 2.f);
  A.set(1, 0, 3.f);  A.set(1, 1, 8.f);  A.set(1, 2, 4.f);
  A.set(2, 0, 5.f);  A.set(2, 1, 2.f);  A.set(2, 2, 9.f);

  // Compute matrix-vector product y = A*x:
  Vec x({1,2,3}), y(N);
  A.product(&x[0], &y[0]);
  res &= y == Vec({15, 31, 36});
  A.iterateProduct(&x[0], &y[0]);
  res &= y == Vec({15, 31, 36});

  // Try to reconstruct x via solving A*x = y:
  //Mat D = A.getDiagonalPart();    // maybe this should be a vector...but maybe not
  //Mat C = A.getNonDiagonalPart();

  //float tol = 1.e-6f;
  float tol = 0.f;
  Vec x2(N);
  //int numIts = Mat::solveGaussSeidel(D, C, &x2[0], &y[0], tol);
  int numIts = Mat::solveGaussSeidel(A, &x2[0], &y[0], tol);
  res &= x == x2;
  res &= numIts == 13;  // converges to machine precision (with 0 tolerance) in 13 iterations


  rsSetZero(x2);
  Vec wrk(N);
  numIts = Mat::solveSOR(A, &x2[0], &y[0], tol, &wrk[0], 1.f); // 1: same as Gauss-Seidel
  res &= x == x2;
  res &= numIts == 13; 
  rsSetZero(x2);
  numIts = Mat::solveSOR(A, &x2[0], &y[0], tol, &wrk[0], 0.f); // 0: same as Jacobi
  res &= x == x2;
  res &= numIts == 40;   // Jacobi takes much longer!

  // try some other values of w:
  rsSetZero(x2);
  numIts = Mat::solveSOR(A, &x2[0], &y[0], tol, &wrk[0], 0.5f); 
  res &= x == x2;
  res &= numIts == 18;   // 0.5 is in between Jacobi and Gauss-Seidel as expected

  // w > 1 gets us into actual *over*relaxation territory - but it doesn't seem to improve 
  // convergence:
  rsSetZero(x2);
  numIts = Mat::solveSOR(A, &x2[0], &y[0], tol, &wrk[0], 1.1f); 
  res &= x == x2;

  // ..so Gauss-Seidel actually seems to converge fastest, which is good because it's also the most
  // efficient algo and uses no additional workspace memory! ..but why does SOR not outperform it? 
  // isn't it supposed to do? maybe there's a bug? ...more tests and experiments needed...

  // ToDo: in certain contexts, it may make sense to use a better initial guess (here, we use the 
  // zero vector)



  // Try to find eigenvalues and vectors:
  float ev;
  x = Vec({1,1,1});
  //A.largestEigenValueAndVector(&ev, &x[0], tol, &wrk[0]);
  // that doesn't work yet


  // Try to find the eigenvalues and -vectors of A = V^T * D * V with
  // 
  //     |1  1  1|                 |5  0  0|          |4  0  10|
  // V = |1  1 -1| / sqrt(3),  D = |0 -3  0|,  -> A = |0  4   2| / 3
  //     |1 -1  1|                 |0  0  2|          |10 2   4|
  //
  // Sage code for the matrices:
  // V = Matrix([[1, 1, 1],[1,1,-1],[1,-1,1]]) / sqrt(3)
  // D = Matrix([[5, 0, 0],[0,-3,0],[0,0,2]])
  // A = V.transpose() * D * V
  // V, D, A

  /*
  A.set(0, 0,  4.f/3);  A.set(0, 1, 0.f/3);  A.set(0, 2, 10.f/3);
  A.set(1, 0,  0.f/3);  A.set(1, 1, 4.f/3);  A.set(1, 2,  2.f/3);
  A.set(2, 0, 10.f/3);  A.set(2, 1, 2.f/3);  A.set(2, 2,  4.f/3);
  x = Vec({1,2,3});
  tol = 1.e-7;
  numIts = A.largestEigenValueAndVector(&ev, &x[0], tol, &wrk[0]);
  */

  // hmm... A.eigenvectors_right()  does not give the eigenvalues and -vectors i expected from
  // the construction isn't V supposed to be the matrix of eigenvectors and D the diagonal matrix
  // with the eignevalues? ...maybe try another example from a book or wikipedia:
  // https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Three-dimensional_matrix_example
  // i think, the matrix v is not orthogonal, so we may need to do A = V * D * inv(V) or something?
  // https://en.wikipedia.org/wiki/Diagonalizable_matrix
  // V = Matrix([[1,1,1],[1,1,-1],[1,-1,1]])
  // V.inverse()


  A.set(0, 0, 2.f);  A.set(0, 1, 0.f);  A.set(0, 2, 0.f);
  A.set(1, 0, 0.f);  A.set(1, 1, 3.f);  A.set(1, 2, 4.f);
  A.set(2, 0, 0.f);  A.set(2, 1, 4.f);  A.set(2, 2, 9.f);
  x = Vec({1,2,3});
  tol = 1.e-7f;
  //numIts = A.largestEigenValueAndVector(&ev, &x[0], tol, &wrk[0]);


  using LA = rsIterativeLinearAlgebra;
  x = Vec({1,2,3});
  numIts = LA::largestEigenValueAndVector(A, &ev, &x[0], tol, &wrk[0]);
  // val: 11, vec: (0,1,2) / sqrt(5)

  Vec vals(3);
  Vec vecs({1,2,3, 1,2,3, 1,2,3});
  numIts = LA::eigenspace(A, &vals[0], &vecs[0], tol, &wrk[0]);

  // todo: 
  // -figure out what happens when we have eigenvalues with multiplicities (algebraic and/or
  //  geometric)


  return res;
}

bool testMatrixConvolution()
{
  bool ok = true;
  using Mat = rsMatrix<double>;
  Mat A(4, 5, {1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20}); // input 1
  Mat B(3, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12});                          // input 2
  Mat T(6, 8, {  1,   4,   10,  20,   30,  34,  31,  20,                // target result
                11,  35,   74, 130,  166, 161, 133,  80,
                50, 133,  252, 410,  488, 441, 346, 200,
               125, 298,  522, 800,  878, 756, 571, 320,
               179, 399,  662, 970, 1038, 857, 625, 340,
               144, 313,  508, 730,  772, 625, 448, 240}); 
  Mat C = A.getConvolutionWith(B);
  ok &= C == T;
  return ok;
}
// for producing reference output, see:
// http://juanreyero.com/article/python/python-convolution.html

bool testMatrix()
{
  bool ok = true;

  // Tests for old matrix class:
  ok &= testMatrixArithmeticOld();
  // maybe eventually, we may bet rid of them or move them into the code attic

  // Tests for array based matrices:
  ok &= testSquareMatrixTranspose();
  ok &= testMatrixVectorMultiply();
  ok &= testMatrixMultiply();

  // Tests for the (new) dense matrix class:
  ok &= testMatrixView();
  ok &= testMatrixOperators();
  ok &= testMatrixAlloc();
  ok &= testKroneckerProduct();
  //ok &= testTransformMatrices();
  ok &= testMatrixConvolution();
  // todo: inverse, pseudo-inverse: (A^T * A)^-1 * A^T

  // Tests for the sparse matrix class:
  ok &= testSparseMatrix(); 
  ok &= testSparseMatrixSolvers();



  return ok;
}

