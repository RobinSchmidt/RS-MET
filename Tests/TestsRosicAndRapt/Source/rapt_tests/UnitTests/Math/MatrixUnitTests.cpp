//=================================================================================================
// Tests for old matrix implementation:

bool testMatrixScalarOperations()
{
  bool ok = true;

  // A = |15 20 26|
  //     |12 18 24|
  rsMatrixDbl A(2, 3);
  A.set(0,0, 15); A.set(0,1, 20); A.set(0,2, 26);
  A.set(1,0, 12); A.set(1,1, 18); A.set(1,2, 24);

  // A = A + 2 = |17 22 28|
  //             |14 20 26|
  A = A + 2.0;
  ok &= A(0,0) == 17 && A(0,1) == 22 && A(0,2) == 28;
  ok &= A(1,0) == 14 && A(1,1) == 20 && A(1,2) == 26;

  // A = A - 2 = |15 20 26|
  //             |12 18 24|
  A = A - 2.0;
  ok &= A(0,0) == 15 && A(0,1) == 20 && A(0,2) == 26;
  ok &= A(1,0) == 12 && A(1,1) == 18 && A(1,2) == 24;

  // A = 2 + A = |17 22 28|
  //             |14 20 26|
  A = 2.0 + A;
  ok &= A(0,0) == 17 && A(0,1) == 22 && A(0,2) == 28;
  ok &= A(1,0) == 14 && A(1,1) == 20 && A(1,2) == 26;

  // A = 30 - A = |13  8  2|
  //              |16 10  4|
  A = 30.0 - A;
  ok &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  ok &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A = 2.0 * A = |26 16  4|
  //               |32 20  8|
  A = 2.0 * A;
  ok &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  ok &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A = A / 2.0 = |13  8  2|
  //               |16 10  4|
  A = A / 2.0;
  ok &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  ok &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A = A * 2.0 = |26 16  4|
  //               |32 20  8|
  A = A * 2.0;
  ok &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  ok &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A /= 2.0 = |13  8  2|
  //            |16 10  4|
  A /= 2.0;
  ok &= A(0,0) == 13 && A(0,1) ==  8 && A(0,2) ==  2;
  ok &= A(1,0) == 16 && A(1,1) == 10 && A(1,2) ==  4;

  // A *= 2.0 = |26 16  4|
  //            |32 20  8|
  A *= 2.0;
  ok &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  ok &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  // A += 2.0 = |28 18  6|
  //            |34 22 10|
  A += 2.0;
  ok &= A(0,0) == 28 && A(0,1) == 18 && A(0,2) ==  6;
  ok &= A(1,0) == 34 && A(1,1) == 22 && A(1,2) == 10;

  // A -= 2.0 = |26 16  4|
  //            |32 20  8|
  A -= 2.0;
  ok &= A(0,0) == 26 && A(0,1) == 16 && A(0,2) ==  4;
  ok &= A(1,0) == 32 && A(1,1) == 20 && A(1,2) ==  8;

  return ok;
}
bool testMatrixAddSubTransEqual()
{
  bool ok = true;

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
  ok &= C.getNumRows()    == 0;
  ok &= C.getNumColumns() == 0;
  //ok &= C.mFlat == nullptr;
  //ok &= C.m     == nullptr;




  // test unary minus:
  // C = -A = |-2  -3  -5|
  //          |-7 -11 -13|
  C = -A;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) == -2 && C(0,1) ==  -3 && C(0,2) ==  -5;
  ok &= C(1,0) == -7 && C(1,1) == -11 && C(1,2) == -13;

  // test transposition:
  // C = A^T = |2  7|
  //           |3 11|
  //           |5 13|
  C = trans(A);
  ok &= C.getNumRows()    == 3;
  ok &= C.getNumColumns() == 2;
  ok &= C(0,0) ==  2 && C(0,1) ==  7;
  ok &= C(1,0) ==  3 && C(1,1) == 11;
  ok &= C(2,0) ==  5 && C(2,1) == 13;

  // test in-place transpostion and equality operator:
  // C = C^T == A
  C.transpose();
  ok &= C == A;

  // test in/equality operators:
  C.transpose();
  ok &= (A ==  B) == false;
  ok &= (A !=  B) == true;
  ok &= (A ==  C) == false;
  ok &= (A !=  C) == true;

  // C = A + B = |19 22 28|
  //             |36 42 50|
  C = A + B;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) == 19 && C(0,1) == 22 && C(0,2) == 28;
  ok &= C(1,0) == 36 && C(1,1) == 42 && C(1,2) == 50;

  // C += A = |21 25 33| 
  //          |43 53 63|
  C += A;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) == 21 && C(0,1) == 25 && C(0,2) == 33;
  ok &= C(1,0) == 43 && C(1,1) == 53 && C(1,2) == 63;

  // C -= A = |19 22 28|
  //          |36 42 50|
  C -= A;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) == 19 && C(0,1) == 22 && C(0,2) == 28;
  ok &= C(1,0) == 36 && C(1,1) == 42 && C(1,2) == 50;

  // C = B - A = |15 16 18|
  //             |22 20 24|
  C = B - A;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) == 15 && C(0,1) == 16 && C(0,2) == 18;
  ok &= C(1,0) == 22 && C(1,1) == 20 && C(1,2) == 24;

  return ok;
}
bool testMatrixMul()
{
  bool ok = true;

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
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 4;
  ok &= C(0,0) ==  362 && C(0,1) ==  414 && C(0,2) ==  464 && C(0,3) ==  492;
  ok &= C(1,0) == 1071 && C(1,1) == 1229 && C(1,2) == 1379 && C(1,3) == 1469;

  // same with the accumulative multiplier:
  C  = A;
  C *= B;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 4;
  ok &= C(0,0) ==  362 && C(0,1) ==  414 && C(0,2) ==  464 && C(0,3) ==  492;
  ok &= C(1,0) == 1071 && C(1,1) == 1229 && C(1,2) == 1379 && C(1,3) == 1469;

  // B2 = B * B^T = |2020 3420  4932|
  //                |3420 5860  8460|
  //                |4932 8460 12220|
  rsMatrixDbl B2 = B * trans(B);
  ok &= B2.getNumRows()    == 3;
  ok &= B2.getNumColumns() == 3;
  ok &= B2(0,0) == 2020 && B2(0,1) == 3420 && B2(0,2) ==  4932;
  ok &= B2(1,0) == 3420 && B2(1,1) == 5860 && B2(1,2) ==  8460;
  ok &= B2(2,0) == 4932 && B2(2,1) == 8460 && B2(2,2) == 12220;

  // C = A * B2 = | 38960  66720  96344|
  //              |115876 198380 286444|
  C = A * B2;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) ==  38960 && C(0,1) ==  66720 && C(0,2) ==  96344;
  ok &= C(1,0) == 115876 && C(1,1) == 198380 && C(1,2) == 286444;
  C  = A;
  C *= B2;
  ok &= C.getNumRows()    == 2;
  ok &= C.getNumColumns() == 3;
  ok &= C(0,0) ==  38960 && C(0,1) ==  66720 && C(0,2) ==  96344;
  ok &= C(1,0) == 115876 && C(1,1) == 198380 && C(1,2) == 286444;

  return ok;
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
  bool ok = true;

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
  ok &= y[0] == 14;
  ok &= y[1] == 32;
  ok &= y[2] == 50;

  // transposed 3x3-matrix times 3-vector:
  // |1 2 3|^T * |1| = |1 4 7| * |1| = |30|
  // |4 5 6|     |2|   |2 5 8|   |2|   |36|
  // |7 8 9|     |3|   |3 6 9|   |3|   |42|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 3, 3);
  ok &= y[0] == 30;
  ok &= y[1] == 36;
  ok &= y[2] == 42;

  // 3x2-matrix times 2-vector:
  // |1 2| * |1| = | 5|
  // |4 5|   |2|   |14|
  // |7 8|         |23|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::matrixVectorMultiply(pA, x, y, 3, 2);
  ok &= y[0] == 5;
  ok &= y[1] == 14;
  ok &= y[2] == 23;

  // transposed 2x3-matrix times 2-vector:
  // |1 2 3|^T * |1| = |1 4| * |1| = | 9|
  // |4 5 6|     |2|   |2 5|   |2|   |12|
  //                   |3 6|         |15|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 2, 3);
  ok &= y[0] == 9;
  ok &= y[1] == 12;
  ok &= y[2] == 15;

  // 2x3-matrix times 3-vector:
  // |1 2 3| * |1| = |14|
  // |4 5 6|   |2|   |32|
  //           |3|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::matrixVectorMultiply(pA, x, y, 2, 3);
  ok &= y[0] == 14;
  ok &= y[1] == 32;
  ok &= y[2] == -1;

  // transposed 3x2-matrix times 3-vector:
  // |1 2|^T * |1| = |1 4 7| * |1| = |30|
  // |4 5|     |2|   |2 5 8|   |2|   |36|
  // |7 8|     |3|             |3|
  rsArrayTools::fillWithValue(y, 3, -1.0);
  rsMatrixTools::transposedMatrixVectorMultiply(pA, x, y, 3, 2);
  ok &= y[0] == 30;
  ok &= y[1] == 36;
  ok &= y[2] == -1;

  return ok;
}

bool testMatrixMultiply3x3()
{
  // This function tests only some special cases where P either equals N or M using matrices where
  // no index-range exceeds 3. It is called from the more general test function
  // testMatrixMultiply3x3 internally.

  bool ok = true;

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
  ok &= C[0][0] == 30 && C[0][1] == 24 && C[0][2] == -1;
  ok &= C[1][0] == 84 && C[1][1] == 69 && C[1][2] == -1;
  ok &= C[2][0] == -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x2-matrix times 2x3 matrix:
  // |1 2| * |9 8 7| = | 21 18 15|
  // |4 5|   |6 5 4|   | 66 57 48|
  // |7 8|             |111 96 81|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiply(pA, pB, pC, 3, 2, 3);
  ok &= C[0][0] ==  21 && C[0][1] == 18 && C[0][2] == 15;
  ok &= C[1][0] ==  66 && C[1][1] == 57 && C[1][2] == 48;
  ok &= C[2][0] == 111 && C[2][1] == 96 && C[2][2] == 81;

  // transposed 3x2-matrix times 3x2 matrix:
  // |1 2|^T * |9 8| = |1 4 7| * |9 8| = |54 42|
  // |4 5|     |6 5|   |2 5 8|   |6 5|   |72 57|
  // |7 8|     |3 2|             |3 2|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 3, 2, 2);
  ok &= C[0][0] == 54 && C[0][1] == 42 && C[0][2] == -1;
  ok &= C[1][0] == 72 && C[1][1] == 57 && C[1][2] == -1;
  ok &= C[2][0] == -1 && C[2][1] == -1 && C[2][2] == -1;

  // transposed 2x3-matrix times 2x3 matrix:
  // |1 2 3|^T * |9 8 7| = |1 4| * |9 8 7| = |33 28 23|
  // |4 5 6|     |6 5 4|   |2 5|   |6 5 4|   |48 41 34|
  //                       |3 6|             |63 54 45|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 2, 3, 3);
  ok &= C[0][0] == 33 && C[0][1] == 28 && C[0][2] == 23;
  ok &= C[1][0] == 48 && C[1][1] == 41 && C[1][2] == 34;
  ok &= C[2][0] == 63 && C[2][1] == 54 && C[2][2] == 45;

  // 2x3-matrix times transposed 2x3 matrix:
  // |1 2 3| * |9 8 7|^T = |1 2 3| * |9 6| = | 46 28|
  // |4 5 6|   |6 5 4|     |4 5 6|   |8 5|   |118 73|
  //                                 |7 4|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 2, 3, 2);
  ok &= C[0][0] ==  46 && C[0][1] == 28 && C[0][2] == -1;
  ok &= C[1][0] == 118 && C[1][1] == 73 && C[1][2] == -1;
  ok &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x2-matrix times transposed 3x2 matrix:
  // |1 2| * |9 8|^T = |1 2| * |9 6 3| = | 25 16  7|
  // |4 5|   |6 5|     |4 5|   |8 5 2|   | 76 49 22|
  // |7 8|   |3 2|     |7 8|             |127 82 37|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 3, 2, 3);
  ok &= C[0][0] ==  25 && C[0][1] == 16 && C[0][2] ==  7;
  ok &= C[1][0] ==  76 && C[1][1] == 49 && C[1][2] == 22;
  ok &= C[2][0] == 127 && C[2][1] == 82 && C[2][2] == 37;

  // transposed 2x3-matrix times transposed 3x2 matrix:
  // |1 2 3|^T * |9 8|^T = |1 4| * |9 6 3| = |41 26 11|
  // |4 5 6|     |6 5|     |2 5|   |8 5 2|   |58 37 16|
  //             |3 2|     |3 6|             |75 48 21|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 2, 3, 3);
  ok &= C[0][0] == 41 && C[0][1] == 26 && C[0][2] == 11;
  ok &= C[1][0] == 58 && C[1][1] == 37 && C[1][2] == 16;
  ok &= C[2][0] == 75 && C[2][1] == 48 && C[2][2] == 21;

  // transposed 3x2-matrix times transposed 2x3 matrix:
  // |1 2|^T * |9 8 7|^T = |1 4 7| * |9 6| = | 90 54|
  // |4 5|     |6 5 4|     |2 5 8|   |8 5|   |114 69|
  // |7 8|                           |7 4|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 3, 2, 2);
  ok &= C[0][0] ==  90 && C[0][1] == 54 && C[0][2] == -1;
  ok &= C[1][0] == 114 && C[1][1] == 69 && C[1][2] == -1;
  ok &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  // 3x3-matrix times 3x3-matrix in place multiplication:
  // |1 2 3| * |9 8 7| = | 30  24 18|
  // |4 5 6|   |6 5 4|   | 84  69 54|
  // |7 8 9|   |3 2 1|   |138 114 90|
  rsMatrixTools::copyMatrix(pA, pC, 3, 3);
  rsMatrixTools::matrixInPlaceMultiply(pC, pB, 3, 3);
  ok &= C[0][0] ==  30 && C[0][1] ==  24 && C[0][2] == 18;
  ok &= C[1][0] ==  84 && C[1][1] ==  69 && C[1][2] == 54;
  ok &= C[2][0] == 138 && C[2][1] == 114 && C[2][2] == 90;

  // 2x3-matrix times 3x3-matrix in place multiplication:
  // |1 2 3| * |9 8 7| = | 30  24 18|
  // |4 5 6|   |6 5 4|   | 84  69 54|
  //           |3 2 1|
  rsMatrixTools::initMatrix(pC, 3, 3, -1.0);
  rsMatrixTools::copyMatrix(pA, pC, 2, 3);
  rsMatrixTools::matrixInPlaceMultiply(pC, pB, 2, 3);
  ok &= C[0][0] ==  30 && C[0][1] ==  24 && C[0][2] == 18;
  ok &= C[1][0] ==  84 && C[1][1] ==  69 && C[1][2] == 54;
  ok &= C[2][0] ==  -1 && C[2][1] == -1 && C[2][2] == -1;

  return ok;
}

bool testMatrixMultiply()
{
  std::string testName = "MatrixMultiply";
  bool ok = testMatrixMultiply3x3();

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
  ok &= C[0][0] ==  822 && C[0][1] ==  864 && C[0][2] ==  898 && C[0][3] ==  944;
  ok &= C[1][0] == 3247 && C[1][1] == 3415 && C[1][2] == 3567 && C[1][3] == 3757;
  ok &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  ok &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // 2x3 matrix times transposed 3x4 matrix:
  // | 2  3  5| * | 59  61  67|^T = | 636  798 1012 1192|
  // |11 13 17|   | 73  79  83|     |2581 3241 4131 4827|
  //              | 97 101 103|
  //              |109 113 127|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplySecondTransposed(pA, pB, pC, 2, 3, 4);
  ok &= C[0][0] ==  636 && C[0][1] ==  798 && C[0][2] == 1012 && C[0][3] == 1192;
  ok &= C[1][0] == 2581 && C[1][1] == 3241 && C[1][2] == 4131 && C[1][3] == 4827;
  ok &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  ok &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // transposed 3x2 matrix times 3x4 matrix:
  // | 2  3|^T * |59  61  67  71| = |3152 3314 3416 3582|
  // |11 13|     |73  79  83  89|   |3939 4139 4267 4473|
  // |23 29|     |97 101 103 107|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplyFirstTransposed(pA, pB, pC, 3, 2, 4);
  ok &= C[0][0] == 3152 && C[0][1] == 3314 && C[0][2] == 3416 && C[0][3] == 3582;
  ok &= C[1][0] == 3939 && C[1][1] == 4139 && C[1][2] == 4267 && C[1][3] == 4473;
  ok &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  ok &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  // transposed 3x2 matrix times transposed 4x3 matrix:
  // | 2  3|^T * | 59  61  67|^T = |2330 2924 3647 4382|
  // |11 13|     | 73  79  83|     |2913 3653 4691 5479|
  // |23 29|     | 97 101 103|
  //             |109 113 127|
  rsMatrixTools::initMatrix(pC, 4, 4, -1.0);
  rsMatrixTools::matrixMultiplyBothTransposed(pA, pB, pC, 3, 2, 4);
  ok &= C[0][0] == 2330 && C[0][1] == 2924 && C[0][2] == 3674 && C[0][3] == 4382;
  ok &= C[1][0] == 2913 && C[1][1] == 3653 && C[1][2] == 4591 && C[1][3] == 5479;
  ok &= C[2][0] ==   -1 && C[2][1] ==   -1 && C[2][2] ==   -1 && C[2][3] ==   -1;
  ok &= C[3][0] ==   -1 && C[3][1] ==   -1 && C[3][2] ==   -1 && C[3][3] ==   -1;

  return ok;
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
  bool ok = true;  // test result
  using Mat = rsMatrix<double>;
  using Vec = std::vector<double>;

  Mat A(3, 3, { 2,1,4, 3,10,3, 1,5,1 });
  Mat X(3, 2, {1,4, 2,5, 3,6});
  Mat B = A * X;
  ok &= B == Mat(3, 2, {16,37, 32,80, 14,35});

  // test matrix-vector multiplication:
  Vec x, y; x = Vec({1,2,3});
  A = Mat(2, 3, {1,2,3, 4,5,6});  y = A*x; ok &= y == Vec({14,32});
  A = Mat(3, 2, {1,2, 3,4, 5,6}); y = x*A; ok &= y == Vec({22,28});

  // test element-wise +=, *=, -=, /=
  double tmp;
  tmp = A(2,1); A(2,1) += 2.0; ok &= A(2,1) == tmp + 2.0;
  tmp = A(2,1); A(2,1) *= 2.0; ok &= A(2,1) == tmp * 2.0;
  tmp = A(2,1); A(2,1) -= 2.0; ok &= A(2,1) == tmp - 2.0;
  tmp = A(2,1); A(2,1) /= 2.0; ok &= A(2,1) == tmp / 2.0;

  return ok;
}

bool testMatrixAlloc() // rename to testMatrixAllocationAndArithmetic
{
  // ToDo: use the new rsLoggingVector as new 2nd template argument for rsMatrix and use that for
  // figuring out the number of allocs. Maybe first make a unit test of rsLoggingVector

  bool ok = true;

  using Vector = std::vector<double>;
  using LV = rsLoggingVector<double>;
  using Matrix = rsMatrix<double, LV>;
  size_t& allocs = LV::numPotentialAllocs; // counts potential re-allocations
  allocs = 0;

  // A = |1 2 3|
  //     |4 5 6|
  Matrix A(2, 3, {1.,2.,3., 4.,5.,6.}); // calls rsMatrix(int, int, std::vector<T>&&)
  ok &= A(0,0) == 1 &&  A(0,1) == 2 && A(0,2) == 3;
  ok &= A(1,0) == 4 &&  A(1,1) == 5 && A(1,2) == 6;
  ok &= allocs == 1;

  // B = |1 2| 
  //     |3 4|
  //     |5 6|
  Matrix B(3, 2, {1.,2.,3., 4.,5.,6.}); // calls rsMatrix(int, int, std::vector<T>&&)
  ok &= B(0,0) == 1 &&  B(0,1) == 2;
  ok &= B(1,0) == 3 &&  B(1,1) == 4;
  ok &= B(2,0) == 5 &&  B(2,1) == 6;
  ok &= allocs == 2;

  // matrix-vector multiplication:
  Vector x({1,2,3});
  Vector y = A * x;
  ok &= y == Vector({14,32});

  // matrix multiplication:
  Matrix C = A*B; 
  // calls: 
  //   operator*(const rsMatrix<T>&)
  //     rsMatrix(int numRows, int numColumns)
  //     rsMatrix(rsMatrix&& B)
  ok &= allocs == 3;
  ok &= C(0,0) == 22 &&  C(0,1) == 28;
  ok &= C(1,0) == 49 &&  C(1,1) == 64;

  C = A;  // calls copy assigment operator
  ok &= allocs == 4;
  ok &= C == A;

  A = A;  // self copy assignment - should not re-allocate
  ok &= allocs == 4;

  Matrix D = A+A;
  ok &= allocs == 5;
  ok &= D == Matrix(2, 3, {2.,4.,6., 8.,10.,12.});
  ok &= allocs == 6;
  ok &= D == A+A;
  ok &= allocs == 7;

  C = B*A;  // calls move assignment operator
  ok &= allocs == 8;   // fails here
  ok &= C(0,0) ==  9 &&  C(0,1) == 12 && C(0,2) == 15;
  ok &= C(1,0) == 19 &&  C(1,1) == 26 && C(1,2) == 33;
  ok &= C(2,0) == 29 &&  C(2,1) == 40 && C(2,2) == 51;


  Matrix E(A+D);  // calls move constructor
  ok &= allocs == 9;
  ok &= E == A+D;  // temporary A+D needs allocation
  ok &= allocs == 10; 

  Matrix F(E);   // calls copy constructor
  ok &= allocs == 11;
  ok &= F == E;   // no allocation here
  ok &= allocs == 11;

  // hmm - these here call move/copy constructors and not move/copy assigment operators - why
  // ..ah - it's because it's not a re-assignment
  Matrix G = A+D; // calls move constructor
  ok &= allocs == 12;
  ok &= G == A+D;  // temporary A+D needs allocation
  ok &= allocs == 13;

  Matrix H = G;   // calls copy constructor
  ok &= allocs == 14;
  ok &= H == G;
  ok &= allocs == 14;


  Matrix I, J;  // no allocations yet
  ok &= allocs == 14;
  I = A+D;      // move assignment
  ok &= allocs == 15;
  J = I;        // copy assignment
  ok &= allocs == 16;


  // multiplication with a scalar:
  I = 2.0 * A;
  ok &= allocs == 17;
  ok &=  I == Matrix(2, 3, {2.,4.,6., 8.,10.,12.});
  ok &= allocs == 18;

  J = A * 2.0;
  ok &= allocs == 19;
  ok &= I == J;
  ok &= allocs == 19;

  // in-place addition and subtraction via +=, -=:
  J += A;
  ok &= allocs == 19;
  ok &= J == 3.0 * A;
  ok &= allocs == 20;
  J -= A;
  ok &= allocs == 20;
  ok &= J == 2.0 * A;
  ok &= allocs == 21;

  // the *= operator can't operate in-place:
  J = A;    // should not re-allocate because J already has the right shape
  ok &= allocs == 21;
  J *= B;   // allocates
  ok &= allocs == 22;
  ok &= J == A*B;
  ok &= allocs == 23;

  // in-place scaling by a scalar:
  J = A;    
  ok &= allocs == 24;
  J *= 2.0;
  ok &= allocs == 24;
  ok &= J == 2.0 * A;
  ok &= allocs == 25;
  J /= 2.0;
  ok &= allocs == 25;
  ok &= J == A;
  ok &= allocs == 25;


  C = Matrix(2, 4, {1,2,3,4, 5,6,7,8});     // 2x4
  ok &= allocs == 26;
  D = Matrix(4, 2, {1,2, 3,4, 5,6, 7,8});   // 4x2
  ok &= allocs == 27;

  E = B*C;                                  // 3x2 * 2x4 = 3x4
  ok &= allocs == 28;
  ok &= E.getNumRows()    == 3;
  ok &= E.getNumColumns() == 4;
  ok &= E(0,0) == 11 && E(0,1) == 14 && E(0,2) == 17 && E(0,3) == 20;
  ok &= E(1,0) == 23 && E(1,1) == 30 && E(1,2) == 37 && E(1,3) == 44;
  ok &= E(2,0) == 35 && E(2,1) == 46 && E(2,2) == 57 && E(2,3) == 68;

  E = D*A;
  ok &= allocs == 29;               // 4x2 * 2x3 = 4x3
  ok &= E.getNumRows()    == 4;
  ok &= E.getNumColumns() == 3;
  ok &= E(0,0) ==  9 && E(0,1) == 12 && E(0,2) == 15;
  ok &= E(1,0) == 19 && E(1,1) == 26 && E(1,2) == 33;
  ok &= E(2,0) == 29 && E(2,1) == 40 && E(2,2) == 51;
  ok &= E(3,0) == 39 && E(3,1) == 54 && E(3,2) == 69;

  // move constructor:
  Matrix K(std::move(E));
  ok &= allocs == 29;
  ok &= K.getNumRows() == 4 && K.getNumColumns() == 3;
  ok &= E.getNumRows() == 0 && E.getNumColumns() == 0;
  ok &= K.getDataVectorConst().size() == 12;
  ok &= E.getDataVectorConst().size() == 0;
  ok &= K.getDataPointerConst() != nullptr;
  ok &= E.getDataPointerConst() == nullptr;

  // move assignment operator:
  E = std::move(K);
  ok &= allocs == 29;
  ok &= E.getNumRows() == 4 && E.getNumColumns() == 3;
  ok &= K.getNumRows() == 0 && K.getNumColumns() == 0;
  ok &= E.getDataVectorConst().size() == 12;
  ok &= K.getDataVectorConst().size() == 0;
  ok &= E.getDataPointerConst() != nullptr;
  ok &= K.getDataPointerConst() == nullptr;

  // transposition:
  E = A;
  ok &= allocs == 30;
  ok &= E.getNumRows() == 2 && E.getNumColumns() == 3;
  E.transpose();
  ok &= allocs == 31;
  ok &= E == Matrix(3, 2, {1,4, 2,5, 3,6});
  ok &= allocs == 32;

  Matrix S4(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});
  ok &= allocs == 33;
  S4.transpose();
  ok &= allocs == 33;
  ok &= S4 == Matrix(4, 4, {1,5,9,13, 2,6,10,14, 3,7,11,15, 4,8,12,16});
  ok &= allocs == 34;

  Matrix V3(3, 1, {1,2,3});     // a 3D column-vector...
  ok &= allocs == 35;
  V3.transpose();               // ...becomes a 3D row-vector
  ok &= allocs == 35;
  ok &= V3 == Matrix(1, 3, {1,2,3});
  ok &= allocs == 36;

  // factory functions:
  E = Matrix::identity(4, 1.0);
  ok &= allocs == 37;
  ok &= E == Matrix(4, 4, {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1});
  ok &= allocs == 38;

  E = Matrix::zero(3, 4, 0.0);
  ok &= allocs == 39;
  ok &= E == Matrix(3, 4, {0,0,0,0, 0,0,0,0, 0,0,0,0});
  ok &= allocs == 40;

  // swapping two matrices:
  RAPT::rsSwapNaive(A, B);
  ok &= allocs == 43;   // naive swapping causes 3 heap allocations...
  RAPT::rsSwapViaMove(A, B);
  ok &= allocs == 43;   // ...whereas swapping via std::move causes none
  RAPT::rsSwap(A, B);
  ok &= allocs == 43;   // rsSwap calls rsSwapViaMove


  // create diagonal matrix



  // todo:
  // -maybe add/subtract scalars (in-place and out-of-place)....but maybe not: mathematically, the
  //  space of MxN matrices is a vector space and addition of a scalar and a vector is not a thing
  //  in vector space

  return ok;
}

bool testKroneckerProduct()
{
  bool ok = true;

  using LV = rsLoggingVector<double>;
  using Matrix = rsMatrix<double, LV>;
  size_t& allocs = LV::numPotentialAllocs; // counts potential re-allocations

  //typedef rsMatrix<double> Matrix;
  //int& allocs  = Matrix::numHeapAllocations;  // to count allocations

  allocs = 0;

  Matrix A(2, 3, {1,2,3,4,5,6});
  ok &= allocs == 1;
  Matrix B(4, 5, {1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20});
  ok &= allocs == 2;

  Matrix C = Matrix::getKroneckerProduct(A, B);
  ok &= allocs == 3;
  ok &= C.getNumRows() == 8 && C.getNumColumns() == 15;

  Matrix T(8, 15,   // target Kronecker product
    { 1, 2, 3, 4, 5, 2, 4, 6, 8, 10, 3,  6,  9, 12, 15,
      6, 7, 8, 9,10,12,14,16,18, 20,18, 21, 24, 27, 30,
     11,12,13,14,15,22,24,26,28, 30,33, 36, 39, 42, 45,
     16,17,18,19,20,32,34,36,38, 40,48, 51, 54, 57, 60,
      4, 8,12,16,20, 5,10,15,20, 25, 6, 12, 18, 24, 30,
     24,28,32,36,40,30,35,40,45, 50,36, 42, 48, 54, 60,
     44,48,52,56,60,55,60,65,70, 75,66, 72, 78, 84, 90,
     64,68,72,76,80,80,85,90,95,100,96,102,108,114,120});
  ok &= allocs == 4;
  ok &= C == T;

  // the following sage code produces the result:
  // A = matrix([[1,2,3],[4,5,6]])
  // B = matrix([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]])
  // C = A.tensor_product(B)
  // A, B, C

  return ok;
}

bool testSparseMatrix()
{
  // We test the class rsSparseMatrix. We first check that the function .toDense() to convert from 
  // sparse to dense format works correctly and then perform many of the subsequent tests for 
  // matrix operations (like matrix addition, multplication, etc.) by doing the same operation in 
  // parallel on dense matrices and after every operation, convert the (converted) result of the 
  // sparse computation to that of the dense computation. Thus, these tests assume that the class 
  // rsMatrix that implements dense matrices is already reliable enough to serve as a reference to 
  // check the sparse computations against. We use letters A,B,C for dense matrices and R,S,T for 
  // sparse matrices.

  bool ok = true;

  using Real = float;
  using Vec  = std::vector<Real>;
  using Mat  = rsMatrix<Real>;
  using MatS = rsSparseMatrix<Real>;


  // We create a matrix of the form:
  // 
  //      |1 1 0 0 0 0 0 0|
  //  A = |0 0 1 1 0 0 0 0|
  //      |0 0 0 0 1 1 0 0|
  //      |0 0 0 0 0 0 1 1|
  //
  // so, in a matrix-vector product y = A*x, it takes an 8-vector x as input and produces a 4-vector 
  // as output y in which each element is the sum of two consecutive elements from the input.


  rsSparseMatrix<Real> R(4, 8);  // still the zero matrix
  Vec x({1,2,3,4,5,6,7,8});      // input vector
  Vec y({1,2,3,4});              // output vector
  R.product(&x[0], &y[0]);  
  ok &= y == Vec({0,0,0,0});    // A was still zero
  // todo: use * operator: y = A*x (it calls the product function)

  // Now build up the actual matrix in a kind of "random" order. The nonzero (i.e. one) elements 
  // are at positions: (0,0),(0,1),(1,2),(1,3),(2,4),(2,5),(3,6),(3,7):
  R.set(1, 2, 1.f);
  R.set(2, 5, 1.f);
  R.set(0, 0, 1.f);
  R.set(0, 1, 1.f);
  R.set(3, 6, 1.f);
  R.set(1, 3, 1.f);
  R.set(2, 4, 1.f);
  R.set(3, 7, 1.f);

  // Create the same matrix as dense matrix:

  Mat A(4, 8, {1, 1, 0, 0, 0, 0, 0, 0,
               0, 0, 1, 1, 0, 0, 0, 0,
               0, 0, 0, 0, 1, 1, 0, 0,
               0, 0, 0, 0, 0, 0, 1, 1});

  // Test conversion of A to dense matrix, structure of A and conversion of a dense matrix
  // to a sparse matrix:
  Mat  dR = MatS::toDense(R);    ok &= dR == A;
  MatS R2 = MatS::fromDense(dR); ok &= R2 == R;

  // Test the multiplication with the new matrix again:
  R.product(&x[0], &y[0]); ok &= y == Vec({3,7,11,15});  // redundant because...
  y = R*x;                 ok &= y == Vec({3,7,11,15});  // ...A*x invokes A.product

  // Test transposition and matrix multiplication:
  MatS RT = R.getTranspose();
  Mat  TT = A.getTranspose();
  ok &= MatS::toDense(RT) == TT;

  // Now we build some slightly more interesting matrices:
  A = Mat(3, 5, { 2, 3,  0,  0,  5,
                  0, 7,  0, 11,  0,
                 13, 0, 17,  0, 19});
  Mat B;
  B = Mat(5, 3, { 23,  0, 29,
                   0, 31,  0,
                  37, 41,  0,
                   0,  0, 53,
                  59, 61, 67 } );

  // Test matrix multiplication:
  Mat C = A*B;

  R.setToZero();
  R.setShape(3, 5);
  R.set(0, 0,  2);
  R.set(0, 1,  3);
  R.set(0, 4,  5);
  R.set(1, 1,  7);
  R.set(1, 3, 11);
  R.set(2, 0, 13);
  R.set(2, 2, 17);
  R.set(2, 4, 19);

  MatS S(5, 3);
  S.set(0, 0, 23);
  S.set(0, 2, 29);
  S.set(1, 1, 31);
  S.set(2, 0, 37);
  S.set(2, 1, 41);
  S.set(3, 2, 53);
  S.set(4, 0, 59);
  S.set(4, 1, 61);
  S.set(4, 2, 67);

  MatS T = R*S;
  Mat  D = MatS::toDense(T);
  ok &= C == D;

  // Test matrix addition, subtraction and negation:
  B.transpose();
  S.transpose();
  C = A + B;
  T = R + S;
  D = MatS::toDense(T);
  ok &= C == D;

  C = A - B;
  T = R - S;
  D = MatS::toDense(T);
  ok &= C == D;

  C = -A;
  T = -R;
  D = MatS::toDense(T);
  ok &= C == D;


  // ...



  // ToDo: 
  // -Implement and test matrix subtraction
  // -Implement matrix negation (unary minus). 
  // -Test replacing elements, also with zero (in which case they should get removed)
  // -Make sure that no element is stored twice or multiple times. Each index pair should be 
  //  unique
  // -test setShape. Verify that elements get removed when they are out of range of the new 
  //  shape in case of decreasing the size. Test also the special case of setting the shape
  //  to 0,0. This should lead to all elements being removed
  // -Implement at test transposition
  // -Implement and test sparse matrix multiplication. The tests for that can be based on 
  //  performing the same products with dense matrices, converting the sparse results to dense
  //  and comparing the results.

  return ok;
}

// move to LinearAlgebraUnitTests:
bool testSparseMatrixSolvers()
{
  bool res = true;

  using Vec = std::vector<float>;
  using Mat = rsSparseMatrix<float>;

  //      |7 1 2|
  //  A = |3 8 4|  ...must be strictly diagonally dominant for Gauss-Seidel to converge
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
  numIts = LA::eigenspace(A, &vals[0], &vecs[0], tol, &wrk[0], 1000);

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


// Convenience function to convert a matrix of integers into a matrix of modular integers with 
// given modulus:
template<class T>
rsMatrix<rsModularInteger<T>> toModular(const rsMatrix<T>& A, T modulus)
{
  int M = A.getNumRows();
  int N = A.getNumColumns();
  rsMatrix<rsModularInteger<T>> B(M, N);
  for(int i = 0; i < M; i++)
    for(int j = 0; j < N; j++)
      B(i, j) = rsModularInteger<T>(A(i,j), modulus);
  return B;
}

bool testModIntMatrix()
{
  // We test matrices of modular integers.

  bool ok = true;

  using Int  = int;
  using Mod  = rsModularInteger<Int>;
  using MatI = rsMatrix<Int>;
  using MatM = rsMatrix<Mod>;

  // Test matrix multiplication with matrices of modular integers with modulus m = 7:
  Int  m = 7;                         // Modulus to use
  MatI Ai(2, 3, {9,4,2, 4,6,8});      // Matrix A as integer matrix
  MatI Bi(3, 2, {7,5, 2,8, 5,3});     // ..same for B
  MatM Am = toModular(Ai, m);         // A converted to modular integers
  MatM Bm = toModular(Bi, m);         // ..same for B
  MatI Ci = Ai * Bi;                  // Product C = A * B as integer matrix
  MatM Ct = toModular(Ci, m);         // Target result: C as modular integers
  MatM Cm = Am * Bm;                  // Actual result
  ok &= Cm == Ct;                     // Test
  Ci = Bi * Ai;                       // Now the same for the product B * A, i.e. with swapped
  Ct = toModular(Ci, m);              // factors
  Cm = Bm * Am;
  ok &= Cm == Ct;





  return ok;

  // ToDo:
  //
  // - Try to solve a linear system of equations.
}


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
  ok &= testModIntMatrix();
  // todo: inverse, pseudo-inverse: (A^T * A)^-1 * A^T

  // Tests for the sparse matrix class:
  ok &= testSparseMatrix(); 
  ok &= testSparseMatrixSolvers();



  return ok;
}

