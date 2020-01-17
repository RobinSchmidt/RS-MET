
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
bool testMatrixArithmetic(std::string &reportString)
{
  std::string testName = "MatrixArithmetic";
  bool testResult = true;

  testResult &= testMatrixScalarOperations();
  testResult &= testMatrixAddSubTransEqual();
  testResult &= testMatrixMul();

  // intitalize some test matrices:
  // A = |2  3  5|
  //     |7 11 13|
  // 
  rsMatrixDbl A(2, 3);
  A.set(0,0, 2); A.set(0,1,  3); A.set(0,2,  5);
  A.set(1,0, 7); A.set(1,1, 11); A.set(1,2, 13);

  A.applyFunction(&square<double>);
  testResult &= A(0,0) ==  4 && A(0,1) ==   9 && A(0,2) ==  25;
  testResult &= A(1,0) == 49 && A(1,1) == 121 && A(1,2) == 169;   

  return testResult;
}


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
  //  0  1  2  3  4  5  6  7  8
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



  // have a function isStorageRowMajor() which determines if it is in row-major format from the 
  // strides - the only case where this is impossible is for 1x1 matrices - which are actually 
  // row-major and column-major at the same time
  // we can determine




  return r;
}

bool testMatrixOperators()
{
  bool r = true;  // test result
  using Matrix = rsMatrix<double>;
  using Vector = std::vector<double>;

  Matrix A(3, 3, { 2,1,4, 3,10,3, 1,5,1 });
  Matrix X(3, 2, {1,4, 2,5, 3,6});
  Matrix B = A * X;
  r &= B == Matrix(3, 2, {16,37, 32,80, 14,35});
  
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

bool testMatrix()
{
  std::string dummy;
  bool testResult = true;

  //rsMatrixDbl A(2, 3, true);
  //A = A + 2.0;


  testResult &= testMatrixArithmetic(dummy);  
  // obsolete - tests the old matrix class - todo: write similar tests for the new class
  //...

  // tests for the nwe matrix class:
  testResult &= testMatrixView();
  testResult &= testMatrixOperators();
  testResult &= testMatrixAlloc();
  testResult &= testKroneckerProduct();


  //testResult &= testTransformMatrices();

  return testResult;
}

