#include "MatrixUnitTests.h"

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

bool testMatrixView()
{
  std::string testName = "MatrixView";
  bool r = true;  // est result

  typedef rsMatrixView<double> MatrixView;

  double A6[6] = { 1,2,3,4,5,6 }; // array of 6 elements

  // 1 2 3
  // 4 5 6
  MatrixView m(2, 3, A6);
  r &= m(0,0) == 1; r &= m(0,1) == 2; r &= m(0,2) == 3;
  r &= m(1,0) == 4; r &= m(1,1) == 5; r &= m(1,2) == 6;

  // 1 2 
  // 3 4 
  // 5 6
  m.reshape(3, 2);
  r &= m(0,0) == 1; r &= m(0,1) == 2; 
  r &= m(1,0) == 3; r &= m(1,1) == 4;
  r &= m(2,0) == 5; r &= m(2,1) == 6; 






  return r;
}

bool testMatrix2()
{
  std::string testName = "Matrix2";
  bool testResult = true;



  return testResult;
}

bool testTransformMatrices()
{
  bool testResult = true;






  return testResult;
}

bool testMatrix()
{
  std::string dummy;
  bool testResult = true;

  //rsMatrixDbl A(2, 3, true);
  //A = A + 2.0;


  testResult &= testMatrixArithmetic(dummy);
  //...


  testResult &= testMatrixView();
  testResult &= testMatrix2();
  testResult &= testTransformMatrices();

  return testResult;
}

