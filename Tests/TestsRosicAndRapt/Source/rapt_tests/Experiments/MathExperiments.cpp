using namespace RAPT;
using namespace std;


// maybe do eigenvalue and eigenvector stuff here and rename the function accordingly
void characteristicPolynomial() 
{
  // We try to figure out how to compute the coeffs (or better: roots) of the characteristic 
  // polynomial of a matrix via Gaussian elimination by creating a matrix of rational functions and
  // see what happens when we unleash the elimination algorithm on it.

  // The characteristic polynomial p(x) of a matrix A is given by p(x) = det(A - x*I) where I is 
  // the identity matrix. The roots of p are the eigenvalues of A. The elimination algorithm does 
  // not change these roots because the elementary row operations can at most multiply the 
  // determinant by a nonzero factor (in case of a regular matrix). But we cannot just run the 
  // elimination algo on a matrix of numbers and then do p(x) = det(A' - x*I) where A' is the 
  // result of the elimination because during the process, the x would have to get mangled into 
  // the matrix itself - we can't introduce it after the elimination. So, what we do instead is to
  // consider a full-blown matrix of rsRationalFunction objects....

  using Poly   = RAPT::rsPolynomial<double>;
  using Matrix = RAPT::rsMatrix<double>;

  // Create matrix A:
  // A = |-4  6|
  //     |-3  5|
  Matrix A(2, 2, {-4,6, -3,5});
  Poly p = getCharacteristicPolynomial(A, 0.0);
  // this is example 39.1 in Karpfinger

  // The matrix plugged into its characteristic polynomial as argument should yield the zero 
  // matrix. Matrices are "roots" of their own characteristic polynomials, see
  // https://en.wikipedia.org/wiki/Cayley%E2%80%93Hamilton_theorem
  Matrix I(2, 2, {1,0, 0,1});      // 2x2 identity matrix
  Matrix pA = p[0]*I + p[1]*A + p[2] * A*A;

  //bool test = pA.isZero();

  // the sum of the eigenvalues should equal the trace of the matrix (Karpf. pg 412)

  //vector<complex<double>> eigenvalues1 = getPolynomialRoots(p);
  //vector<complex<double>> eigenvalues2 = getEigenvalues(A);

  // this is wrong:
  //vector<vector<complex<double>>> eigenVectors = getEigenvectors(A);
  // they are both computed as zero - can that be? -> ask sage!
  // well, the zero-vector actually is a solution - seems we compute the trivial solution
  // -> how do we compute the nontrivial solution? maybe, we need a different right-hand side?
  // correct eigenvectors are (1,1) and (1,1/2) ~= (2,1)


  // What is actually the set of matrices that gives zero when plugged into the polynomial? Are
  // these the matrices that are "similar" to A? I think, they have all the same eigenvalues 
  // because they must have the same characteristic polynomial (and the eigenvalues are its roots).
  // What about the eigenvectors - can they be different? What about the "companion matrix" - it 
  // should also be a member of the set - is this in some way a "special" member? Maybe another
  // special member is the one whose eigenvectors are the canoncial basis vectors.

  //int dummy = 0;

  // This Sage code:
  // A = matrix([[-4,6],[-3,5]])
  // E = matrix([[ 1,0],[ 0,1]])
  // R.<x>=QQ[]   # to prevent x from vanshing
  // B = A - x*E  # for computing eigenvalues x of A
  // R = B.echelon_form()
  // A,E,B,R
  //
  // produces:
  // 1  1/(3x) - 5/3
  // 0  x^2 - x - 2
  //
  // see here, why the R.<x>=QQ[] is needed:
  // https://ask.sagemath.org/question/8386/row-echelon-form-of-a-matrix-containing-symbolic-expresssions/
  // http://www.cfm.brown.edu/people/dobrush/am34/sage/echelon.html

  // Conclusion:
  // There seems to be no practical way to compute the characteristic polynomial via Gaussian 
  // elimination for a matrix of numbers. One really has to work with a matrix of rational 
  // functions in order to keep track of the mangling of the x-variable in the course of the 
  // elimination process. One cannot introduce the x-variable after the elimination :-(
  // Or is there some other way?

  // Maybe try this with a larger matrix - 3x3 or 4x4

  // Soo - how do we diagonalize a matrix, i.e. find eigenvalues and -vectors? I think, this may be 
  // possible only numerically because the eigenvalues are in general irrational (algebraic?) 
  // numbers. This is in contrast to finding an inverse or solving a linear system - in this case,
  // the solution can be computed exactly - if the matrix elements are rational numbers, the 
  // inverse matrix elements are also rational because they can be computed via Gauss elimination.
  // But what about the eigenvectors? Can't they be rational anyway? Probably not but maybe?

  /*
  // copied - maybe delete later:
  // should go to experiment:
  // check diagonalization:
  //A = Matrix(3, 3, { 6,1,9, 0,3,6, 0,0,2 });
  A = Matrix(3, 3, { -15,17,-3, 17,24,-14, -3,-14,48 });
  Matrix E(3, 3); E.setToIdentity();
  LA::makeSystemDiagonal(A, E);
  // after that E contains the inverse of A with its rows scaled by the diagonal elements that are
  // still in A - that doesn't seem to be useful - remove from library - move to experiments

  solveDiagonalSystem(A, E, E);

  // is this actually true diagonalization? nope doesn't seem to be - the results that are now in A
  // and E seem to have nothing to do with the eigenvectors and eigenvalues - but what have we done
  // what does the result represent? is it the inverse with rows scaled by the diagonal elements
  // of the resulting A? -> it's the reduced row echelon form. its determinant is the same (up to a 
  // factor) to that of the original matrix - but the characteristic polynomial has nothing to do
  // with that of the original matrix-

  // Sage code:
  // D  = matrix([[2,0,0],[0,3,0],[0,0,-5]])
  // M  = matrix([[-1,1,-5],[1,3,-1],[2,-1,1]])
  // MT = M.transpose()
  // A  = MT*D*M
  // Ai = A.inverse()
  // MT,D,M,A,Ai

  // figure out how to diagonalize a matrix using the Gaussian elemination - maybe we can figure 
  // out the coeffs of the characteristic polynomial by making it triangular (in which case the 
  // determinant is the product of the diagonal elements) - but each elementary transformation
  // modifies the determinant in some way - we need to keep track of these changes - maybe the
  // Gaussian elimination can be extended in a suitable way? -nope we are fine: swapping rows
  // inverts the sign, multiplying a row by a factor multiplies the determinant by the same factor
  // and adding a row to another doesn't change the determinant - in the end, all it does is to
  // scale our characteristic polynomial by a factor which is irrelevant for finding the roots
  // algorithm: make matrix diagonal - then we can compute the polynomial coeffs from the diagonal
  // elements...wait - oculd it be that the diagonal elements are already the roots of the
  // characteristic polynomial? after calling makeSystemDiagonal? i think so - figure it out!
  // ...hmm - it probably doesn't work because the "lambda" in det(A-lambda*E) = 0 is actually 
  // introduced before the elimination starts, so in every step, terms involving lambda get mangled
  // in - we can't just introduce the lambda after the elimination...i think - but maybe it's 
  // possible to keep track of all the mangling and arrive at an algo for finding the 
  // characteristic polynomial anyway....
  */
}
// todo: try some crazy stuff: use a matrix of rational functions of (maybe complex) rational 
// numbers made of big integers - maybe for this, it makes sense to use my own rsComplex class 
// again. Nest the different math types in variuos ways: complexes of matrices, matrices of 
// complexes, matrices of ratfuncs, ratfuncs of matrices, polynomials, vectors of matrices,
// matrices of vectors - maybe out this into a unit test

void determinant()
{
  // Computation of determinant by the (totaly inefficient) textbook method.

  using Matrix = RAPT::rsMatrix<double>;

  // todo: turn into unit test...
  bool r = true;

  Matrix A(4, 4, { 1,0,2,0, 1,0,3,0, -1,2,3,4, 2,0,5,1});
  double det;
  det = getDeterminantColumnWise(A, 0); r &= det == -2;
  det = getDeterminantColumnWise(A, 1); r &= det == -2;
  det = getDeterminantColumnWise(A, 2); r &= det == -2;
  det = getDeterminantColumnWise(A, 3); r &= det == -2;
  det = getDeterminantRowWise(   A, 0); r &= det == -2;
  det = getDeterminantRowWise(   A, 1); r &= det == -2;
  det = getDeterminantRowWise(   A, 2); r &= det == -2;
  det = getDeterminantRowWise(   A, 3); r &= det == -2;

  // ...use such a function to find the characteristic polynomial without resorting to rational 
  // functions - only class rsPolynomial itself is needed -> compare both ways. For this way, we 
  // must promote the matrix to a matrix of polynomials instead of rational functions.

  //int dummy = 0;
}




// move this code to unit tests:
template<class T>
bool testNullSpace(RAPT::rsMatrix<T> A)
{
  T tol = T(1.e-12);
  RAPT::rsMatrix<T> B = getNullSpace(A, tol); // B is basis for the nullspace
  RAPT::rsMatrix<T> null = A*B; 
  return null.isZero(tol) && !B.containsZeroColumn(tol);
}

template<class T>
bool testNullSpaceTailParams(RAPT::rsMatrix<T> A)
{
  T tol = T(1.e-12);  
  RAPT::rsMatrix<T> B = getNullSpaceTailParams(A, tol); 
  RAPT::rsMatrix<T> null = A*B; 
  return null.isZero(tol) && !B.containsZeroColumn(tol);
}

template<class T>
bool testNullSpaceBoth(RAPT::rsMatrix<T> A)
{
  return testNullSpace(A) && testNullSpaceTailParams(A);
}

bool testSubSpaces()
{
  // todo: move the getNullSpace functions into the library and turn this into unit-test
  // maybe add tests for checking eigenspaces
  bool r = true;

  using Matrix = RAPT::rsMatrix<double>;
  using LA     = RAPT::rsLinearAlgebraNew;

  // Examples from Höhere Mathematik in Rezepten (3.Aufl.), page 142 - our found nullspace basis 
  // vectors aggree with the results there only up to multiplying the basis vectors by a factor. 
  // That's ok - the basis of a vector space is not unique - we apparently have sometimes made
  // different choices for the free parameters.

  Matrix A, A2, B, B2, C, null, BB, z, x;
  int rank;
  double tol = 1.e-12;
  bool test;
   
  // Example 15.4 in Karpfinger, pg 141. When we don't manually put it in row-echelon-form, our 
  // result is different from the book - but that doesn't mean it's wrong. Apparently, in the book,
  // different swaps were performed in the row elimination process.
  A  = Matrix(3, 3, {1,1,1, 1,2,4, 2,3,5});   // original matrix
  A2 = Matrix(3, 3, {1,1,1, 0,1,3, 0,0,0});   // row echelon form of A
  B  = getRowSpace(A,  tol);                  // returns <(2,3,5),(0,0.5,1.5)>
  B2 = getRowSpace(A2, tol);                  // returns <(1,1,1),(0,1,3)> -  as in the book
  test = spanSameSpace(B, B2, tol);
  test = spanSameSpace(B.getTranspose(), B2.getTranspose(), tol);
  // Both tests evaluate to true - that's interesting - i expected only the test with the transpose
  // to return true -> figure out, how row-space and column-space are related - they have the same 
  // dimension but are they the same space?


  A  = Matrix(3, 3, {1,1,1, 1,2,4, 2,3,5});   // book says, col-space is <(1,1,2),(0,1,1)>
  A2 = Matrix(3, 3, {1,0,0, 1,1,0, 2,1,0});   // again, A in row echelon form
  B  = getColumnSpace(A,  tol);
  B2 = getColumnSpace(A2, tol);
  test = spanSameSpace(B, B2, tol);
  test = spanSameSpace(B.getTranspose(), B2.getTranspose(), tol);
  // These tests also both evaluate to true


  A = Matrix(3, 5, {1,-2,0,-1,3, 0,0,1,2,-2, 0,0,0,0,0}); r &= testNullSpace(A);  
  B = getNullSpace(A, tol); 
  // https://www.wikihow.com/Find-the-Null-Space-of-a-Matrix

  // test rank computation:
  A = Matrix(3, 3, {0,0,0, 0,0,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 0;
  A = Matrix(3, 3, {1,0,0, 0,0,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 1;
  A = Matrix(3, 3, {0,1,0, 0,0,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 1;
  A = Matrix(3, 3, {0,0,1, 0,0,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 1;
  A = Matrix(3, 3, {1,0,0, 0,1,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 2;
  A = Matrix(3, 3, {0,1,0, 0,0,1, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 2;
  A = Matrix(3, 3, {0,0,1, 0,1,0, 0,0,0}); rank = getRankRowEchelon(A, 0.0); r &= rank == 2;// can this occur when producing the ref?
  A = Matrix(3, 3, {0,0,1, 0,1,0, 0,0,2}); rank = getRankRowEchelon(A, 0.0); r &= rank == 3;

  // test nullspace computation:
  A = Matrix(2, 2, {1,0,  0,1}); r &= testNullSpaceBoth(A);
  A = Matrix(2, 2, {0,1, -1,0}); r &= testNullSpaceBoth(A);
  A = Matrix(2, 2, {0,1,  0,0}); r &= testNullSpace(A);

  A = Matrix(3, 3, {-1,-1,2, 1, 2,3, -1,0, 7}); r &= testNullSpaceBoth(A); // B = {(7,-5,1)}
  A = Matrix(3, 3, { 4, 2,2, 2, 1,1,  2,1, 1}); r &= testNullSpaceBoth(A);
  A = Matrix(3, 3, {-2, 2,2, 2,-5,1,  2,1,-5}); r &= testNullSpaceBoth(A);
  A = Matrix(3, 3, { 1, 2,3, 4, 5,6,  7,8, 9}); r &= testNullSpaceBoth(A); // B = {(1,-2,1)}
  A = Matrix(3, 3, { 1, 0,0, 0, 2,0,  0,0, 3}); r &= testNullSpaceBoth(A); // B = empty
  A = Matrix(3, 3, { 0, 1,2, 0, 0,4,  0,0, 0}); r &= testNullSpace(A); 

  A = Matrix(4, 4, {1,2,3,4, 2,4,6,8, 3,6,9,12, 4,8,12,16}); r &= testNullSpaceBoth(A); // B = {(2,-1,0,0), (3,0,-1,0), (4,0,0,-1)} // sign is different
  A = Matrix(4, 4, {1,5,6,7, 0,1,2,3, 0,0,0, 4, 0,0, 0, 0}); r &= testNullSpace(A); 
 
  A = Matrix(5, 4, {1,5,6,7, 0,1,2,3, 0,0,0,4, 0,0,0,0, 0,0,0,0}); r &= testNullSpace(A); 
 
  A = Matrix(5, 5, {1,2,3,4,5, 0,1,6,7,8, 0,0,1,9,1, 0,0,0,0,0, 0,0,0,0,0}); r &= testNullSpace(A);
  A = Matrix(5, 5, {0,0,3,4,5, 0,0,6,7,8, 0,0,1,9,1, 0,0,0,0,0, 0,0,0,0,0}); r &= testNullSpace(A);  
  A = Matrix(5, 5, {1,2,3,0,0, 0,1,6,0,0, 0,0,1,0,0, 0,0,0,0,0, 0,0,0,0,0}); r &= testNullSpace(A);  
  A = Matrix(5, 5, {1,2,0,4,5, 0,1,0,7,8, 0,0,0,9,1, 0,0,0,0,0, 0,0,0,0,0}); r &= testNullSpace(A);
  A = Matrix(5, 5, {0,2,3,4,5, 0,1,6,7,8, 0,0,1,9,1, 0,0,0,0,0, 0,0,0,0,0}); r &= testNullSpace(A); 



  // Function that takes a matrix A and its nullspace and checks, if our nullspace computation 
  // function produces a result that agrees with the target nullspace
  using Vec = std::vector<double>;
  auto checkNullSpace = [&](int M, int N, int R, Vec vecA, Vec vecB)->bool // maybe pass matrices
  { 
    Matrix A(M, N,   vecA);
    Matrix T(M, N-R, vecB);      // target nullspace
    Matrix B = getNullSpace(A, tol);
    //rsAssert(B == T);
    Matrix null = A*B; 
    C = B-T;
    return C.isZero(tol) && null.isZero(tol);
  };

  r &= checkNullSpace(3, 3, 0, {0,0,0, 0,0,0, 0,0,0}, {1,0,0, 0,1,0, 0,0,1}); // rank 0
  r &= checkNullSpace(3, 3, 1, {0,1,0, 0,0,0, 0,0,0}, {1,0, 0,0, 0,1});       // rank 1
  r &= checkNullSpace(3, 3, 1, {1,0,0, 0,0,0, 0,0,0}, {0,0, 1,0, 0,1});       // rank 1
  r &= checkNullSpace(3, 3, 2, {1,0,0, 0,1,0, 0,0,0}, {0,0,1});               // rank 2
  r &= checkNullSpace(3, 3, 2, {1,0,0, 0,0,0, 0,1,0}, {0,0,1});               // rank 2
  r &= checkNullSpace(3, 3, 3, {1,0,0, 0,1,0, 0,0,1}, {});                    // rank 3
  // create some more examples - also using 3x2 and 2x3 matrices

  // ahh - wait - the nullspace is the space that is spanned by the rows! not by the columns that 
  // wouldn't make any sense since it's the rows that live in the embedding space



  // check a couple of examples for which we have produced the correct nullspaces with sage via
  // commands like:
  // A = matrix(QQ, [[ 1,  4,  0, -1,  0,   7, -9],
  //                 [ 2,  8, -1,  3,  9, -13,  7],
  //                 [ 0,  0,  2, -3, -4,  12, -8],
  //                 [-1, -4,  2,  4,  8, -31, 37]])
  // A.right_kernel()
  //
  // which gives as answer:
  //
  // Vector space of degree 7 and dimension 4 over Rational Field
  // Basis matrix:
  // [    1     0     0     0  -3/7  -1/7     0]
  // [    0     1     0     0 -12/7  -4/7     0]
  // [    0     0     1     0  -3/7 -9/14  -1/2]
  // [    0     0     0     1   1/7 13/28   1/4]
  A = Matrix(4, 7, { 1,  4,  0, -1,  0,   7, -9,
                     2,  8, -1,  3,  9, -13,  7,
                     0,  0,  2, -3, -4,  12, -8,
                    -1, -4,  2,  4,  8, -31, 37});
  B  = getNullSpace(A, tol);
  B2 = Matrix(4, 7, {1, 0, 0, 0,  -3./7, -1./7,  0,
                     0, 1, 0, 0, -12./7, -4./7,  0,
                     0, 0, 1, 0,  -3./7, -9./14,-1./2,
                     0, 0, 0, 1,   1./7, 13./28, 1./4 }).getTranspose();
  test = spanSameSpace(B, B2, tol);
  // why do we need to transpose? what convention uses sage to return the nullspace? are the rows 
  // shown there supposed to be columns? maybe we should also transpose our result to be consistent 
  // with sage? ...figure out which convention is used in math for representing baes of nullspaces
  // and adapt code, if necessarry


  B = getOrthogonalComplement(A.getTranspose(), tol);
  C = A*B;
  // C is zero - that's good - but A^T has 4 cols, B also - shouldn't B have only 3? the number of 
  // dimensions of both bases should add up to the dimension of the embedding space. Or maybe A has
  // one vector that is superfluous in it? maybe we should have a function 
  // removeSuperfluousColumns ..or maybe take A^T -> turn into row-echelon -> 
  // remove bottom zero lines -> transpose - that should get rid of supefluous vectors


  Matrix rsp = getRowSpace( A, tol);
  Matrix nsp = getNullSpace(A, tol);
  Matrix cmp = getOrthogonalComplement(rsp.getTranspose(), tol); 
  r &= spanSameSpace(cmp, nsp, tol);
  // we need to transpose because a rowspace consists of rows - is that the common convention?
  // figure out! probabyl we may use these relations to figure out if getNullSpace has actually 
  // produced the whole embedding space - we may create the union of the nullspace and its
  // complement and check, if it spans the whole space

  // here:
  // https://en.wikipedia.org/wiki/Row_and_column_spaces#Dimension it says:
  // The dimension of the column space is called the rank of the matrix. 
  // we may check that

  //rsEigenSpaceSet<complex<double>> eigen = getEigenSpaces(A, tol);
  // we should notneed to wrap it into complex here - getEigenSpaces should automatically 
  // complexify for eal matrices - or it should take complex matrices - yes - that would perhaps
  // be the best


  // todo:
  // add a test that the returned matrix does not contain the zero vector - figure out, how it 
  // comes about that the produced nullspaces sometimes contain it
  // try various MxN matrices - maybe with random elements ...although, for random matrices, the 
  // nullspaces will probably always come out as the whole embedding space - it's unlikely that
  // random matrices are singular - maybe programmatically construct singular matrices (use a set
  // of random vectors and produce random linear combinations of them)
  // how can we check that we actually get the whole nullspace and not just a subspace thereof?
  // maybe the only way is to use matrices with knowm nullspaces - again we can construct matrices
  // randomly from a known set of basis vectors - if we write some random basis vectors into the 
  // rows and construct the other rows as linear combinations of the them, we have a known 
  // nullspace - we can check with spanSameSpace if we got the full nullspace
  // try various "mean" matrices (with zero columns in various places

  // todo: orthogonalize the obtained nullspace bases and compute BB = B.getTranspose() * B; 
  // this should give the identity matrix, iff B is orthonormal - maybe make such testt later when we have
  // Gram-Schmidt orthogonalizations

  if(r)
    std::cout << "nullspace works";
  rsAssert(r);
  return r;
}

RAPT::rsMatrix<double> getRandomMatrix(int numRows, int numCols, int seed,
  double min = -1, double max = +1)
{
  RAPT::rsMatrix<double> A(numRows, numCols);
  RAPT::rsArrayTools::fillWithRandomValues(A.getDataPointer(), A.getSize(), min, max, seed);
  return A;
}

RAPT::rsMatrix<double> getRandomOrthogonalMatrix(int size, int seed)
{
  RAPT::rsMatrix<double> A = getRandomMatrix(size, size, seed);
  orthonormalizeColumns1(A); 
  return A;
}


bool testSigularValueDecomp()
{
  using Vec = std::vector<double>;
  using Matrix = RAPT::rsMatrix<double>;
  double tol = 1.e-12;
  //double tol = 1.e-8;
  bool r = true;




  
  // this does not yet work:
  Matrix A, U, S, V;  // the original matrices
  Matrix a, u, s, v;  // the computed matrices

  int M = 7;
  int N = 5;
  V = getRandomOrthogonalMatrix(N, 2);
  U = getRandomOrthogonalMatrix(M, 1);
  S.setShape(7, 5); S.setToZero();
  S(0,0) = 5; S(1,1) = S(2,2) = S(3,3) = 3; S(4,4) = 2; // the middle sv mas multiplicity 3
  A = U * S * V.getTranspose();
  decomposeRealUSV(A, u, s, v, tol);
  //decomposeRealUSV(A, u, s, v, 1.e-6); // higher tolerance doesn't help
  a = u * s * v.getTranspose();
  // doesn't work - characteristic polynomial computation fails - maybe we need higher tolerance?
  // ...could be that the reduction of the rational functions doesn't work - increasing tol 
  // doesn't seem to help - oh - wait the rsRationalFunction class doesn't get its tolerance from 
  // here ...hmm rsRationalFunction<T>::reduce is not even called
  // but ratReduce is called - with tol == 0 - the operators of RatFunc use ratMul, ratAdd, etc. - 
  // they do not pass anything for the tolerance, so it defaults to zero - maybe RatFunc should
  // have a member reductionTolerance
  
  // ..maybe try to use matrices with integer elements - but how to we produce
  // orthogonal matrices with integer elements - if we use random values and Gram-Schmidt, the 
  // elements will often be irrational due to the sqrt in the normalization
  





  auto checkSVD = [&](int M, int N, Vec vecA)->bool
  { 
    Matrix A(M, N, vecA);
    Matrix U, S, V;
    decomposeRealUSV(A, U, S, V, tol);
    Matrix T = U * S * V.getTranspose();
    bool result = (A-T).isZero(tol);
    result &= isOrthogonal(U, tol);
    result &= isOrthogonal(V, tol);
    //result &= isDiagonal(  S, tol);
    return result;
  };

  // examples/excercises from Karpfinger - Höhere Mathematik in Rezepten:
  r &= checkSVD(2, 3, {-1,1,0, -1,-1, 1});              // pg. 448.
  r &= checkSVD(2, 3, { 1,1,3,  1, 1,-3});              // pg. 450, ex 42.3 (a)
  r &= checkSVD(3, 1, { 2,2,1 });                       // (b)
  r &= checkSVD(3, 2, { 1,1, 1,1, 3,-3});               // (c)
  r &= checkSVD(1, 3, { 2,2,1 });                       // (d)
  r &= checkSVD(3, 4, {8,-4,0,0, -1,-7,0,0, 0,0,1,-1}); // (e)
  r &= checkSVD(2, 3, { 2,1,2,  1,-2,1});               // (f) 



  r &= checkSVD(2, 4, {1,0,1,0, 0,1,0,1});
  // https://mysite.science.uottawa.ca/phofstra/MAT2342/SVDproblems.pdf - has multiplicity
  // ....also uses A * A^T ...why? how?


  // for the matrices where the number of rows is larger than the number of columns, the assertion 
  // triggers but the test returns true anyway - maybe it doesn't actually matter, if we fill up
  // U? ...i think, whether we fill up U or not, the condition A = U * S * V^T will always hold - 
  // but if we don't fill U up, it won't be orthogonal (not even invertible) - which may be 
  // disadvantageous, if we want to do some manipulation, such as left-multiplying and equation by
  // U^-1 (= U^T, iff U is orthogonal)
  // try more examples with m > n, m < n and different multiplicities for the eigenvalues of 
  // A^T * A, i.e. cases where there is no one-to-one correspondence between eigenvalues and 
  // eigenvectors
  // n: dimensionality of input space
  // m: dimensionality of output space
  // r: dimensionality of image, rank of A, r <= m







  // see also:
  // https://en.wikipedia.org/wiki/Singular_value_decomposition#Example
  // https://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm
  // https://www.d.umn.edu/~mhampton/m4326svd_example.pdf


  // maybe construct examples by specifying U,S,V - when specifying S, make sure to cover cases 
  // where we have singular values with multiplicities




  return r;
}


void linearCombinations()
{
  // We figure out, if a given vector is a linear combination of a set of given vectors, see:
  // https://www.mathbootcamps.com/determine-vector-linear-combination-vectors/
  // https://math.stackexchange.com/questions/851766/determining-whether-or-not-a-vector-is-a-linear-combination-of-a-given-matrix


  using Matrix = RAPT::rsMatrix<double>;
  using LA     = RAPT::rsLinearAlgebraNew;

  Matrix B, x;
  bool linComb;
  double tol = 1.e-13;

  //      b1 b2
  //     |2  7|                    |25|
  // B = |5  3|, x = 2*b1 + 3*b2 = |19|
  //     |3  5|                    |15|
  B = Matrix(3, 2, {2,7, 5,3, 3,5});  
  x = Matrix(3, 1, {25,19,21});
  // our basis vectors (2,5,3),(7,3,5) spanning a 2D space within R^3, as colums of matrix B, x is
  // a linear combination of b1,b2 with coeffs 2,3

  linComb = isInSpanOf(B, x, tol);
  //LA::makeTriangular(B, x);
  // x is a linear combination of the columns of B, if this results in no nonzero line in x for 
  // which there is a zero line in B - or the other way around: whereever B has a zero line, x must
  // also have a zero entry. I think, this means the rank(B|x) <= rank(B)

  // try the same procedure now with some other vector x:
  B = Matrix(3, 2, {2,7, 5,3, 3,5});  // bcs the elimination has destryed it
  x = Matrix(3, 1, {20,15,-10});
  linComb = isInSpanOf(B, x, tol);
  //LA::makeTriangular(B, x);
  // this is not e linear combination of b1,b2

  // make some more linear combinations and random vectors
  B = Matrix(3, 2, {2,7, 5,3, 3,5});
  x = Matrix(3, 4, {25,20,20,3, 19,21,15,4, 21,19,-10,1});
  // the 2nd is 3*b1 + 2*b2, the first as above and 3,4 are random vectors
  linComb = isInSpanOf(B, x, tol);


  // how can we find the coefficients of the linear combinations - that must also be a linear 
  // system -> figure out, how to set it up

  //int dummy = 0;
}

void linearIndependence()
{
  // Tests, whether a set of vectors (given as columns of a matrix) is linearly independent. A set
  // of vectors a1,a2,..,aN is linearly independent, if k1*a1 + k2*a2 + ... + kN*aN = 0 implies 
  // that k1 = k2 = ... = kN = 0. That means, we have to find the nullspace of the matrix 
  // A = (a1 a2 ... aN), where the ai are the columns. The vectors are independent if that 
  // nullspace is trivial (i.e. contains only the zero vector).


  linearCombinations();

  using Matrix = RAPT::rsMatrix<double>;
  using LA     = RAPT::rsLinearAlgebraNew;

  double tol = 1.e-12;

  //      a1 a2 a3     we have a3 = a1-a2, so the columns are not linearly independent
  //     |1  3  -2|
  // A = |2  1   1|
  //     |1  2  -1|
  Matrix A(3, 3, {1,3,-2, 2,1,1, 1,2,-1});
  Matrix nullspace = getNullSpaceTailParams(A, tol);
  // A basis for the nullspace is {(-1,1,1)} - that means: -1*a1 + 1*a2 + a3 = 0 which can be 
  // simplified to a3 = a1-a2 - our free parameter is a3 - we may choose a1,a2 in any way that 
  // satifies this equation ...i think

  // A = |1 2 1|
  //     |0 0 1|
  A = Matrix(2, 3, {1,2,1, 0,0,1});
  nullspace = getNullSpaceTailParams(A, tol);  
  // Returns (-2,1,0) - that means: -2*a1 + 1*a2 + 0*a3 = 0 - we can choose a3 freely, the 
  // equation is always satisfied. Or we can solve this equation - for example - for a1:
  // -2*a1 = -1*a2 - 0*a3 -> a1 = (1/2)*a2 - 0*a3 - so we have expressed a1 as linear combination
  // of the other two, so a1 is linearly dependent of a2,a3. However, it does not mean that every
  // vector in the set is dependent on some others - a3 cannot be expressed as linear combination
  // of a1,a2, for example ...is that because it has coefficient zero?

  //int dummy = 0;
}


void orthogonalizedPowerIteration()
{
  // Experiments with an algorithm that i came up with that is a simple extension of the power 
  // iteration method to compute the largest eigenvalue and it's corresponding eigenvector. The 
  // algo is extended to find more (potentially all) eigenvalues and -vectors by forcing the 
  // iterates to be orthogonal to the subspace spanned by all the already found eigenvectors. It 
  // does so by doing a Gram-Schmidt like orthogonalization step after the multiplication of the
  // iterate by the matrix. It turns out (as far as i can tell, so far), that the resulting algo 
  // is indeed able to find all eigenvalues but the eigenvectors can only be found if they happen
  // to be orthogonal (which is the case for eigenvectors corresponding to distinct eigenvalues of
  // a symmetric matrix). For the time being, i would suggest to call the algorithm "orthogonalized
  // power iteration" (OPI).

  using Real = double;
  using Mat  = rsMatrix<Real>;
  using Vec  = std::vector<Real>;
  using ILA  = rsIterativeLinearAlgebra;
  using AT   = rsArrayTools;

  // Create 2x2 matrix with eigenvalues 3,2 and (orthogonal) eigenvectors (6,8),(-8,6):
  int N = 2;                      // dimensionality
  Vec s({3.0,  2.0});             // eigenvalues* (distinct, positive real)
  Vec v({6,8, -8,6});             // eigenvectors (orthogonal, both have norm 100)
  Mat A = fromEigenSystem(s, v);  // A = (2.36, 0.48,  0.48, 2.64) (symmetric)
  int maxIts = 500;               // maximum number of iterations
  // (*) I use s for the eigenvalues to indicate that they are scalars ("eigenscalars")

  // Output variables for algo:
  Vec wrk(N), t(N), w(N*N);       // wrk: workspace, t: recovered eigenvalues, w: orthovectors*
  int its;                        // number of iteratiosn taken by the algo
  bool ok = true;
  // (*) "orthovector" is a term i just made up myself to call the orthogonalized eigenvectors that
  // the algo produces.

  // Retrieve eigensystem:
  AT::fillWithRandomValues(&w[0], N*N, -1.0, +1.0, 0);             // initial guess
  its = ILA::eigenspace(A, &t[0], &w[0], 1.e-13, &wrk[0], maxIts); // its = 81
  ok  = checkEigensystem(t, w, s, v, 1.e-12);                      // yep, works!

  // Now try the same thing with eigenvectors (.6,.8),(.8,.6). These have both have norm 1 and are
  // not orthogonal anymore:
  v = Vec({.6,.8, .8,.6});                                         // not orthogonal anymore
  A    = fromEigenSystem(s, v);                                    // A is now antisymmetric
  AT::fillWithRandomValues(&w[0], N*N, -1.0, +1.0, 0);             // initial guess
  its = ILA::eigenspace(A, &t[0], &w[0], 1.e-13, &wrk[0], maxIts); // its = 75
  ok  = rsIsPermutation(s, t, 1.e-12);                             // eigenvalues are correct
  ok  = checkEigensystem(t, w, s, v, 1.e-12);                      // eigenvectors are wrong

  // Test if the found vectors are orthonormal (they should be):
  Mat W = rsToMatrixColumnWise(w, N, N);
  ok = areColumnsOrthogonal( W, 1.e-12);                   // yes
  ok = areColumnsOrthonormal(W, 1.e-12);                   // yes

  // Test if the found vectors are compatible with the Gram-Schmidt orthogonalized set of the 
  // actual, true eigenvectors:
  Mat V   = rsToMatrixColumnWise(v, N, N);
  Mat V_o = V; orthonormalizeColumns1(V_o);       // V_o is (Gram-Schmidt) orthogonalized V
  Mat D   = W - V_o;
  // D is the zero matrix, so the computed vectors are indeed equal to the Gram-Schmidt 
  // orthogonalized set of the original eigenvectors. I think, in general, we would expect them to
  // equal to but only compatible with them, i.e. they may have different signs. (...something to
  // figure out...)

  // Now try to recover the orginal eigenvector v2 from w1,w2 by undoing the Gram-schmidt process
  // ...but how can this be done? Can it be done at all?
  Vec v1(N); V.copyColumn(0, &v1[0]); // todo: implement convenience function: v1 = V.getColumn(0)
  Vec v2(N); V.copyColumn(1, &v2[0]);
  Vec w1(N); W.copyColumn(0, &w1[0]);
  Vec w2(N); W.copyColumn(1, &w2[0]);
  Vec Aw2 = A*w2;
  Vec u2p = w1 + Aw2;  // nah - does not equal v2
  Vec u2m = w1 - Aw2;  // ditto
  // The Gram-Schmidt process applied to v1,v2 does (verify!):
  //
  //   v1_o = v1 / |v1|
  //
  //            (v2 - <v1_o, v2> * v1_o)
  //   v2_o = ----------------------------
  //           |(v2 - <v1_o, v2> * v1_o)|
  //
  // Can this be undone? The normalization steps obviously can't be undone because the length 
  // information is indeed lost in such a step. But we are not interested in reconstructing the 
  // length information anyway. Only the direction information is relevant. Can we hope to 
  // reconstruct from v1_o, v2_o a vector v2_n that is a normalized version of v2 not necessarily
  // orthogonal to v1_o but instead pointing in the same direction as v2 originally did? Maybe not
  // without additional information. But we do have additional information: the matrix A, its
  // eigenvalues s1,s2 and we know that v1 and v1_o are eigenvectors to s1 and we know that
  // v2_o is the perpendicular part of v2 with respect to v1 (or v1_o). Maybe writing:
  //   v2_n = par(v2_n, {v1}) + perp(v2_n, {v1})
  // could help? Here, par and perp denote the parallel and perpendicular part of a vector with 
  // respect to a set of vectors (here, a 1 element set). The term perp(v2_n, {v1}) has the same 
  // direction as v2_o (but possibly different length?)
  // See:
  // https://www.physicsforums.com/threads/matrix-which-reverses-gram-schmidt-linear-algebra.981879/
  // https://arxiv.org/ftp/arxiv/papers/1607/1607.04759.pdf - this apaper discusses an inversion
  // algo for Gram-Schmidt, but that makes use of an additional matrix r(M,N) that is computed 
  // during the forward orthogonalization along with the orthogonalized set of vectors. We don't
  // have such a matrix available here. Maybe define the (unknown) denominator as:
  //   k2 := |(v2 - <v1_o, v2> * v1_o)|. 
  // Then we may write 
  //   v2_o*k2 = v2 - <v1_o, v2> * v1_o
  // here v2_o,v1_o are the knwons and k2,v2 are the unknowns.
  // I think, in this special case here, we may be able to reconstruct the eigenvector by flipping
  // the sign of x or y. Try it! but this may work only when A is antisymmetric? ...maybe if it 
  // works, the algo can be applied to symmetric and antisymmetic parts (A_s, A_a) separately? But 
  // that would only help, if we could find the eigenvectors of A from those of A_s, A_a

  // Other idea:
  // Form a vector:
  //   u2 := w1 + w2
  // and multiply it by the matrix A:
  //   p2 : = A*u2
  // then find the projections of p2 onto w1 and w2:
  //   q12 := <w1, p2>
  //   q22 := <w2, p2>
  // and their ratio:
  //   r12 := q12/q22
  // i think, this ratio should be the same as the ratio of the corresponding eigenvalues s1/s2? 
  // Let's try it:
  Vec  u2  = w1 + w2;
  Vec  p2  = A*u2;
  Real q12 = rsDot(w1, p2);  // -0.42857142857126318
  Real q22 = rsDot(w2, p2);  //  2 == 2nd eigenvalue
  Real r12 = q12/q22;        // -0.21428571428564130
  // ...hmmm - nope! but can it somehow help us to reconstruct the projection coeff of v2 onto v1?
  // -q12 and therfore r12 grows with growing largest eigenvalue s1
  // -todo: 
  //  -plot the function r12(s1) for a given s2
  //  -plot also the dependency on the angle bewteen the eigenvectors - i expect r12 to be +- s1/s2
  //   when the angle is +-90°. make a function that takes a reference vector and angle and two 
  //   eigenvalues that creates the plot

  auto getProjectionRatio = [&](Real angle, Real s1, Real s2)
  {
    // Function to compute the ratio r12 (as above) as function of the angle between the 
    // eigenvectors
    Real c = cos(angle);
    Real s = sin(angle);
    Vec e({s1,  s2 });              // eigenvalues
    Vec v({1,0, c,s});              // eigenvectors
    Mat A = fromEigenSystem(e, v);  // our matrix
    Vec wrk(2), t(2), w(2*2);       // wrk: workspace, t: recovered eigenvalues, w: orthovectors
    int its;                        // number of iterations taken by the algo
    AT::fillWithRandomValues(&w[0], 2*2, -1.0, +1.0, 0);     // initial guess
    its = ILA::eigenspace(A, &t[0], &w[0], 1.e-8, &wrk[0], maxIts);  // that's a high tolerance!
    Mat W = rsToMatrixColumnWise(w, 2, 2);
    Vec w1(2); W.copyColumn(0, &w1[0]);
    Vec w2(2); W.copyColumn(1, &w2[0]);
    Vec  u2  = w1 + w2;             // or maybe s1*w1 + s2*w2?
    //Vec  u2  = s1*w1 + s2*w2;       // test
    Vec  p2  = A*u2;
    Real q12 = rsDot(w1, p2);
    Real q22 = rsDot(w2, p2);       // always(?) equals the smaller eigenvalue
    Real r12 = q12/q22;
    Real r21 = q22/q12;
    //return r21;
    //return r12;
    return q12;
  };

  // Test:
  //Real Q12 = rsDot(v1, p2);  // -0.42857142857126318
  //Real Q22 = rsDot(v2, p2);  //  0.14857142857163841
  //Real tmp = acos(q12);      //  2.0137073708683526
  //tmp = cos(q12);            //  0.90956035167423543

  // This is our target value, that we want to find - the projection of the 2nd eigenvector onto 
  // the 1st. If we know that projection coeff, we should be able to reconstriuct v2 from the 
  // available information (i hope):
  Real t12 = rsDot(v1, v2);  // 0.96


  Real minAngle = 0.05;  // avoid singularity at 0
  Real maxAngle = 3.10;  // avoid singularity at pi
  int numAngles = 500;
  Vec angles = rsLinearRangeVector(numAngles, minAngle, maxAngle);
  Vec ratios(numAngles);
  Vec cot_p(numAngles);  // positive cotangent
  Vec cot_m(numAngles);  // negative cotangent
  for(int i = 0; i < numAngles; i++)
  {
    ratios[i]  = getProjectionRatio(angles[i], s[0], s[1]);
    ratios[i] -= s[0];
    cot_m[i]   = tan(angles[i] - PI/2);
    cot_p[i]   = tan(PI/2 - angles[i]);
  }
  rsPlotVectorsXY(angles, ratios, cot_p, cot_m); 
  // The black curve (representing q12 shifted by 3 which is the larger eigenvalue) agrees 
  // partially with the positive and partially with the negative cotangent function. ...but this
  // seems to work only, if the eigenvalues are 3 and 2 - WTF?




  int dummy = 0;


  // Now try it with a 3x3 matrix with orthogonal eigenvectors (6,8,0),(-8,6,0),(0,0,10):
  // ...

  // Observations:
  // -For the symmetric 2x2 matrix A = (2.36, 0.48,  0.48, 2.64) with eigenvalues 3,2 and 
  //  orthogonal eigenvectors (6,8),(-8,6), the algorithm can indeed recover both eigenpairs from
  //  the matrix A. The algo takes 81 iterations.
  // -When we use the non-orthogonal eigenvectors (6,8),(8,6), the matrix A becomes antisymmetric.
  //  The algo will still produce the correct eigenvalues, but the returned vectors are not the
  //  eigenvectors. The algo takes 75 iterations - that's a bit less than before.

  // Conclusion (preliminary):
  // -I think, the set of vectors that the algorithm produces is forced to be orthogonal and 
  //  therefore it can match the actual eigenvectors only when they are orthogonal themselves.
  // -I think in any step, (i.e. for any n), the set of the n already found vectors spans the same
  //  space as the first n eigenvectors. Indeed, i think the set of vectors is the Gram-Schmidt
  //  orthogonalized set of the first n eigenvectors.
  // -Maybe the algorithm could at least be useful for symmetric matrices, because their 
  //  eigenvectors are indeed orthogonal. Oh - no - that holds only for eigenvectors corresponding 
  //  to distinct eigenvalues:
  //  https://math.stackexchange.com/questions/82467/eigenvectors-of-real-symmetric-matrices-are-orthogonal

  // ToDo:
  // -Try to apply Gram-Schmidt othogonalization to the actual, true eigenvectors and see if the
  //  result matches the output of the algo.
  // -Try to reconstruct the actual eigenvectors from the returned orthonormal set by a 
  //  postprocessing step. Maybe in the 2x2 case try w2' = w1 + A*w2 or w2' = w1 - A*w2 where
  //  w2' is supposed to be the reconstructed eigenvector v2 and w1 is the first (found) 
  //  eigenvector. I think the computed w2 is perp(v2; w1), i.e. the perpendicular componenent of
  //  the actual eigenvector v2 with respect to the 1st eigenvector w1 (which itself is compatible
  //  with v1, i.e. equals v1 up to scaling)
  // -If this doesn't work out, we may use the known eigenvalues s_i to compute the corresponding 
  //  eigenvectors v_i in a 2nd step, one at a time, by solving the linear system 
  //  (A - s_i*I) * v_i = 0. If this system is also to be solved iteratively, maybe it makes sense
  //  to use the computed "orthovector" as initial guess? Maybe it would make sense to have a 
  //  function shiftedProduct that can be used instead of the regular call to product in the 
  //  solvers. shiftedProduct(const rsMatrix<T>& A, const T& s, const T* x, T* y) should work
  //  like product(const rsMatrix<T>& A, const T* x, T* y) but instead use the shifted matrix
  //  (not yet sure, if the shift should be applied with positive or negative sign). 
  //  To do this, we should first move the solvers from rsSparseMatrix into 
  //  rsIterativeLinearAlgebra
  // -Rayleigh iteration may also be a candidate to find the eigenvectors
}

void eigenstuff()
{
  orthogonalizedPowerIteration();


  // create a matrix and find its eigenvalues and eigenvectors

  using Matrix  = RAPT::rsMatrix<double>;
  using MatrixC = RAPT::rsMatrix<complex<double>>;

  double tol = 1.e-12;
  bool r = true;
  Matrix  A, T, z;                       // matrix and dummy vector (get rid of the dummy)
  MatrixC B;                             // correct basis of the eigenspace (target)
  std::vector<rsEigenSpace<double>> eig; // hold the computed eigenspaces

  A = Matrix(2,2, {0,1, 0,0});
  Matrix nullspace = getNullSpace(A, tol);
  // returns |1 0| - the canonical basis of R^2 - but this is wrong (i think)
  //         |0 1|
  // maybe it fails because the rank is not computed correctly

  z = Matrix(2, 1); // dummy
  int rank = RAPT::rsLinearAlgebraNew::makeTriangular(A, z);
  rank = getRankRowEchelon(A, 0.0);
  // yup - rank is returned as 0 - but actual rank is 1 - the number of steps taken is actually not
  // the rank in all cases - often it is, but not always - try this with various matrices with 
  // leading zeros in the rows - these make problems - the actual rank is the numer of nonzero rows
  // in row echelon form - after i steps, we have zeroed out i columns - but what does this tell us 
  // about the rows? not so much - hmm - here:
  // https://en.wikipedia.org/wiki/Rank_(linear_algebra)
  // it says: "Once in row echelon form, the rank is clearly the same for both row rank and column 
  // rank, and equals the number of pivots" ...hmm - maybe it's because we use partial pivoting?

  // examples from https://www.youtube.com/watch?v=lyXwcXjJdYM
  // todo: add checks (r &= ...)

  A = Matrix(2,2, {1,1, 0,1});
  eig = getEigenSpaces(A, tol);
  //findEigenSpacesReal(A);
  // correct [1,(1,0)] - or maybe weitz didn't compute the other?, found: [1,{(1,0),(0,1)}] twice
  // we can actually check, if an eigenvector v is indeed eigenvector - just compute 
  // A * v and x_i * v

  A = Matrix(2,2, {1,1, 1,1});
  eig = getEigenSpaces(A, tol);
  //findEigenSpacesReal(A);
  // found: [0,(-1,1)],[2,(1,1)] -> correct


  A = Matrix(2,2, {0,-1, 1,0});
  eig = getEigenSpaces(A, tol);
  //findEigenSpacesReal(A);
  // found: [-i,{(-i,1)}],[i,{(i,1)}] - todo: check this manually

  // Examples from Ahrens,pg.658:
  A = Matrix(2,2, {3,-1, 1,1});
  eig = getEigenSpaces(A, tol);
  //findEigenSpacesReal(A);
  // [2,(1,1)], [2,(1,1)] -> correct

  A = Matrix(3,3, {1,2,2, 2,-2,1, 2,1,-2});
  //eig = getEigenSpaces(A, tol);    // doesn't work - we need higher tolerance
  eig = getEigenSpaces(A, 1.e-7);  // works
  //findEigenSpacesReal(A);
  // [3,(2,1,1)],[-3,(-1,0,2),(-1,2,0)] - (2,1,1) is found correctly, the other eigenspaces are 
  // empty - numerical issues?, also, the -3 eigenvalue is found twice - once as -3.000...x and 
  // once as -2.999...x - in each case with an empty eigenspace - it smells like a precision issue
  // ...with tol = 1.e-7, they are collapsed into one - with a 2D eigenspace - this seems correct
  // but why are the eigenspaces empty in case of too small tolerance - shouldn't they just be 
  // duplicated, too?
  // maybe revisit the polynomial root-finder, maybe the root-polishing can be refined - maybe
  // using Newton iteration? but what if it's a double (or multiple) root? - in this case, Newton
  // iteration may not converge - maybe somehow work with deflated polynomials? the characteristic
  // polynomial here is 27 + 9*x - 3*x^2 - 1*x^3, with a simple root at +3 and a double root at -3

  // Example from Ahrens,pg.659 - has a single eigenvalue of -2 (with multiplicity 5) with a 2D 
  // eigenspace spanned by {(1,1,1,1,1),(0,0,1,0,0)}:
  A = Matrix( 5, 5, {-3,1,0,0,0, -1,-1,0,0,0, -3,1,-2,1,1, -2,1,0,-2,1, -1,1,0,0,-2});
  B = MatrixC(5, 2, {1,0, 1,0, 1,1, 1,0, 1,0});  // correct basis for eigenspace
  eig = getEigenSpaces(A, tol);
  r &= eig.size() == 1;
  r &= eig[0].eigenValue == -2.0;
  r &= eig[0].getAlgebraicMultiplicity() == 5;
  r &= eig[0].getGeometricMultiplicity() == 2;
  r &= spanSameSpace(eig[0].eigenSpace, B, complex<double>(tol)); // get rid of complexification of tol
  // we get different basis vectors from what the book says - but that doesn't mean, they are 
  // wrong - they are just a different basis that spans the same space


  // todo: try some more examples with sage, clean up the eigenvalues - use a function like
  // cleanUpIntegers that does it for real and imaginary parts


  // try orthonormalization - todo: make extra function for this
  A = Matrix( 3, 3, {1,1,1, 0,1,1, 0,0,1});       // Karpf. pg160
  orthonormalizeColumns1(A);
  r &= A == Matrix(3, 3, {1,0,0, 0,1,0, 0,0,1});  // should give the standard-basis
  r &= areColumnsOrthonormal(A, tol);

  // other algorithm:
  A = Matrix( 3, 3, {1,1,1, 0,1,1, 0,0,1}); 
  orthonormalizeColumns2(A, tol);
  r &= A == Matrix(3, 3, {1,0,0, 0,1,0, 0,0,1});


  // example from:
  // https://www.khanacademy.org/math/linear-algebra/alternate-bases/orthonormal-basis/v/linear-algebra-gram-schmidt-example-with-3-basis-vectors
  A = Matrix( 4, 3, {0,0,1, 0,1,1, 1,1,0, 1,0,0});
  T = A;
  orthonormalizeColumns1(A);
  r &= areColumnsOrthonormal(A, tol);
  r &= spanSameSpace(A, T, tol);
  double s1 = 1. / sqrt(2), s2 = sqrt(2./3.), s3 = 1./(2*sqrt(3.));
  T = Matrix( 4, 3, {0,0,3*s3, 0,s2,s3, s1,0.5*s2,-1*s3, s1,-0.5*s2,s3});   // target
  A = T - A;          // A is now the error (target - computed) - should be the zero-matrix
  r &= A.isZero(tol);
  // we could also make tests base on spanSameSpace and areColumnsOrthonormal

  // QR decomposition:
  A = Matrix(4, 3, {2,0,2, 1,0,0, 0,2,-1, 2,0,0}); // Karpf. pg.185
  Matrix Q, R;
  decomposeQR(A, Q, R); 
  T = Q*R;                 // should be equal to A
  r &= (A-T).isZero(tol);



  //int dummy = 0;

  // Some Notes
  // -x_i is an eigenvalue of matrix A, iff (A - x_i * I) * v = 0 has solutions v different from 
  //  the zero vector, the equation can also be written as A * v = x_i * v
  // -a matrix A is invertible, iff 0 is not an eigenvalue (Ahrens, 655)
  // -the dimensionality of an eigenspace to an eigenvalue x_i has a dimensionality of at least
  //  1 and at most the multiplicity (as root of the characteristic polynomial) of x_i

  // ToDo:
  // -Figure out use cases and implement 
  //  https://en.wikipedia.org/wiki/Proper_orthogonal_decomposition
  //  ...Google says, it's the same as PCA (principle component analysis) but for time-series data
  //  as opposed to random vector data. Figure out relation to SVD (singular value decomposition). 
  //  See also:
  //  https://www.youtube.com/watch?v=axfUYYNd-4Y
  //  https://www.youtube.com/watch?v=OhyksL-1vew
  // -Maybe apply these techniques to results of simualtions of the 2D wave-equation and/or maybe
  //  Navier-Stokes
};




template<class T>
int rsSolveRichardson2(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b, T alpha,
  T tolR, int maxIts, T smooth = T(0))
{
  // This version uses 1st order FIR smoothing of the error.
  int its = 0;
  int N = (int) x.size();
  std::vector<T> err, errOld(N), dx;
  smooth *= 0.5;
  while(its < maxIts) 
  {
    err = A*x - b;
    //dx  = -alpha * err;          // this is the basic Richardson iteration rule
    dx  = -alpha * ((1-smooth)*err + smooth*errOld);
    if(rsStaysFixed(x, dx, T(1), tolR)) 
      return its;  // x has converged
    x = x + dx;
    errOld = err;
    its++; 
  }
  return its;
  // maybe we should smooth only when the errors err and errOld have different signs...or maybe 
  // reduce the alpha? ...and increase it, when they have the same size. maybe the decrease should
  // large (like a factor of 0.5) but the increase smaller (like 1.1)
}
template<class T>
int rsSolveRichardson3(const rsMatrix<T>& A, std::vector<T>& x, const std::vector<T>& b, T alpha0,
  T tolR, int maxIts, T grow, T shrink, T mom = T(0))
{
  // This version uses an adaptive update rate per coordinate that grows and shrinks according to 
  // the signs of the errors in two successive iterations. If the error (for some coordinate i) has
  // different signs in successive iterations, it means that we have overshot the optimum and 
  // therefore, the update step was too large, so the update rate for coordinate i is reduced by 
  // the shrink factor which should be <= 1. On the other hand, if the error has the same sign, the
  // update rate is multiplied by the grow factor >= 1.
  int its = 0;
  int N = (int) x.size();
  std::vector<T> err, errOld(N), alpha(N), dx(N);
  for(int i = 0; i < N; i++)       // initialize update rates
    alpha[i] = alpha0;
  while(its < maxIts) 
  {
    err = A*x - b;

    for(int i = 0; i < N; i++) {
      if(err[i]*errOld[i] < T(0))          // adapt update rate...
      {
        alpha[i] *= shrink;                // ...shrink, if error has alternating sign

        // maybe take back previous step for this coordinate?
        //x[i] -= dx[i]; // nope - not good - slows down convergence

      }
      else
        alpha[i] *= grow;                  // ...grow, if not
      dx[i] = mom*dx[i] + (T(1)-mom)*(-alpha[i]*err[i]);  // ...do the actual update
    }

    if(rsStaysFixed(x, dx, T(1), tolR))    // x has converged
      return its;

    x = x + dx;
    errOld = err;
    its++; 
  }
  return its;
}


void iterativeLinearSolvers()
{
  using Real = double;
  using Mat  = rsMatrix<Real>;
  using Vec  = std::vector<Real>;
  //using ILA  = rsIterativeLinearAlgebra;
  //using AT   = rsArrayTools;


  bool ok = true;

  int maxIts = 1000;

  int N = 3;
  Mat A(3, 3, {5,-1,2, -1,7,3, 2,3,6}); // is symmetric and positive definite (SPD) (verify!)
  Vec x({1,2,3});
  Vec b = A*x;
  Vec x2(N);
  int its = rsSolveCG(A, x2, b, 1.e-12, maxIts);
  ok &= rsIsCloseTo(x, x2, 1.e-12);
  ok &= its == 3; // conjugate gradient is supposed to find the solution after at most N steps

  rsFill(x2, 0.0);
  its = rsSolveRichardson(A, x2, b, 0.16, 1.e-13, maxIts); // around 0.16 seems best, 75 iterations
  Real err = rsMaxDeviation(x2, x);
  ok &= err <= 1.e-12;

  rsFill(x2, 0.0); 
  its = rsSolveRichardson2(A, x2, b, 0.16, 1.e-13, maxIts, 0.5);         // 66
  err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;

  // no momentum, only grow/shrink:
  rsFill(x2, 0.0); 
  its = rsSolveRichardson3(A, x2, b, 0.16, 1.e-13, maxIts, 1.25, 0.75);  // 51
  err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;

  // no grow/shrink, only momentum:
  rsFill(x2, 0.0); 
  its = rsSolveRichardson3(A, x2, b, 0.16, 1.e-13, maxIts, 1.0, 1.0, 0.25); // 58
  err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;

  // momentum and grow/shrink:
  rsFill(x2, 0.0); 
  its = rsSolveRichardson3(A, x2, b, 0.16, 1.e-13, maxIts, 1.25, 0.75, 0.125); // 45
  err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;

  // todo: 
  // -try adaptive momentum per coordinate
  // -try to combine smoothing with grow/shrink


  // Now, let's use a non-symmetric matrix:
  A = Mat(3, 3, {-1,5,2, 2,3,7, 6,3,2});
  b = A*x;

  // Try Richardson iteration:
  rsFill(x2, 0.0); 
  its = rsSolveRichardson3(A, x2, b, 0.1, 1.e-13, maxIts, 1.25, 0.75, 0.125);
  //err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;
  // result is infinite -> algo has diverged

  // Try conjugate gradient algo:
  rsFill(x2, 0.0); 
  its = rsSolveCG(A, x2, b, 1.e-13, maxIts);
  // its = maxIts, x2 is wrong -> algo has failed

  // Try least squares conjugate gradient algo:
  rsFill(x2, 0.0); 
  its = rsSolveLSCG(A, x2, b, 1.e-13, maxIts);
  err = rsMaxDeviation(x2, x); 
  ok &= err <= 1.e-12 && its == 3;

  // Can we do some sort of iterative improvement to polish the solution?
  its = rsSolveLSCG(A, x2, b, 0.0, maxIts);
  err = rsMaxDeviation(x2, x); 
  // ...hmm - this takes 30 iterations and the result is still inexact. ToDo: look up the 
  // iterative improvement in Numerical Recipies and see, if it can be applied here


  // soon obsolete:
  // Idea: The conjugate gradient algorithm needs a symmetric and positive definite matrix. Let's
  // construct a related problem, that features such a matrix and has the same solution as A*x = b.
  // Consider the equation: dot(A*x-b, A*x-b) = (A*x-b)^T * (A*x-b) = min. Multiplying it out, 
  // taking the derivative and setting it to zero gives the equation: P*x = q, where P := A^T*A,
  // q = A^T*b. So, what we do now instead of solving A*x = b is:
  //   Define: P := A^T*A, q := A^T*b
  //   Solve:  P*x = q
  Mat A_T = A.getTranspose();
  Mat P   = A_T * A;
  Vec q   = A_T * b;
  rsFill(x2, 0.0); 
  its = rsSolveCG(P, x2, q, 1.e-13, maxIts);
  err = rsMaxDeviation(x2, x); ok &= err <= 1.e-12;
  Vec t = P*x2; // test - should equal q - seems ok
  // ...yes, that seems to work indeed. I think, if the problem has no solution, the algo will 
  // produce a least squares approximation to a solution and if it has multiple solutions, it will
  // produce the minimum norm solution (ToDo: verify and document). For sparse matrices A, it does 
  // not make sense to explicitly create the matrix A^T*A because it may not be so sparse anymore. 
  // Its number of nonzero elements can in worst case be the product of the numbers of nonzero 
  // elements of the two factor matrices (I think - verify that!).
  
  
  // Instead, we need an adapted CG algorithm
  // that computes the product (A^T * A) * x by first computing y = A*x and then computing 
  // z = A^T*y where in both products, the sparsity can be exploited. So the algo will need to 
  // compute twice as many matrix-vector products as the original CG algo. More specifically, it 
  // will have to compute 2 matrix-vector products per iteration. But that is more than compensated 
  // for by the fast convergence. The vector q = A^T * b on the other hand can be precomputed. See:
  // https://en.wikipedia.org/wiki/Gramian_matrix
  // https://math.stackexchange.com/questions/158219/is-a-matrix-multiplied-with-its-transpose-something-special
  // https://www.quora.com/Whats-the-meaning-of-matrixs-transpose-multiplied-by-the-matrix-itself
  // The algo can be derived more simply by just premultiplying both sides of A*x = b by 
  // A^T

  // ToDo:
  // -implement a function transProduct which computes A^T * x for rsMatrix, rsSparseMatrix

  int dummy = 0;
}


void linearSolverPrecision()
{
  // We construct matrices of different sizes by shuffling diagonal matrices for which we know
  // the exact solution. Then we compare the exact solution to the computed solution and find the
  // maximum error. When the numbers that occur in the computation cannot be represented exactly 
  // (which we may ensure by scaling all matrix entries by pi, say), the solver will introduce 
  // roundoff error. We plot this error as function of the size. This may help to find sensible
  // values for the singularity detection threshold in the solver.
  // This is not yet finished


  using Matrix = RAPT::rsMatrix<double>;
  using LA     = RAPT::rsLinearAlgebraNew;

  int N = 100;  // size

  Matrix A(N, N), B(N, N), X(N, N);
  A.setToIdentity();   // coefficient matrix
  B.setToIdentity();   // solution matrix

  // ...in this case, the solution is expected to be the identity matrix

  shuffle(A, B, N, 1);
  plotMatrix(A, true);
  shuffle(A, B, N, 2);
  plotMatrix(A, true);
  shuffle(A, B, N, 3);
  plotMatrix(A, true);
  shuffle(A, B, N, 4);
  plotMatrix(A, true);

  // this doesn't look good - the shuffling algo is still bad
  // the first pass of shuffling leads to L-shaped regions of similar values, the second pass to 
  // blocks - why

  LA::solve(A, X, B);

  // ok - but we can see the numerical error increase with the size

  //int dummy = 0;
}






void ellipseLineIntersections()
{
  // create and set up ellipse:
  //rsEllipseF ellipse;
  rsEllipse<float> ellipse;
  ellipse.setParameters(1.5f, 3.f, float(PI/4), 0.5f, 1.f); 

  // create line parameters:
  float x, y, dx, dy;
  x  = 0.5f;
  y  = 0.5f;
  dx = 0.5f;
  dy = 0.2f;

  // create data for drawing the ellipse:
  static const int Ne = 500;
  float xe[Ne], ye[Ne];
  float err;
  for(int n = 0; n < Ne; n++)
  {
    float phi = float(2*PI*n) / (Ne-1);
    ellipse.getPointOnEllipse(phi, &xe[n], &ye[n]);
    err = ellipse.evaluate(xe[n], ye[n]);
  }

  // create data for drawing the line:
  float tMin = -3.5;
  float tMax = +3.0;
  float xl[2], yl[2];
  xl[0] = x + tMin*dx;
  xl[1] = x + tMax*dx;
  yl[0] = y + tMin*dy;
  yl[1] = y + tMax*dy;

  // find intersection points between line and ellipse:
  float ti1, ti2, xi1, yi1, xi2, yi2;
  ellipse.lineIntersectionParameter(x, dx, y, dy, &ti1, &ti2);
  xi1 = x + ti1*dx;
  yi1 = y + ti1*dy;
  xi2 = x + ti2*dx;
  yi2 = y + ti2*dy;

  // find tangent line to intersection point:
  float A, B, C, a, b;
  float xt1[2], yt1[2], xt2[2], yt2[2];
  ellipse.getTangentCoeffs(xi1, yi1, &A, &B, &C); // implicit  A*x + B*y + C = 0
  a = -A/B;     // explicit y = a*x + b
  b = -C/B;
  xt1[0] = -1.2f;   // use points where x=0, x=2 to draw the tangent:
  yt1[0] = a*xt1[0] + b;   
  xt1[1] = -0.5f;
  yt1[1] = a*xt1[1] + b; 

  // same for the 2nd intersection:
  ellipse.getTangentCoeffs(xi2, yi2, &A, &B, &C);
  a = -A/B;
  b = -C/B;
  xt2[0] = 0;
  yt2[0] = a*xt2[0] + b;   
  xt2[1] = 2;
  yt2[1] = a*xt2[1] + b; 

  GNUPlotter plt;
  plt.setRange(-1.5, +2.5, -1, +3);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");   // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(Ne, xe,  ye);   // ellipse
  plt.addDataArrays(2,  xl,  yl);   // line
  plt.addDataArrays(2,  xt1, yt1);  // tangent at 1st intersection
  plt.addDataArrays(2,  xt2, yt2);  // tangent at 2nd intersection
  plt.plot();
}

// move after void numericDiffAndInt();
void iteratedNumDiff()
{
  // under construction

  // todo: test taking numeric derivatives of numeric derivatives to get 2nd derivatives - maybe go 
  // up to the 4th derivative. figure out, what choices of the stepsize give the most accurate 
  // results - maybe the outer derivatives need smaller or larger stepsizes than the inner ones?
  // plot error between numeric and analytic result in some range
  // example functions: sin, exp, log, 1/x, x^a, 1/(1+x^2), 

  static const int N = 500;
  double xMin = -10;
  double xMax = +10;

  double x[N], y[N];
  double y1[N], y2[N], y3[N], y4[N];  // computed derivatives
  double t1[N], t2[N], t3[N], t4[N];  // target derivatives
  double e1[N], e2[N], e3[N], e4[N];  // errors

  rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);


  // Set up approximation stepsizes - in order to have the x+h, x+2h inputs to be machine numbers, 
  // we use a stepsize h that is an (inverse) power of 2 - this should get rid of at least one 
  // source of inaccuracy (right? ...verify conditions under which this makes sense)
  //double h1 = pow(2, -15);      // 1.5 e-10, sine
  //double h1 = pow(2, -16);    // 4.e-11, noisy sine
  //double h1 = pow(2, -17);  // 1.5e-11, more noisy sine
  //double h1 = pow(2, -17.5);  // 
  double h1 = pow(2, -18);    // 1.5e-11, noise with a little bit of sine
  //double h1 = pow(2, -19);     // 3.e-11, noise with some modulation
  //double h1 = pow(2, -21); 

                            // obtained using h1 = 2^-(-18):
  double h2 = pow(2, -11);        // 6.0e-8, noisy sine
  //double h2 = pow(2, -12);      // 6.0e-8, modulated noise
  //double h2 = pow(2, -13);      // 1.0e-7, noise
  //double h2 = pow(2, -14);      // 2.0e-7, noise
  //double h2 = pow(2, -15);      // 3.5e-7, noise
  //double h2 = pow(2, -16);      // 8.0e-7, noise
  //double h2 = pow(2, -17);      // 1.5e-6, noise
  //double h2 = pow(2, -18);      // 3.0e-6, noise
  //double h2 = pow(2, -19);      // 4.0e-6, noise


  h1 = pow(2, -17); h2 = pow(2, -12);  // 3.0e-8
  //h1 = pow(2, -12); h2 = pow(2, -17);  // 3.0e-8


  //h1 = pow(2, -16); h2 = pow(2, -12);  // 2.0e-8
  //h1 = pow(2, -16); h2 = pow(2, -13);  // 2.0e-8
  //h1 = pow(2, -15); h2 = pow(2, -12);    // 1.5e-8, noisy sine
  //h1 = pow(2, -14); h2 = pow(2, -12);    // 1.2e-8, less noisy sine
  //h1 = pow(2, -13); h2 = pow(2, -12);     // 1.3e-8, even less noisy sine
  //h1 = pow(2, -13); h2 = pow(2, -13);     // 8.0e-9
  //h1 = pow(2, -12); h2 = pow(2, -13);   // 1.3e-8, same as using h1 = 2^(-13), h2 = 2^(-12) - 
                                          // the swap makes no difference - is this generally true?

  double h3 = pow(2, -15);
  double h4 = pow(2, -15);

  // function f0(x) := f(x) = sin(x) and its first 4 derivatives:
  auto f0 = [=](double x)->double{ return  sin(x); };
  auto f1 = [=](double x)->double{ return  cos(x); };
  auto f2 = [=](double x)->double{ return -sin(x); };
  auto f3 = [=](double x)->double{ return -cos(x); };
  auto f4 = [=](double x)->double{ return  sin(x); };

  // g0 = f0 and then g1,... are numerical derivatives:
  using Func = std::function<double(double)>;  // function from double to double
  Func g0 = [=](double x)->double{ return  sin(x); };
  Func g1 = RAPT::rsDerivative(g0, h1);
  Func g2 = RAPT::rsDerivative(g1, h2);
  Func g3 = RAPT::rsDerivative(g2, h3);
  Func g4 = RAPT::rsDerivative(g3, h4);

  // compute functions values and numerical and analytical derivatives:
  for(int n = 0; n < N; n++)
  {
    y[n]  = f0(x[n]);

    y1[n] = g1(x[n]);
    y2[n] = g2(x[n]);
    y3[n] = g3(x[n]);
    y4[n] = g4(x[n]);

    t1[n] = f1(x[n]);
    t2[n] = f2(x[n]);
    t3[n] = f3(x[n]);
    t4[n] = f4(x[n]);

    e1[n] = y1[n] - t1[n];
    e2[n] = y2[n] - t2[n];
    e3[n] = y3[n] - t3[n];
    e4[n] = y4[n] - t4[n];
  }

  //using  Vec = std::vector<double>;


  //rsPlotArraysXY(N, x, y, t1, t2, t3);  // needs more inputs


  //rsPlotArraysXY(N, x, e1);

  rsPlotArraysXY(N, x, e2);

  //rsPlotArraysXY(N, x, e1, e2);

  //rsPlotArraysXY(N, x, e1, e2, e3, e4); // e4 is large

  // Observations:
  // -for the sine function, the optimal stepsize h for the first derivative seems to be 2^(-18)
  //  -in this case, the error is noise-like with some sort of sinusoidal modulation with a maximum
  //   error of around 1.5e-11
  //  -going lower, e.g. h = 2^(-19), the character of the error stays modulated-noise-like but 
  //   with higher amplitude
  //  -going higher, e.g. h = 2^(-17), the error looks like a noisy sine-wave, where at 2^-(15), 
  //   the sinusoidal shape clearly dominates - there's not much noise anymore at this setting and
  //   the maximum amplitude of the error also increases
  //  -using 2^(-17.5) actually further reduces the maximum error to 1.2e-11, 2^(-17.4) or 
  //   2^(-17.6) is really strange
  // -if h1 is chosen to be h1 = 2^(-18), i.e. the optimal value, then the optimal choice for h2 
  //  seems to be h2 = 2^(-12), giving a noisy error with a maximum around 6.e-8
  // -using h1 = 2^(-13), h2 = 2^(-12) seems to give the same error as h1 = 2^(-12), h2 = 2^(-13) 
  //  so swapping h1 and h2 seems to make no difference -> figure out, if this generally true
  //  -> at least with -12 and -17, the swap also makes no difference - maybe derive the final
  //  formula analytically and see, if it's symmetrical in h1 and h2 -> done (below): yes, it is

  // f=f0: function, f1: 1st numeric derivative, f2: 2nd numeric derivative:
  //   f1(x) = (f0(x+h1) + f0(x-h1)) / (2*h1)
  //   f2(x) = (f1(x+h2) + f1(x-h2)) / (2*h2)
  //         = (f(x+h2+h1) - f(x+h2-h1) - f(x-h2+h1) + f(x-h2-h1)) / (4*h1*h2)
  // defining:
  //   hs = h1+h2, hd = h1-h2  (sum and difference)
  // i arrived at:
  //   f2(x) = (f(x+hs) + f(x-hs) - f(x+hd) - f(x-hd)) / (4*h1*h2)
  // this is a formula which is indeed symmetric in h1,h2. It's interesting to note that it uses 4
  // evaluation points - the 2nd derivative via 5-point stencil would use 5 (and yes, the center 
  // coeff would be nonzero in this case) - so chaining two 1st order central differences gives not
  // the same formula as using a 2nd order central difference directly, even when h1==h2. The 2nd 
  // order, 5-point formula is:
  //   f''(x) ~= (-f(x-2h) + 16*f(x-h) - 30*f(x) + 16*f(x+h) -f(x+2h) ) / (12*h^2)
  // with our chained 1st central differences, we get for h1==h2:
  //   f''(x) ~= (f(x-2h) + f(x-2h) - 2*f(x)) / (4*h^2)
  // which is actually just the 2nd order accurate 3-point stencil formula with stepsize 2h. This 
  // means, when h1==h2, we get effectively only a 3-point formula, even though we do 4 function
  // evaluations. What when h1 != h2? Do we get better accuracy in this case? This may be relevant
  // for how we compute the higher order partial derivatives in rsManifold...but for partial 
  // derivatives, it may be inconvenient to use direct 2nd order formulas - how would that work 
  // anyway? Would the final division be replaced by a matrix inversion?

  // when using this 4-point stencil: -2,-1,1,2 for the 2nd derivative here:
  // http://web.media.mit.edu/~crtaylor/calculator.html
  // we get this:
  //   f''(x) ~= (f(x-2h) - f(x-h) - f(x+h) + f(x+2h)) / (3*h^2)
  // which is yet another formula...but maybe it becomes the same when h2 = 1.5h, h1 = 0.5*h, yes:
  //   f2(x) = (f(x+h2+h1) - f(x+h2-h1) - f(x-h2+h1) + f(x-h2-h1)) / (4*h1*h2)
  //         = (f(x+2h)    - f(x+h)     - f(x-h)     + f(x-2h))    / (4*1.5*h*0.5*h)


  // todo: 
  // -figure out the optimal stepsize for the 2nd derivative, given the one for the first has 
  //  been chosen optimally, i.e h1 = 2^(-18)
  // -try to tweak the stepsize for the 1st derivative - maybe values that are suboptimal for 
  //  computing the 1st derivative as such could be better, wehn the 1st derivative is computed 
  //  only as intermediate result for the 2nd derivative? (well - that would be strange, but who
  //  knows - floating point arithmetic sometimes works in mysterious ways)
  //  -> plot the maximum error of the 2nd derivative as 2D function of h1 and h2


  // (preliminary) conclusions:
  // -for computing the 2nd derivative as central difference of first central differences, the 
  //  optimal stepsize h1 in the inner approximation may be different from the optimal stepsize
  //  h1 used for computing the 1st derivative itself. It seems weird, that we should use a 
  //  suboptimal inner derivative computation - maybe inaccuracies in the two neighbouring 
  //  1st derivatives somehow cancel in the computaion of the 2nd derivative?
}

/*
double getMax2ndDerivativeErrorSin(double h1, double h2,
  double xMin = -10, double xMax = +10, int numSamples = 500)
{

  return 0;  // preliminary
}
*/

//-------------------------------------------------------------------------------------------------





/** Map composed of two maps of the form given by rsRationalMap, one applied to the domain 
0..0.5 and the other to 0.5..1 where inputs and outputs are appropriately scaled and shifted. The 
maps are combined in such a way to give a range of shapes between sigmoidal (s-shaped) and a 
saddle-ish, i.e. shape similar to a cubic like x^3, i.e. flatter in te middle. */
template<class T>
T rsBiRationalMap_01(T x, T a)
{
  if(x < 0.5)
    return T(0.5) * rsRationalMap_01(T(2)*x, a);
  else
    return T(0.5) * rsRationalMap_01((x-T(0.5))*T(2), -a) + T(0.5);
}
// Maybe rename this to rsBiMoebiusMap_01 and rsRationalMap_01 to rsMoebiusMap_01
// Or rsSplitMoebiusMap or rsSymmetrizedMoebiusMap, rsMirroredMoebiusMap, rsReflectedMoebiusMap,
// rsSigmoidizedMoebiusMap, rsSplitLinFracMap

/** Map composed of a rational map, a birational map and another rational map. ...TBC.. */
template<class T>
T rsTetraRationalMap_01(T x, T a, T b, T c)
{
  x = rsRationalMap_01(  x, a);  // pre convexity or concavity
  x = rsBiRationalMap_01(x, b);  // sigmoidity or saddleness
  x = rsRationalMap_01(  x, c);  // post convexity or concavity
  return x;
}
// Maybe rename this to rsComposedMoebiusMap_01. Maybe implement a variant that uses the birational
// map as outer functions and the normal rational map as inner function.
// Maybe intead of tetra use tri. The split-moebius map should perhaps not count as two because it
// just switches between two maps rather than applying two maps in sequence.
// The full transformation is a sequence of:
//   Moebius -> affine -> Moebius -> affine -> Moebius
// where in the x <= 0.5 case, the affine map is actually even linear (i.e. no offset is involved). 
// But wait, the switch between the cases happens not on the original x but on the intermediate 
// value after the first map has been applied. But be that as it may, I think, in both cases, the
// full transformation is again Moebius map because they form a closed set (and affine maps are 
// also in the set - they are very simple Moebius maps). So it seems, it would be OK to call the 
// full map also just a Moebius map. But, it's actually a segmented map built from two Moebius 
// maps stitched together. Maybe call it rsBiMoebiusMap
// https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
// https://de.wikipedia.org/wiki/M%C3%B6biustransformation


bool moebiusMapTest()  // rename to linFracTest
{
  bool ok  = true;

  using Real = double;
  using Vec  = std::vector<Real>;

  // Setup for the tests:
  int  N   = 257;     // Number of samples to generate between 0 and 1
  Real tol = 1.e-13;  // Tolerance for floating point rounding errors

  // Coefficients:
  Real a = 0.5, b = -0.6, c = -0.3; // moderately flat at0, very flat at 1
  //a = 0.5;  b = 0.75; c = -0.25; // very steep at 0, moderately steep at 1
  //a = 0.75; b = 0.75; c = -0.5;  

  // Allcoate some arrays to write the functions into:
  Vec x = rsLinearRangeVector(N, 0, 1);
  Vec y(N), z(N);
  Vec err;                          // Error between expected and actual

  // Check inversion of the rational map via negating the paraneter:
  for(int n = 0; n < N; n++)
  {
    y[n] = rsRationalMap_01(x[n],  a);
    z[n] = rsRationalMap_01(y[n], -a);  // Should give back x
  }
  err = z-x; ok &= rsIsAllZeros(err, tol);
  //rsPlotVectorsXY(x, y, z);

  // Check inversion of the birational map via negating the paraneter:
  for(int n = 0; n < N; n++)
  {
    y[n] = rsBiRationalMap_01(x[n],  a);
    z[n] = rsBiRationalMap_01(y[n], -a);  // Should give back x
  }
  err = z-x; ok &= rsIsAllZeros(err, tol);
  //rsPlotVectorsXY(x, y, z);

  // Check inversion of the tetrarational map via negating the paraneters and reversing their 
  // order. The tetrarational map is what we want to use for the invertible interpolation scheme,
  // so if that works, we are good.
  for(int n = 0; n < N; n++)
  {
    y[n] = rsTetraRationalMap_01(x[n],  a,  b,  c);
    z[n] = rsTetraRationalMap_01(y[n], -c, -b, -a);  // Should give back x
  }
  err = z-x; ok &= rsIsAllZeros(err, tol);
  //rsPlotVectorsXY(x, y, z);

  // Test combining two rational functions into a single one. This is just to verify the formula
  // for c.
  c = (a + b) / (a*b + 1); // Parameter of the resulting function
  for(int n = 0; n < N; n++)
  {
    // Compute y by applying two rational maps in sequence:
    y[n] = rsRationalMap_01(x[n], a);
    y[n] = rsRationalMap_01(y[n], b);

    // Compute z by applying a single rational map:
    z[n] = rsRationalMap_01(x[n], c);  // should be equal to y[n]
  }
  err = z-y; ok &= rsIsAllZeros(err, tol);
  //rsPlotVectorsXY(x, y, z);
  // This same way of combination works also for the birational maps. 
  // ToDo: write code that verifies this! Extract the formula for c into a function or at least
  // document the formula in some easier to find place.

  // Compare rational and birational map for the same a:
  for(int n = 0; n < N; n++)
  {
    y[n] = rsRationalMap_01(  x[n], a);
    z[n] = rsBiRationalMap_01(x[n], a);
  }
  //rsPlotVectorsXY(x, y, z); // have same slope at x = 0 and reciprocal slopes at x = 1
  // This is just a plot for inspection. No programmatic unit test to be made here.


  // Now amalgamate the whole compose map into a single piecewise Moebius map:
  Real splitX = rsRationalMap_01(0.5, -a);    // point at which to switch between lower and upper piece
  Real test   = rsRationalMap_01(splitX, a);  // should be 0.5
  ok &= rsIsCloseTo(test, 0.5, tol);
  Real ab  = a*b, ac = a*c, bc = b*c;
  Real abc = ab*c;
  Real A1  = -(abc + ab + ac + bc + a + b + c + 1);
  Real B1  = 0;
  Real C1  = -2*(ab - bc + a + 2*b + c);
  Real D1  = abc - ab - ac - bc + a + b + c - 1;
  Real A2  = -2*(abc + ab + ac - 3*(bc + b) + a + c + 1);
  Real B2  = 4*(abc + ab - bc - b);
  Real C2  = - 4*(ab - bc + a - 2*b + c);
  Real D2  = 2*(abc + 3*(ab - b) - ac - bc + a + c - 1);
  for(int n = 0; n < N; n++) {
    y[n] = rsTetraRationalMap_01(x[n],  a,  b,  c);
    if(x[n] <= splitX)
      z[n] = (A1*x[n] + B1) / (C1*x[n] + D1);
    else
      z[n] = (A2*x[n] + B2) / (C2*x[n] + D2);   }
  err = z-y; ok &= rsIsAllZeros(err, tol);
  //rsPlotVectorsXY(x, y, z);


  // ToDo:
  // -Maybe instead of implementing these big formulas for A1, B1, ...etc., write a general 
  //  function to compute coeffs of a composed Moebius trafo and call that. This might be simpler
  //  and more efficient. I even have code for that somewhere already. It's in
  //  rsMoebiusTransform<T>::followedBy. Maybe implement that as a static function. But maybe the 
  //  way we do it here is actually not so inefficient. The code above does 13 multiplications 
  //  (and a whole lot of additions). Usingthe code from rsMoebiusTransform<T>::followedBy, we 
  //  would require 8 multiplications for each composition of 2 Moebius transforms. The code for
  //  that is actually the same as for 2x2 matrix multiplication. For the lower piece of the map,
  //  we would have to call that twice (giving 16 muls already) and then we would have to do the 
  //  same for the upper piece (given yet another 16 muls). In between we would have to do a few 
  //  more muls to take care of the affine trafos pre and post the inner map. So, even though it 
  //  looks ugly, maybe the code above is actually not so bad.
  // -Wrap the amalgamation/composition code into a function that can be used as library function
  //  and test it for some more cases, i.e. more choices for a,b,c. Maybe let a,b,c each run 
  //  through -0.9...+0.9 in 0.1 steps. That should give 19^3 = 6859 test cases. Maybe use steps 
  //  of 0.3 giving 7^3 = 343 test cases.
  // -Maybe for this purpose here, it is inconvenient to parameterize the linfrac map via the
  //  parameter a. Maybe we should parametrize it directly via its slope parameter. Maybe that can 
  //  get rid of all the complicated equations for A1, B1, ... and replace them by simpler 
  //  formulas. On the other hand, it's convenient that the inverse map just requires negation
  //  of the parameter. But: when using the slope, inversion would just reuqire to take the 
  //  reciprocal which would be equally convenient and intuitive. Maybe the conversion between 
  //  slope and parameter is generally useful for conveting parameter ranges back and forth between
  //  (-1,0,+1) <--> (0,1,inf) where in the middle, we have the neutral value. ...DONE!!...

  return ok;
}

void linearFractionalInterpolation()
{
  bool ok = moebiusMapTest();  // Remnant from the development process
  rsAssert(ok);
  // A function that performs some unit tests on some stuff used in the old implementation. The 
  // whole function may be deleted someday. But it contains some code implementing some rather 
  // complicated formulas that allow to combine 3 maps into 1 when they are parametrized via the 
  // old p in [-1,+1]. ...But should we ever need such a thing, it's perhaps better to 
  // reparametrize -> combine -> reparametrize back. But this will require two divisions, so I'm 
  // not sure yet. That's why the code is still there.

  // Experiments with linear fractional interpolation - an interpolation method that I have 
  // constructed around the so called linear fractional transformations. It is, as far as I know,
  // new. If not, then I have yet again re-invented the wheel. We create a plot with several graphs 
  // showing the normalized interpolants for a fixed choice of the slope at the origin an a bunch 
  // of exponentially spaced slopes at 1,1. These normalized interpolants with adjustable slopes
  // at the endpoints can later be combined with a numerical differentiation algo (to produce
  // target values for such slopes) and a bit of back-and-forth de/normalization to construct an
  // interpolation scheme. That is still to be done - but that's supposedly grunt work because
  // it will work exactly analogous to the cubic Hermite interpolation scheme. The hard math work 
  // is done and the results are presented here.

  using Real = double;
  using Vec  = std::vector<Real>;
  using LFI  = rsLinearFractionalInterpolator<Real>;

  // User parameters for the plots:
  int  N           = 257;       // Number of samples
  Real shape       = 0.0;       // 0.0: symmetric (default)
  Real slopeAt0    = 1.0/1.0;   // Slope of all graphs at x,y = 0,0
  Real minSlopeAt1 = 1.0/128.0; // Minimum slope at x,y = 1,1
  Real maxSlopeAt1 = 128.0;     // Maximum slope at x,y = 1,1
  //Real tol         = 1.e-13;    // Tolerance for some tests that we do along the way
  //bool ok          = true;      // Flag to indicate a failed test

  Vec x = rsLinearRangeVector(N, 0, 1);
  Vec y(N);
  GNUPlotter plt;
  //setToDarkMode(&plt);
  plt.setToDarkMode();
  Real slopeAt1 = minSlopeAt1;           // Variable slope at x,y = 1,1. Goes up in the loop
  while(slopeAt1 <= maxSlopeAt1) {
    for(int n = 0; n < N; n++)
      y[n] = LFI::getNormalizedY(x[n], slopeAt0, slopeAt1, shape);
    plt.addDataArrays(N, &x[0], &y[0]);
    slopeAt1 *= 2;  }                    // Double the slope for the next graph
  plt.addCommand("set size square");
  plt.setPixelSize(600, 600);
  plt.plot();

  // Observations:
  // -When shape = 0.5, the interpolants of some slopes K and 1/K at 1 have corresponding 
  //  interpolants that are symmetric to the y = x line, i.e. are inverses of each other.
  // -The shapes look beautifully onion shaped or maybe like christmas tree decoration. Maybe they
  //  could be useful for graphics applications.

  // Notes:
  // -Check literature about rational splines. We are doing something similar here, I think.
  //  https://www.alglib.net/interpolation/rational.php  ...but it's not quite the same thing.
  // -My initial idea was to use a linear combination of self-inverse functions for interpolation.
  //  However, the linear combination of self-inverse functions does not lead to a function that is
  //  of the same type. The forward combination stacks the functions on top of each other along the
  //  y-axis and the backward combination would stack them to the right of each other along the
  //  x-axis. This is not the same thing, so the approach was a dead end road. What we actually 
  //  need is a set of functions, all of whose inverses are also members of the set. Our set of
  //  functions should have the structure of a group in the sense of abstract algebra.
  // -When thinking about cubic splines vs cubic Hermite interpolation, we actually note that 
  //  Hermite with numerical derivatives may make more sense than spline with respect to how well
  //  the data is modeled because we also model the *measured* derivative whereas in spline 
  //  interpolation, we *adjust* the derivatives to make the spline smooth to 2nd order. But if the
  //  underlying data does not happen to be described by a cubic polynomial, this smoothness 
  //  constraint may be meaningless with respect to data modeling. Maybe try what happens when we 
  //  use numeric derivatives that are at least 3rd order accurate in cubic Hermite interpolation. 
  //  That accuracy seems to be a good fit. When the input data actually *is* a cubic polynomial, 
  //  the Hermite interpolant should be able to reconstruct it exactly
  // -For a pdf documentation, it would make sense to produce 3 plots with slopeAt0 = 1/4, 1, 4. 
  //  The plot for 4 contains graphs which inverted versions of those in the plot for 1/4.
  // -Maybe for a given set of coeffs a,b,c, plot the 3 individual maps together with the combined 
  //  map to develop intuitions about what happens in the nesting process. Or maybe make a diagram 
  //  with 7 plots and arrows between them like this:
  //
  //            A          B              C
  //            |          |              |
  //            v          v              v
  //    x ---> A(x) ---> B(A(x)) ---> C(B(A(x)))
  //
  //  where x plots the identity function, A and A(x) plot the first map (controlled by a), B plots
  //  second map applied to x, B(A(x)) plots the second map applied to A(x), C plots the third map
  //  applied to x, C(B(A(x))) plot the third map applied to B(A(x)) which is the final output.
  //  Maybe scrap x and A because A is redundant with A(x) and x is just the identity. But maybe it
  //  is clearer when the are included in the diagram nonetheless.
  // -Plot the cubic Hermite interpolant for comparison.
  // -The shape controls, how the graphs "fan out" at 0,0. With lower values, the fan is denser at
  //  the top-left and with higher values, it's denser at the bottom-right. With 0.5, the graphs 
  //  fan out 
  //  symmetrically and evenly. Well...at least, that holds for s0 = 1. If we set s0 = 4, then 
  //  shape = 0.75 leads to a more evenly spaced fan out.
  //  It seems like the shape parameter can even go beyond the range 0..1.
  //  2 or -1 produce rather extreme results but they still seem to satisfy the slope conditions at
  //  the endpoints. Maybe if we provide this parameter to the user, we shopuld rescale it from 
  //  0..1 to -1..+1 such that the user gets a symmetric shape when the parameter is 0. ...but then 
  //  check, if inverse interpolation can use the negative shape. If not, then maybe this rescaling 
  //  is not a good idea. We'll see...

  // -ToDo:
  // -Explain, how this can be used for invertible interpolation - how to obtain the inverse
  //  interpolant etc. I think, it will just use the exact same algorithm just with slope data 
  //  reciprocated (due to inverse function rule of differentiation). If the shape parameter is 
  //  used (i.e. != 0.5), I think for the inverse, we need to use 1 - shape
  // -Use this map for 1st order smooth monotonic and invertible interpolation for 
  //  data arrays x[n], y[n], s[n] containing the x,y values and desired slopes. The API should
  //  be like: interpolateMoebius(const T* x, const T* y, const T* s, N, x0). It should take the 
  //  data arrays of length N and a desired interpolation point x0 and return the interpolated y0 
  //  value at x0. Make the API consistent with other interpolation functions. Note that for the 
  //  inverse interpolation, the s array should be reciprocated because of the inverse function 
  //  rule (I think). See https://en.wikipedia.org/wiki/Inverse_function_rule. Doing:
  //    y0   = interpolateMoebius(x, y, s,  N, x0);
  //    x0_r = interpolateMoebius(y, x, sr, N, y0);
  //  should be an identity operation, i.e. we should have x0_r == x0 up to roundoff error. Here, 
  //  sr are the reciprocals of s and_x0_r is the reconstructed x0.
  // -Use that interpolation scheme in rsTimeWarper (let the user switch between linear and 
  //  rational/Moebius)
  // -Maybe target values for the slopes can be obtained from the data by numerical dervatives or
  //  maybe we can apply a constraint that the curvatures should match at the nodes similar to
  //  what is done in cubic spline interpolation. But maybe matching to 2nd order is impossible?
  //  Dunno - figure out!

  // See also:
  // https://math.stackexchange.com/questions/4162828/interpolation-with-exact-inverse
  // https://en.wikipedia.org/wiki/Linear_fractional_transformation
  // https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
  // https://en.wikipedia.org/wiki/Laguerre_transformations
  // https://de.wikipedia.org/wiki/M%C3%B6biustransformation

  // Questions:
  // Is there a 2D and nD generalization of that group? Like, a function R^2 -> R^2 given by
  //   (z,w) = f(x,y) = ( (a*x+b*y+c)/(d*x+e*y+f), (g*x+h*y+i)/(j*x+k*y+l) )
  // I think, for a fixed y, each of these 2 functions just boil down to a 1D linfrac - the 
  // y-dependent terms just become constants. How does this relate to linfracs C -> C? Is it the 
  // same thing? Could 2D linfracs be useful for 2D interpolation? In general, for R^n -> R^n:
  //   f(x) = ( (a_1*x + b_1)/(c_1*x + d_1), ..., (a_n*x + b_n)/(c_n*x + d_n)
  // where a_i,c_i are vectors, x is a vector, b_i,d_i are scalars and * is the dot product. 
  // Is this a group? I think, the neutral element should have a_i = 1 in the i-th component 
  // function and all other coeffs zero. But what about inverses? Try to find them with Sage in 
  // the 2D. Solve both equations for x and for y in terms of z,w

  // Other ideas (spin offs - move to some experiment involving waveshaping):
  // -The function -cbrt(1-x^3) might be interesting for waveshaping. It does something strange
  //  around zero but away from zero, it is the identity.
  // -Implement a smoothQuantize function based on the birational map:
  //    i = floor(x);
  //    f = x - i;
  //    f = rsBiRationalMap(f, a);
  //    y = i + f;
  //  This scheme quantizes to integers for a = -1 or a = +1, is neutral for a = 0. We could 
  //  quantize to other intervals q, by pre-scaling x by 1/q and post-scaling y by q (or the other
  //  way around). But this scheme quatizes to floor or ceil (depending on whether a = -1 or +1). 
  //  To quantize to the nearest integer, we may shift f by 0.5 before applying the map and shift y
  //  back by -0.5 (or something like that).
  // -Our group seems to have finite subgroups:
  //  https://math.stackexchange.com/questions/381066/showing-that-this-set-of-functions-is-a-group
  //  https://math.stackexchange.com/questions/2758562/what-causes-a-set-of-functions-to-form-a-group-under-composition
  //  Actually, the 1-parametric family is also an infinite subgroup of the group of all linfracs.
  //  I think, it is a Lie-group...and the mainfold is parametrized by the parameter and maybe we 
  //  have to pick a fixed x to actually evaluate points on the manifold?
  //  Countable subgroups:
  //  https://mathworld.wolfram.com/ModularGroupGamma.html
  //  https://mathworld.wolfram.com/ModularGroupLambda.html
  //  https://en.wikipedia.org/wiki/Modular_group
  //  https://mathworld.wolfram.com/ModularForm.html
  //  https://en.wikipedia.org/wiki/Modular_form
}

//-------------------------------------------------------------------------------------------------

void monotonicInterpolation1()
{
  // We compare linear, cubic and linear fractional ("linfrac") interpolation methods applied to 
  // monotonic data. We create a small set of 5 data points and interpolate that up to a couple
  // of hundreds to see the continuous function that our interpolation schemes create from our
  // data set. For the linfrac method, we will also plot the inverese interpolation function, i.e.
  // the function that results when swapping the roles of y and y.

  using Real = double;
  using Vec  = std::vector<Real>;
  using AT   = RAPT::rsArrayTools;

  // Setup:
  bool decrease = false;  // Switch between monotonically increasing and decreasing
  bool linExtra = false;  // Switch linear extrapolation on (on is the default)
  Real shape    = 0.0;

  // Define datapoints:
  static const int N = 5;
  Real x[N] = { 0, 2, 7, 8,  9 };  // x-values always increase monotonicallly anyway
  Real y[N] = { 0, 5, 6, 9, 10 };  // y-values also increase monotonically here

  // Convert to monotonically decreasing data, if desired:
  if(decrease) {
    for(int n = 0; n < N; n++)
      y[n] = -y[n]; }

  // Allocate arrays for the interpolated data:
  static const int Ni = 501;    // Number of interpolated values
  Real xi[Ni];  
  Real yL[Ni];                  // The L stands for linear
  Real xiMin = -2;              // If < 0, we'll get front-extrapolation
  Real xiMax =  11.0;           // If > 9, we'll get tail-extrapolation
  AT::fillWithRangeLinear(xi, Ni, xiMin, xiMax);

  // Do linear inter-/extrapolation:
  rsInterpolateLinear(x, y, N, xi, yL, Ni);
  // This is our simplemost interpolation scheme and serves as baseline reference. It is actually
  // monotonic by nature, i.e. monotonic data will give rise to monotonic interpolating functions.

  // Compute slopes via numeric differentiation:
  Real s[N];
  //rsNumericDifferentiator<Real>::derivative(x, y, s, N, true);
  rsNumericDifferentiator<Real>::derivative(x, y, s, N, false);
  // These slopes will be used for cubic Hermite and linear fractional interpolation.

  // Do cubic Hermite interpolation:
  Real yH[Ni];                  // The H stands for Hermite
  {
    // Sub block to not litter outer scope with ps, pps.
    Real* ps = s; Real** pps = &ps; // Needed for technical reasons
    rsInterpolateSpline(x, y, pps, N, 1, xi, yH, Ni);
  }

  // Do linear fractional interpolation:
  using LFI = rsLinearFractionalInterpolator<Real>;
  Real yF[Ni];                  // The F stands for fractional
  LFI::interpolate(x, y, s, N, xi, yF, Ni, linExtra, shape);

  // Now let's attempt an inverse interpolation using reversed roles of x- and y-values. According
  // to the inverse function theorem for differentiation, we need the inverse slopes.
  Real yi[Ni];
  Real xF[Ni];
  AT::fillWithRangeLinear(yi, Ni, -4.0, 12.0);
  for(int n = 0; n < N; n++)
    s[n] = 1/s[n];
  LFI::interpolate(y, x, s, N, yi, xF, Ni, linExtra, -shape);


  // Set up the plotter an plot the data along with the interpolants:
  {
    GNUPlotter plt;
    setToDarkMode(&plt);
    plt.setRange(-2, 12, -2, 12);
    plt.setPixelSize(600, 600);
    plt.addCommand("set size square");
    plt.addCommand("set size ratio -1");
    plt.addCommand("set key top left");    // Legend appears top-left
    plt.addCommand("set xtics 1.0");       // x-gridlines at integers
    plt.addCommand("set ytics 1.0");       // y-gridlines at integers
    plt.addDataArrays(N, x, y);            // index 0: samples
    plt.addDataArrays(Ni, xi, yL, yH, yF); // index 1: linear, cubic, linfrac
    plt.addDataArrays(Ni, yi, xF);         // index 2: inverse lin frac
    plt.addDataArrays(Ni, yi, yi);         // index 3: identity for orientation
    plt.addGraph("index 3 using 1:2 with lines lw 1 lc rgb \"#444444\" notitle ");  // identity
    plt.addGraph("index 1 using 1:2 with lines lw 1.5 lc rgb \"#409090\" title \"Linear\"");
    plt.addGraph("index 1 using 1:3 with lines lw 1.5 lc rgb \"#9060A0\" title \"Cubic Hermite\"");
    plt.addGraph("index 2 using 1:2 with lines lw 2 lc rgb \"#666666\" title \"LinFrac Inverse\"");
    plt.addGraph("index 1 using 1:4 with lines lw 2.5 lc rgb \"#BBBBBB\" title \"Linear Fractional\"");
    plt.addGraph("index 0 using 1:2 with points pt 7 ps 1.25 lc rgb \"#FFFFFF\" title \"Samples\"");
    plt.plot();
  }

  // Now try different shapes. We repurpose the array yL, yH (linear and Hermite). Now the L,H will
  // mean low and high values for the shape:
  shape = 8.0;
  linExtra = true;
  for(int n = 0; n < N; n++)
    s[n] = 1/s[n];  // Invert it yet again to get it back to original
  AT::fillWithRangeLinear(xi, Ni, xiMin, xiMax);  // is this needed?
  LFI::interpolate(x, y, s, N, xi, yL, Ni, linExtra, -shape);
  LFI::interpolate(x, y, s, N, xi, yH, Ni, linExtra, +shape);
  {
    GNUPlotter plt;
    setToDarkMode(&plt);
    plt.setRange(-2, 12, -2, 12);
    plt.setPixelSize(600, 600);
    plt.addCommand("set size square");
    plt.addCommand("set size ratio -1");
    plt.addCommand("set key top left");    // Legend appears top-left
    plt.addCommand("set xtics 1.0");       // x-gridlines at integers
    plt.addCommand("set ytics 1.0");       // y-gridlines at integers
    plt.addDataArrays(N, x, y);            // index 0: samples
    plt.addDataArrays(Ni, xi, yF, yL, yH); // index 1: interpolants with different shapes
    plt.addGraph("index 1 using 1:2 with lines lw 2.5 lc rgb \"#BBBBBB\" title \"Shape = 0\"");
    plt.addGraph("index 1 using 1:3 with lines lw 1.5 lc rgb \"#409090\" title \"Shape = -8\"");
    plt.addGraph("index 1 using 1:4 with lines lw 1.5 lc rgb \"#9060A0\" title \"Shape = +8\"");
    plt.addGraph("index 0 using 1:2 with points pt 7 ps 1.25 lc rgb \"#FFFFFF\" title \"Samples\"");
    plt.plot();
  }







  // Observations:
  // -The cubic interpolant clearly wiggles and produces a nonmontonic function.
  // -The linfrac interpolant nicely interpolates our data monotonically.
  // -The inverse interpolant is symmetric to the forward interpolant with respect to the y = x 
  //  line, as would be expected for an inverse function.
  // -Highly negative values (like -8) of shape drag the plateaus toward the left node, letting 
  //  more of the curve lie convexely below the linear connection, highly positive values 
  //  (like +8) will drag the plateaus toward the right node, making a greater portion of the 
  //  curve concavely lying above the linear connection.
  // -With shape = -8, we get a pole in the right extrapolation zone
  // -Positive shape values ten to make the spline "stiffer" at the left data point and negative
  //  values make it stiffer at the right data point

  // -Using xMax = 12 produces garbage when trying to extrapolate using the last computed a,b,c,d
  //  coeffs rather than linear. 11.5..1..7 shows the shooting off behavior. At around 11.6, the
  //  linfrac passes the cubic. For the lower boundary, no  such shooting off is observed. I guess,
  //  at the right side, the pole of the final segment happens to be in the extrapolated zone 
  //  whereas at the left side, this is not the case. we could perhaps let the code detect such 
  //  situations (we can easily calculate the location of the pole) and automatically switch to 
  //  linear extrapolation, if the pole happens to be in the extrapolated zone and otherwise use
  //  linear fractional extrapolation.
  //  

  //
  // ToDo:
  // -Maybe make a plot comparing different values for the shape parameter
  // -[Done] Maybe make the shape parameter have zero as neutral value
  // -Maybe obtain and plot numerical derivatives of the interpolants. We want to see the 
  //  smoothness. We expect that the 1st derivative is continuous but has corners at the nodes and
  //  somewhere in between the nodes where the two linfracs are switched over.
  // -Try to do a reverse interpolation via the linfrac method - yi shoudl go form -4 to 16
  // -Add natural spline interpolation
  // -Maybe use different types for x and y (like float and double) to make it more interesting.
  //  For this, the templates need to be adapted to accept two different template parameters for
  //  x and y. It will make the library code more general.
  // -Test what happens, if xiMin < 0 - we should get extrapolation to the left.
  // -Maybe try different interpolation schemes on some monotonic mathematical functions, for 
  //  example 1/(1+x^n), pow(exp(x), 1/n), log(x), asinh(x), exp(-x^2), tanh(x^3), ...
  //  Some of them are monotonic only when we restrict the domain to the nonegative reals which is
  //  what we should do.
  // -Figure out what happens, if a derivative/slope value of zero occurs. In linear fractional
  //  interpolations, all slopes are expected to be either strictly positive or stricly negative, 
  //  i.e. it expects a *strictly* monotonic function. If that's not the case, it will fail. But I 
  //  guess, when the data is strictly monotonic, so will be the numeric estimates for the desired
  //  slopes. But if we use analytically computed slopes, some of the example functions will have
  //  a zero slope at zero, for example 1/(1+x^2).
}

void monotonicInterpolation2()
{
  using Real = double;
  using Vec  = std::vector<Real>;
  using LFI = rsLinearFractionalInterpolator<Real>;

  int  N     = 11;     // Number of input data points
  int  Ni    = 1001;   // Number of output data points
  Real xMin  =  0.0; 
  Real xMax  = 10.0;
  Real shape =  0.0;

  // Some functions that compute function value and derivative of some example functions

  // The Runge function f(x) = 1 / (1 + x^2) is often used to demonstrate the shortcomings of high
  // order polynomial interpolation. Let's see, how linfrac deals with it. It is not monotonic 
  // though, but if we only take the positive wing, it is - but not strictly so, which is also 
  // already problematic for linfrac interpolation. That's why we shift it to enter a strictly 
  // montonic section:
  auto runge = [](Real x, Real* y, Real* s) 
  { 
    Real shift = 1.0;
    x += shift;

    Real d = 1 + x*x;
    *y = 1 / d;
    *s = -2*x / (d*d);
  };

  // Arc sinh looks a lot like the logarithm but without the pesky pole at zero:
  auto arcsinh = [](Real x, Real* y, Real* s)
  {
    *y = asinh(x);
    *s = 1 / sqrt(1 + x*x);
  };

  // Exponential functions are of course also very important, so we want to see, how our linfrac
  // handles theses as well:
  auto expon = [](Real x, Real* y, Real* s)
  {
    Real k = -0.3;
    *y = exp(k*x);
    *s = k*exp(k*x);
  };

  // The hyperbolic tangent is one of the most famous sigmoid shapes:
  auto hyptan = [](Real x, Real* y, Real* s)
  {
    Real shift = +5; 
    x -= shift;
    Real a = 0.5;
    *y = tanh(a*x);
    Real c = cosh(a*x);
    *s = a / (c*c);               // (tanh(x))' = 1 / (cosh(x))^2
    //*s = a * (1 - (a * *y * *y));  // (tanh(x))' = 1 - (tanh(x))^2, then use chain rule ...seems to be wrong!
  };
  // ToDo: shift it to see the sigmoid behavior

  // The linear fractional (a x + b) / (c x + d) should be perfectly interpolated by the linfrac 
  // scheme. We need to shift the pole outside our range. 
  auto linFrac = [](Real x, Real* y, Real* s)
  {
    Real a = 0, b = 1, c = 1, d = 1;  // d must be > 0 to shift the pole out to the left
    Real D = c*x + d;
    *y = (a*x + b) / D;
    *s = (a*d - b*c) / (D*D);
  };

  // A cubic polynomial should be perfectly interpolated by the cubic Hermite scheme:
  auto cubic = [](Real x, Real* y, Real* s)
  {
    Real a0 = 0, a1 = 3, a2 = 2, a3 = 1;  // coeffs
    Real x2 = x*x, x3 = x2*x;
    *y =  a0 + a1*x + a2*x2 + a3*x3;
    *s =  a1 + 2*a2*x + 3*a3*x2;
  };


  // Select the function to be interpolated:
  auto func = runge;
  //auto func = arcsinh;
  //auto func = expon;
  //auto func = hyptan;
  //auto func = linFrac;
  //auto func = cubic;


  // Generate input data:
  Vec x = rsLinearRangeVector(N, xMin, xMax);
  Vec y(N), s(N);
  for(int n = 0; n < N; n++)
  {
    func(x[n], &y[n], &s[n]);
    //s[n] *= 0.25;  // Test: Intentionally mess up the slope values
  }
  //rsPlotVectorsXY(x, y, s);

  // Compute the correct reference data for y(x) and its slope:
  Vec xi = rsLinearRangeVector(Ni, xMin, xMax);
  Vec yR(Ni), sR(Ni);
  for(int n = 0; n < Ni; n++)
    func(xi[n], &yR[n], &sR[n]);

  // Do cubic Hermite interpolation:
  Vec yH(Ni);                  // The H stands for Hermite
  {
    // Sub block to not litter outer scope with ps, pps.
    Real* ps = &s[0]; Real** pps = &ps; // Needed for technical reasons
    rsInterpolateSpline(&x[0], &y[0], pps, N, 1, &xi[0], &yH[0], Ni);
  }

  // Generate interpolated data:
  //s[0] = -0.01; // fudge to make it nonzero for Runge Function
  Vec yF(Ni);
  LFI::interpolate(&x[0], &y[0], &s[0], N, &xi[0], &yF[0], Ni, true, shape);

  { // Plot:
    GNUPlotter plt;
    setToDarkMode(&plt);
    plt.addDataArrays(N, &x[0], &y[0]);  // index 0: samples
    plt.addDataArrays(Ni, &xi[0], &yF[0], &yH[0], &yR[0]); // index 1: dense, pseudo-continuous
    plt.addGraph("index 1 using 1:4 with lines lw 1.5 lc rgb \"#BBBB33\" title \"Reference\"");
    plt.addGraph("index 1 using 1:3 with lines lw 1.5 lc rgb \"#BB77BB\" title \"Cubic Hermite\"");
    plt.addGraph("index 1 using 1:2 with lines lw 1.5 lc rgb \"#55BBBB\" title \"Linear Fractional\"");
    plt.addGraph("index 0 using 1:2 with points pt 7 ps 1.25 lc rgb \"#FFFFFF\" title \"Samples\"");
    plt.plot();
  }

  // Plot errors:
  rsPlotVectorsXY(xi, yR - yH, yR - yF);


  // Numerically differentiate the interpolants and plot these estimates together with the correct
  // derivative:
  using ND = rsNumericDifferentiator<Real>;
  Vec sH(Ni), sF(Ni);
  ND::derivative(&xi[0], &yH[0], &sH[0], Ni, true);
  ND::derivative(&xi[0], &yF[0], &sF[0], Ni, true);
  rsPlotVectorsXY(xi, sR, sH, sF);



  int dummy = 0;

  // Observations:
  // -For the Runge function with shift = 0, which has a derivative of zero at zero, we need to 
  //  manually fudge the derivative to some small negative value to ensure strict monotonicity, 
  //  i.e. the derivative must be nonzero at all data points. If we don't, we hit an assertion 
  //  and the first segment will contain NaNs and therefore not be plotted by the plotter.
  // -If we do this kind of fudging for the Runge function, the first segment of the linfrac 
  //  interpolant looks really weird! It tries to create a plateau around the peak. The other 
  //  segments look good and are visually indistinguishible from the Hermite interpolant.
  // -Numerically differentiating the interpolants shows a sharp spike of the numerical derivative
  //  of the linfrac interpolant at around x = 0.73. That's within the first, i.e. within the 
  //  problematic segment and is likely the point at which the linfrac interpolator switches 
  //  between the two partial maps.
  // -By tweaking the shape parameter into the negtaive like -1, we can make the first segment look 
  //  almost linear until shortly before it reaches the point 0,1 - we can kind of squeeze the 
  //  weird shape into a smaller region
  // -The numerical derivative of the linfrac interpolant of arcsinh has some corner at 0.48 and is
  //  a bit off compared to the cubic and reference. Also, the sections left and right to the 
  //  corner look suspiciously linear. 
  // -When we intentionally mess up the slope values by scaling them down by 0.25, the linfrac
  //  interpolant tends to become a smoothed stairstep function, i.e. with sigmoidish segments 
  //  between the nodes. The cubic interpolant shows this behavior as well but to a much lesser
  //  extent - it stays closer to the correct function. Might this indicate that linfrac 
  //  interpolation is more sensitive to errors in the target slope values?
  // -The error between reference and interpolated functions for the exponential with k = -0.3 is 
  //  around 4 * 10^-5 for linfrac and around 2 * 10^-5 for cubic - so cubic is twice as good 
  //  (actually even a bit better). For k = +0.3, the situation is similar. For sinh, linfrac seem
  //  to slightly better than cubic in terms of maximum error. May it have to do with concave vs 
  //  convex? May cubic be better for convex functions and linfrac (very slightly) better for
  //  concave functions?
  // -For tanh with a = 0.3, we also see a corner in the numeric derviative. That corner does 
  //  indeed occur at the splitting point between the two half-segments. that can be verified in 
  //  the debugger by looking at xSplit after the line:
  //    T xSplit = getSplitPoint(s1); 
  //  rsLinearFractionalInterpolator<T>::interpolate()
  // -The linfrac function (a x + b) / (c x + d) is perfectly interpolated by the linfrac scheme 
  //  as expected. It's the same thing like the cubic interpolant applied to the cubic polynomial. 
  //  Maybe this could make for a nice unit test.
  // -The error in the derivative seems to be greatest at the split point? That would make seem to
  //  be plausible.

  // Conclusion:
  // -For data for which the derivative tends to zero at one of the end point, linfrac 
  //  interpolation does not seem to give good results.
  // -The additional discontinuities in the 2nd derivative in the interpolant at the split points
  //  between the half segments can be quite visible when looking at the numeric derivative of
  //  the interpolated data. It's not the end of the world, but they seem to be bigger than
  //  the ones at the nodes themselves. 
  // -In terms of smoothness and error, cubic seems to be better than linfrac. Maybe that's the 
  //  price we have to pay for the easy invertibility - a somewhat less smooth interpolant with 
  //  larger maximum error.

  // ToDo:
  // -Maybe it could be useful to let the user specify the shape parameter per segment. Maybe
  //  it could even be possible to devise some algorithm that produces somehow "optimal" shape
  //  values? What would that even mean? Maybe, if we have the true function available, they 
  //  should minimize some error criterion - like minimax error.
  // -Maybe try finding optimal shape values for some common functions like x^2, sqrt, exp, log,
  //  etc. Maybe the optimal value depends of concave vs convex, increasing vs decreasing (nah! 
  //  implausible) or some other features of the function.
  // -Check if the numerical differentiation routine that is used to produce the target values for
  //  the derivatives is guaranteed to produce values > 0 for increasing and < 0 for decreasing 
  //  data at all datapoints. ...I think it should. It doesn't use any negative coeffs...or does 
  //  it?
  // -Plot error between reference and interpolants
  // -Try what happens when we use numerical derivatives instead of analytic ones for the target
  //  values of the slopes. Could these perhaps be better in terms of the resulting overall shape?
}

void monotonicInterpolation()
{
  monotonicInterpolation1(); // 5 datapoints, interpolated via Hermite and linfrac
  monotonicInterpolation2(); // Some interesting functions
}


void interpolatingFunction()
{
  typedef RAPT::rsInterpolatingFunction<float, double> IF;

  // create data to interpolate:
  int N = 5;                          // why not static const?
  float  x[5] = { 2, 4, 5, 7, 8 };    // use x[N] instead of x[5]
  double y[5] = { 1, 3, 1, 2, 3 };    // dito
  // ToDo: 
  // -Explain why x and y are of different type. Probably to make the template instantiations
  //  more "interesting", i.e. to cover a more general case.

  // create and set up interpolating function object:
  IF intFunc;
  //intFunc.setMode(IF::LINEAR);
  intFunc.setMode(IF::CUBIC_HERMITE);
  intFunc.setPreMap( &log);
  intFunc.setPostMap(&exp);
  intFunc.setPreMap( nullptr);
  intFunc.setPostMap(nullptr);


  // do extra/interpolation:
  static const int M = 500; // number of interpolated values, rename to Ni
  float  xi[M];  
  double yi[M];
  float  xiMin = 0;
  float  xiMax = 10;
  RAPT::rsArrayTools::fillWithRangeLinear(xi, M, xiMin, xiMax);
  intFunc.interpolate(x, y, N, xi, yi, M);

  // convert xi to double for plotter and plot:
  double xid[M];
  RAPT::rsArrayTools::convert(xi, xid, M);
  GNUPlotter plt;
  plt.addDataArrays(M, xid, yi);
  plt.setRange(xiMin, xiMax, 0.0, 4.0);
  plt.plot();
  // todo: plot the original data as points
  // todo: plot all the different possible interpolants in one plot for comparison
}

void linearRegression()
{
  static const int N = 500;  // number of samples
  float minDist = 0.1f;       // minimum distance between successive samples
  float maxDist = 1.0f;      // maximum ...
  float a       = 0.8f;
  float b       = 20.0f;
  float noise   = 10.0f;
  int   seed    = 0;

  // create input data:
  float x[N], y[N];
  float dx;      // delta between samples
  x[0] = 10.f;
  y[0] = a*x[0] + b + noise * (float)round(rsRandomUniform(-1, +1, seed)); 
  int n;
  for(n = 1; n < N; n++){
    dx = (float)rsRandomUniform(minDist, maxDist);
    x[n] = x[n-1] + dx;
    y[n] = a*x[n] + b + noise * (float)rsRandomUniform(-1, +1); 
  }

  // retrieve a,b from the data:
  float ar, br = 0;
  rsStatistics::linearRegression(N, x, y, ar, br);
  //ar = Statistics::proportionalRegression(N, x, y);

  // create lines for correct and retrieved line parameters:
  float yc[N], yr[N];
  for(n = 0; n < N; n++){
    yc[n] = a *x[n] + b;
    yr[n] = ar*x[n] + br; }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0], &yc[0], &yr[0]);
  plt.plot();
}

void multipleRegression()
{
  // under construction - we use this example:
  // https://www.youtube.com/watch?v=vvay3t19HU8&list=PLb0zKSynM2PCmp5J5LWM3PcZXBaCoQkXj&index=40
  // The price of pizza is assumed to have a linear relationship with the diameter and the number
  // of extra ingredients, so the model is y = b0 + b1*x1 + b2*x2 + noise.

  static const int N = 7;                                   // number of datapoints (pizza orders)
  double x0[N] = {   1,   1,   1,    1,    1,    1,    1 }; // 1st regressor: constant term 
  double x1[N] = {  24,  25,  28,   30,   40,   36,   32 }; // 2nd regressor: pizza diameters
  double x2[N] = {   0,   1,   0,    1,    2,    1,    3 }; // 3rd regressor: number of extras
  double y[N]  = { 5.5, 8.9, 7.4, 10.9, 19.4, 13.8, 14.0 }; // regressand:    pizza prices

  // create input data matrix:
  double* X[3];
  X[0] = x0;
  X[1] = x1;
  X[2] = x2;

  // Compute right hand side for linear equation system (X * X^T) * b = X * Y:
  double Y[3];

  // on wikipedia, the formula reads (X^T * X) * b = X^T * Y because there, the data arrays with 
  // the independent variables are stored in the columns of the X-matrix whereas we store them in 
  // the rows

  //RAPT::rsMatrixTools::transposedMatrixVectorMultiply(X, y, Y, N, 3);
  // something is wrong - maybe we don't need the transposed version bcs our X matrix is defined 
  // differently? - yes - we have the regressors in the rows, wikipedia has them in the columns

  //RAPT::rsMatrixTools::matrixVectorMultiply(X, y, Y, N, 3); // access violation
  RAPT::rsMatrixTools::matrixVectorMultiply(X, y, Y, 3, N);   // works but shouldn't according to doc

  // establish matrix on left-hand-side:
  double flatXX[9];     // flat storage array
  double* XX[3];        // array of pointers to rows
  XX[0] = &flatXX[0];
  XX[1] = &flatXX[3];
  XX[2] = &flatXX[6];
  RAPT::rsMatrixTools::matrixMultiplySecondTransposed(X, X, XX, 3, 7, 3); // is this correct?

  // estimate parameter vector by solving the linear system:
  double b[3];  
  RAPT::rsLinearAlgebra::rsSolveLinearSystemInPlace(XX, b, Y, 3); // get rid of rs prefix
  // verify, if the result is correct!

  // Use the model to predict the pizza prices from their diameters and extras:
  double yp[N], ea[N], er[N];
  for(int n = 0; n < N; n++) {
    yp[n] = b[0]*x0[n] + b[1]*x1[n] + b[2]*x2[n]; // predicted price
    ea[n] = y[n]  - yp[n];                        // absolute error
    er[n] = ea[n] / y[n];                         // relative error
  }

  // plot - draw true datapoints and prdicted ones:
  GNUPlotter plt;
  plt.addDataArrays(N, x1, x2, y);  // the triplets should be interpreted as 3D points
  plt.addDataArrays(N, x1, x2, yp);
  plt.setGraphStyles("points", "points");  // maybe select style - or use better default
  plt.plot3D(); 
  // that doesn't look too bad actually :-) ...but verify and clean up
  // todo: plot the plane (semitransparently?) that is defined by our model
}

// https://en.wikipedia.org/wiki/Linear_regression#Least-squares_estimation_and_related_techniques
// https://en.wikipedia.org/wiki/Weighted_least_squares

// ToDo: generalize and factor out the code for multiple regression and move it to rapt 
// (Numerics/DataFitting.h/cpp), then 
// implement polynomial regression on top of that (use 1,x,x^2,x^3,... as regressors for a given
// x-array) - this may e generalized to use a set of arbitrary functions of x - the model is a 
// weighted sum of these functions of x - for example, we could model a signal with sines/cosines
// of given frequencies - maybe create experiments for that (polynomialRegression, 
// sinusoidalRegression, exponentialRegression) - problem: the frequnecies and decays must be known
// in advance - can these be estimated too? ...for exponential decays, there's code somewhere but 
// what about sine frequencies? maybe the exponential can be made complex? -> figure out

void polynomialRegression()
{
  // We create data from a polynomial function f(x) = 1 - 2x + 3x^2 - 4x^3 + 5x^4 with added noise 
  // and try to estimate the polynomial coefficients from the data. For this estimation, we once 
  // use the clean data and once the noisy data. When the clean data is used and the degree of the 
  // model is >= 4, it is expected to get a perfect fit.

  // User parameters:
  double noise         = 2.0;   // amount of noise
  int    numDataPoints = 100;   // number of data points
  int    modelDegree   = 4;     // degree of our model polynomial
  double xMin          = -1.0;
  double xMax          = +2.0;

  typedef std::vector<double> Vec;
  typedef RAPT::rsPolynomial<double> Poly;
  typedef RAPT::rsArrayTools AT;
  typedef RAPT::rsCurveFitter CF;

  // The polynomial to generate our data is 1 - 2x + 3x^2 - 4x^3 + 5x^4:
  Poly p(Vec({ 5,-4,3,-2,1 }));  // todo: have a constructor that takes an inititalizer list

  // Generate clean and noisy data:
  int N = numDataPoints;
  Vec x(N), yc(N), yn(N);  // x-values, clean and noisy y-values
  AT::fillWithRangeLinear(&x[0], N, xMin, xMax); // todo: maybe use noisy x-data, too
  RAPT::rsNoiseGenerator<double> prng;
  for(int n = 0; n < N; n++) {
    yc[n] = p(x[n]);
    yn[n] = yc[n] + noise * prng.getSample(); }

  // Fit polynomials to clean and noisy data:
  Poly qc = CF::fitPolynomial(&x[0], &yc[0], N, modelDegree); // fit to the clean data
  Poly qn = CF::fitPolynomial(&x[0], &yn[0], N, modelDegree); // fit to the noisy data

  // Create prediction data from the two models:
  Vec zc(N), zn(N);
  for(int n = 0; n < N; n++) {
    zc[n] = qc(x[n]);
    zn[n] = qn(x[n]); }

  // Plot:
  rsPlotVectorsXY(x, yc, zc);      // plot clean true data and model-predicted data
  rsPlotVectorsXY(x, yc, yn, zn);  // plot noisy
  //GNUPlotter plt;

  // Observations:
  // -When modelDegree >= 4, we get a perfect fit of the model to the clean data - which is 
  //  expected because the clean data was generated by a 4th degree polynomial.
  // -The fit to the noisy data is a bit below the true generation function in the range x = 0..1, 
  //  presumably because the noise has a local mean < 0 in that region. This leads to the strange
  //  circumstance that our model fits the data actually better than the original generator 
  //  function - at least when modelDegree >= 4. Is this already overfitting? Kind of, i guess - 
  //  although, for modelDegree == 4, we are actually not able to overfit...soo...welll...
  // -When modelDegree < 4, we don't get perfect fits for the clean data anymore - but still 
  //  reasonable ones.
  // -Conclusion: it works! :-)

  // todo: 
  // -make plots that show several model-orders in a single plot to see, how the fit gets 
  //  better as the model degree goes up
  // -maybe the library should provide a function that suggests a model order based on the drop-off
  //  of the approximation error as function of the model degree - when the gains diminish, we have
  //  found the "right" model degree
  
  // -compare approximating a sinewave with polynomials - once using a Taylor series and once using
  //  a least squares fit over one cycle (maybe from -pi...pi)
  //  -yet another way to approximate a sinewave is to match a number of derivatives at the start 
  //   and at the end of a cycle (Taylor series matches N derivatives at t = 0, this variant should 
  //   match N/2 derivatives at t = 0 and N/2 at t = 2*pi
  //  -what about minimax fits? how would they be computed? maybe iteratively, using a 
  //   least-squares fit as initial guess?
  // -make a similar function for sinusoidalRegression - what about spline-regression? that 
  //  could be very useful, too
  // -the API should probably allow the user to provide a set of arbitrary basis functions to use
  //  -for polynomial regression, these basis functions are the powers of x
  //  -we could also use chebychev polynomials as basis functions
  //   -this is probably numerically better (the Vandermonde matrix in polynomial regression is 
  //    notoriously ill conditioned)
  //   -it requires the data to be normalized, such that x and y goes from -1..+1 
  //    -we need an affine trafo before and one after the actual polynomial
  //   -we may want to convert the chebychev regression coeffs to power-coeffs after fitting

  // -ToDo: 
  //  -implement fitting a Gaussian bell cirve to data, see:
  //     https://www.youtube.com/watch?v=jezAWd6GFRg
  //   it fits a 2nd order polynomial to the log of the data. the polynomial coeffs can the diractly
  //   be converted to mean, variance and scale-factor
  //  -adds function to the library that does regression in terms of arbitrary basis functions, 
  //   passed as std::vector<std::function> or something like that...oh...seems like we have 
  //   something like that already done below in gaussianRegression
}

void gaussianRegression()
{
  // Under construction
  // We try to fit a bunch of cosines with gaussian envelopes to the "multipass Butterworth" 
  // function 
  //    f(x) = 1 / (1 + x^2N)^M
  //  -> nice test for fitting an arbitrary set of basis functions to data
  //  -> may be useful for implementing the IIR lens blur effect


  static const int N = 200;  // number of data points
  double xMax       = 3.0;
  double freqFactor = 1.0;   // scales all frequencies, 1.0 seems most suitable
  int order  = 6;            // Butterworth order
  int passes = 1;            // number of passes for Butterworth filter
  int numGaussians = 7;      // number of Gaussians in the approximation

  using Vec = std::vector<double>;
  using AT  = rsArrayTools;

  Vec f = rsRangeLinear(0.0, numGaussians-1.0, numGaussians); // freqs of the cosines
  Vec w = freqFactor*0.5*PI*f;                                // omegas
  // todo: variances (sigma)

  // create the target data:
  double x[N], y[N], z[N];
  AT::fillWithRangeLinear(x, N, 0.0, xMax);
  for(int n = 0; n < N; n++)
    y[n] = 1.0 / pow((1 + pow(x[n], 2*order)), passes);

  // create the regressands:
  rsMatrix<double> R(numGaussians, N);
  for(int i = 0; i < numGaussians; i++)
    for(int n = 0; n < N; n++)
      R(i, n) = exp(-x[n]*x[n]) * cos(w[i]*x[n]); // todo: use adjustable variances


  plotMatrixRows(R, x);   // plot the regressors
  //rsPlotArraysXY(N, x, y); 

  Vec A = rsCurveFitter::multipleRegression(R, y); // amplitudes of the Gaussians
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  AT::fillWithZeros(z, N);
  for(int i = 0; i < numGaussians; i++)
    for(int n = 0; n < N; n++)
      z[n] += A[i] * R(i, n);
  plt.addDataArrays(N, x, z);
  plt.plot();

  // Observations:
  // -using w = 0.5*PI*freqs seems to give the best results
  // -with order = 5, passes = 1 numGaussians = 7, the approximation result is quite good
  // -when making the transition steeper, the Gaussians overshoot/ripple more - stronger Gibbs-like
  //  effect
  // -we seem to need more Gaussians to approximate a steeper transition well
  // -a higher number of passes also seems to increase the ripple

  // ToDo:
  // -optimize the frequencies and variances of the Gaussians by gradient descent, starting with an
  //  initial guess obtained like above - the scale-factors should also be refined in this 
  //  optimization procedure
  // -define the parameter vector as p = (a0 w0 s0 a1 w1 s1 a2 w2 s2 ...)

  // initial guess for the parameter vector using 7 Gaussians:
  Vec p = Vec({0.768,0,1, 0.759,1,1, -0.427,2,1, -0.258,3,1, 0.0536,4,1, 0.2297,5,1, -0.1325,6,1});
  // we need to define an error function and optimize the amplitudes, frequencies and widths with 
  // that initial guess - so we need to do an optimization in and 21-dimensional space

  //int dummy = 0;
}

/** Error function for the approximation of a Butterworth function by a sum of Gaussians.. */
class ButterworthGaussError : public RAPT::MultivariateErrorFunction<double>
{

public:

  double getGaussSum(const rsVectorDbl& p, double x)
  {
    rsAssert(p.dim % 3 == 0, "Length of p should be a multiple of 3");
    int numGaussians = p.dim / 3;
    double result = 0.0;
    for(int i = 0; i < numGaussians; i++)
    {
      double a = p[3*i];    // amplitude, weight
      double f = p[3*i+1];  // frequency of cosine
      double s = p[3*i+2];  // sigma (variance, width of envelope)
      double w =  0.5*PI*f; // omega (radian frequency of cosine)
      result += a * exp(-x*x/(s*s)) * cos(w*x);
    }
    return result;
  }

  double getTarget(double x)
  {
    return 1.0 / pow((1 + pow(x, 2*order)), passes);
  }

  double getX(int i)
  {
    return i*xMax / double(numSamples-1);
  }

  double getError(rsVectorDbl p, double x)
  {
    //double target = 1.0 / pow((1 + pow(x, 2*order)), passes);
    double target = getTarget(x);
    double actual = getGaussSum(p, x);
    return target - actual;
  }

  double getValueMaxAbs(rsVectorDbl p)
  {
    double maxAbs = 0.0;
    for(int i = 0; i < numSamples; i++)
    {
      //double x = i*xMax / double(numSamples-1);
      double x = getX(i);
      double delta  = rsAbs(getError(p, x));
      if(delta > maxAbs)
        maxAbs = delta;
    }
    return maxAbs;
  }

  double getValueMeanSquare(rsVectorDbl p)
  {
    double sum = 0.0;
    for(int i = 0; i < numSamples; i++)
    {
      //double x   = i*xMax / double(numSamples-1);
      double x   = getX(i);
      double err = getError(p, x);
      sum += err*err;
    }
    return sum / numSamples;
  }

  virtual double getValue(rsVectorDbl p) override
  {
    //return getValueMaxAbs(p);
    return getValueMeanSquare(p);
    //return sqrt(getValueMeanSquare(p));
  }

  void plot(rsVectorDbl p)
  {
    GNUPlotter plt;
    std::vector<double> x(numSamples), t(numSamples), y(numSamples), e(numSamples);
    for(int i = 0; i < numSamples; i++)
    {
      x[i] = getX(i);
      t[i] = getTarget(x[i]);
      y[i] = getGaussSum(p, x[i]);
      e[i] = t[i] - y[i];
    }
    plt.addDataArrays(numSamples, &x[0], &t[0], &y[0], &e[0]);
    plt.plot();
  }

protected:

  int numSamples = 200;
  double xMax = 3.0;

  int order  = 6; // lower orders are easier to approximate
  int passes = 1;

};

void butterworthViaGaussians()
{
  // Under construction
  // Optimizes the parameters to approximate a Butterworth function by a sum of cosines with 
  // Gaussian envelopes...

  using Minimizer = GradientBasedMinimizer<double>;

  ButterworthGaussError errFunc;
  Minimizer minimizer;
  //minimizer.setAlgorithm(minimizer.GRADIENT_DESCENT);
  //minimizer.setAlgorithm(minimizer.BOLD_DRIVER_WITH_MOMENTUM);
  minimizer.setAlgorithm(minimizer.SCALED_CONJUGATE_GRADIENT);
  //minimizer.setAlgorithm(minimizer.CONJUGATE_GRADIENT);
  //minimizer.setMaxNumSteps(300);
  minimizer.setStepSize(1.0);
  minimizer.setConvergenceThreshold(1.e-7);  // with 1.e-8, it doesn't converge

  minimizer.setPrintInfo(true);
  // when using the max-abs error, the info-printing messes up the screen and the result is not 
  // good

  using Vec = std::vector<double>;

  Vec ps({0.768,0,1, 0.759,1,1, -0.427,2,1, -0.258,3,1, 0.0536,4,1, 0.2297,5,1, -0.1325,6,1});
  // maybe make a function to compute an initial guess based on linear regression

  // test: use a perturbed start vector
  //Vec p({0.768,0,1, 0.759,1.1,1, -0.427,2.1,1, -0.258,3,1, 0.0536,4,1, 0.2297,5,1, -0.1325,6,1});

  rsVectorDbl pInitial(ps);  // the requirement to convert is unelegant!

  errFunc.plot(pInitial);

  rsVectorDbl pFinal = minimizer.minimizeFunction(&errFunc, pInitial);
  Vec pe = pFinal.toStdVector();  
  // almost unchanged - something is wrong! it uses the maximum number of iterations and the 
  // stepsize is denormal at the end
  // OK - using minimizer.setAlgorithm(minimizer.GRADIENT_DESCENT); fixes this - but then, the
  // error oscillates - no convergence

  // todo: plot the sum of Gaussians after optimization - hmm - the improvement is not very 
  // impressive (from 0.027 to 0.011) i think, something is still wrong
  // -> optimize also with scipy and compare results

  errFunc.plot(pFinal);

  //int dummy = 0;

  // ToDo: move the old minimization code into the rs_testing module - it's not good enough for the
  // library yet
}




bool testNumericMinimization()
{
  // Under construction - we compare several approaches to minimize a nonlinear multidimensional 
  // function

  bool r = true;

  using Vec = std::vector<double>;


  // Our example is the bivariate scalar function:
  //   f(x,y) = 4*x^2 + y^4 + x*y
  // where v = (x y) is the 2D input vector.
  static const int N = 2;  // dimensionality of the input space
  std::function<double(double*)> f;
  double hx = pow(2, -18);
  double hy = hx/2;
  double h[N] = {hx, hy};     // h as array
  Vec v({5,3});               // initial guess
  Vec x;
  int evals;

  double tol = 1.e-12;
  double y;

  // try simpler functions: f(x,y) = x^2 + y^2, f(x,y) = x^2, ... 


  f = [=](double* v)->double
  { 
    double x = v[0], y = v[1];
    return 2*x*x + 16*y*y;
    // a*x^2 + b*y^2 - should converge in the 1st iteration, irrespective of a,b (but both must
    // be positive)
  };
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, tol);  // 16 evals
  x = v; evals = minimizeNewton(f, &x[0], N, h, tol);


  //x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1.0, tol); // diverges

  // the minimum is at (0,0)
  // 16 evaluations: it converges in the first iteration, but a 2nd iteration is needed to detect 
  // the convergence, so we get 2 iterations, each taking  4*N = 4*2 = 8 evaluations - so this 
  // seems ok
  x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1./16, tol); 
  // oscillates, "converges" after 276 evals but to the wrong location

  x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1./32,  tol);  // converges, 571 evals
  x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1./64,  tol);  // 1151 evals
  x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1./128, tol);  // 2271 evals

  x = v; evals = minimizeGradientAutoStep( f, &x[0], N, h, tol); // 287 evals

  //minimizeGradientAutoStep(const F& f, T* v, int N, const T* h, T tol = 1.e-8)


  // try a function like f(x,y) = a*x^2 + b*y^2 + c*x + d*y


  // A function that has a general quadratic form:
  //   f(v) = v*A*v + b*v + v with matrix A, vector b and scalar c
  // with
  //   A = |1 2|, b = |5|, c = 7
  //       |3 4|      |6|
  // gives:
  //   f(x,y) = x^2 + 5*x*y + 4*y^2 + 5*x + 6*y + 7
  // the coeff 5 for the x*y term comes from the off-diagonal elements: 2+3. The function as is has 
  // its minimum at (0,0) but if we want it to be at (xMin,yMin), we simply replace x,y by x-xMin 
  // and y-yMin in the formula
  f = [=](double* v)->double
  { 
    double xMin = 0; 
    double yMin = 0;
    double x = v[0] - xMin;
    double y = v[1] - yMin;
    return x*x + 5*x*y + 4*y*y + 5*x + 6*y + 7;
  };
  x = v; evals = minimizeNewton(f, &x[0], N, h, tol);
  // wait: this function has a saddle:
  // https://www.wolframalpha.com/input/?i=+x*x+%2B+5*x*y+%2B+4*y*y+%2B+5*x+%2B+6*y+%2B+7%3B

  // 3*x*x + 2*x*y + 4*y*y + 5*x + 6*y + 7;
  f = [=](double* v)->double
  { 
    double x = v[0]; double y = v[1];
    return 3*x*x + 2*x*y + 4*y*y + 5*x + 6*y + 7; 
  };
  x = v; evals = minimizeNewton(f, &x[0], N, h, tol); y = f(&x[0]);
  // result looks correct, see here:
  // https://www.wolframalpha.com/input/?i=Minimize%283*x*x+%2B+2*x*y+%2B+4*y*y+%2B+5*x+%2B+6*y+%2B+7%29
  // minimum is at (-7/11,-13/22) = (-0.6363..,-0.59090..) and has a value of 40/11 = 3.6363...

  f = [=](double* v)->double
  { 
    double x = v[0], y = v[1];
    return 4*x*x + y*y*y*y + x*y;
  };
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-8);  // 14 iterations, 112 evals (112/14 = 8)
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-9);  // 15 its, 120 evals
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-10);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-11);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-12); // 18 its, 144 evals
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-13);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-14);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-15);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-16);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-17);
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-18);
  // 14 iterations, 112 evaluations with tol = 1.e-8    -> 112/14 = 8 
  // with each power of 10, we need 8 evaluations (i.e. one iteration) more - this looks like
  // linear convergence :-( ...i hoped that with a parabolic approximation, we could get quadratic
  // convergence - ToDo: plot trajectories
  // 18 iterations, 144 evaluations with tol = 1.e-12
  // do these number suggest quadratic convergence? todo: compare to gradient descent - this is
  // the baseline against which all other algos are measured - maybe use an optimal constant 
  // stepsize for gradient descent

  x = v; evals = minimizeGradientDescent(f, &x[0], N, h, 1./16,  1.e-12); //    31 - but wrong
  x = v; evals = minimizeGradientDescent(f, &x[0], N, h, 1./32,  1.e-12); //  5356 - to min2
  x = v; evals = minimizeGradientDescent(f, &x[0], N, h, 1./64,  1.e-12); //  9841 - to min1
  x = v; evals = minimizeGradientDescent(f, &x[0], N, h, 1./128, 1.e-12); // 19246 - to min1
  x = v; evals = minimizeGradientDescent(f, &x[0], N, h, 1./256, 1.e-12); // 36816 - to min1 

  // to converge to the same minimum as the partial parabolic algo, we nee a stepSize of 1/64
  // and 9841 evals compared to 144 evals with the partial parabolic algo

  x = v; evals = minimizeGradientAutoStep( f, &x[0], N, h, 1.e-12);  // 599



  f = [=](double* v)->double
  { 
    double x = v[0], y = v[1];
    return sin(x*y);    // 
    //return sin(3*x*x + y*y*y + x*y);    // 
  };
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h,          1.e-12);  //    16
  x = v; evals = minimizeGradientDescent( f, &x[0], N, h, 1./1024, 1.e-12);  // 216701
  x = v; evals = minimizeGradientAutoStep(f, &x[0], N, h,          1.e-12);
  // the results are slightly different - how can this be - they should converge to the same
  // minimum - oh - the minmima are valleys - that can be the reason
  // maybe the convergence criterion should not be the change of function value


  // https://www.wolframalpha.com/input/?i=4*x*x+%2B+y*y*y*y+%2B+x*y
  // roots: x = -(3 sqrt(3))/128, y = sqrt(3)/8; x = (3 sqrt(3))/128, y = -sqrt(3)/8
  // minima: min1 = -0.000976563 at (x, y)=(-0.0220971,  0.176777)
  //         min2 = -0.000976563 at (x, y)=( 0.0220971, -0.176777)


  // todo: try the minimizePartialParabolic on the rosenborck function - see here:
  // https://en.wikipedia.org/wiki/Adaptive_coordinate_descent
  //

  // https://en.wikipedia.org/wiki/Rosenbrock_function
  f = [=](double* v)->double
  { 
    double x = v[0], y = v[1];
    double a = 1;
    double b = 100;
    double ax  = a-x;
    double yx2 = y-x*x;
    return ax*ax + b*yx2*yx2;
    // minimum: (a,a^2)
  };
  v = Vec({-3,-4});
  x = v; evals = minimizePartialParabolic(f, &x[0], N, h, 1.e-12);  y = f(&x[0]); // 32344

  //x = v; evals = minimizeNewton(f, &x[0], N, h, tol); y = f(&x[0]);
  // doesn't converge - wildly jumps around in the parameter space


  // other test functions to try:
  // https://en.wikipedia.org/wiki/Test_functions_for_optimization
  // https://en.wikipedia.org/wiki/Himmelblau%27s_function
  // https://en.wikipedia.org/wiki/Rosenbrock_function
  // https://en.wikipedia.org/wiki/Shekel_function
  // https://en.wikipedia.org/wiki/Rastrigin_function
  // https://en.wikipedia.org/wiki/Ackley_function


  // Algorithms:
  // https://en.wikipedia.org/wiki/Coordinate_descent
  // https://en.wikipedia.org/wiki/Adaptive_coordinate_descent
  // https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
  // https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  // https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
  // https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm
  // https://en.wikipedia.org/wiki/Limited-memory_BFGS

  // https://en.wikipedia.org/wiki/Gradient_descent
  // https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization

  return r;
}

void numericOptimization()
{
  testNumericMinimization();

  //int dummy = 0;
}

void numericMinimization1D()
{
  // We test the numeric minimization routines on the following problem: find the parameter value a
  // which minimizes the following error function: E(a) = |f(0.25,a)-0.5| + |f(0.75,a)-2| where 
  // f(x,a) = (x / (1-x))^a. This is a real world problem that occured in the discord chat on
  // "The Audio Programmer" server. The problem was to find a function f(x) that goes through the 
  // following points: (0,0), (0.25,0.5), (0.5,1), (0.75,2), (1,inf). It was solved by using
  // (x/(1-x))^a  where  a = log_3(2)  does the job exactly, i.e. the function matches all the 
  // desired points exactly. It matches (0,0), (0.5,1), (1,inf) by construction for any value of a
  // but to get a match also at (0.25,0.5) and (0.75,2), we need to pick the specific value of 
  // a = log_3(2). The value can be found analytically but here we do it numerically. The function
  // is meant to map a normalized parameter in the range 0..1 to the range 0..inf in some "natural"
  // way. One potential disadvantage of that mapping function is that it has an infinite derivative
  // at x=0 which is kinda bad for dialing in values near zero. 

  using Real = double;
  std::function<Real(Real, Real)> f;  // Parametric mapping function f(x,a) = (x/(1-x))^a
  std::function<Real(Real)>           // Error functions of a
    E1,                               // E1(a) = |f(0.25,a)-0.5| + |f(0.75,a)-2|
    E2,                               // E2(a) = max(|f(0.25,a)-0.5|, |f(0.75,a)-2|)
    E3;                               // E3(a) = (f(0.25,a)-0.5)^2 + (f(0.75,a)-2)^2

  // The mapping function with input x and parameter a is defined as  f(x,a) = (x/(1-x))^a
  f = [](Real x, Real a)
  {
    Real y = x / (1-x);
    return pow(y, a);
  };

  int numCalls = 0;              // We keep track of the number of calls in the optimization
  Real aTgt = rsLogB(2.0, 3.0);  // Target value a = log_3(2) - the analytical solution

  // The error function is a function of the parameter a:
  E1 = [&](Real a)
  {
    numCalls++;
    Real err1 = fabs(f(0.25,a) - 0.5);
    Real err2 = fabs(f(0.75,a) - 2.0);
    return err1 + err2;
  };

  E2 = [&](Real a)
  {
    numCalls++;
    Real err1 = fabs(f(0.25,a) - 0.5);
    Real err2 = fabs(f(0.75,a) - 2.0);
    return max(err1, err2);
  };

  E3 = [&](Real a)
  {
    numCalls++;
    Real err1 = f(0.25,a) - 0.5;
    Real err2 = f(0.75,a) - 2.0;
    return err1*err1 + err2*err2;
  };


  // The minimum of the error function E(a) should occur at a = log3(2) at which point the error is
  // actually zero. Let's verify this fact using a numeric minimization:
  numCalls  = 0;
  Real aOpt = rsMinimizer1D<Real>::goldenSectionMin(E1, 0.0, 1.0);
  bool ok   = aOpt == aTgt;  // They are indeed exactly equal - there's no roundoff error.
  // numCalls is now 77, so the golden section algo needed to eavluate E1(a) 77 times.

  // Test it with the alternative error function basd on maximum deviation:
  numCalls = 0;
  aOpt = rsMinimizer1D<Real>::goldenSectionMin(E2, 0.0, 1.0);
  ok   = aOpt == aTgt;
  // numCalls is again 77

  // Test it with the error function based on the sum of squares:
  numCalls = 0;
  aOpt = rsMinimizer1D<Real>::goldenSectionMin(E3, 0.0, 1.0);
  ok   = aOpt == aTgt;
  // numCalls is again 77

  rsAssert(ok);

  // Observations:
  // -The number of iterations of the golden section algo is the same for all the different 
  //  variants of the error function. It doesn't matter if we take the sum-of-squares, sum-of-abs 
  //  or max-abs error. We always get 77 iterations. The algorithm is insensitive to the exact 
  //  shape of the function as long as the minimum stays in place - which is actually quite 
  //  plausible when thinking about it. That may actually be a nice feature in certain contexts.

  // ToDo:
  // -Implement Brent's method and compare the number of evaluations. I guess, for Brent's method, 
  //  the number will depend on the exact shape of the function.
  // -Maybe define an alternative error function based on the sum of squares or on max(err1, err2)
}

void numericRootFinding1D()
{
  bool ok = true;

  using Real = double;
  std::function<Real(Real)> f;

  int numCalls;
  f = [&](Real x)
  { 
    numCalls++;
    return (x+1)*(x-1)*(x-2); 
  }; // Polynomial with 3 roots at -1,+1,+2

  Real x;
  Real tol = std::numeric_limits<Real>::epsilon();
  tol *= 0.5;

  int maxIts = 1000;
  numCalls = 0;
  x = brents_fun(f, -1.3f, -0.8f, tol, maxIts);
  ok &= x == -1;
  // numCalls = 45

  // Notes:
  // -When using epsilon as tolerance, the root at -1 in not found exactly. The last 2 decimal 
  //  digits are wrong. When using 0.5*epsilon, we find the root at -1 exactly.

  // ToDo:
  // -Compare the number of calls to other methods like bisection, falsePosition


  rsAssert(ok);
}


void polynomialSinc()
{
  // We want to find a polynomial that resembles a windowed sinc to be used for high order 
  // interpolation that is easy to vectorize for simd. The first feature of the sinc that we want
  // to replicate is to have zero crossings at all integers except 0. Such a polynomial has the
  // general form:
  //
  //   p(x) = c*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)*...
  //
  // For some constant c that can be determined by evaluating the product with c=1 and then taking
  // the reciprocal, such that the scaled polynomial has an amplitude of 1 at x = 0. This 
  // polynomial by itself does not have the right behavior with regard to its envelope. The 
  // amplitudes of the oscillations grow as one goes further outward:
  //   https://www.desmos.com/calculator/bh2kkuj2t4
  // So, the next step is to produce a suitable window polynomial w(x) that can be multiplied by 
  // p(x) to make the product look like a windowed sinc. If given a fractional position into the 
  // sample with integer index n and fractional part f, we would evaluate p(x) and w(x) at 
  // f,f+1,f-1,f+2,f-2,f+3,f-3,f+4 and sum the values (maybe dividing by the sum of the weights, 
  // such that the sum of all weights is always 1). The evaluation of p(x) can be nicely optimized
  // by observing that: (x+n)(x-n) = x^2 - n^2 := x2 - n2, so our polynomial above simplifies to:
  //   p(x) = (x2-1)*(x2-4)*(x2-9)*(x2-16) 
  // where x^2 can be precomputed. By further observing that (x2-n)*(x2-m) = x4 -(n+m)*x2 + n*m
  // = x2*(x2-(n+m)) + n*m...hmm...no...i don't think, that this leads to further simplification...


  // User parameters:
  int N = 512;              // number of datapoints
  int windowExponent = 10;  // determines shape of the window and therefore frequency response
  int halfNumZeros   =  8;  // number of positive zeros (== 1/2 the total number)


  using Real = double;
  using Vec  = std::vector<Real>;

  Real xMin  = -halfNumZeros;
  Real xMax  = +halfNumZeros;

  // The raw function with zeros at all integers except 0, evaluated naively:
  auto p1raw = [&](Real x)
  {
    Real y = 1;
    for(int i = 1; i <= halfNumZeros; i++)
      y *= (x+i)*(x-i);
    return y;
  };

  // The constant by which we need to scale the raw polynomial. This can be known at compile time
  // when a fixed halfNumZeros is chosen:
  Real c = 1.0 / p1raw(0.0);

  // The scaled function p(x), evaluated naively: 8 mul, 4 add, 4 sub for halfNumZeros == 4:
  auto p1 = [&](Real x) { return c*p1raw(x);  };

  // p(x) with optimized evaluation using x^2: 5 mul, 4 sub for halfNumZeros == 4:
  auto p2 = [&](Real x)  
  { 
    Real x2 = x*x; 
    Real y  = 1;
    for(int i = 1; i <= halfNumZeros; i++)
      y *= (x2-i*i);  // i*i is a compile-time constant
    return c*y;
  };

  // The window function:
  auto w1 = [&](Real x)
  { 
    x /= halfNumZeros;
    Real t = -(x-1)*(x+1);  // zero at +-1
    return pow(t, windowExponent);
  };
  // ToDo: form a linear combination of t, t^2, t^3...or more generally, a polynomial in t, with 
  // the goal to achieve a nice frequency response of the resulting interpolator. 


  Vec x = RAPT::rsRangeLinear(xMin, xMax, N);
  Vec p(N), w(N), pw(N);
  for(int n = 0; n < N; n++)
  {
    p[n]  = p1(x[n]);
    w[n]  = w1(x[n]);
    pw[n] = p[n] * w[n];

    //Real err = p2(x[n]) - p[n];
    //Real tol = 1.e-11;   // we need high tolerance for high halfNumZeros
    //rsAssert(rsAbs(err) <= tol);
  }
  //rsPlotVectorsXY(x, p, w, pw);
  rsPlotVectorsXY(x, pw, w);
  //rsPlotVectorsXY(x, pw);

  int fftSize = 32*N;
  int plotMax = 1024;
  //int plotMax = 4*N;
  Vec mag(fftSize), phs(fftSize); 
  rosic::fftMagnitudesAndPhases(&pw[0], N, &mag[0], &phs[0], fftSize);
  mag = mag * (1.0/mag[0]);  // normalize at DC
  for(int n = 0; n < fftSize; n++)
    mag[n] = rsMax(rsAmp2dB(mag[n]), -200.0);  // convert to dB with floor
  rsPlotArray(&mag[0], plotMax);

  //rsPlotSpectrum(mag, 0.0, -100.0, false);


  // Observations:
  // -The frequency response is generally lowpassish in nature, with passband and stopband ripple
  // -With increasing windowExponent, the cutoff frequency goes down (bad), the passband ripple
  //  decreases (good) and the spectral rolloff increases (good)
  // -When increasing the number of zeros, keeping everything else constant, the amount of 
  //  passband ripple increases. It can be counteracted by choosing a higher window exponent. 
  //  Seems like the windoExponent should equal halfNumZeros for a good looking freq resp.
  // -For higher num zeros, the optimized evaluation seems to produce progressively less 
  //  accurate results (or is it the naive version that gets less accurate? anyway, they diverge)
  // -A choice of 8,8 seems to give a reasonably good looking response. Maybe, if a datatype
  //  float32x16 is available, that should the the first choice - we'll see.

  // ToDo: 
  // -maybe also plot the spectrum of the linear interpolator and of a windowed sinc interpolator
  //  for comparison
  // -try a linear combination of different window exponents with weigths summing to 1. Maybe try
  //  to use t, t^2, t^4, t^8 (easy to compute)
  // -try to expand the optimized form and evaluate it via Horner's rule -> check numerical and 
  //  performance properties

  // Altenatively, we may try a polynomial that has a zero also at zero and then divied the whole
  // thing by x...but we have that expensive division and we also need to take care about avoiding
  // division by zero, so that should be considred only, when the above approach fails...
}


/*
Idea: we model the signal x(t) by a polynomial such that:
x(t) = a0 + a1*t + a2*t^2 + a3*t^3 + ... + aN*t^N
and we have signal values x(0), x(-1), x(-2), x(-3),... available. We want to find the 
polynomial coeffs and then evaluate x(1) to predict/extrapolate the signal value at t=1 
(conveniently, at t=1, the evaluation of the polynomial boils down to a sum over the coeffs)
We get a linear system of equations:
x( 0) = a0
x(-1) = a0 -   a1 +    a2 -    a3 +     a4 -      a5 + ...
x -2) = a0 - 2*a1 +  4*a2 -  8*a3 +  16*a4 -   32*a5 + ... 
x(-3) = a0 - 3*a1 +  9*a2 - 27*a3 +  81*a4 -  243*a5 + ...
x(-4) = a0 - 4*a1 + 16*a2 - 64*a3 + 256*a4 - 1024*a5 + ...
...
hmm...the matrix entries grow exponentially - that's probably not good for numerical 
conditioning - maybe the prediction can be made without explictly solving the system by 
somehow constructing the Netwon polynomial?....figure that out

but maybe a rational function is better than a polynomial?
        a0 + a1*t + a2*t^2 + a3*t^3 + ... + aN*t^N
x(t) = --------------------------------------------
        b0 + b1*t + b2*t^2 + b3*t^3 + ... + bM*t^M
this also leads to a linear system of equations:
x( 0) *  b0                             = a0
x(-1) * (b0 -   b1 +   b2 -   b3 + ...) = a0 -   a1 +    a2 -    a3 + ...
x(-2) * (b0 - 2*b1 + 4*b2 - 8*b3 + ...) = a0 - 2*a1 +  4*a2 -  8*a3 + ...
...
the advantage of the rational function might be that it doesn't necessarily grow so much, in 
fact, if N==M, it approaches aN/bM as t -> inf and if M > N, it approaches zero
maybe try a biquadratic -> 6 coeffs a0,a1,a2,b0,b1,b2 but actually we have only 5 degrees of 
freedom due to being able to normalize the polynomials (make them monic and drag a "gain" 
factor before the fraction or something) so we would use 5 samples

Acually, it is well known that polynomials are not really well suited for data extrapolation, not 
even interpolation, actually....hmm.
*/

void polynomialPrediction()
{

}


/*
        a0 + a1*t + a2*t^2
x(t) = --------------------
        b0 + b1*t + b2*t^2

x( 0) *  b0                 = a0
x(-1) * (b0 -   b1 +    b2) = a0 -   a1 +    a2
x(-2) * (b0 - 2*b1 +  4*b2) = a0 - 2*a1 +  4*a2
x(-3) * (b0 - 3*b1 +  9*b2) = a0 - 3*a1 +  9*a2
x(-4) * (b0 - 4*b1 + 16*b2) = a0 - 4*a1 + 16*a2

let's normalize, such that b2 = 1 (make denominator monic, overall scale factor is a2), we get 5
equations for 5 unknowns and solve them with sage:

var("x0 x1 x2 x3 x4 a0 a1 a2 b0 b1 b2")
e1 = x0 *  b0              == a0
e2 = x1 * (b0 -   b1 +  1) == a0 -   a1 +    a2
e3 = x2 * (b0 - 2*b1 +  4) == a0 - 2*a1 +  4*a2
e4 = x3 * (b0 - 3*b1 +  9) == a0 - 3*a1 +  9*a2
e5 = x4 * (b0 - 4*b1 + 16) == a0 - 4*a1 + 16*a2
solve([e1,e2,e3,e4,e5],[a0,a1,a2,b0,b1])

which gives the result:

a0 = 12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)*x0
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
a1 = (x2*x3*x4 + (x1*(7*x2 - 36*x3 + 30*x4) + x2*(45*x3 - 64*x4) + 18*x3*x4)*x0 - (x2*(16*x3 - 27*x4) + 12*x3*x4)*x1)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
a2 = (x2*x3*x4 + (x1*(x2 - 6*x3 + 6*x4) + x2*(9*x3 - 16*x4) + 6*x3*x4)*x0 - (x2*(4*x3 - 9*x4) + 6*x3*x4)*x1)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
b0 = 12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
b1 = (x0*(x1 - 12*x2 + 27*x3 - 16*x4) + x1*(18*x2 - 64*x3 + 45*x4) + 6*x2*(5*x3 - 6*x4) + 7*x3*x4)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4)

*/
double biquadraticPrediction(double x0, double x1, double x2, double x3, double x4)
{
  double s, a0, a1, a2, b0, b1;
  s  = 1.0 / (x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4);
  a0 = s * (12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)*x0);
  a1 = s * ((x2*x3*x4 + (x1*(7*x2 - 36*x3 + 30*x4) + x2*(45*x3 - 64*x4) + 18*x3*x4)*x0 
            - (x2*(16*x3 - 27*x4) + 12*x3*x4)*x1));
  a2 = s * ((x2*x3*x4 + (x1*(x2 - 6*x3 + 6*x4) + x2*(9*x3 - 16*x4) + 6*x3*x4)*x0 
            - (x2*(4*x3 - 9*x4) + 6*x3*x4)*x1));
  b0 = s * (12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4));
  b1 = s * ((x0*(x1 - 12*x2 + 27*x3 - 16*x4) + x1*(18*x2 - 64*x3 + 45*x4) 
            + 6*x2*(5*x3 - 6*x4) + 7*x3*x4));
  // maybe optimize ..the product x3*x4 appears often

  return (a0 + a1 + a2) / (b0 + b1 + 1.0); // evaluate rational function at x=1
}
//....but: what if the function has a pole at or near x=1? or what if s = inf ...maybe rational 
// extrapolation is not such a good idea after all...hmmm


// maybe move to RAPT into the Statistics section
double variance(double *x, int N)
{
  double mx  = RAPT::rsArrayTools::mean(x, N); // sample mean of x
  double sum = 0;
  for(int n = 0; n < N; n++) {
    double d = x[n] - mx;
    sum += d*d;
  }
  return sum / (N-1);

  // ToDo: 
  // -Make an optimized version that takes the mean as argument. But then we may want to 
  //  distinguish two cases: (1) the mean is exact, (2) the mean is itself an estimate. The correct
  //  (unbiased) formula is different for these two cases. If the exact mean is known from theory, 
  //  we get an unbiased estimator by dividing by N. Usually, the mean is not known exactly but 
  //  also estimated from the data (as we do here). In such a case, we must divide by N-1 (as we 
  //  do here) to get an unbiased estimator for the variance. See:
  //    https://www.statlect.com/fundamentals-of-statistics/variance-estimation
  //  for an explanation. I don't fully understand it, though. I wonder why they call it
  //  "degrees of freedom" here:
  //    "...the number of degrees of freedom is equal to the number of sample points (n) 
  //    minus the number of other parameters to be estimated."
  //  By "degrees of freedom", I usually think of a bunch of variables (degrees of freedom) and a 
  //  bunch of constraints (equations that the variables must satisfy). But I don't see, how this
  //  relates to this context here. Are the datapoints seen as variables and the estimated values
  //  somehow seen as "constraints"? Like:
  //    Variables:   x1, x2, x3, ... 
  //    constraints: mean(x1,x2,x3,...) = M; var(x1,x2,x3,...) = V
  //  or something? If so, why? Why do variables map to the data and constraints map to estimates?
  //  What has estimation to do with constraining? I have never seen that explained. They always 
  //  just say: "Hey, we call this number "degrees of freedom"" without explaining why. But even 
  //  without fully understanding it, it seems clear that dividing by N-1 is the apparently the 
  //  right thing to do, if the mean is itself just an estimate or sample-mean.
  // -See also: https://mathworld.wolfram.com/SampleVariance.html
  //  Interesting: Even if we use an unbiased estimator for the variance sigma^2, the estimated 
  //  standard variation sigma by taking the square-root of sigma^2 is biased again. Statistics is
  //  complicated and confusing. If this code goes into production code, these issues should be
  //  carefully documented (and ideally also explained why things are the way they are).
}
double standardDeviation(double *x, int N)
{
  return sqrt(variance(x, N));
}
double covariance(double *x, double *y, int N)
{
  double mx = RAPT::rsArrayTools::mean(x, N); // mean of x
  double my = RAPT::rsArrayTools::mean(y, N); // mean of y
  double sum = 0;
  for(int n = 0; n < N; n++)
    sum += (x[n]-mx) * (y[n]-my);
  return sum / (N-1);
}
double correlation(double *x, double *y, int N)
{
  double vxy = covariance(x, y, N);
  double vx  = variance(x, N);
  double vy  = variance(y, N);
  return vxy / sqrt(vx*vy);
}
double conditionalProbability(double* a, double* b, int N)
{
  // conditional probability of "a given b"
  // a and b must be arrays of boolean values (0 or 1) but as double data type
  int na = 0, nb = 0;
  for(int n = 0; n < N; n++) {
    if(b[n] == 1) {
      nb++;
      if(a[n] == 1)
        na++;
    }
  }
  return (double) na / (double) nb;

  // can this somehow be generalized...like sum(prod(a, b)) / sum(b) ...sum over the elementwise 
  // product of the realizations of a,b divided by sum of realizations of b - should allow
  // realizations between 0..1
}
// or maybe it's called conditional relative frequency?
// https://mathbitsnotebook.com/Algebra1/StatisticsReg/ST2TwoWayTable.html

double jointProbability(double* a, double* b, int N)
{
  int sum = 0;
  for(int n = 0; n < N; n++)
    if(a[n] == 1 && b[n] == 1)
      sum += 1;
  return sum / (double) N;
  // maybe generalize by using a sum-of-products (divided by N)
}
// maybe rename to andProbability and write also an orProbability


void probabilityLogic()
{
  // consider an example: produce random number x in 0..1, 

  // let: A:  x < 0.5, B: 0.3 < x < 0.6,
  // -> P(A)=0.5, P(B)=0.3, P(A|B)=2/3, P(B|A)=2/5=0.4 -> what's the correlation C(A,B)?


  // setup:
  int N = 100000;   // number of realizations to produce
  double aL = 0.0;  // lower limit for x, such that event A is considered to have occured
  double aU = 0.5;  // upper limit ....
  double bL = 0.3;  // ...same for event B
  double bU = 0.6;


  typedef std::vector<double> Vec;
  Vec A(N), B(N);   // realizations of events A and B (true/false, represented as 0/1)
  Vec x(N);
  int n;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0, 1);
  for(n = 0; n < N; n++)
  {
    double xn = prng.getSample();
    x[n] = xn;

    if(xn > aL && xn < aU)  
      A[n] = 1;
    else         
      A[n] = 0;

    if(xn > bL && xn < bU)
      B[n] = 1;
    else
      B[n] = 0;
  }

  // compute relative frequencies of events A and B (should approximate their probabilities):
  //double fA = RAPT::rsArrayTools::sum(&A[0], N) / N;
  //double fB = RAPT::rsArrayTools::sum(&B[0], N) / N;
  // are actually the mean values

  // compute sample mean values for event A and B:
  double mA = RAPT::rsArrayTools::mean(&A[0], N);
  double mB = RAPT::rsArrayTools::mean(&B[0], N);

  // compute sample variances:
  double vA = variance(&A[0], N);
  double vB = variance(&B[0], N);

  // compute sample covariance and correlation:
  double cov = covariance( &A[0], &B[0], N);
  double cor = correlation(&A[0], &B[0], N);

  // compute empirical probabilities (by relative frequencies):
  double pa  = RAPT::rsArrayTools::sum(&A[0], N) / N;        // P(A), empirical prob of event A
  double pb  = RAPT::rsArrayTools::sum(&B[0], N) / N;        // P(B), empricial prob of event B
  double cab = conditionalProbability(&A[0], &B[0], N); // P(A|B), empirical prob of A given B
  double cba = conditionalProbability(&B[0], &A[0], N); // P(B|A), empirical prob of B given A
  double jab = jointProbability(      &A[0], &B[0], N); // P(A,B), empirical prob of A and B

  // compute joint probability by formulas via conditional probability:
  //double jab1 = cab * pb;   // P(A,B) = P(A|B) * P(B)
  //double jab2 = cba * pa;   // P(A,B) = P(B|A) * P(A)
  // ok, this works - both are equal to jab

  // soo, let's say P(A) and P(B) are known and let's assume that the correlation (or maybe 
  // covariance) between a and b is also known - can we compute P(A,B)?
  // how does the correlations relate to the conditional probabilities P(A|B) or P(B|A)?

  // https://math.stackexchange.com/questions/1751950/from-correlation-coefficient-to-conditional-probability
  // says:
  // Corr(A,B) = (P(A,B) - P(A)*P(B)) / sqrt( P(A)*(1-P(A)) * P(B)*(1-P(B)) )

  // let's try it:
  double cor2 = (jab-pa*pb) / sqrt(pa*(1-pa)*pb*(1-pb));
  // yes, looks good

  // from this, we may build a general continuous "and" formula that incorporates correlation
  // ...and from that, we can also create an "or" formula

  //int dummy = 0;
}

double productLog(const double z) 
{
  // Evaluates the product-log function W(z) defined as the inverse of the function f(W) = W * e^W.
  // This function is also known as Lambert-W or Omega function.

  // adapted from http://keithbriggs.info/software/LambertW.c
  // here's more: http://keithbriggs.info/software.html

  if(0.0 == z) 
    return 0.0;

  double eps = 4.0e-16;
  double em1 = 0.3678794411714423215955237701614608; // 1/Euler
  double p, e, t, w;

  if(z < -em1+1e-4) // series near -em1 in sqrt(q)
  { 
    double q = z+em1, r = sqrt(q), q2 = q*q, q3 = q2*q;
    return 
      -1.0
      +2.331643981597124203363536062168*r
      -1.812187885639363490240191647568*q
      +1.936631114492359755363277457668*r*q
      -2.353551201881614516821543561516*q2
      +3.066858901050631912893148922704*r*q2
      -4.175335600258177138854984177460*q3
      +5.858023729874774148815053846119*r*q3
      -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }

  // Initial approximation for Halley iteration:
  if(z < 1.0)     // series near 0
  { 
    p = sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));                 // euler-number
    w = -1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); // -1/3, 11/72
  } 
  else 
    w = log(z);   // asymptotic

  if(z > 3.0) 
    w -= log(w);  // useful?

  // Halley iteration:
  for(int i = 0; i < 10; i++) 
  { 
    e  = exp(w); 
    t  = w*e-z;
    p  = w+1.0;
    t /= e*p-0.5*(p+1.0)*t/p; 
    w -= t;
    if(fabs(t) < eps*(1.0+fabs(w)))  // rel-abs error
      return w;
  }

  //rsAssertFalse(); // no convergence
  return 0.0;
}
void productLogPlot()
{
  int N = 1000;        // number of values
  float xMin = -3.0;
  float xMax = +15.0;

  // evaluate function:
  vector<float> x(N), y(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
    y[n] = (float) productLog(x[n]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}






// try something based on recursive subdivision of an interval at a ratio r: Take the interval
// 0..1, split it into to sections of length r and 1-r, apply the same subdivision to both
// sub-intervals (starting with the larger one), maybe alternate which side get r and 1-r
// ...but what if the number of desired splits is not a power of 2?...always subdivide the
// currently largest subinterval...and the decision which side gets r and which gets 1-r is based
// on which of the neighbouring sections is larger
// example: r=0.8
// 0...1 -> 0....0.8..1 ...that's good for the supersaw because when changing the number of
// saws (intervals/splits), the old splits stay where they are and it is just a new split/saw
// introduced - that means we can just fade it in when implementing a continuous density
// parameter (we may have to do something special for the transition between 1 and 2, though - like
// moving the first saw in frequency while fading in the second - but as soon as the 2 outer freqs
// are established at their places at the interval limits, more saws are just faded in when 
// increasing density continuously)
std::vector<double> intervalSplitting(int numSegments, double midpoint)
{
  int N = numSegments;
  double r = midpoint;
  std::vector<double> a(N+1);
  a[0] = 0.0;
  a[N] = 1.0;
  int n = 2;    // number of points
  int kl = 0;   // lower index
  int ku = N;   // upper index
  while(n < N)
  {
    int km = (int)round(r*(ku-kl));
    a[km] = a[kl] + r * (a[ku] - a[kl]);

    // update kl, ku
    // ...

    n++;
  }
  return a;
}


bool rangeStartLess(const rsRange<double>& r1, const rsRange<double>& r2)
{
  return r1.getMin() < r2.getMin();
}
void splitRange(const rsRange<double>& r, rsRange<double>& rl, rsRange<double>& ru, double ratio)
{
  rl = rsRange<double>(r.getMin(),  r.getMin() + ratio*r.getSize()); // lower part of range r
  ru = rsRange<double>(rl.getMax(), r.getMax());                     // upper part of range r
}
std::vector<double> intervalSplittingProto(int numSegments, double midpoint, int splitStrategy = 1)
{
  // prototype implementation - may not be optimal in terms of efficiency, but is algorithmically
  // easier to understand

  int N = numSegments;
  double r = midpoint;
  typedef rsRange<double> Range;
  std::vector<Range> s(N);  // array/set of segments/intervals/ranges (initially empty)
  s[0] = Range(0.0, 1.0);   // seed set-of-ranges is { [0,1) }
  int n = 1;                // current number of intervals
  while(n < N) 
  {
    int k = rsArrayTools::maxIndex(&s[0], n);    // index of largest range in the current set

    // split s[k] into two ranges:
    Range rl, ru;
    if(splitStrategy == 0) splitRange(s[k], rl, ru, r);  // always use r
    else if(splitStrategy == 1)                          // alternating, ...
      if(rsIsOdd(n))  splitRange(s[k], rl, ru, r);       // ...odd n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...even n uses 1-r
    else if(splitStrategy == 2)                          // alternating, ...
      if(rsIsEven(n)) splitRange(s[k], rl, ru, r);       // ...even n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...odd n uses 1-r
    // is it somehow possible to continuously morph between the odd and even version of this
    // strategy? ...maybe just "crossfade" the result split-point arrays? that would add another
    // potentially interesting dimension for tweaking ...maybe even vector-crossfade between
    // skew-right/skew-left/alternate-odd/alternate-even?

    // the lower part of s[k] replaces sk, the upper part gets appended to the array:
    s[k] = rl;
    s[n] = ru;
    n++;
  }
  // the cost of this algorithm is O(N^2) because in each iteration, we search for the maximum 
  // which is itself an O(N) operation - can we avoid this search by somehow keeping the array of
  // ranges sorted (by size), thereby turning it into an O(N) algorithm?

  // sort the array ranges by their start-point:
  rsHeapSort(&s[0], N, &rangeStartLess);
  // this is an O(N*log(N)) operation - so if we can turn the above into an O(N) operation, the 
  // overall comlexity of range splitting would still be O(N*log(N)) - certainly much better than
  // O(N^2) - but to achieve O(N), we would have to avoid the final sorting too - maybe by always
  // keeping a version sorted by size and another version sorted by start around?

  // turn the ranges into a vector of split-points:
  std::vector<double> a(N+1); // split-points
  for(n = 0; n < N; n++)
    a[n] = s[n].getMin();
  a[n] = 1.0;
  return a;
}

void ratioGenerator()
{
  int numRatios = 6;     // number of ratios (i.e. "density")
  int numParams = 200;   // number of sampe values for the parameter of the ratio algo
  int numPrimes = numRatios+1;

  std::vector<rsUint32> primes(numPrimes);
  rsFillPrimeTable(&primes[0], numPrimes);

  typedef rsRatioGenerator<double>::RatioKind RK;
  rsRatioGenerator<double> ratGen;
  ratGen.setPrimeTable(&primes);
  ratGen.setRatioKind(RK::linToExp);
  ratGen.setParameter1(0.5);

  std::vector<double> ratios(numRatios);
  ratGen.fillRatioTable(&ratios[0], numRatios);

  std::vector<double> p(numParams);  // rename to params
  std::vector<double> r(numRatios);  // rename to tmp
  rsArrayTools::fillWithRangeLinear(&p[0], numParams, 0.0, 1.0);
  double** y; // rename to r(atios)
  rsMatrixTools::allocateMatrix(y, numRatios, numParams);

  for(int i = 0; i < numParams; i++) {
    ratGen.setParameter1(p[i]);
    ratGen.fillRatioTable(&r[0], numRatios);
    rsArrayTools::transformRange(&r[0], &r[0], numRatios, 1., 2.);  // ratios in 1...2
    for(int j = 0; j < numRatios; j++)
      y[j][i] = r[j];
  }



  GNUPlotter plt;
  plt.addDataArrays(numParams, &p[0], numRatios, y);
  plt.setRange(-0.1, 1.1, 0.9, 2.1);
  plt.plot();


  rsMatrixTools::deallocateMatrix(y, numRatios, numParams);

  // primePowerDiff has graphs that cross each other - that means that there are values for the
  // parameter for which not all ratios are distinct - but we want them to be distict

  // make an LinToExp mapping where p=0 means linear spacing, p=1 means exponential spacing 
  // (lin-spacing of the logs) and in between, there's some intermediate setting (crossfade)

  // f = i / (numRatios-1); // fraction 
  // linVal = 1 + i / (numRatios-1);
  // expVal = exp(linVal) / exp(2);
  // mixVal = (1-p)*linVal + p*expVal;




  // try the self-similar interval splitting algorithm:
  //std::vector<double> a = intervalSplittingProto(100, 1.0/GOLDEN_RATIO, 1);
  //std::vector<double> a = intervalSplittingProto(100, 0.9, 2);
  std::vector<double> a = intervalSplittingProto(16, 0.75, 1);
  double mean = rsMean(a);;
  //rsPlotMarkers(&a[0], (int) a.size());
  // on the gui for the supersaw audio module, we should visualize the spread in a similar way
  // maybe with vertical lines instead of + markers - the visualization may or may not use a 
  // normalized interval - if it doesn't, it will also visualize the amount of spread ..maybe
  // also show the sawtooth amplitudes by the height of the lines - especially useful when we later
  // introduce amp envelopes (bell-shaped curves, but mayby also inverse bell curves emphasizing 
  // the outer freqs - how about: 0: flat, 1: bell, -1: inverted bell)

  // todo: plot the ratios as function of the parameter, let the parameter go from 0 to 1
}

void ratiosLargeLcm()
{
  // We want to produce rational numbers r between 1 and 2 that have the property that when their
  // decimal expansion (or any other expansion) is truncated, the truncated numbers have a large
  // lowest common multiple

  // We want to see all LCMs of a all pairs in a range of numbers in order to pick out those, for
  // which we have a resonably large LCM for any pair of numbers between some nMin and nMax. When 
  // nMin is fixed, we actually only need to consider numbers up to nMax = 2*nMin-1 - but we plot 
  // more...
  //  ...maybe it makes sense to choose prime numbers from the upper half
  // ..to make them irrational, square them, add 1 and take the square root

  unsigned int nMin = 100;
  unsigned int nMax = 2*nMin-1;
  unsigned int N = nMax - nMin;

  // fill LCM matrix:
  std::vector<double> axes(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&axes[0], N, double(nMin), double(nMax));
  double** lcmMatrix;
  RAPT::rsMatrixTools::allocateMatrix(lcmMatrix, N, N);
  for(unsigned int i = 0; i < N; i++)
    for(unsigned int j = 0; j < N; j++)
      lcmMatrix[i][j] = (double) RAPT::rsLcm(nMin+i, nMin+j);


  GNUPlotter plt;
  plt.addDataMatrix(N, N, &axes[0], &axes[0], lcmMatrix);
  plt.setPixelSize(800, 800);
  plt.addCommand("set size square");

  // factor out into fucntion plotHeatMap
  plt.addGraph("i 0 nonuniform matrix w image notitle");
  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  plt.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale
  plt.plot();
  //plt.plot3D();

  RAPT::rsMatrixTools::deallocateMatrix(lcmMatrix, N, N);
}

void ratiosEquidistantPowers()
{
  // We produce ratios in the following way: the user gives a minimum and maximum value a,b. Then,
  // we produce equidistant values between A := a^p and B := b^p with some user parameter p. Then,
  // we take the p-th root of the so produced values.

  using Real = float;   // to switch between float and double at compile time
  using Vec  = std::vector<Real>;
  using AT   = rsArrayTools;

  int numRatios = 11;     // number of ratios (i.e. "density")
  int numParams = 201;    // number of sample values for the parameter of the ratio algo
  Real pMin     = -10.0;  // lower value of the exponent parameter
  Real pMax     = +12.0;  // upper value of the exponent parameter
  Real a        = 10.0;   // lower frequency or period
  Real b        = 20.0;   // upper frequency or period

  // for testing a suitable threshold/tolerance for the absolute value of p for treating it as 
  // zero:
  //Real smalll = 2000 * RS_EPS(Real);  // should be around twice the tolerance for good plot
  //pMax = smalll; pMin = -pMax;      // uncomment to test the edge case of p ~= 0


  Vec col(numRatios);  // holds one matrix column at a time
  Vec p(numParams);    // hold the array of parameter values 
  AT::fillWithRangeLinear(&p[0], numParams, pMin, pMax);

  // create the matrix for plotting:
  rsMatrix<Real> R(numRatios, numParams);
  for(int j = 0; j < numParams; j++) {
    AT::fillWithRange(&col[0], numRatios, a, b, p[j]);
    R.setColumn(j, &col[0]); } 
  plotMatrixRows(R, &p[0]);

  // plot the generalized mean of a,b and of the produced columns:
  Vec gmc(numParams);   // generalized mean of current column
  Vec gmab(numParams);  // generalized mean of a and b
  for(int i = 0; i < numParams; i++) {
    R.copyColumn(i, &col[0]);
    gmc[i] = AT::generalizedMean(&col[0], numRatios, p[i]);
    col[0] = a; col[1] = b;
    gmab[i] = AT::generalizedMean(&col[0], 2, p[i]); }
  rsPlotVectorsXY(p, gmc, gmab);
  // -OK - gmc and gmab do match indeed
  // -the curve looks a bit like atan - could it be atan indeed?
  // -maybe plot also the regular mean - maybe we should somehow "fix" the mean such that the 
  //  perceived pitch always stays the same, regardless of p
  // -figure out, how it behaves when the the minimum is neagtive...it probably will not work in
  //  any resonable way...what would we actually want/expect in such a case?

  // Notes:
  // Maybe we should somehow make a bipolar version of the shape - use different (but related) 
  // formulas for the lower and upper half.
  // ToDo: figure out which values of lead to the least phasiness when used for frequencies in the
  // supersaw...maybe p = 0.5 is nice because it's "halfway" between linear and exponential 
  // spacing, both of which are musically meaningful in their own way (but both of which lead to
  // heavy phasing when sued for supersaw) ...maybe -0.5 is also interesting - we should really 
  // make this a user parameter...
  // -figure out, if we can use the AGM (arithmetic-geometric mean) in this context in a 
  //  meaningful way
}

void ratiosMetallic()
{
  // The metallic ratio rn for integer n is given by (n + sqrt(n^2+4))/2 - so the ratio of two
  // metallic ratios for integers m, n is: rn/rm = (n + sqrt(n^2+4)) / (m + sqrt(m^2+4))
  // ...maybe obtain the continued fraction expansion of these numerically and see, if we get small
  //  coefficients


  //int dummy = 0;
}


void sinCosTable()
{
  // A test for the rsSinCosTable class.

  rsSinCosTable<float> table(8); // parameter is the table size 

  // create data:
  int N = 2000;  // number of values to plot
  float xMin = -15.0;
  //float xMin =   0.0;
  float xMax = +15.0;
  vector<float> x(N), ySin(N), yCos(N), ySinTbl(N), yCosTbl(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    ySin[n] = sin(x[n]);
    yCos[n] = cos(x[n]);
    //table.getValuesRounded(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesTruncated(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesLinear(x[n], &ySinTbl[n], &yCosTbl[n]);
    table.getValuesCubic(x[n], &ySinTbl[n], &yCosTbl[n]);
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &ySin[0], &yCos[0], &ySinTbl[0], &yCosTbl[0]);
  plt.plot();
}

void expBipolar()
{
  // find coeffs for a parametrized exponential function y = f(x) = a * e^(b*x) + c such that 
  // f(0) = y0, f(1) = y1, f'(0) = s. We have f'(x) = a*b*e^(b*x). The equation system is:
  // y0 = a + c         
  // y1 = a * e^b + c
  // s  = a * b
  // solving 1 and 3 for c and b and plugging into 2 gives:
  // d = y1 - y0 == a e^(s/a) - a
  // which can be given to wolfram alpha as:
  // Solve[d == a e^(s/a) - a, a ]

  double y0 = 0.2;
  double y1 = 8.0;
  double s  = 1.0;

  // compute coeffs:
  //double dy  = y0 - y1;
  //double tmp = (exp(s/dy)*s)/dy;
  //tmp = productLog(tmp);
  //double a = s*dy / (-dy * tmp + s); // nope - formula must be wrong 

  double d = y1 - y0;
  double w = productLog(- (exp(-s/d)*s) / d);
  double a = -d*s / (d*w+s);
  double b = s  / a;
  double c = y0 - a;

  // compute values at the endpoints for test:
  double f0  = a * exp(b * 0) + c; // f0 is wrong
  double f1  = a * exp(b * 1) + c;
  double fp0 = a * b;

  //int dummy = 0;
}

void expGaussBell()
{
  // This is still wrong. The idea is to find a,b,c,d parameters for the function
  // f(x) = a * exp(b*x) + c * exp(d*x^2)
  // and impose the constraints f(0) = 1, f'(0) = s0, f(1) = 0, f'(1) = s1.
  // we have:
  // f (x) = a * exp(b*x) +   c  *  exp(d*x^2)
  // f'(x) = a*b*exp(b*x) + 2*c*d*x*exp(d*x^2)
  // so:
  // (1) f (0) = 1  = a + c
  // (2) f'(0) = s0 = a * b
  // (3) f (1) = 0  = a*exp(b) + c*exp(d)
  // (4) f'(1) = s1 = a*b*exp(b) + 2*c*d*exp(d)
  // solve(1=a+c,s_0=a*b,0=a*exp(b)+c*exp(d),s_1=a*b*exp(b)+2*c*d*exp(d);a,b,c,d) -> nope
  // we still need to solve the system - maybe an expression that contains only one unknown 
  // - say a - can be derived and we can st up a Newton iteration to the 1D root finding problem
  // we have: (1) c = (1-a), (2) b = (s0/a), 
  // (3) d = ln(-a*exp(b)/c) = ln(-a*exp(s0/a)/(1-a))
  // (4) s1 = a*b*exp(b) + 2*c*d*exp(d)
  // s1 = s0*exp(s0/a) + 2*(1-a)*d*exp(d)
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a))*exp(ln(-a*exp(s0/a)/(1-a))) - s1
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a)) * -a*exp(s0/a)/(1-a) - s1
  // ...more math to do....


  // settings:
  int N = 300;        // number of values
  float xMin = 0.0;
  float xMax = 1.0;
  double s0   = -0.1f; // slope at 0
  double s1   = -0.2f; // slope at 1

                       // compute coefficients:
  double a, b, c, d;
  a = -s0 / productLog(0.5*s0*s1);
  b = -s0 / a;
  c = 1 - a;
  d = b; // wrong!

         // evaluate function:
  vector<float> x(N), y(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    double t = x[n]; // temporary
    t = a*exp(b*t) + c*exp(d*t);
    y[n] = (float)t;
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}


void spreadAmounts(double amount, double spread, double* amount1, double *amount2)
{
  double absAmt = amount;
  double target = 1;
  if(amount < 0) {
    absAmt = -amount;
    target = -target;
  }
  *amount1 = RAPT::rsLinToLin(absAmt, 0.0, 1.0,  spread, target);
  *amount2 = RAPT::rsLinToLin(absAmt, 0.0, 1.0, -spread, target);
}

double spreadAmounts1(double x, double y)
{ 
  double out1, out2; spreadAmounts(x, y, &out1, &out2); 
  return out1;
}
double spreadAmounts2(double x, double y)
{
  double out1, out2; spreadAmounts(x, y, &out1, &out2);
  return out2;
}
void twoParamRemap()
{
  // ...figure out, if for every output pair (z1,z2) there is an input pair (x,y) that leads to
  // the desired output pair...maybe try to find the inverse function 
  // we have (z1,z2) = f(x,y) -> find (x,y) = g(z1,z2) where g is the inverse of z

  // maybe plot contour lines on both functions z1(x,y), z2(x,y) and see if every possible pair
  // of such contour lines has at least one intersection - done - yes - that seems to work

  int N = 21;
  GNUPlotter plt;
  plt.addDataBivariateFunction(N, -1.0, +1.0, N, -1.0, +1.0, &spreadAmounts1); 
  plt.addDataBivariateFunction(N, -1.0, +1.0, N, -1.0, +1.0, &spreadAmounts2);
  plt.addCommand("set hidden3d");
  //plt.addCommand("set pm3d");
  plt.addCommand("set palette gray");
  plt.addCommand("set lmargin 0");         // margin between plot and left border
  plt.addCommand("set tmargin 0");         // margin between plot and top border
  plt.addCommand("set contour");       // draws contours on the base plane
  //plt.addCommand("set contour surface"); // draws contours on the surface
  //plt.addCommand("set contour both"); // draws contours on the surface and base plane
  plt.addCommand("set cntrparam levels incr -1,0.2,1");

  //plt.addCommand("set palette rgbformulae 8, 9, 7");
  //plt.addCommand("set style fill solid 1.0 noborder");
  //plt.addCommand("set pm3d depthorder noborder");
  //plt.addCommand("set pm3d lighting specular 0.6");

  //plt.addCommand("set view 0, 0");
  plt.addCommand("set view 60, 45");
  plt.plot3D();



  //http://gnuplot.sourceforge.net/demo/contours.html
}

void numberTheoreticTrafo()
{
  // We implement a fast convolution by means of a number theorectic transform (NTT). An NTT is 
  // similar to an FFT, but instead of working in the field of complex numbers, it works in a 
  // finite field Z_p of modular integers. Here, we work with p = 97, which is a prime, which 
  // guarantees that all the required "magic numbers" exist in this field. We need the N-th roots
  // of unity and their inverses where N = 2^k for all N less than maxN/2, where maxN is the
  // maximum supported transform size and we also need the inveres of all these N (verify - i may
  // have the factor 2 wrong...). ...I'm also not sure anymore that p being prime ensures the 
  // existence of all desired N-th roots of unity (where N = 2^k for all k < p). Maybe it only 
  // ensures, that they are invertible, if they happen to exist? I think, typically, there is some 
  // highest k, for which the root exists, but it may be much less than p/2...figure out..

  // Sage code to find the primitive n-th roots of unity for n = 2^k, k = 1,..,5
  //
  // p = 97
  // Zp = Integers(p)
  // root_list = Zp.zeta(16, all=True);  # replace the 16 by 2,4,8,32 for the other roots
  // root_list
  //
  // gives the following lists of roots:
  //
  //  2: 96
  //  4: 75, 22
  //  8: 64, 50, 47, 33
  // 16: 89, 85, 79, 70, 27, 18, 12, 8
  // 32: 78, 77, 69, 67, 63, 55, 52, 51, 46, 45, 42, 34, 30, 28, 20, 19
  //
  // Interestingly, for n = 64, there are no solutions anymore. Is it generally true that n-th 
  // roots exist only for n < p/2?. For any n, we can just pick any of the n-th roots. Right? Or 
  // should we pick the smallest? For the time being, i have just picked the first from the list.
  // The inverses of these roots can be found by, e.g.:
  //
  // root_list = [96, 75, 64, 89, 78]
  // [inverse_mod(x, 97) for x in root_list]
  //
  // which gives [96, 22, 47, 12, 51] Similarly, the inverses of the possible FFT lengths can be 
  // found by:
  //
  // length_list = [2, 4, 8, 16, 32]
  // [inverse_mod(x, 97) for x in length_list]
  //
  // which gives [49, 73, 85, 91, 94]. These are all the magic numbers, we need to do NTTs in Z_97.
  // We'll keep them around in a table. For practical work, a much larger prime number should be 
  // chosen. I think, ideally, it should be as large as possible to give maximum headroom for the
  // computations, so perhaps the largest prime less than 2^64 is the best choice?

  using uint = unsigned int;

  // Our magic numbers:
  //uint numRoots     = 5;
  //uint maxN         = 32;                      // = 2^5 = 2^numRoots, maximum supported trafo size
  //uint modulus      = 97;
  //uint roots[5]     = { 96, 75, 64, 89, 78 };  // We can pick any, i think,
  //uint rootsInv[5]  = { 96, 22, 47, 12, 51 }; 
  //uint lengthInv[5] = { 49, 73, 85, 91, 94 };

  // Try it with the (first five) numbers from rsModularIntegerNTT_64
  uint numRoots     = 5;
  uint maxN         = 32;                      // = 2^5 = 2^numRoots, maximum supported trafo size
  uint modulus      = 3221225473;
  uint roots[5]     = { 3221225472, 2207278994, 2607818977, 2831384513, 3154294145 };
  uint rootsInv[5]  = { 3221225472, 1013946479, 2190011530, 1626607911,  122509875 }; 
  uint lengthInv[5] = { 1610612737, 2415919105, 2818572289, 3019898881, 3120562177 };


  // Check, if the magic numbers satisfy the requirements:
  using ModInt = RAPT::rsModularInteger<rsUint64>;       // Why do we use rsUnt64 internally? To
  ModInt zero = ModInt(rsUint64(0), rsUint64(modulus));  // increase headroom befor overflow 
  ModInt one  = ModInt(rsUint64(1), rsUint64(modulus));  // occurs? Document!
  ModInt a, b, c;
  bool ok = true;
  for(int i = 0; i < (int)numRoots; i++) {
    int n = rsPowInt(2, i+1);
    a = ModInt(n,            modulus);
    b = ModInt(lengthInv[i], modulus);
    c = a * b;
    ok &= c == one;

    a = ModInt(roots[i],    modulus);
    b = ModInt(rootsInv[i], modulus);
    c = a * b;
    ok &= c == one;

    c = a;
    for(int j = 1; j < n; j++) {
      ok &= c != one;              // Root should be primitive, we should not get 1 for any power
      c *= a;  }                   // ...less than n (i.e. n/2, n/3, etc.)
    ok &= c == one; 
  }
  rsAssert(ok);


  // Create two sequences of integers and convolve them using naive convolution and reduce to the 
  // modulus afterwards. It doesn't matter whether we do the reduction afterwards or after each 
  // operation:
  using VecI = std::vector<int>;
  VecI x  = VecI({11,32,15,75,51});       // input signal
  VecI h  = VecI({65,72,42,28,91,14,32}); // impulse response
  int  Lx = (int)rsSize(x);               // length of x
  int  Lh = (int)rsSize(h);               // length of h
  int  Ly = Lx+Lh-1;                      // length of result
  VecI y(Ly);                             // result (of convolution)
  rsArrayTools::convolve(&x[0], Lx, &h[0], Lh, &y[0]);
  for(int i = 0; i < Ly; i++) 
    y[i] %= modulus;

  // Prepare NTT buffers:
  using VecM = std::vector<ModInt>;
  using LT = rsLinearTransforms;
  int N = 32;                       // FFT/NTT size - todo: determine from Lx,Lh the desired length
  int k = 5-1;                      // log2(N) - 1
  VecM buf1(N), buf2(N);
  for(int i = 0;  i < Lx; i++) buf1[i] = ModInt(x[i], modulus);
  for(int i = Lx; i < N;  i++) buf1[i] = ModInt(0,    modulus);
  for(int i = 0;  i < Lh; i++) buf2[i] = ModInt(h[i], modulus);
  for(int i = Lh; i < N;  i++) buf2[i] = ModInt(0,    modulus);

  // Do NTT-based convolution:
  ModInt WN = ModInt(rootsInv[k], modulus);  // twiddle factor for forward NTT
  LT::fourierRadix2DIF(&buf1[0], N, WN);     // transform buffer 1
  LT::fourierRadix2DIF(&buf2[0], N, WN);     // transform buffer 2
  for(int i = 0; i < N; i++)
    buf1[i] = buf1[i] * buf2[i];             // spectral multiplication, result goes to buffer 1
  WN = ModInt(roots[k], modulus);            // twiddle factor for inverse NTT
  LT::fourierRadix2DIF(&buf1[0], N, WN);     // unscaled inverse NTT
  ModInt S(lengthInv[k], modulus);           // scaler = 1/N
  for(int i = 0; i < N; i++)
    buf1[i] = S * buf1[i];                   // do the scaling
  for(int i = 0;  i < Ly; i++) ok &= (buf1[i] == ModInt(y[i], modulus));  // check result
  for(int i = Ly; i < N;  i++) ok &= (buf1[i] == ModInt(0,    modulus));
  rsAssert(ok);

  VecI y2 = rsConvolveNTT(x, h);
  ok &= y2 == y;
  rsAssert(ok);

  // Observations:
  // -With the sequence lengths Lx = 5, Lh = 7, we need an NTT size of 32 to make it work. 
  //  Apparently, with shorter length, we run into circular convolution. The length of the result
  //  is 7+5-1 = 11. I think, the NTT size must be at least twice the length of the result, so 
  //  16 < 22 is not enough

  // Questions:
  // -Does the modulus really need to be a prime or is it enough, when all desired N-th roots of 
  //  unity, their inverses and the inverses of all desired lengths exist? Being prime is a 
  //  sufficient condition for all of this, but maybe it's not a necessary one?
  // -Do we actually need to store all roots of unity? isn't the highest order root enough, i.e.
  //  the 16th, and the others can be obtained by repeated squaring?

  // Notes:
  // -I think, the interpretation of the NTT coeffs is that they show the correlations of the 
  //  input sequences with different white noise sequences, where the noise is generated by a 
  //  linear congruential generator. The factor in this generator is different for each coeff.
  //  -> verify this: actually create the noise sequences and compute the correlations. Maybe
  //  one noise sequence can be computed from the twice-as-long one by decimating and repeating?
}

void numberTheoreticTrafoModuli()
{
  // We try to find moduli (and related magic numbers) that are suitable for NTT convolution. It 
  // would be nice, if we could use a power-of-2 modulus, because then the modular arithmetic 
  // simplifies (we can use bitmasks), ideally 2^64 (then we don't even need a mask). But that 
  // immediately rules out radix-2 NTTs (i think, the length must be coprime with the modulus or 
  // something). But maybe radix-3, radix-5, radix-7 etc. can be made to work? Their sizes are 
  // coprime with power-of-2 moduli. We investigate some powers of 2 for the modulus and try to 
  // figure out which N = 3^k th roots of unity exist (for radix-3), if their inverses also exist
  // and if the inverses of the 3^k values themsleves exist. That's all that is needed for there to
  // be an NTT of the desired kind, if i'm not mistaken.


  //rsUint64 B = 2;               // base
  //rsUint64 P = 16;              // power
  //rsUint64 M = rsPowInt(B, P);  // candidate modulus, 2^64 wraps around to 0

  rsUint64 M = 3221225473;        // from https://www.programmersought.com/article/8906403181/
  //rsUint64 M = 7681;
  //rsUint64 M = 97;
  rsUint64 R = 2;               // candidate radix

  // Sage can not help, this code:
  //
  // R = 3
  // k = 1
  // M = 2^16
  // Zp = Integers(M)
  // root_list = Zp.zeta(R^k);  # replace the k
  // root_list
  //
  // produces an error message: "NotImplementedError: factorization of polynomials over rings with
  // composite characteristic is not implemented". Maybe, for a preliminary investigation using 
  // small numbers (small P), we can just do an exhaustive search. Later, for bigger and more 
  // realistic P, we may have to look for some software that can handle the job. Maybe Mathematica?

  using ModInt = RAPT::rsModularInteger<rsUint64>;

  auto findRoot = [](rsUint64 N, rsUint64 M, int k) // N: N-th root of unity, M: modulus
  {
    for(rsUint64 x = M-1; x > 0; x--) 
    {
      ModInt yM(x, M);
      int i;
      for(i = 1; i <= k; i++) {
        yM = yM * yM;      // i think, repeated squaring works only for radix R = 2 - in general,
        if(yM.value == 1)  // we must repeatedly do: yM = yM^R
          break;    }
      if(yM.value == 1 && i == k)
        return x;



      //ModInt xM(x, M);
      //ModInt yM(1, M);
      //int i;
      //for(i = 1; i <= N; i++) {
      //  yM = yM * xM;
      //  if(yM.value == 1)
      //    break;    }
      //if(yM.value == 1 && i == N)
      //  return xM.value;
    }
    return 0ULL;
  };

  auto findInverse = [](rsUint64 x, rsUint64 M) // x: number, M: modulus
  {
    rsInt64 tmp = rsModularInverse(rsInt64(x), rsInt64(M));
    if(tmp >= 0) return rsUint64(tmp);
    else         return M - rsUint64(tmp);
  };

  using Vec = std::vector<rsUint64>;
  Vec roots, rootsInv, lengths, lengthsInv;
  rsUint64 N = R;  // N =  R^k

  RSLib::rsFileStream file("MagicNumbers.txt");
  int k = 1;
  while(N < M)
  {
    rsUint64 r = findRoot(N, M, k);
    if(r == 0)
      break;
    if(r != 0)
    {
      rsUint64 ri = findInverse(r, M);
      rsUint64 Ni = findInverse(N, M);

      lengths.push_back(N);
      roots.push_back(r);
      lengthsInv.push_back(Ni);
      rootsInv.push_back(ri);

      bool openSuccess = file.openForAppend();
      if(openSuccess) {
        std::string str;
        str  = "N = "   + to_string(N)  + "   ";
        str += "1/N = " + to_string(Ni) + "   ";
        str += "r = "   + to_string(r)  + "   ";
        str += "1/r = " + to_string(ri) + "\n";
        file.appendText(str);
        file.close(); }
    }
    N *= R;
    rsPrintLine(std::string("Number of roots found: ") + to_string(k));
    k++;
  }
  // This loops finds all the roots but after that, it doesn't stop searching! I think we must 
  // break, when N *= R would lead to overflow.


  int dummy = 0;
  // I think, it would make more sense to run the loop from high to low values, too, because if 
  // a root for a high power is found, i think, we can be sure that the other roots also exist and 
  // we may obtain them by multiplying by some power of R

  // Observations:
  // -For B=2,R=3, it seems that the only root of unity is the trivial one: 1 itself. We don't find
  //  any others. Maybe that means the whole idea is not workable, because the desired roots of
  //  unity do not exist?

  // ToDo:
  // -Print some multiplication tables for general small composite numbers like 4,6,8,9,10,12,...
  //  and try to find rules for when our desired magic numbers may exist. That may help in the 
  //  search for suitable moduli. To this end, maybe implement a general 
  //  rsToString(const rsMatrix&) function (that could be useful for other stuff as well) and 
  //  collect the multiplication tables into rsMatrix
  // -Maybe we can use Mersenne primes of the form 2^n - 1 or 2^n + 1 for the modulus? I think, the
  //  modular arithmetic can still be implemented via bit masking. Then we can perhaps again use 
  //  the radix-2 NTT. The Mersenne prime 2^31 - 1 and 2^61 - 1 could be both useful (for 32 and 64
  //  bit based NTTs https://oeis.org/A000668
  // -Or a Fermat prime of the form 2^n + 1, like 2^32 + 1 = 4294967297 or 
  //  2^64 + 1 = 18446744073709551617 ...but 2^64 + 1 does not fit into 64 bit, even 2^64 doesn't.
  //  But maybe with soem trickery with the overflow flag, we can nevertheless use them? Maybe we
  //  can detect the overflow and repair the result in this case

  // ToDo:
  // -Implement a class template similar to rsModularInteger but with the modulus as template 
  //  parameter, ...also the number of lengths, i.e. numRoots.
  // -Use that in the partinioned convolver - maybe ConvolverFFT and ConvolverNTT (to be written)
  //  should have a common baseclass with common interface
  // -Implement a class for that, storing the modulus and the magic numbers.
  // -

  // The largest prime less than 2^64 is p = 18446744073709551557 and can be found with sage via:
  //   random_prime(2^64, proof=true, lbound=2^64-59)
  // Now, let's find our magic numbers for p ...damn! Z_p seems to have only 2nd and 4th roots of 
  // unity. I thought, primality is a sufficient condition for all N-th roots to exist. Is this not
  // the case?
  // See:
  // http://www.apfloat.org/ntt.html
  // https://www.nayuki.io/page/number-theoretic-transform-integer-dft
  // https://ieeexplore.ieee.org/document/1162555
  // https://ieeexplore.ieee.org/document/1672090
  // https://web.archive.org/web/20130425232048/http://www.cse.psu.edu/~furer/Papers/mult.pdf
  // https://en.wikipedia.org/wiki/Talk%3ANumber-theoretic_transform

  // http://www.faginfamily.net/barry/Papers/Discrete%20Weighted%20Transforms.pdf
  // says: p = 2^61 - 1 is a Mersenne prime and 
  // h = 2147483648 + 1033321771269002680i is primitive root of unity in Z_p...oh! it has a 
  // imaginary part? Is this some sort of complex modular arithmetic?

  // https://www.programmersought.com/article/13432387902/
  // 998244353,1004535809, 469762049, 998244353

  // https://www.programmersought.com/article/8906403181/ has list of primes fot NTT, the largest
  // is 4179340454199820289, which seems to correspond to around 61 bit: giving wolfram alpha 
  // log2(4179340454199820289.0), gives back 61.85798... OK, that seems good - let's try to find 
  // the roots of unity:
  //
  // p = 4179340454199820289 
  // Zp = Integers(p)
  // root_list = Zp.zeta(2^10, all=False);  # replace the 2^10 by 2^11, etc
  // root_list
  //
  //  1: 4179340454199820288
  //  2: 3360066027580426122
  //  3: 3324705732702508476
  //  4: 4093416561646622555
  //  5: 4129893269131444668
  //  6: 4086220048833014884
  //  7: 4075462463776626479
  //  8: 4128470768322469725
  //  9: 4161514758139463127
  // 10: 4178097261067721820
  // 11: 4176218026832679610
  // 12: 4178374021926307362
  // 13: 4177450540047517585
  // 14: 4179324170293557359
  // 15: 4179136626770643812  takes about 2 minutes
  // 16: SageMathCell gave up
  //
  // Oh, oh! The algo to compute them seems to scale exponentially in k. With k=15, it already 
  // takes 2 minutes. And with k=16, SageMathCell did not return anything anymore. There was no 
  // result, but the busy indicator disappeared. Probably, a maximum computation time was 
  // exceeded. Maybe try it with wolfram alpha - perhaps they have a faster algo? Or use a local
  // installation of Sage.

  // Why do so many start with 417? The distribution is interesting - very skewed towards the upper
  // limit - why is that the case? Maybe the algo does just a linear search starting at the highest
  // possible number? Maybe when starting at the high end, the algo is more likely to find a root 
  // more quickly?


  // Table of suitable prime moduli for NTT from here
  // https://www.programmersought.com/article/8906403181/:
  //
  //  g is the original root of mod(r*2^k+1)
  //  Prime number        r   k   g
  //  3                   1   1   2
  //  5                   1   2   2
  //  17                  1   4   3
  //  97                  3   5   5
  //  193                 3   6   5
  //  257                 1   8   3
  //  7681                15  9   17
  //  12289               3   12  11
  //  40961               5   13  3
  //  65537               1   16  3
  //  786433              3   18  10
  //  5767169             11  19  3
  //  7340033             7   20  3
  //  23068673            11  21  3
  //  104857601           25  22  3
  //  167772161           5   25  3
  //  469762049           7   26  3
  //  1004535809          479 21  3
  //  2013265921          15  27  31
  //  2281701377          17  27  3
  //  3221225473          3   30  5
  //  75161927681         35  31  3
  //  77309411329         9   33  7
  //  206158430209        3   36  22
  //  2061584302081       15  37  7
  //  2748779069441       5   39  3
  //  6597069766657       3   41  5
  //  39582418599937      9   42  5
  //  79164837199873      9   43  5
  //  263882790666241     15  44  7
  //  1231453023109121    35  45  3
  //  1337006139375617    19  46  3
  //  3799912185593857    27  47  5
  //  4222124650659841    15  48  19
  //  7881299347898369    7   50  6
  //  31525197391593473   7   52  3
  //  180143985094819841  5   55  6
  //  1945555039024054273 27  56  5
  //  4179340454199820289 29  57  3
  //
  // I'm not sure, what the r,k,g values mean. I think, g may be the primitive (2^k)th root of 
  // unity. ...hmm...dunno...the webpage surrounding that table is a total mess. I think, we want
  // modulus with large k, because length 2^k is the maximum supported NTT length?`..i think, 
  // that's why I picked 3221225473, IIRC. Because we can do NTTs up to length 2^30 with it?
  // Perhaps 65537 = 2^16+1 could be useful as modulus? Are r and g relevant for this application?

  // Some potentially interesting facts given by wolfram alpha:
  //   7340033    = 1063^2     + 2492^2
  //   7340033^2  = 5080095^2  + 5297992^2
  //   23068673   = 2288^2     + 4223^2
  //   23068673^2 = 12598785^2 + 19324448^2
  // Hey! I tried to enter a handful of numbers from the list into wolfram alpha and they all are 
  // the sum of two two squares and the hypothenuse of a Pythagorean triple. Is that a coincidence 
  // or are these features so common that we should expect them? Or do these features have anything 
  // to do with the suitability for NTT? If so, what?
}


void powerIterator()
{
  // We want to iteratively compute an approximatiosn of y(x) = x^p and z(x) = 1/x. The z = 1/x 
  // function is needed in the computation of y = x^p.

  int    N  = 1000;   // number of  data points to generate
  double x0 = 0.2;    // start value for the x-values
  double h  = 0.01;   // step size
  double p  = 0.7;    // the power


  using Vec = std::vector<double>;

  // Compute abscissa values and true values for y(x) and z(x):
  Vec xt(N), yt(N), zt(N);
  for(int n = 0; n < N; n++)
  {
    xt[n] = x0 + n*h;
    zt[n] = 1.0 / xt[n];
    yt[n] = pow(xt[n], p);
  }

  // Compute approximation using a 1st order iterative formula. Computing xn iteratively is 
  // actually trivial and useless, but for completeness the computation has been included here:
  Vec x1(N), y1(N), z1(N);
  double xn, yn, zn; // x[n], y[n], z[n]
  xn = xt[0];
  yn = yt[0];
  zn = zt[0];
  for(int n = 0; n < N; n++) {
    x1[n] = xn; z1[n] = zn; y1[n] = yn;
    xn += h;
    zn -= h*zn*zn;
    yn += h*zn*p*yn; 
  }

  // Alternative way to compute it, that more clearly resembles the textbook formula:
  y1[0] = yt[0];
  z1[0] = zt[0];
  for(int n = 0; n < N-1; n++)
  {
    z1[n+1] = z1[n] + h*(-z1[n]*z1[n]);

    //y1[n+1] = y1[n] + h*(p*z1[n+1]*y1[n]);   // use the new z: works not so well
    y1[n+1] = y1[n] + h*(p*z1[n]*y1[n]);       // use the old z: works actually better

    // Try using a linear combination of old and new z:
    //double c = -0.7;     // coeff for the new z
    //double z = c*z1[n+1] + (1-c)*z1[n];
    //y1[n+1] = y1[n] + h*(p*z*y1[n]);         // yes! this work even better!
  }
  // ToDo: maybe use different arrays for this - maybe z1c, y1c where the c stands for the
  // usage of the c-coefficient


  //rsPlotVectorsXY(x, y, z);

  rsPlotVectorsXY(xt, yt, zt, y1, z1);


  // Observations:
  // -For the 1st order method, we can improve the accuracy of the x^p approximation by using a 
  //  linear combination of the old and new approximate z value, like: zUsed = c*zNew + (1-c)*zOld.
  //  Empirically, a value of c = -0.7 seems to give good results. Maybe c should be -p?

  // ToDo:
  // -Figure out, if using a linear combination of old and new z can also improve the higher order
  //  methods.

  // also interesting:
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=567563&sid=c1cfc136c438d3e75d69d39d3912b84f&start=15
  // https://www.desmos.com/calculator/iceqafamt8

}

void gaussianIterator()
{
  // Note: The technique demonstrated here is now superseded by RAPT::rsExpPolyIterator. The 
  // code here actually just implements what this class does for the special case of N=2.

  // We iteratively produce a Gaussian bell function y(x) = exp(a*x^2) where a < 0 and compare it
  // to a directly computed function. To derive the iterative/recursive rule, we consider:
  //   y(x+h) = exp(a*(x+h)^2) = exp(a*(x^2+2*h*x+h^2)) = exp(a*x^2)*exp(2*a*h*x)*exp(a*h^2)
  //          = y(x) * exp(2*a*h*x) * exp(a*h^2)
  // where the 2nd factor is just a regular decaying exponential, which we shall define as z(x). 
  // We already routinely compute such things recursively and in this case, it can be done by 
  // considering:
  //   z(x)   = exp(2*a*h*x)
  //   z(x+h) = exp(2*a*h*(x+h)) = exp(2*a*h*x) * exp(2*a*h*h) = z(x) * exp(2*a*h*h) = z(x) * d
  // where we have defined d as the decay factor. The 3rd factor is just a constant scale factor 
  // given by:
  //   s = exp(a*h^2)
  // from which we also see that d = s^2, which is also quite convenient.
  // Alternatively to consider y(x+h), we could have considered the derivative y'(x) = 2*a*x*y and
  // solve that approximately via an initial value solver. But since this will only give an 
  // approximate solution, the technique above is better. However, other functions may not admit 
  // such a computationally efficient exact recursion. In these cases, the approach via an ODE may 
  // be useful.


  int    N  = 1000;    // number of data points to generate
  double x0 = -5.0;    // start value for the x-values
  double h  =  0.01;   // step size
  double a  = -0.2;    // factor in front of the x in the exponent

  // Compute abscissa values and true values for y(x):
  using Vec = std::vector<double>;
  Vec x(N), yt(N), zt(N);
  for(int n = 0; n < N; n++)
  {
    x[n]  = x0 + n*h;
    yt[n] = exp(a*x[n]*x[n]);
    zt[n] = exp(2*a*h*x[n]);
  }

  // Compute z(x) and y(x) by exact recursion:
  Vec yr(N), zr(N);
  yr[0] = yt[0]; zr[0] = zt[0];  // initial values copied from directly computed values
  double s = exp(a*h*h);         // scaler
  double d = s*s;                // == exp(2*a*h*h), decay for exponential z(x)
  Vec errAbs(N), errRel(N);      // absolute and relative error
  for(int n = 1; n < N; n++)
  {
    //yr[n] = yr[n-1] * zt[n-1] * s;  // using the exact z(x)
    yr[n] = yr[n-1] * zr[n-1] * s;    // using the recursive z(x)
    zr[n] = zr[n-1] * d;
    errAbs[n] = yr[n] - yt[n];
    errRel[n] = errAbs[n] / yt[n];
  }

  // Compute y(x) by an approximate recursion using the ODE, solved via the forward Euler method:
  Vec ye(N);
  ye[0] = yt[0];
  for(int n = 0; n < N-1; n++)
    ye[n+1] = ye[n] + h * (2*a*x[n]*ye[n]);  // y[n+1] = y[n] + h * y'[n]

  // Plot results:
  rsPlotVectorsXY(x, yt, yr, ye);
  rsPlotVectorsXY(x, errAbs, errRel);

  // Observations:
  // -When using the exact z(x) together with the recursion for y(x), the (accumulated) relative 
  //  error in the last datapoint is around 2.7e-15. When we use the recursively computed z(x), 
  //  it's around 1.7e-11.
  // -With the exact z(x) the errors look quite noisy but with the recursive z(x), the absolute 
  //  error has also a bell shape (with maximum a bit to the right of y(x)'s peak) and the 
  //  relative error increases monotonically in a superlinear way (looks like a parabola).
  // -The solution of the inexact recursion via the forward Euler ODE solver lies entirely below 
  //  the correct function. The error at the peak is around 2.2%. The shape looks OK, though.

  // ToDo:
  // -Figure out, how the error of exact and approximate recursions behave when we choose a smaller
  //  stepsize h. Maybe the error of the exact recursion increases due to more error accumulation 
  //  and that of the approximate recursion decreases due to lower approximation error? ...could 
  //  there even be a break-even point? ...but that would weird. And what about better ODE solver 
  //  formulas such as Adams-Bashforth?
  // -Maybe introduce another parameter "mu" for the center of the peak, so we can trigger it at
  //  t = 0 like an envelope. In this case, we may also want to shift and scale the resulting 
  //  envelope to exactly hit zero at both ends and 1 at the peak.
  // -Can this recursion some be modified to admit for an input as in a filter? Maybe something 
  //  like: zr[n] = zr[n-1] * d + c * in; where c is some coeff that makes the Gaussian the unit 
  //  impulse response? But this is not a  (linear) filter - we from a product of past outputs: 
  //  y[n-1] * z[n-1] - so the response to a sum of inputs would not be the sum of the responses
  //  to the separate inputs...hmmm...
  // -Use it as attack or release shape for an envelope. It think, tacking it to a flat sustain
  //  gives infinite smoothness at the junction.

  // -Try to derive a recursion for the bump function https://en.wikipedia.org/wiki/Bump_function
  //    y(x) = exp(-1/(1-x^2))
}

void expPolyIterator()
{
  // Tests and demostartes the use RAPT::rsExpPolyIterator to iteratively evaluate the exponential
  // function of a polynomial y(x) = exp(p(x))

  // As a generalization of the recursion for the Gaussian y(x) = exp(a*x^2), we now consider
  //   y(x) = exp(p(x)) 
  // where p(x) is a polynomial:
  //   p(x) = a0 + a1*x + a2*x^2 + a3*x^3 + ...
  // we can write y(x) as:
  //   y(x) = exp(a0) * exp(a1*x) * exp(a2*x2) * exp(a3*x^3) * ...
  // Let's consider the exp(a3*x^3) factor, which we shall denote by z3(x)
  //   z3(x)   = exp(a3*x^3) 
  //   z3(x+h) = exp(a3*(x+h)^3) = exp( a3 * (x^3 + 3*x^2*h + 3*x*h^2 + h^3) )
  //           = exp(a3*x^3) * exp(a3*3*h*x^2) * exp(a3*3*h^2*x) * exp(a3*h^3)
  //           = z3(x) * ...
  // we see that the second factor exp(a3*3*h*x^2) also produces an x^2 factor which would have to
  // multiplied by the x^2 factor that already arises from z2(a2*x^2). I think, in general, we get
  // 4 factors and the coefficient in each factor is obtained by a sum of the a-coeffs weighted by
  // binomial coefficients -> todo: work out the details....

  using Real = double;
  using Vec  = std::vector<Real>;
  using Poly = rsPolynomial<Real>;

  int  N  = 1000;    // number of data points to generate
  Real x0 = -5.0;    // start value for the x-values
  Real h  =  0.01;   // step size
  Real k  =  0.1;    // overall scaler for p(x), controls width
  Real a[5];         // polynomial coeffs excluding the k factor
  a[0] =  0.0; 
  a[1] = -0.5;
  a[2] = -0.4;
  a[3] = -0.3;
  a[4] = -0.2;

  // Bake the k-factor into the polynomial coeffs:
  for(int i = 0; i < 5; i++)
    a[i] *= k;

  // Compute x-axis and ground truth pt, yt for polynomial p(x) and for exp(p(x)):
  using Vec = std::vector<double>;
  Vec x(N), pt(N), yt(N);
  //Real p;
  for(int n = 0; n < N; n++)
  {
    x[n]  = x0 + n*h;
    pt[n] = Poly::evaluate(x[n], a, 4);  // p(x)
    yt[n] = exp(pt[n]);                  // y(x), true/target value
  }

  // Compute p(x) and y(x) iteratively:
  rsPolynomialIterator<Real, 4> polyIter;  polyIter.setup(a, h, x0);
  rsExpPolyIterator<   Real, 4> expIter;   expIter.setup( a, h, x0);
  Vec yr(N), pr(N);
  for(int n = 0; n < N; n++)
  {
    pr[n] = polyIter.getValue();
    yr[n] = expIter.getValue();
  }

  // Compute absolute and relative error and plot results:
  Vec errAbs = yr - yt;
  Vec errRel = errAbs / yt;
  rsPlotVectorsXY(x, pt, pr);
  rsPlotVectorsXY(x, yt, yr);
  rsPlotVectorsXY(x, errAbs, errRel);

  // Observations:
  // -The absolute error again resembles the original bump.
  // -The relative error looks like an x^4 function. Maybe the general rule is that the accumulated 
  //  relative error always increases with the same degree as the polynomial?

  // ToDo:
  // -Try to use a cubic rsExpPolyIterator for std::complex to generate a complex sinusoid with 
  //  cubic envelope for the instantaneous phase. Compare result to directly calculated values, in
  //  particular, pay attention to roundoff accumulation

  int dummy = 0;
}


void reciprocalIterator() // rename to odeMultiStep or odeSolvers
{
  // We want to iteratively compute an approximation of y(x) = 1/x. To do this, we find the ODE
  // whose solution is given by our target function y(x). This can be found by differentiating:
  // y' = (x^-1)' = -x^(-2) = -y^2, so our ODE is y' = -y^2. We define f(y) := -y^2 to be 
  // consistent with much of the ODE literature.
  // Over time, this experiment has turned more into one of testing different multistep ODE solver 
  // methods (applied to the given ODE above). Maybe rename the function to multiStepSolvers. Or 
  // maybe reciprocalMultiStep or multistepReciprocal. The code here can serve as reference 
  // prototype for eventually implementing a multistep ODE solver for RAPT.
  
  int    N  = 1000;   // number of data points to generate
  double x0 = 0.2;    // start value for the x-values
  double h  = 0.01;   // step size

  // Compute abscissa values and true values for y(x):
  using Vec = std::vector<double>;
  Vec x(N), yt(N);
  for(int n = 0; n < N; n++) {
    x[n]  = x0 + n*h;
    yt[n] = 1.0 / x[n]; }

  // Some local variables for convenience:
  using Func = std::function<double(double)>;
  Func f  = [](double y) { return -y*y; }; // y' = f(y) = -y^2
  Func fp = [](double y) { return -2*y; }; // derivative of f(y) for Newton iteration
  double* y;                               // shorthand for currently computed approximation of y


  // Now, we compute solution iteratively by various methods...

  // The explicit methods:

  // 1st order forward Euler:
  Vec y1f(N); y = &y1f[0];
  y[0] = yt[0];
  for(int n = 0; n < N-1; n++)
    y[n+1] = y[n] + h * f(y[n]);
  Vec e1f = y1f - yt;

  // 2nd order Nyström:
  Vec y2n(N); y = &y2n[0];
  y[0] = yt[0]; y[1] = yt[1];
  for(int n = 0; n < N-2; n++)
    y[n+2] = y[n] + 2*h * f(y[n+1]);
  Vec e2n = y2n - yt;

  // 3rd order Nyström:
  Vec y3n(N); y = &y3n[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2];
  for(int n = 0; n < N-3; n++)
    y[n+3] = y[n+1] + (h/3) * (7*f(y[n+2]) - 2*f(y[n+1]) + f(y[n]));
  Vec e3n = y3n - yt;

  // 4th order Nyström:
  Vec y4n(N); y = &y4n[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3];
  for(int n = 0; n < N-4; n++)
    y[n+4] = y[n+2] + (h/3) * (8*f(y[n+3]) - 5*f(y[n+2]) + 4*f(y[n+1]) - f(y[n]));
  Vec e4n = y4n - yt;

  // 2nd order Adams-Bashforth:
  Vec y2ab(N); y = &y2ab[0];
  y[0] = yt[0]; y[1] = yt[1];
  for(int n = 0; n < N-2; n++)
    y[n+2] = y[n+1] + (h/2) * (3*f(y[n+1]) - f(y[n]));
  Vec e2ab = y2ab - yt;

  // Form a linear combination of 1st order and 2nd order Adams-Bashforth. This reduces the error 
  // further:
  double c = 0.92; // weight for 2nd order solution
  Vec yMix = c*y2ab + (1-c)*y1f;
  Vec eMix = yMix - yt;
  // todo: move to the plotting code below

  // 3rd order Adams-Bashforth:
  Vec y3ab(N); y = &y3ab[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2];
  for(int n = 0; n < N-3; n++)
    y[n+3] = y[n+2] + (h/12) * (23*f(y[n+2]) - 16*f(y[n+1]) + 5*f(y[n]));
  Vec e3ab = y3ab - yt;

  // 4th order Adams-Bashforth:
  Vec y4ab(N); y = &y4ab[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3];
  for(int n = 0; n < N-4; n++)
    y[n+4] = y[n+3] + (h/24) * (55*f(y[n+3]) - 59*f(y[n+2]) + 37*f(y[n+1]) - 9*f(y[n]));
  Vec e4ab = y4ab - yt;

  // 5th order Adams-Bashforth:
  Vec y5ab(N); y = &y5ab[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3]; y[4] = yt[4];
  for(int n = 0; n < N-5; n++)
    y[n+5] = y[n+4] + (h/720) * (  1901*f(y[n+4]) - 2774*f(y[n+3]) + 2616*f(y[n+2]) 
                                 - 1274*f(y[n+1]) +  251*f(y[n]));
  Vec e5ab = y5ab - yt;

  /*
  // Test the newton iteration (todo: move into a unit test):
  double a = 9;  // number from which we want to extract the sqrt
  Func testFunc  = [&](double x) { return x*x - a; };
  Func testFuncD = [&](double x) { return 2*x; };
  double b = newton(testFunc, testFuncD, 5.0, 0.0); // should compute sqrt(9) = 3, uses 5 as guess
  */


  // Now for the implicit methods. They are more complicated because the value that we want to 
  // compute appears on the left and right hand side of the equation (in general nonlinearly, so we 
  // can't just isolate it on one side). This is solved by setting up the equation as a 
  // root-finding problem at each step and solving it via Newton iteration, taking as initial guess 
  // the previous value.

  // The backward Euler method. It's defined by: y[n+1] = y[n] + h * f(y[n+1]). For convenience, 
  // let's define z := y[n+1]. Now we need to solve y[n] + h*f(z) - z = 0 for z at each step. 
  Func fn, fnp;
  Vec y1b(N); y = &y1b[0];
  y[0] = yt[0];
  fnp = [&](double z) { return h*fp(z) - 1; };
  for(int n = 0; n < N-1; n++)
  {
    fn  = [&](double z) { return y[n] + h*f(z) - z; };
    y[n+1] = newton(fn, fnp, y[n]);
  }
  Vec e1b = y1b - yt;

  // The trapezoidal rule: y[n+1] = y[n] + (h/2) * (f(y[n+1]) + f(y[n])) is an implicit 1-step 
  // method of order 2 (i think - verify!):
  Vec y1t(N); y = &y1t[0];
  y[0] = yt[0];
  fnp = [&](double z) { return (h/2)*fp(z) - 1; };
  for(int n = 0; n < N-1; n++)
  {
    fn = [&](double z) { return y[n] + (h/2)*(f(z)+f(y[n])) - z; };
    y[n+1] = newton(fn, fnp, y[n]);
  }
  Vec e1t = y1t - yt;

  // 2-step Adams-Moulton: y[n+2] = y[n+1] + (h/12) * (5*f(y[n+2]) + 8*f(y[n+1]) - f(y[n]))
  Vec y2am(N); y = &y2am[0];
  y[0] = yt[0]; y[1] = yt[1];
  fnp = [&](double z) { return (5*h/12)*fp(z) - 1; };
  for(int n = 0; n < N-2; n++)
  {
    fn = [&](double z) { return y[n+1] + (h/12)*(5*f(z) + 8*f(y[n+1]) - f(y[n])) - z; };
    y[n+2] = newton(fn, fnp, y[n+1]);
  }
  Vec e2am = y2am - yt;

  // 3-step Adams-Moulton: y[n+3] = y[n+2] + (h/24) * ...
  // ... (9*f(y[n+3]) + 19*f(y[n+2]) - 5*f(y[n+1]) + f(y[n]) )
  Vec y3am(N); y = &y3am[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2];
  fnp = [&](double z) { return (9*h/24)*fp(z) - 1; };
  for(int n = 0; n < N-3; n++)
  {
    fn = [&](double z) { return 
      y[n+2] + (h/24)*(9*f(z) + 19*f(y[n+2]) - 5*f(y[n+1]) + f(y[n])) - z; };
    y[n+3] = newton(fn, fnp, y[n+2]);
  }
  Vec e3am = y3am - yt;

  // 4-step Adams-Moulton: y[n+4] = y[n+3] + (h/720) * ...
  // ... ( 251*f(y[n+4]) + 646*f(y[n+3]) - 264*f(y[n+2]) + 106*f(y[n+1]) - 19*f(y[n]) ):
  Vec y4am(N); y = &y4am[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3];
  fnp = [&](double z) { return (251*h/720)*fp(z) - 1; };
  for(int n = 0; n < N-4; n++)
  {
    fn = [&](double z) { return y[n+3] 
      + (h/720)*(251*f(z) + 646*f(y[n+3]) - 264*f(y[n+2]) + 106*f(y[n+1]) - 19*f(y[n])) - z; };
    y[n+4] = newton(fn, fnp, y[n+3]);
  }
  Vec e4am = y4am - yt;

  // 2-step Milne-Simpson: y[n+2] = y[n] + (h/3) * (f(y[n+2]) + 4*f(y[n+1]) + f(y[n])):
  Vec y2ms(N); y = &y2ms[0];
  y[0] = yt[0]; y[1] = yt[1];
  fnp = [&](double z) { return (1*h/3)*fp(z) - 1; };
  for(int n = 0; n < N-2; n++)
  {
    fn = [&](double z) { return y[n] + (h/3)*(f(z) + 4*f(y[n+1]) + f(y[n])) - z; };
    y[n+2] = newton(fn, fnp, y[n+1]);
  }
  Vec e2ms = y2ms - yt;

  // 4-step Milne-Simpson: y[n+4] = y[n+2] + (h/90) * ... 
  // ... ( 29*f(y[n+4]) + 124*f(y[n+3]) + 24*f(y[n+2]) + 4*f(y[n+1]) - f(y[n]) ):
  Vec y4ms(N); y = &y4ms[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3];
  fnp = [&](double z) { return (29*h/90)*fp(z) - 1; };
  for(int n = 0; n < N-4; n++)
  {
    fn = [&](double z) { return y[n+2] 
      + (h/90)*(29*f(z) + 124*f(y[n+3]) + 24*f(y[n+2]) + 4*f(y[n+1]) - f(y[n])) - z; };
    y[n+4] = newton(fn, fnp, y[n+3]);
  }
  Vec e4ms = y4ms - yt;

  // 2-step BDF: y[n+2] = -(1/3)*y[n] + (4/3)*y[n+1] + (2/3)*h*f(y[n+2]), 1-step BDF would be 
  // backward Euler again, so we don't include it here.
  Vec y2bdf(N); y = &y2bdf[0];
  y[0] = yt[0]; y[1] = yt[1];
  fnp = [&](double z) { return (2*h/3)*fp(z) - 1; };
  for(int n = 0; n < N-2; n++)
  {
    fn = [&](double z) { return (-1./3)*y[n] + (4./3)*y[n+1] + (2./3)*h*f(z) - z; };
    y[n+2] = newton(fn, fnp, y[n+1]);
  }
  Vec e2bdf = y2bdf - yt;

  // 3-step BDF: y[n+3] = (18/11)*y[n+2] - (9/11)*y[n+1] + (2/11)*y[n] + (6/11)*h*f(y[n+3]):
  Vec y3bdf(N); y = &y3bdf[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2];
  fnp = [&](double z) { return (6*h/11)*fp(z) - 1; };
  for(int n = 0; n < N-3; n++)
  {
    fn = [&](double z) { return (18./11)*y[n+2] - (9./11)*y[n+1] + (2./11)*y[n] 
      + (6./11)*h*f(z) - z; };
    y[n+3] = newton(fn, fnp, y[n+2]);
  }
  Vec e3bdf = y3bdf - yt;
  // todo: factor out the 1/11

  // 4-step BDF: y[n+4] =   (48/25)*y[n+3] - (36/25)*y[n+2] + (16/25)*y[n+1] - (3/25)*y[n] 
  //                      + (12/25)*h*f(y[n+4])
  Vec y4bdf(N); y = &y4bdf[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3];
  fnp = [&](double z) { return (12*h/25)*fp(z) - 1; };
  for(int n = 0; n < N-4; n++)
  {
    fn = [&](double z) { return (48./25)*y[n+3] - (36./25)*y[n+2] + (16./25)*y[n+1] 
      - (3./25)*y[n] + (12./25)*h*f(z) - z; };
    y[n+4] = newton(fn, fnp, y[n+3]);
  }
  Vec e4bdf = y4bdf - yt;

  // 5-step BDF: y[n+5] =  (300/137)*y[n+4] - (300/137)*y[n+3] + (200/137)*y[n+2] - (75/137)*y[n+1] 
  //                     + (12/137)*y[n] + (60/137)*h*f(y[n+4])
  Vec y5bdf(N); y = &y5bdf[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3]; y[4] = yt[4];
  fnp = [&](double z) { return (60*h/137)*fp(z) - 1; };
  for(int n = 0; n < N-5; n++)
  {
    fn = [&](double z) { return (300./137)*y[n+4] - (300./137)*y[n+3] + (200./137)*y[n+2] 
      - (75./137)*y[n+1] + (12./137)*y[n] + (60./137)*h*f(z) - z; };
    y[n+5] = newton(fn, fnp, y[n+4]);
  }
  Vec e5bdf = y5bdf - yt;


  // todo: use simplified notation: y5 instead of y[n+5], etc. apply this to the comments above 

  // 6-step BDF: y6 = (1/147) * (360*y5 - 450*y4 + 400*y3 - 225*y2 + 72*y1 - 10*y0 + 60*h*f(y6))
  Vec y6bdf(N); y = &y6bdf[0];
  y[0] = yt[0]; y[1] = yt[1]; y[2] = yt[2]; y[3] = yt[3]; y[4] = yt[4]; y[5] = yt[5];
  fnp = [&](double z) { return (60*h/147)*fp(z) - 1; };
  for(int n = 0; n < N-6; n++)
  {
    fn = [&](double z) { return (1./147) * (360*y[n+5] - 450*y[n+4] + 400*y[n+3] - 225*y[n+2] 
      + 72*y[n+1] - 10*y[n+0] + 60*h*f(z)) - z; };
    y[n+6] = newton(fn, fnp, y[n+5]); 
  }
  Vec e6bdf = y6bdf - yt;



  // Results and errors of 1st order forward and backward Euler and trapezoidal methods:
  rsPlotVectorsXY(x, yt, y1f, y1b, y1t);
  rsPlotVectorsXY(x, e1f, e1b, e1t);

  // Results and error of 1st and 2nd order Adams-Bashforth method and error of the mix between 1st
  // and 2nd order method:
  rsPlotVectorsXY(x, yt, y1f, y2ab);
  rsPlotVectorsXY(x, e1f, e2ab, eMix);

  // Errors of Adams-Bashforth methods of various orders:
  rsPlotVectorsXY(x, e1f, e2ab, e3ab, e4ab, e5ab);

  // Errors of Adams-Moulton methods of various orders:
  rsPlotVectorsXY(x, e1b, e1t, e2am, e3am, e4am);

  // Errors of BDF methods:
  rsPlotVectorsXY(x, e1b, e2bdf, e3bdf, e4bdf, e5bdf, e6bdf);

  // Errors of Adams-Bashforth and Adams-Moulton methods of various orders:
  rsPlotVectorsXY(x, e4ab, e2am, e5ab, e3am);

  // Errors of Adams-Moulton and BDF methods of various orders:
  rsPlotVectorsXY(x, e2am, e4bdf, e3am, e5bdf, e4am, e6bdf);

  // Errors of Adams-Moulton and Milne-Simpson for 2-step and 4-step
  rsPlotVectorsXY(x, e2am, e2ms, e4am, e4ms);

  // Results and error of 1st and 2nd order Nyström method:
  rsPlotVectorsXY(x, yt, y1f, y2n);
  rsPlotVectorsXY(x, e1f, e2n);

  // Results and error of 2nd and 3rd order Nyström method:
  rsPlotArraysXY(800, &x[0], &yt[0], &y2n[0], &y3n[0]);
  rsPlotArraysXY(800, &x[0], &e2n[0], &e3n[0]);

  // Results and error of 3rd and 4th order Nyström method:
  rsPlotArraysXY(350, &x[0], &yt[0], &y3n[0], &y4n[0]);
  rsPlotArraysXY(350, &x[0], &e3n[0], &e4n[0]);


  // Observations:
  // -The 2nd order Nyström method for z(x) = 1/x is initially more accurate than the 1st order 
  //  method but later builds up unstable oscillations. Also, towards the later section, even the
  //  mean value of the alternative Nyström solution detoriates more from the true solution than 
  //  the 1st order (forward Euler) solution.
  // -The 3rd order Nyström method explodes even faster, even though the oscillations get dampened
  //  out over time. They increase initially, but slower than in the 2nd order method, but later
  //  an non-alternating explosion takes over.
  // -The 2nd order Adams-Bashforth method seems to be stable and reduces the by an order of 
  //  magnitude compared to the 1st order method. However, the error can be further reduced, by 
  //  taking a linear combination of the results of 1st and 2nd order solutions because the error
  //  has a similar shape but opposite signs. Taking about 0.92 of the 2nd order and 0.08 = 1-0.92
  //  of the 1st order solution gives a very accurate result.
  // -Higher order Adams-Bashforth methods do also seem to be stable for this problem and with each
  //  order increase, the error goes down by some factor of around 5.
  // -The 1st order forward and backward methods produce errors of similar shape but opposite sign
  //  which suggests to use a mix of both to increase accuracy. But maybe such a mix yields the 
  //  trapezoidal method?
  // -The accuracy of the n-step Adams-Moulton method is better than an (n+2)-step Adams-Bashforth
  //  method.
  // -The 2- and 4-step Milne-Simpson methods also oscillate, the 2-step method more so than the
  //  4-step method. This is different from the Nyström methods, where higher order methods 
  //  oscillate more than the lower order ones.
  // -What Nyström and Milne-Simpson have in common is that the use the value 2 samples past as the
  //  starting point: to compute y[n], they do somthing like y[n-2] + ... rather than y[n-1] + ...
  //  as the other methods do. Both show oscillations, so oscillation my be related to that 
  //  feature. Q: Is this like the "leapfrog" method? -> figure out relation to that
  // -Seems like the BDF-n method have an accuracy comparable to Adams-Moulton methods of order 
  //  n-2. Apparently the better stability is payed for by less accuracy, so maybe we should resort
  //  to BDF only in case of stiff equations
  // -Generally, the absolute error follows an up-down shape not quite unlike the attack/decay 
  //  envelope. I guess, the initial increase is what should be expected (the error starts at zero 
  //  and then builds up) and maybe the later decay can be explained by the fact that the solution 
  //  itself has decaying behavior. Maybe we should plot the relative error instead. This will 
  //  probably show a steady increase...

  // Conclusions:
  // -Adams-Bashforth methods seem to be well suited for this problem. 
  // -The implicit Adams-Moulton method is more costly per computed value than the explicit 
  //  Adams-Bashforth method but make up for it by being more accurate so one may be able to choose
  //  a larger step size. I think, they also have better stability (verify!).
  // -Nyström and Milne-Simpson methods are not suitable for this problem.

  // ToDo:
  // -Implement implicit solver schemes (Adams-Moulton (done), Milne-Simpson, BDF)
  // -compare errors of Adams-Bashforth and Adams-Moulton
  // -Maybe plot the relative error.
  // -Figure out, if it's possible to stabilize the Nyströms methods. Maybe by some sort of
  //  2-point averaging filter?
  // -Implement Adams-Bashforth of orders 2..5 (formulas on wikipedia:
  //  https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods)
  // -Try to figure out general formulas for the coefficients of arbitrary order Nyström and 
  //  Adams-Bashforth methods. Maybe implement the formula given on wikipedia in Sage.
  // -Implement a multistep ODE solver. It should take a couple of initial values that the client
  //  code supplies. Client can choose to compute them directly of via 1-step method, such as 
  //  Runge-Kutta. Maybe implement also an implicit solver using Adams-Moulton and Milne-Simplson
  //  (see Numerik, p 322) and BDF methods (backward differentiation formula), see:
  //  https://en.wikipedia.org/wiki/Backward_differentiation_formula
  //  The explicit solvers can be implemented with just one evaluation of f per generated value by
  //  storing the past f-values from previous calls. Implicit solvers may use Newton iteration 
  //  using an initial guess produced by an explicit method. In a practical implementation, we 
  //  would not re-evaluate f so many times as in, for example:
  //    y[n+4] = y[n+3] + (h/24) * (55*f(y[n+3]) - 59*f(y[n+2]) + 37*f(y[n+1]) - 9*f(y[n]));
  //  (taken from 4-step Adams-Bashforth). Instead, we would store the previously computed f values
  //  together with the previously computed y values (the f values are also known as y' values).
  //  Maybe we can use a lower accuracy explicit Adams-Bashforth predictor together with a higher 
  //  accuracy implicit Adams-Moulton or BDF corrector and use the absolute difference (divided by 
  //  the corrected value) as error estimate to drive the step-size adaption. For implementing a 
  //  predictor/corrector scheme, see:
  //  https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods
  //  https://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method
  // -Implement a function that computes the coefficients for arbitrary order Adams-Bashforth
  //  methods (done)
  // -It seems Adams methods use a lot of past y' (or f) values and only one past y-value whereas
  //  BDF method use a lot of past y-values but only one past f value. What about methods that use
  //  a bunch of both types? wouldn't that give yet better accuray? ah - see 
  //  "A First Course in the Numerical Analysis of Differential Equations" pg 25 ff. order of 
  //  accuracy is not the only relevant criterion
  // -Try to apply the technique to other interesting functions such as: Gaussian, 1/(1+x^2n), tan,
  //  tanh, atan, exp(-0.1*x) + exp(-10*x) (stiff?)
  //    y(x) = 1/(1+x^(2*n))  ->  y' = -y^2 * 2*n * x^(2*n-1)
  //    y(x) = exp(-a*x^2)    ->  y' = -2*a*x*y
  //  these require our function f(y) = y' to also take x as argument. Maybe the general ODE 
  //  solver should have an API that lets the user specify the function as f(x,y). Strictly 
  //  speaking, if y can a vector, that would not be needed because we could add a 0th component
  //  that just uses the identity function and has derivative 1. But i think, the API is more
  //  convenient, if we treat the independent variable separately.


  // Ideas:
  // -Maybe an adaptive stepsize control can be implemented by letting the stepsize grow or 
  //  shrink by a factor of 2. Let's say, we have a 4-step method. In the case of shrinking the 
  //  stepsize, we need a memory of y[n],..,y[n-4]. We pass a 4th degree polynomial to our most 
  //  recent 4 datapoints and evaluate it at values halfway between our knwon values to produce
  //  the denser history, i.e. we produce y[n-0.5], y[n-1.5], etc. and use these (interleaved with
  //  our original history data) as new most recent 4 points - thus we have halved the stepsize. In
  //  case of growing (doubling) the stepsize, we just discard y[n-1],y[n-3],y[n-5],y[n-7]...yes, 
  //  we need to keep a history twice as long as the current stepsize dictates and if we have grown 
  //  the stepsize in one step, we can't grow it again in the very next step because we don't have 
  //  enough history. So, growth is more limited than shrinkage, which is better than the other way 
  //  around, Maybe in case of growth we should also filter the past outputs witha 2-point average.
  //  ..and maybe do the growth only if the filtered values are close enough to the unfiltered ones
  // -Maybe instead of modeling the solution as a polynomial that passes through the past points, 
//    model it as a sum of exponentials: y(x) = sum_i a_i * exp(b_i * x) where i = 1,..,N. The 
//    initial values y_0, y_0' can be used to determine the coeffs if N=1, for N=2, we'll also need
//    y_-1, y_-1', etc.


  // See also:
  // https://www.researchgate.net/publication/257143339_Construction_of_Improved_Runge-Kutta_Nystrom_Method_for_Solving_Second-Order_Ordinary_Differential_Equations

  // Sundials ODE library:
  // https://computing.llnl.gov/projects/sundials/sundials-software
  // https://github.com/LLNL/sundials
}


// fun stuff:

//=================================================================================================

rsGroupString add(rsGroupString A, rsGroupString B)
{
  return A + B;
  // associative, not distributive, also, it would be weird to have addition and multipication 
  // doing the same thing (although, it's not explicitly forbidden, i think)
}
rsGroupString mul1(rsGroupString A, rsGroupString B)
{
  return (-A) + (-B); // reverse and add
  // neither associative nor distributive :-(
}
rsGroupString mul2(rsGroupString A, rsGroupString B)
{
  return -((-A) + (-B));  // reverse, add, reverse result
  // maybe associative, not distributive
}
// other ideas
// -ab*cde = (a+c)(a+d)(a+e)(b+c)(b+d)(b+e) where + is modular addition of the char-numbers, so we 
//  have 6=2*3 letters and each is given by modular addition
// -compute the (modular) sum of the letters of the first input and mod-add them to each letter in
//  the 2nd input
// -mod-add element--wise...but what if inputs have different length?
// -some sort of "convolution"

rsGroupString mul3(rsGroupString A, rsGroupString B)
{
  rsGroupString C = A;
  for(int i = 0; i < B.length(); i++) {
    if(C.length() > 0 && C[C.length()-1] != B[i])
      C.append(B[i]);
    else if(C.length() > 0)     // avoid popping on empty vector
      C.removeLast();
  }
  return C;
}
// is the same as the current "+" operator of the class - but is probably better to use as mul
// maybe we should give the operations neutral names here, such as op1, op2, etc. or name them
// what they do - concatDel, modAdd


rsGroupString add2(rsGroupString A, rsGroupString B, unsigned int m) // m: modulus
{
  // do modular addition of of the characters where the shorter string is appropriately zero-padded
  // (zero standing for a "blank" character). if the result contains trailing 0s, these will be 
  // removed. the inverse element is obtained by taking charcter-wise m-x where x is the input char
  // and m is the modulus
  // could this be distributive with this concatenation and delete pairs thing that we use as the
  // other addition? if so, we could switch to using concatenation as multiplication

  //std::vector<unsigned int> s1 = A.get(), s1 = B.get();

  int LA = A.length();
  int LB = B.length();
  //int m  = 5;              // modulus

  int i;                   // loop index 
  rsGroupString C;
  if(LA >= LB) {           // A is longer than B or has same length
    C.resize(LA);
    for(i = 0; i < LB; i++)   C[i] = (A[i] + B[i]) % m;
    for(i = LB; i < LA; i++)  C[i] =  A[i];
  }
  else {                   // A is shorter than B
    C.resize(LB);
    for(i = 0; i < LA; i++)   C[i] = (A[i] + B[i]) % m;
    for(i = LA; i < LB; i++)  C[i] =  B[i];
  }
  while(C.last() == 0)  C.removeLast();  // remove trailing zeros
  return C;
}
rsGroupString add2(rsGroupString A, rsGroupString B) // to make compatible with isAsso..., etc.
{
  return add2(A, B, 7);  // uses fixed modulus
}


// tests associativity of the given operation for the given triple of arguments
bool isAssociative(rsGroupString (*op) (rsGroupString A, rsGroupString B),
  rsGroupString A, rsGroupString B, rsGroupString C)
{
  rsGroupString s1 = op(op(A, B), C); // (A+B)+C
  rsGroupString s2 = op(A, op(B, C)); // A+(B+C)
  return s1 == s2; 
}

// tests, if the given "mul" operation is distributive over the given "add" operation for the given
// triple of input arguments
bool isDistributive(
  rsGroupString (*add) (rsGroupString A, rsGroupString B),
  rsGroupString (*mul) (rsGroupString A, rsGroupString B),
  rsGroupString A, rsGroupString B, rsGroupString C)
{
  rsGroupString s1 = mul(A,add(B,C));         // A*(B+C)
  rsGroupString s2 = add(mul(A,B), mul(A,C)); // A*B + A*C
  return s1 == s2; 
}

bool testStringMultiplication(rsGroupString (*mul) (rsGroupString A, rsGroupString B))
{
  //bool r = true;
  typedef rsGroupString2 GS;

  // test associativity:
  GS abc("abc"), cde("cde"), efg("efg");
  GS A = abc, B = cde, C = efg;

  std::string t1, t2;

  // do this in a loop with various A, B, C (maybe random or by systematically checking all 
  // possible strings up to a given length) - it should be a function that automatically
  // checks a lot of strings and also takes the operations as inputs (as pointers), so we can try
  // various things


  bool asso = true;
  GS s1 = mul(mul(A, B), C); // (A*B)*C
  GS s2 = mul(A, mul(B, C)); // A*(B*C)
  t1 = s1.toString();        // for inspection/debugging
  t2 = s2.toString();
  asso &= s1 == s2; 

  // distributivity:
  bool dist = true;
  s1 = mul(A,(B+C));        // A*(B+C)
  s2 = mul(A,B) + mul(A,C); // A*B + A*C
  t1 = s1.toString();       // for inspection/debugging
  t2 = s2.toString();
  dist &= s1 == s2;
  // should also use a loop over many strings

  // maybe try all strings of length up to 4 from the alphabet a,b,c,d - each of the 3 inputs
  // A,B,C should take on all possible values - so we need a doubly nested loop


  return asso && dist;
}

bool groupString()
{
  bool r = true;  // test result - todo: turn into unit test
  typedef rsGroupString2 GS;
  GS abc("abc"), cde("cde"), efg("efg");
  GS abde = abc + cde;    r &= abde == "abde"; 
  GS edba = -abde;        r &= edba == "edba";
  GS empty = abde - abde; r &= empty == "";

  // test associativity of addition:
  GS s1 = (abc + cde) + efg;
  GS s2 = abc + (cde  + efg);
  r &= s1 == s2;  // of course, this is only an example - do more random tests

  // test neutral element of addition:
  s1 = abc + GS(""); r &= s1 == "abc";
  s1 = GS("") + abc; r &= s1 == "abc";


  // test multiplication (with the various candidate rules):
  //r &= testStringMultiplication(&add); // asso, not distri
  //r &= testStringMultiplication(&mul1); // not asso, not distri - bad!
  r &= testStringMultiplication(&mul2);


  typedef std::vector<unsigned int> Vec;

  unsigned int m = 7; // the modulus
  rsGroupString s2314 = (Vec({2, 3, 1, 4}));
  rsGroupString s546  = (Vec({5, 4, 6}));
  rsGroupString s3613 = (Vec({3, 6, 1, 3}));
  rsGroupString t1, t2, t3;     // temporaries
  t1 = add2(s2314, s546,  m);   // 0004
  t2 = add2(s546, s2314,  m);   // 0004 -> commutative
  t1 = add2(add2(s2314, s546, m), s3613, m);  // 361
  t2 = add2(s2314, add2(s546, s3613, m), m);  // 361 -> associative
  // OK - this looks good - addition is commutative and associative - at least for the tried 
  // examples - todo: try more examples - if it works out, try inverse elements and then 
  // distributivity with the concat-delete operation ...hmm - this addition actually allows
  // pairs of equal characters - which is a good thing

  bool asso = isAssociative(&add2, s2314, s3613, s546);
  bool distri = isDistributive(&add2, &mul3, s2314, s3613, s546);



  // test multiplication: define various candidate multiplication functions
  // rsGroupString mul1(rsGroupString A, rsGroupString B), mul2, etc. and make a function 
  // testMultiplication passing a function pointer to one of these candidates. Inside this 
  // function, test associativity, distributivity

  // -maybe first try to get a ring with a couple of candidate multiplication rules and then try to 
  //  find among the candiates a rule that turns the ring into a field
  // -if no suitable multiplication rule between strings can be found, mybe try a multiplication 
  //  rule between strings and single characters - some sort of "scalar-multiplication" just like 
  //  in a vector space - for example modular-addition of the scalar value to all chars in the 
  //  string ('a' is the the neutral element and the is no null element)

  return r;
}

//=================================================================================================

void primeAlternatingSums()
{
  // Take the array of prime numbers with alternating signs
  // 2 -3 5 -7 11 -13 17 -19 ...
  // and take running sums of various order ...just for fun to see what happens
  // ...do they also all have alternating signs? what happens, if we add them to one of the
  // previous arrays...just mess around a little - may interesting patterns emerge...
  // -what changes, if we give the even-indexed primes a negative sign?
  // -what happens, if we flip every k-th sign instead of every 2nd?
  // -what happens, if we filter the arrays with integer-coefficient IIR filters? this is a 
  //  generalization of (iterated) running sums (i think)
  // -what about nonlinear and/or time-variant recursive filters with coeffs depending on index?
  // -what about trying to use non/linear prediction methods on the time series?
  //  -maybe the equations (correlation coeffs, etc) should use rational numbers?
  //  -test on number sequences with known structure first (like fibonacci numbers), see, if they
  //   pick it up correctly
  // -what if we sign invert not every second number but every second pair of numbers? ..and then 
  //  combine the resulting sequences with those obtained from inverting every other number ...then
  //  generalize: invert every triple, quadruple, etc...

  int N = 200; // number of primes

  //typedef std::vector<int> IntVec;
  typedef RAPT::rsArrayTools AR;

  // create array of primes:
  //std::vector<int> primes(N);
  std::vector<int> primes(N);
  RAPT::rsFillPrimeTable(&primes[0], N);  // make convenience function getPrimes

  // create arrays of primes with alternating signs:
  int n;
  std::vector<int> ev, od; // ev: have even-numberd signs flipped, od: odd numbered signs flipped
  ev = primes; od = primes;  
  for(n = 0; n < N; n += 2) ev[n] = -ev[n];
  for(n = 1; n < N; n += 2) od[n] = -od[n];

  // 1st order running sums:
  std::vector<int> ev1(N), od1(N);
  AR::cumulativeSum(&ev[0], &ev1[0], N);
  AR::cumulativeSum(&od[0], &od1[0], N);

  // 2nd order running sums:
  std::vector<int> ev2(N), od2(N);
  AR::cumulativeSum(&ev1[0], &ev2[0], N);
  AR::cumulativeSum(&od1[0], &od2[0], N);

  // 3rd order running sums:
  std::vector<int> ev3(N), od3(N);
  AR::cumulativeSum(&ev2[0], &ev3[0], N);
  AR::cumulativeSum(&od2[0], &od3[0], N);

  // ...at this point, we really should use a loop



  // plot:
  GNUPlotter plt;
  plt.setPixelSize(1000, 500);
  //plt.addDataArrays(N, &primes[0]);
  //plt.addDataArrays(N, &ev[0]);
  plt.addDataArrays(N, &od[0]);
  //plt.addDataArrays(N, &ev1[0]);
  plt.addDataArrays(N, &od1[0]);
  //plt.addDataArrays(N, &ev2[0]);
  plt.addDataArrays(N, &od2[0]);
  //plt.addDataArrays(N, &ev3[0]);
  plt.addDataArrays(N, &od3[0]);
  plt.plot();
  // todo: maybe plot with stems or at least, get rid of the linear interpolation (draw steps)

  // the 3rd order sum is very smooth ... subtract even and odd 3rd order sums...

  // ..but maybe it's actually faster to noodle around with that stuff in python or sage
}

// naive implementation - can probably be optimized, if we know the prime factorization of n
// see: https://www.quora.com/What-is-the-fastest-way-to-find-the-divisors-of-a-number#
void rsFindNonTrivialDivisors(rsUint32 n, std::vector<rsUint32>& d)
{
  d.clear();
  if(n == 0) return; // or should we say that 0 is divisible by any number?
  for(rsUint32 i = 2; i <= rsIntSqrt(n); i++)
    if(n % i == 0) {
      d.push_back(i);
      rsUint32 j = n/i;
      if(j != i)
        d.push_back(j);
      // todo: optimize by using a divmod operation: divmod(n, j, i)
    }
  //RAPT::rsHeapSort(&d[0], (int)d.size());
  std::sort(d.begin(), d.end());
  //d.push_back(n);
}
// we do not add the trivial divisors 1 and n to the array...or should we?
// the definition of divisors would include them: https://en.wikipedia.org/wiki/Divisor
// ...unless we qualify them as nontrivial

// computes the number of divisors of a number from the exponents of its prime-factorization
// which must be passed as argument.
rsUint32 rsNumDivisors(std::vector<rsUint32>& exponents)
{
  rsUint32 nd = 1;
  for(size_t i = 0; i < exponents.size(); i++)
    nd *= exponents[i] + 1;
  return nd;
}

class rsNumberDivisibilityInfo
{
public:
  rsNumberDivisibilityInfo(rsUint32 number)
  {
    this->number = number;
    RAPT::rsPrimeFactors(number, factors, exponents);
    rsFindNonTrivialDivisors(number, divisors);
  }

  rsUint32 number;
  std::vector<rsUint32> factors;
  std::vector<rsUint32> exponents;
  std::vector<rsUint32> divisors;

protected:

};

void divisibility()
{
  rsUint32 max = 5040;
  std::vector<rsNumberDivisibilityInfo> numInfos;
  numInfos.reserve(max+1);
  for(rsUint32 i = 0; i <= max; i++)
    numInfos.push_back(rsNumberDivisibilityInfo(i));

  // todo: find highly composite and largely composite numbers...they can be useful for GUI sizes

  // plot the number of non-trivial divisors as function of the number itself:
  std::vector<int> numDivisors(numInfos.size());
  for(size_t i = 0; i < numInfos.size(); i++)
    numDivisors[i] = (int) numInfos[i].divisors.size();

  GNUPlotter plt;
  plt.addDataArrays((int)numDivisors.size(), &numDivisors[0]);
  plt.plot();
  // todo: mark highly composite numbers with a big mark and largely composite numbers with a 
  // smaller mark
}


// computes the "arithmetic derivative" of given natural number (todo: generalize to integers and
// rationals)
rsUint32 numDeriv(rsUint32 n)  // number derivative
{
  if(n == 0 || n == 1)
    return 0;

  std::vector<rsUint32> p, e; // prime factors and exponents
  RAPT::rsPrimeFactors(n, p, e);

  double s = 0; // sum of fractions - use a rational number class later
  for(size_t i = 0; i < p.size(); i++)
    s += (double)e[i] / (double)p[i];

  double ns = n*s;
  rsUint32 d = (rsUint32) round(ns); // the derivative

  return d;
}
void arithmeticDerivative()
{
  int N = 100;
  std::vector<rsUint32> x(N), d(N); // numbers and their derivatives
  for(int i = 0; i < N; i++) {
    x[i] = i;
    d[i] = numDeriv(i);
  }
  //int dummy = 0;
}
/*
https://oeis.org/A003415/list
[0,0,1,1,4,1,5,1,12,6,7,1,16,1,9,8,32,1,21,1,24,
10,13,1,44,10,15,27,32,1,31,1,80,14,19,12,60,1,21,
16,68,1,41,1,48,39,25,1,112,14,45,20,56,1,81,16,
92,22,31,1,92,1,33,51,192,18,61,1,72,26,59,1,156,
1,39,55,80,18,71]
*/
// https://en.wikipedia.org/wiki/Arithmetic_derivative
// https://web.archive.org/web/20050426071741/http://web.mit.edu/lwest/www/intmain.pdf
// http://oeis.org/wiki/Arithmetic_derivative
// https://www.youtube.com/watch?v=lh1qbp7LpFU
//
// How to Differentiate a Number:
// https://cs.uwaterloo.ca/journals/JIS/VOL6/Ufnarovski/ufnarovski.pdf - has generalization to irrationals
//
// The Arithmetic Derivative and Antiderivative:
// https://cs.uwaterloo.ca/journals/JIS/VOL15/Kovic/kovic4.pdf
