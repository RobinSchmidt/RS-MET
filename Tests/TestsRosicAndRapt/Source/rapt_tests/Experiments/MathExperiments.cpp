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







void eigenstuff()
{
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
};


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



void interpolatingFunction()
{
  typedef RAPT::rsInterpolatingFunction<float, double> IF;

  // create data to interpolate:
  int N = 5;
  float  x[5] = { 2, 4, 5, 7, 8 };
  double y[5] = { 1, 3, 1, 2, 3 };

  // create and set up interpolating function object:
  IF intFunc;
  //intFunc.setMode(IF::LINEAR);
  intFunc.setMode(IF::CUBIC_HERMITE);
  intFunc.setPreMap( &log);
  intFunc.setPostMap(&exp);
  intFunc.setPreMap( nullptr);
  intFunc.setPostMap(nullptr);


  // do extra/interpolation:
  static const int M = 500; // number of interpolated values
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
  // impressive...i think, something is still wrong
  // optimize also with scipy and compare results

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
  double mx  = RAPT::rsArrayTools::mean(x, N); // mean of x
  double sum = 0;
  for(int n = 0; n < N; n++) {
    double d = x[n] - mx;
    sum += d*d;
  }
  return sum / (N-1);
  // todo: make an optimized version that takes the mean as argument
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

  // initial approximation for Halley iteration:
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

