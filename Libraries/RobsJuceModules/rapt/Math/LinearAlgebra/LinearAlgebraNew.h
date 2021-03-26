#ifndef RAPT_LINEARALGEBRA_NEW_H_INCLUDED
#define RAPT_LINEARALGEBRA_NEW_H_INCLUDED

/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc. 

Many algorithms destroy the original input matrices and/or right-hand-sides in the course of their 
operation. If a parameter like a pointer to an array or a reference to an rsMatrixView is passed as
non-const, client code should always assume that the argument will contain either a desired result 
or garbage when the function returns. If you need to keep the inputs, pass a copy. */

class rsLinearAlgebraNew
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Convenience functions 
  These functions are convenient to use but not recommendable for production code - especially, 
  when that code should be run in realtime (they allocate memory). But for prototyping, they are 
  fine. */

  /** Returns the vector x that solves the linear system of equations A * x = b. */
  template<class T>
  static std::vector<T> solve(const rsMatrixView<T>& A, const std::vector<T>& b);
  // allocates
  // maybe rename to solveRegular because this function can't deal with singular matrices. It 
  // assumes A to be regular and will fail otherwise. A more general "solve" function should be 
  // able to produce solutions for soluble singular systems, i.e. those singular systems which have 
  // infinitely many solutions - in this case return a particular solution (the set of solutions
  // is then given by the particular solution plus any linear combination of vectors that span the
  // nullspace of A...right?)

  /** Returns the inverse of the given matrix A which is assumed to be a square matrix. */
  template<class T>
  static rsMatrix<T> inverse(const RAPT::rsMatrixView<T>& A);
  // allocates, todo: pseudoInverse

  /** Computes the determinant of A by Gaussian elemination. */
  template<class T>
  static T determinant(const rsMatrixView<T>& A);
  // allocates


  //-----------------------------------------------------------------------------------------------
  /** \name Solvers */

  /** Solves the system(s) of linear equations A * X = B. Often X and B are vectors, i.e. Nx1 
  matrices in which case lowercase letters are typically used: A * x = b. However, the function 
  works also for matrices X,B - this means, the system A * x = b is solved simultaneously for a 
  bunch of right-hand-side vectors b collected into a matrix whose columns are the desired solution
  vectors. A and B are the inputs, X is the output. X and B must have the same shape and both must 
  have a number of rows equal to the size A, which assumed to be square. The function works in 
  place, i.e. it does not allocate any extra memory. The process destroys the contents of A and B. 
  In general, X and B must be distinct but in certain cases, one may use it in place, i.e. X == B. 
  I'm not sure yet, but i think, it works, iff B is diagonal -> figure out.  */
  template<class T>
  static bool solve(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate, todo: document, when this may be used in place
  // maybe pass arguments as A,B,X - outputs should come last - check, how the old versions does 
  // it
  // document return value: it returns true, iff the linear system has a unique solution - it will
  // fail whenever there are no solutions or a continuum of solutions
  // maybe it should also return the rank?

  //-----------------------------------------------------------------------------------------------
  /** \name Subspaces */


  template<class T>
  static rsMatrix<T> rowSpace(rsMatrix<T> A, T tol);
  // allocates, needs tests

  /*
  template<class T>
  static rsMatrix<T> columnSpace(const RAPT::rsMatrixView<T>& A);

  template<class T>
  static rsMatrix<T> nullSpace(const RAPT::rsMatrixView<T>& A);

  template<class T>
  static rsMatrix<T> eigenSpace(const RAPT::rsMatrixView<T>& A, const T& eigenvalue);
  */




  //-----------------------------------------------------------------------------------------------
  /** \name Subroutines */

  /** Transforms the augmented coefficient matrix A|B for a set of linear systems of equations of 
  the form A * X = B into an upper triangular form, a.k.a. "row echelon form". In this form, the 
  corresponding linear system is easy to solve. It is a *set* of linear systems because B may 
  contain several right hand side vectors and the matrix X should contain just as many solutions. 
  The process used is Gaussian elimination with partial pivoting. */
  template<class T>
  static int makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // todo: document return value: it returns the ietration number, i.e. the number of row 
  // reductions that have been performed. for a NxN matrix, if this number is < N it means, the 
  // function encountered a zero pivot and returns early and indicates that A is singular, i think
  // the returned value is the rank of the matrix
  // doesn't allocate, maybe rename to rowElimination or rowEchelonForm, or rowEchelonRegular
  // return value is the number of iterations taken until no pivot could be found - is this the 
  // rank? i think so -> figure out, if yes, add the info to the documentation

  /** Like makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B) but with additional the output 
  parameter numSwaps that returns the number of row-swaps that occurred due to the pivoting. This
  information is required, when the function is used as subroutine in a function for computing
  determinants via Gaussian elimination. */
  template<class T>
  static int makeTriangular(rsMatrixView<T>& A, rsMatrixView<T>& B, int* numSwaps);

  /** Simplified version that doesn't use pivoting - this may fail even for non-singular 
  matrices, so it's not recommended for general use, but if you know that the elimination will 
  never encounter a zero diagonal element, you may use this version for optimization purposes. It 
  also has less requirements on the datatype T, so it may be used in cases where the datatype 
  doesn't meet all reuqirements needed for pivoting. It doesn't require comparison operators <,>
  or an rsAbs function defined for it - it does, however, require an == operator. */
  //template<class T>
  //static void makeTriangularNoPivot(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate - commented because it's useless for production code

  /** This produces the so called *reduced* row echelon form of a linear system of equations. In 
  this form, the coefficient matrix is diagonal and it is the simplemost representation of a 
  linear system of equations. In this form, the system is trivial to solve. Note that this has 
  *nothing* to do with another common process called "diagonalization". The latter refers to 
  finding the eigenvalues and -vectors and representing the coefficient matrix as product: 
  A = inv(M) * D * M where M is a matrix whose columns are the eigenvectors and D is a diagonal 
  matrix with the eigenvalues on the main diagonal - this is *NOT* what this function does. */
  template<class T>
  static int makeDiagonal(rsMatrixView<T>& A, rsMatrixView<T>& B);
  // doesn't allocate, todo: needs test, maybe renmae to reducedRowEchelonForm

  /** Solves the system(s) of linear equations A * X = B for the special case where A is an upper
  triangular matrix. */
  template<class T>
  static void solveTriangular(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // doesn't allocate, todo: document, when this may be used in place



  //-----------------------------------------------------------------------------------------------
  /** \name Deprecated */

  template<class T>
  static std::vector<T> solveOld(rsMatrix<T> A, std::vector<T> b);


};

#endif