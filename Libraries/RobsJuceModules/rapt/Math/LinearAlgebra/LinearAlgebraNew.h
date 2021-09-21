#ifndef RAPT_LINEARALGEBRA_NEW_H_INCLUDED
#define RAPT_LINEARALGEBRA_NEW_H_INCLUDED

/** Collection of functions for linear algebra such as solving systems of linear equations, matrix 
inversion, etc. 

Many algorithms destroy the original input matrices and/or right-hand-sides in the course of their 
operation. If a parameter like a pointer to an array or a reference to an rsMatrixView is passed as
non-const, client code should always assume that the argument will contain either a desired result 
or garbage when the function returns. If you need to keep the inputs, pass a copy. 


ToDo:
-Maybe have independent types for the matrix elements and the vector elements. In some computer 
 graphics algorithms, the solution- and right-hand-side vectors are themselves composed of vectors 
 as their elements (whereas the matrix elements are still just scalars).

*/

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
  I'm not sure yet, but i think, it works, iff B is diagonal ...or maybe triangular is enough? 
  and/or if X and B are vectors?  ->  figure out and document.  */
  template<class T>
  static bool solve(rsMatrixView<T>& A, rsMatrixView<T>& X, rsMatrixView<T>& B);
  // -doesn't allocate, todo: document, when this may be used in place
  // -maybe pass arguments as A,B,X - outputs should come last - check, how the old versions does 
  //  it
  // -document return value: it returns true, iff the linear system has a unique solution - it will
  //  fail whenever there are no solutions or a continuum of solutions
  // -maybe it should also return the rank?
  // -maybe rename to solveGaussPartialPivot or ..RowPivot, implement also solveGaussFullPivot, 
  //  maybe keep the unqualified "solve" as alias for convenience..or maybe the unqualified solve 
  //  should take an additional parameter that selects the algorithm (which may default to Gaussian
  //  elimination with partial pivoting). maybe the algo parameter should be a bitfield: one 
  //  segment selects the core algo, another the preconditioner, yet another the incremental 
  //  refinement, etc.

  /** Solves the tridiagonal system of equations defined by a NxN matrix having the 3 nonzero 
  diagonals 'lowerDiag', 'mainDiag' and 'upperDiag' where the 'mainDiag' array should have N 
  elements and 'lowerDiag' and 'upperDiag' represent the the diagonals below  and above the main 
  diagonal respectively and should have N-1 elements. They are read from top-left to bottom-right. 
  The 'rightHandSide' argument represents the right hand side of the equation (vector with N 
  elements) and 'solution' ís where the solution vector will be stored (also N elements). The 
  mainDiag and rightHandSide vectors are destroyed in the process, the other diags are only for 
  read access. It can work in place, i.e. the solution may point to the same array as the RHS 
  vector, such that the RHS will be overwritten by the solution.
  
  In some sources, the algorithm used here is called "Thomas algorithm":
    https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
  while in others, the "Thomas algorithm" is defined in a slightly different way:
    https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  They are apparently both just variants of the same idea. We use the variant here that modifies 
  the main diagonal and leaves the upper and lower diagonals as is. (The other variant is 
  implemented somewhere in the prototypes. ToDo: maybe investigate the differences with regard to 
  numerical accuracy, stability and efficiency.). */
  template<class T>
  static void solveTridiagonal(int N, const T* lowerDiag, T* mainDiag, const T* upperDiag, 
    T* solution, T* rightHandSide);

  /** Solves a trididagonal system for 2 right hand side vectors b1,b2 in one go, storing the 
  corresponding solution vetors in x1,x2. It's intended to facilitate the application of the 
  Sherman-Morrison-Woodbury formula to problems with periodic boundary conditions. It's used by
  solveWrappedTridiagonal. It can work in place, i.e. x1 may be the same array as b1 and x2 may be 
  the same as b2. */
  template<class T>
  static void solveTridiagonal2Rhs(int N, const T* lowerDiag, T* mainDiag, const T* upperDiag, 
    T* x1, T* b1, T* x2, T* b2);

  /** Solves a perturbed tridiagonal system where the matrix is of the form (6x6 example):

     D0 U0 0  0  0  L0
     L1 D1 U1 0  0  0
     0  L2 D2 U2 0  0
     0  0  L3 D3 U3 0
     0  0  0  L4 D4 U4
     U5 0  0  0  L5 D5

  so, in the first and last row, some elements on the lower and upper diagonal are sort of wrapped
  around horizontally: to fit the pattern, L0 would actually belong to the left of D0 and U5 to the
  right of D5, but that would be beyond the matrix boundaries, so their positions are wrapped 
  around to the rightmost and leftmost position in the same row. Such systems occur in cubic spline
  interpolation problems with periodic boundary conditions (maybe also in ODEs?). 

  So that means that all 3 diagonals U,D,L (for upper,diag,lower) must have N elements (in contrast
  to the non-wrapped solveTridiagonal where U and L only have N-1 elements). The algorithm is based 
  on combining the regular tridiagonal ("Thomas") solver algorithm with the 
  Sherman-Morrison-Woodbury formula. In this case here, the solution and right-hand-side vectors 
  may not point to the same array, i.e. it cannot work in place. */
  template<class T>
  static void solveWrappedTridiagonal(int N, const T* lowerDiag, T* mainDiag, const T* upperDiag, 
    T* solution, T* rightHandSide);

  /** Like solveTridiagonal but solves the system for multiple right hand sides at ones, collected 
  as the columns of B. The solutions vectors are the corresponding columns of X. The size N of the
  system is inferred from the number of rows in B and X (which must match). */
  template<class T>
  static void solveTridiagonal(const T* lowerDiag, T* mainDiag, const T* upperDiag, 
    rsMatrixView<T>& X, rsMatrixView<T>& B);


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