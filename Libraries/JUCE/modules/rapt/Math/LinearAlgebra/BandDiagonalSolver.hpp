#pragma once

/** A convenience wrapper class for LAPACK's routines for solving band-diagonal systems of linear 
equations. It wraps the driver routines _gbsv, _gbsvx and _gbsvxx where the underscore is to be 
understood as a placeholdder for the datatype, i.e. s,d,c,z for single, double, single-complex or 
double-complex numbers. So far, it has been tested only for double precision real numbers...
todo: check, if it actually works for all these 4 datatypes (and maybe more)... */

template<class T>
class rsBandDiagonalSolver
{

public:

  /** Constructor. You have to pass the size of the system, i.e. the number of rows and columns of
  the coefficient matrix and the number of sub- and superdiagonals. Allocates memory for the matrix 
  and its factored form. Memory for all the other workspace arrays will be allocated in a call to 
  solve. */
  rsBandDiagonalSolver(int matrixSize, int numSubDiagonals, int numSuperDiagonals) 
  { setSystemSize(matrixSize, numSubDiagonals, numSuperDiagonals); }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Enumeration of the available algorithms. */
  enum class Algorithm
  {
    gbsv,    // LAPACK's simplest driver routine for band diagonal systems
    gbsvx,   // LAPACK's expert driver routine
    gbsvxx   // LAPACK's extended expert driver 
  };

  /** Selects the algorithm to be used for computing the solution. In effect, this switches between
  different lapack driver routines that have different levels of algorithmic sophistication. The 
  simplest choice is gbsv (simple driver), gbsvx is intermediate (expert driver) and the most 
  sophisticated is gbsvxx (extended expert driver). gbsvxx has the highest numerical precision. */
  void setAlgorithm(Algorithm algorithm) { algo = algorithm; }

  /** Sets up the various size parameters of the system to be solved. The matrixSize is the number
  of rows and columns of the matrix, the others are self-explanatory */
  void setSystemSize(int matrixSize, int numSubDiagonals, int numSuperDiagonals);

  /** Sets one of the values in one of the diagonals. The diagIndex indicates which diagonal is 
  meant where 0 is the main diagonal, -1 is the first subdiagonal, +1 the first superdiagonal, 
  -2 the second subdiagonal and so on. Within each diagonal, the elemIndex counts from the 
  top-left to the bottom right, where the main-diagonal has N elements, the first sub- and 
  superdiagonals N-1 elements and so on. */
  void setDiagonalElement(int diagIndex, int elemIndex, const T& value)
  {
    int i = diagElemIndex(diagIndex, elemIndex);
    if(i < 0 || i >= (int)A.size())
      error("invalid index in rsBandDiagonalSolver::setDiagonalElement");
    A[i] = value;
  }

  /** Sets matrix elements via regular dense-matrix row- and column indices. */
  void setElement(int rowIndex, int columnIndex, const T& value)
  {
    int i = rowColToArrayIndex(rowIndex, columnIndex);
    if(i < 0 || i >= (int)A.size())
      error("invalid index in rsBandDiagonalSolver::setElement");
    A[i] = value;
  }

  /** Sets the elements of the coefficient matrix in LAPACK's band-sotrage format. The passed 
  "newMatrix" array must be of size (kl+ku+1)*N where N is the size of the linear system, kl is the
  number of subdiagonals and ku the number of superdiagonals. These numbers are given by a call to
  the class constructor or a call to setSystemSize. For details on the band-storage  format, see 
  the comment at the bottom of the .cpp file. */
  void setCoefficientMatrix(T* newMatrix)
  {
    for(int i = 0; i < N*lda; i++)
      A[i] = newMatrix[i];
  }

  /** Selects whether or not equilibration should be used (if necessarry). Affects computations 
  only when gbsvx or gbsvxx is selected for the algorithm - gbsv never equilibrates. */
  void setEquilibration(bool shouldEquilibrate);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns information about the computations in the most recent call to solve(). A zero value 
  indicates successful exit. For the meaning of other values, refer to the documentation of "INFO" 
  parameter of the lapack routines gbvsv, gbsvx and gbsvxx. */
  long getInfo() const { return info; }

  /** Returns reciprocal condition number of the matrix in the most recent call to solve() where 
  gbsvx or gbsvxx has been used. */
  T getReciprocalConditionNumber() const { return rcond; }

  /** Returns reciprocal pivot growth factor in the most recent call to solve() where gbsvxx has 
  been used. */
  T getReciprocalPivotGrowth() const { return rpvgrw; }

  /** Returns the row scaling factors that have been used to equilibrate the matrix in the most 
  recent call to solve() where gbsvx or gbsvxx has been used. */
  std::vector<T> getRowScaleFactors() const { return R; }

  /** Returns the column scaling factors that have been used to equilibrate the matrix in the most 
  recent call to solve() where gbsvx or gbsvxx has been used. */
  std::vector<T> getColumnScaleFactors() const { return C; }

  //T getMaxErrorBound()
  // todo: error-bounds, etc.
  // getDiagonalElement, getElement, operator()(int i, int j)

  //-----------------------------------------------------------------------------------------------
  /** \name Computation */

  /** After the matrix has been set up via a call to setSystemSize and a bunch of calls to 
  setDiagonalElement, a call to solve will actually solve the system for a given number of right 
  hand sides and produce an equal number of solution vectors. */
  void solve(T* solutions, T* rightHandSides, int numRightHandSides);
  // todo: try to make rightHandSides const - figure out, if it's allowed that rightHandSides may
  // point to the same array as solutions
  // check, how we handle the ordering of solution and rhs in rapt - make it consistent

  //-----------------------------------------------------------------------------------------------
  /** \name Index conversion */

  /** Conversion for index of the diagonal (ranging from -kl to ku) and index of element within the 
  diagonal (ranging from 0 to N-abs(diagIndex)) to the flat array index for LAPACK's band storage
  format. Mainly for internal use in setDiagonalElement, but you can also use it yourself - but 
  only after setting up setSystemSize appropriately, of course. */
  inline int diagElemIndex(int diagIndex, int elemIndex)
  {
    int row = ku - diagIndex; 
    int col = elemIndex;
    if(diagIndex > 0) 
      col += diagIndex;
    return bandedRowColToIndex(row, col, kl, ku);
  }

  /** Converts from regular (dense) matrix indices for the row and column to the actual index in 
  the storage array. */
  inline int rowColToArrayIndex(int row, int col)
  {
    return bandedRowColToIndex(row+ku-col, col, kl, ku);
  }

  /** Converts from a row- and column index given in banded storage format to a flat array 
  index for the given number of subdiagonals (kl) and superdiagonals (ku). */
  static inline int bandedRowColToIndex(int row, int col, int numSubDiags, int numSuperDiags)
  {
    return (numSubDiags+numSuperDiags+1)*col + row;
  }


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Buffering */

  /** Allocates the storage memory for the matrix. Called from setSystemSize(). */
  void allocateMatrix();

  /** Allocates the temporary workspace buffers. Called from solve(). */
  void allocateBuffers();

  /** Prepares for a call to the simple driver routine gbsv by copying the matrix A into AF (with 
  appropriate offsets to leave space for the factored form) and B into X because gbsv replaces the
  right-hand side B with the solution X. */
  void prepareForGbsv(T* B, T* X);

  /** Error function - delegates to the blas function xerbla - but we don't want to call xerbla 
  here in the hpp file because calling blas/lapack/etc. functions here in the header screws up the 
  build-system for the RS-MET codebase. */
  void error(char* message);

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  // system size parameters:
  long N;       // matrix size is N x N
  long kl;      // number of subdiagonals
  long ku;      // number of superdiagonals
  long lda;     // leading dimension (ld) of matrix A in band storage format
  long ldab;    // leading dimension af matrix AB in gbsv, same as ld of AF in gbsvx and gbsvxx

  // right-hand-side size parameters:
  long nrhs;    // number of right hand sides
  long ldb;     // = N, redundant - maybe get rid

  // workspace buffers:
  std::vector<T> A;        // coefficient matrix in band-storage, size N*(kl+ku+1)
  std::vector<T> AF;       // factored form of A or matrix AB in gbsv, size N*(2*kl+ku+1)
  std::vector<T> work;     // workspace, size 4*N
  std::vector<T> R, C;     // row and column scale factors for equilibration, size N
  std::vector<long> iwork; // integer workspace, size N
  std::vector<long> ipiv;  // pivot indices, size N

  // options:
  Algorithm algo = Algorithm::gbsvxx; // extended expert driver - most precise
  char trans = 'N';        // input matrix is not transposed
  char fact  = 'E';        // input matrix is not factored and shall be equilibrated
  long nparams = 0;        // number of additional parameters for gbsvxx (not yet used)
  T params[1];             // dummy - not referenced, if nparams == 0

  // computed outputs:
  long info;                     // 0 on exit, if successful, see doc of gbsv/x/x for other values
  char equed;                    // returns the form of equlibration that was done
  T rcond;                       // reciprocal condition number
  T rpvgrw;                      // reciprocal pivot growth
  std::vector<T> ferr;           // componentwise relative forward error, size nrhs
  std::vector<T> berr;           // componentwise relative backward error, size nrhs
  long n_err_bnds = 3;           // number of error bounds 
  std::vector<T> err_bnds_norm;  // norm-wise error bounds, size n_err_bnds*nrhs
  std::vector<T> err_bnds_comp;  // component-wise error bounds, size n_err_bnds*nrhs
};