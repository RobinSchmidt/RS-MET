// setup:

template<class T>
void rsBandDiagonalSolver<T>::setSystemSize(int matrixSize, int numSubDiagonals,
  int numSuperDiagonals)
{
  N    = matrixSize;
  kl   = numSubDiagonals;
  ku   = numSuperDiagonals;
  lda  = kl+ku+1;    // number of rows in A, maybe rename to M
  ldab = 2*kl+ku+1;  // number of rows in AF (factored version of A)
  ldb  = N;          // redundant - maybe get rid
  allocateMatrix();
}

template<class T>
void rsBandDiagonalSolver<T>::setEquilibration(bool shouldEquilibrate)
{
  if(shouldEquilibrate)
    fact = 'E';
  else
    fact = 'N';
}

// computation:

template<class T>
void rsBandDiagonalSolver<T>::solve(T* B, T* X, int numRightHandSides)
{
  nrhs = numRightHandSides;
  allocateBuffers();
  if(algo == Algorithm::gbsv) {
    prepareForGbsv(B, X);
    gbsv(&N, &kl, &ku, &nrhs, &AF[0], &ldab, &ipiv[0], &X[0], &ldb, &info);
  }
  else if(algo == Algorithm::gbsvx) {
    gbsvx(&fact, &trans, &N, &kl, &ku, &nrhs, &A[0], &lda, &AF[0], &ldab, &ipiv[0], &equed, &R[0],
      &C[0], &B[0], &ldb, &X[0], &ldb, &rcond, &ferr[0], &berr[0], &work[0], &iwork[0], &info,
      0, 0, 0); 
  }
  else if(algo == Algorithm::gbsvxx) {
    gbsvxx(&fact, &trans, &N, &kl, &ku, &nrhs, &A[0], &lda, &AF[0], &ldab, &ipiv[0], &equed, &R[0], 
      &C[0], &B[0], &ldb, &X[0], &ldb, &rcond, &rpvgrw, &berr[0], &n_err_bnds, &err_bnds_norm[0], 
      &err_bnds_comp[0], &nparams, &params[0], &work[0], &iwork[0], &info, 0, 0, 0); 
  }
  else
    xerbla("invalid algo in rsBandDiagonalSolver::solve", 0, 0);
}

// memory allocation:

template<class T>
void rsBandDiagonalSolver<T>::allocateMatrix()
{
  A.resize(lda*N);
  AF.resize(ldab*N);
}

template<class T>
void rsBandDiagonalSolver<T>::allocateBuffers()
{
  // todo: maybe re-allocate only those buffers that are needed for the selected routine - gbsv 
  // needs less buffers than gbsvxx, for example
                                          // used in routines (verify, if this is complete):
  ipiv.resize(N);                         // gbsv, gbsvx, gbsvxx 
  R.resize(N);                            // gbsvx, gbsvxx 
  C.resize(N);                            // gbsvx, gbsvxx 
  ferr.resize(nrhs);                      // gbsvx
  berr.resize(nrhs);                      // gbsvx, gbsvxx
  err_bnds_norm.resize(n_err_bnds*nrhs);  // gbsvxx
  err_bnds_comp.resize(n_err_bnds*nrhs);  // gbsvxx
  work.resize(4*N);                       // gbsvxx: 4*N, gbsvx: 3*N
  iwork.resize(N);                        // gbsvx, gbsvxx
}

template<class T>
void rsBandDiagonalSolver<T>::prepareForGbsv(T* B, T* X)
{
  for(int i = 0; i < N*nrhs; i++)        // copy B into X for being replaced by gbsv
    X[i] = B[i];
  for(int c = 0; c < N; c++)             // copy A into AF (leaving appropriate free spaces)...
    for(int r = 0; r < lda; r++)         // ...for being replaced by gbsv
      AF[c*ldab + kl+r] = A[c*lda + r];
}

/* Notes on the storage format for band-diagonal matrices in LAPACK:

Lapack generally uses flat arrays with column-major ordering for storing matrices, so if the 
row-index is "i" and the column index is "j", then the flat-array index for a matrix element
A(i,j) is given by k = j*numRows + i. Matrices in band storage format use the same number of 
columns as the original matrix, say N, and the number of rows is given by kl+ku+1 where kl is the
number of sub-diagonals, ku is the number of super-diagonals and the +1 is for the main diagonal.
For example, the 9x9 matrix (N=9) with kl=3 sub- and ku=2 superdiagonals:

  11 12 13 00 00 00 00 00 00
  21 22 23 24 00 00 00 00 00
  31 32 33 34 35 00 00 00 00
  41 42 43 44 45 46 00 00 00
  00 52 53 54 55 56 57 00 00
  00 00 63 64 65 66 67 68 00
  00 00 00 74 75 76 77 78 79
  00 00 00 00 85 86 87 88 89
  00 00 00 00 00 96 97 98 99

would look in band storage format:

  ** ** 13 24 35 46 57 68 79   upper diagonal 2
  ** 12 23 34 45 56 67 78 89   upper diagonal 1
  11 22 33 44 55 66 77 88 99   main diagonal
  21 32 43 54 65 76 87 98 **   lower diagonal 1
  31 42 53 64 75 86 97 ** **   lower diagonal 2
  41 52 63 74 85 96 ** ** **   lower diagonal 3

where the ** indicate unused memory locations that are not referenced. Taking into account the
column-major ordering, an array to represent that matrix in band storage format could be created
like:

  double _             = 0.0/sqrt(0.0);   // we init unused storage cells with NaN 
  static const int N   = 9;
  static const int kl  = 3;
  static const int ku  = 2;
  static const int lda = kl+ku+1;         // leading dimension of A (i.e. number of rows)
  double A[lda*N] =
  {  _, _,11,21,31,41,     // 1st column of banded format
     _,12,22,32,42,52,     // 2nd column
    13,23,33,43,53,63,
    24,34,44,54,64,74,
    35,45,55,65,75,85,
    46,56,66,76,86,96,
    57,67,77,87,97, _,
    68,78,88,98, _, _,
    79,89,99, _, _, _ };   // 9th column

This is the format that is used, for example, by gbmv (general banded matrix-vector multiply) and 
also for the input matrix in gbsvx and gbsvxx. The simple driver gbsv, on the other hand, expects
the matrix passed in a format that can hold the additional entries that are required to represent
the LU-factored form of the matrix which means, it needs kl additional free rows on the top. In 
this case, we need to store the matrix in the form:

  ** ** ** ** ** ++ ++ ++ ++   fields with ++ are used internally by the algorithm, fields
  ** ** ** ** ++ ++ ++ ++ ++   with ** are not accessed
  ** ** ** ++ ++ ++ ++ ++ ++
  ** ** 13 24 35 46 57 68 79   upper diagonal 2
  ** 12 23 34 45 56 67 78 89   upper diagonal 1
  11 22 33 44 55 66 77 88 99   main diagonal
  21 32 43 54 65 76 87 98 **   lower diagonal 1
  31 42 53 64 75 86 97 ** **   lower diagonal 2
  41 52 63 74 85 96 ** ** **   lower diagonal 3

which, again taking column-major ordering into account, looks for an array like:

  static const int ldab = 2*kl+ku+1;
  double AB[ldab*N] = 
  {  _, _, _, _, _,11,21,31,41,     // 1st column of banded format
     _, _, _, _,12,22,32,42,52,     // 2nd column
     _, _, _,13,23,33,43,53,63,
     _, _, _,24,34,44,54,64,74,
     _, _, _,35,45,55,65,75,85,
     _, _, _,46,56,66,76,86,96,
     _, _, _,57,67,77,87,97, _,
     _, _, _,68,78,88,98, _, _,
     _, _, _,79,89,99, _, _, _ };  // 9th column

This is the "AB" matrix that has to be passed to gbsv and on return will contain the factored form.
In gbsvx and gbsvxx, we pass the A matrix as above and additionally pass an uninitialized matrix
of the shape of AB which will on return also contain the factored form of A. The user of this class
has normally not to deal with this - all this stuff is taken care of by this wrapper class. */