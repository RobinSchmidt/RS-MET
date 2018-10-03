//#include "rosic_LinearAlgebra.h"
//using namespace rosic;

bool rosic::solveLinearSystem(double** A, double* x, double* b, int N)
{
  int i, j;

  // allocate memory for temporary copies of the coefficient matrix A and target vector b:
  double*  tmpB  = new double[N];
  double*  tmpA  = new double[N*N];
  double** tmpAP = new double*[N];

  // assign the pointer array for the matrix (tmpAP) to to the beginnings of the rows:
  for(i=0; i<N; i++)
    tmpAP[i] = &(tmpA[i*N]);

  // copy the data into the temporary arrays:
  for(i=0; i<N; i++)
  {
    tmpB[i] = b[i];
    for(j=0; j<N; j++)
      tmpAP[i][j] = A[i][j];
  }

  // solve the linear system in place with the temporary arrays:
  bool success = solveLinearSystemInPlace(tmpAP, x, tmpB, N);

  // free allocated memory:
  delete[] tmpB;
  delete[] tmpA;
  delete[] tmpAP;

  return success;
}

bool rosic::solveLinearSystemInPlace(double** A, double* x, double* b, int N)
{
  bool        matrixIsSingular = false;
  int         i, j, k, p;
  double      biggest, multiplier;
  long double tmpSum;

  // outermost loop over the rows to be scaled and subtracted from the rows below them:
  for(i=0; i<N; i++) 
  {

    // search for largest pivot in the i-th column from the i-th row downward:
    p       = i;
    biggest = 0.0;
    for(j=i; j<N; j++)
    {
      if( fabs(A[j][i]) > biggest )
      {
        biggest = fabs(A[j][i]);
        p       = j;
      }
    }
    if( RAPT::rsIsCloseTo(biggest, 0.0, 1.e-12) )
    {
      matrixIsSingular = true;
      DEBUG_BREAK; 
      // uh oh: the coefficient-matrix seems to be (close to) singular - Gaussian elimination
      // is likely to return meaningless results.
    }

    // swap the current row with pivot row (a pointer switch in the first index of A and an 
    // exchange of two values in b):
    if( p != i )
    {
      RAPT::rsSwap(A[i], A[p]);
      RAPT::rsSwap(b[i], b[p]);
      p = i;
    }

    // subtract a scaled version of the pivot-row p (==i) from all rows below it (j=i+1...N-1), 
    // the multiplier being the i-th column element of the current row divided by the i-th column
    // element of the pivot-row - but we do the subtraction only for those columns which were not 
    // already zeroed out by a previous iteration, that is: only from the i-th column rightward 
    // (k=i...N-1):
    for(j=i+1; j<N; j++)
    {
      multiplier = A[j][i] / A[p][i];
      b[j] -= multiplier * b[p];
      for(k=i; k<N; k++)
        A[j][k] -= multiplier * A[p][k];
    }

  }

  // the matrix A is now in upper triangular form - we now solve for the unknowns x[0...N-1] by 
  // backsubstitution:  
  x[N-1] = b[N-1] / A[N-1][N-1];
  for(i=N-2; i>=0; i--)
  {
    // multiply the already obtained x-values by their coefficients from the current row and 
    // accumulate them:
    tmpSum = 0.0;
    for(j=i+1; j<N; j++)
      tmpSum += A[i][j] * x[j];

    // this accumulated sum is to be subtracted from the target value, and the result of that 
    // subtraction must be divided diagonal-element which corresponds to the index of the new 
    // x-value which is to be computed:
    x[i] = (b[i] - tmpSum) / A[i][i];
  }

  // done: the vector x now contains the solution to the system A*x=b (unless the matrix was 
  // singular in which case it contains meaningless numbers or not-a-numbers)

  return !matrixIsSingular;
}

bool rosic::invertMatrix(double** A, int N)
{
  bool   matrixIsSingular = false;

  int    i, j, k, p;
  double biggest, multiplier;

  double*  tmpA  = new double[N*N];
  double** tmpAP = new double*[N];

  // assign the pointer array for the matrix (tmpAP) to to the beginnings of the rows:
  for(i=0; i<N; i++)
    tmpAP[i] = &(tmpA[i*N]);

  // copy the data from A into the temporary matrix and initialize the matrix A with the unit 
  // matrix:
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
    {
      tmpAP[i][j] = A[i][j];
      A[i][j]     = 0.0;
    }
    A[i][i] = 1.0;
  }

  // outermost loop over the rows to be scaled and subtracted from the rows below them:
  for(i=0; i<N; i++) 
  {

    // search for largest pivot in the i-th column from the i-th row downward:
    p       = i;
    biggest = 0.0;
    for(j=i; j<N; j++)
    {
      if( fabs(tmpAP[j][i]) > biggest )
      {
        biggest = fabs(tmpAP[j][i]);
        p       = j;
      }
    }
    if( RAPT::rsIsCloseTo(biggest, 0.0, 1.e-12) )
    {
      matrixIsSingular = true;
      DEBUG_BREAK; 
      // uh oh: the coefficient-matrix seems to be (close to) singular - Gaussian elimination
      // is likely to return meaningless results.
    }

    // swap the current row with pivot row (a pointer switch in the first index of A and an 
    // exchange of two values in b):
    if( p != i )
    {
      RAPT::rsSwap(tmpAP[i], tmpAP[p]);
      RAPT::rsSwap(A[i],     A[p]);
      p = i;
    }

    // divide the pivot-row row by the pivot element to get a 1 on the diagonal:
    multiplier = 1.0 / tmpAP[i][i];
    for(k=0; k<N; k++)
    {
      tmpAP[p][k] *= multiplier;
      A[p][k]     *= multiplier;
    }

    // subtract a properly scaled version of the pivot-row from all other rows to get zeros in this
    // column (this is different from Gaussian eleimination where it is subtracted only from the 
    // rows below):
    for(j=0; j<N; j++)
    {
      multiplier = tmpAP[j][i];  
      if( j != i )
      {
        for(k=0; k<N; k++)
        {
          tmpAP[j][k] -= multiplier * tmpAP[p][k];
          A[j][k]     -= multiplier * A[p][k];
        }
      }
    }

  }

  // free temporarily allocated memory:
  delete[] tmpA;
  delete[] tmpAP;

  return !matrixIsSingular;
}

bool rosic::solveTriDiagonalSystem(double *lower, double *main, double *upper, double *rhs, 
                                   double *solution, int N)
{
	if( main[0] == 0.0 )
  {
    DEBUG_BREAK;    // division by zero
    return false;
  }

	double divisor = main[0];   
  double *tmp    = new double[N];

	solution[0] = rhs[0] / divisor;
	for(int n=1; n<N; n++) 
  {
		tmp[n]  = upper[n-1] / divisor;
		divisor = main[n] - lower[n-1]*tmp[n];
		if( divisor == 0.0 )	
    {
      DEBUG_BREAK;    // division by zero
      delete[] tmp;
      return false;
    }
		solution[n] = (rhs[n]-lower[n-1]*solution[n-1]) / divisor;
	}

	for(int n=N-2; n>=0; n--)
		solution[n] -= tmp[n+1] * solution[n+1];

  delete[] tmp;
  return true;
}