// third party includes:
#include "../_third_party/tnt/tnt.h"
#include "../_third_party/tnt/jama_svd.h"

//#include "rosic_Matrix.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

rosic::Matrix::Matrix()
{
  numRows    = 0;
  numColumns = 0;
  mFlat      = NULL;
  m          = NULL;
}

rosic::Matrix::Matrix(int numRows, int numColumns, bool initElementsWithZeros)
{
  this->numRows    = numRows;
  this->numColumns = numColumns;
  mFlat            = new double[numRows*numColumns];
  m                = new double*[numRows];
  for(int r = 0; r < numRows; r++)
    m[r] = &mFlat[r*numColumns];
  if( initElementsWithZeros == true )
    initWithZeros();
}

rosic::Matrix::Matrix(int numRows, int numColumns, double **values)
{
  this->numRows    = numRows;
  this->numColumns = numColumns;
  mFlat            = new double[numRows*numColumns];
  m                = new double*[numRows];
  for(int r=0; r<numRows; r++)
    m[r] = &mFlat[r*numColumns];

  for(int r=0; r<numRows; r++)
  {
    for(int c=0; c<numColumns; c++)
      m[r][c] = values[r][c];
  }
}

rosic::Matrix::Matrix(const Matrix& other)
{
  numRows    = other.numRows;
  numColumns = other.numColumns;
  mFlat      = new double[numRows*numColumns];
  m          = new double*[numRows];
  for(int r=0; r<numRows; r++)
    m[r] = &mFlat[r*numColumns];

  // copy the values:
  memcpy(mFlat, other.mFlat, numRows*numColumns*sizeof(double));
}

rosic::Matrix::~Matrix()
{
  delete[] m;
  delete[] mFlat;
}

//-------------------------------------------------------------------------------------------------
// setup:

//-------------------------------------------------------------------------------------------------
// inquiry:

//-------------------------------------------------------------------------------------------------
// matrix computations:

// copies the content of a rosic::Matrix into a corresponding object of class TNT::Array2D, the 
// latter of which must have the correct size (number of rows and columns) beforehand
void copyRosicMatrixToTntArray2D(const rosic::Matrix &rm, TNT::Array2D<double> &tm)
{
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++)
      tm[r][c] = rm.m[r][c];
  }
}

// copies the content of a TNT::Array2D into a corresponding object of class rosic::Matrix, the 
// size of the rosic::Matrix will be adjusted if it doesn't match
void copyTntArray2DToRosicMatrix(const TNT::Array2D<double> &tm, rosic::Matrix &rm)
{
  rm.setSize(tm.dim1(), tm.dim2());
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++) 
      rm.m[r][c] = tm[r][c];
  }
}

void rosic::Matrix::getSingularValueDecomposition(rosic::Matrix *U, rosic::Matrix *S, 
                                                  rosic::Matrix *V)
{
  TNT::Array2D<double> tntA(numRows, numColumns); copyRosicMatrixToTntArray2D(*this, tntA);
  JAMA::SVD<double> svd(tntA);
  svd.getU(tntA); copyTntArray2DToRosicMatrix(tntA, *U);
  svd.getS(tntA); copyTntArray2DToRosicMatrix(tntA, *S);
  svd.getV(tntA); copyTntArray2DToRosicMatrix(tntA, *V);
}

//-------------------------------------------------------------------------------------------------    
// others:

void rosic::Matrix::print()
{
  printf("%s %d %s %d %s", "Matrix - rows: ", numRows, " columns:", numColumns, "\n");
  for(int r=0; r<numRows; r++)
  {
    for(int c=0; c<numColumns; c++)
      printf("%.4f %s", m[r][c], "  ");
    printf("%s", "\n");
  }
}