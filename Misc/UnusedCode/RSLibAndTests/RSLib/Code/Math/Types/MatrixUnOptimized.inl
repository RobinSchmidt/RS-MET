#ifndef RS_MATRIX_INL
#define RS_MATRIX_INL

using namespace RSLib;

// construction/destruction:

template<class T>
rsMatrix<T>::rsMatrix()
{
  numRows    = 0;
  numColumns = 0;
  mFlat      = NULL;
  m          = NULL;
}
// not tested

template<class T>
rsMatrix<T>::rsMatrix(int numRows, int numColumns, bool initElementsWithZeros)
{
  this->numRows    = numRows;
  this->numColumns = numColumns;
  mFlat            = new T[numRows*numColumns];
  m                = new T*[numRows];
  for(int r = 0; r < numRows; r++)
    m[r] = &mFlat[r*numColumns];
  if( initElementsWithZeros == true )
    initWithZeros();
}
// not tested

template<class T>
rsMatrix<T>::rsMatrix(int numRows, int numColumns, T **values)
{
  this->numRows    = numRows;
  this->numColumns = numColumns;
  mFlat            = new T[numRows*numColumns];
  m                = new T*[numRows];
  for(int r=0; r<numRows; r++)
    m[r] = &mFlat[r*numColumns];

  for(int r=0; r<numRows; r++)
  {
    for(int c=0; c<numColumns; c++)
      m[r][c] = values[r][c];
  }
}
// not tested

template<class T>
rsMatrix<T>::rsMatrix(const rsMatrix<T>& other)
{
  numRows    = other.numRows;
  numColumns = other.numColumns;
  mFlat      = new T[numRows*numColumns];
  m          = new T*[numRows];
  for(int r=0; r<numRows; r++)
    m[r] = &mFlat[r*numColumns];

  // copy the values:
  memcpy(mFlat, other.mFlat, numRows*numColumns*sizeof(T));
}
// not tested

template<class T>
rsMatrix<T>::~rsMatrix()
{
  delete[] m;
  delete[] mFlat;
}
// not tested

//-------------------------------------------------------------------------------------------------
// setup:

//-------------------------------------------------------------------------------------------------
// inquiry:



/*
//-------------------------------------------------------------------------------------------------
// TNT-based computations:

// copies the content of a rsMatrix into a corresponding object of class TNT::Array2D, the 
// latter of which must have the correct size (number of rows and columns) beforehand
void copyRosicMatrixToTntArray2D(const rsMatrix<T> &rm, TNT::Array2D<double> &tm)
{
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++)
      tm[r][c] = rm.m[r][c];
  }
}
// not tested

// copies the content of a TNT::Array2D into a corresponding object of class rsMatrix<T>, the 
// size of the rsMatrix<T> will be adjusted if it doesn't match
void copyTntArray2DToRosicMatrix(const TNT::Array2D<double> &tm, rsMatrix<T> &rm)
{
  rm.setSize(tm.dim1(), tm.dim2());
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++) 
      rm.m[r][c] = tm[r][c];
  }
}
// not tested

void rsMatrix<T>::getSingularValueDecomposition(rsMatrix<T> *U, rsMatrix<T> *S, rsMatrix<T> *V)
{
  TNT::Array2D<double> tntA(numRows, numColumns); copyRosicMatrixToTntArray2D(*this, tntA);
  JAMA::SVD<double> svd(tntA);
  svd.getU(tntA); copyTntArray2DToRosicMatrix(tntA, *U);
  svd.getS(tntA); copyTntArray2DToRosicMatrix(tntA, *S);
  svd.getV(tntA); copyTntArray2DToRosicMatrix(tntA, *V);
}
// not tested
*/

//-------------------------------------------------------------------------------------------------    
// others:

template<class T>
void rsMatrix<T>::print()
{
  printf("%s %d %s %d %s", "rsMatrix - rows: ", numRows, " columns:", numColumns, "\n");
  for(int r=0; r<numRows; r++)
  {
    for(int c=0; c<numColumns; c++)
      printf("%.4f %s", m[r][c], "  ");
    printf("%s", "\n");
  }
}

#endif