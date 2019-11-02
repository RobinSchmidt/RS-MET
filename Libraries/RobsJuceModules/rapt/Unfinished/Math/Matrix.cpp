#    
/*
// doesn't compile in VC - why? ...i have dragged it into the header file, for the time being
template<class T>
rsMatrixOld<T>::rsMatrixData<T>::rsMatrixData(int numRows, int numColumns)
{
  this->numRows    = numRows;
  this->numColumns = numColumns;
  allocateMemory();
  numReferences = 0;
}

template<class T>
rsMatrixOld<T>::rsMatrixData<T>::~rsMatrixData()
{
  if(numReferences == 0)
    freeMemory();
}

template<class T>
void rsMatrixOld<T>::rsMatrixData<T>::allocateMemory()
{
  mFlat = new T[numRows*numColumns];
  m     = new T*[numRows];
  for(int i = 0; i < numRows; i++)
    m[i] = &mFlat[i*numColumns];
}

template<class T>
void rsMatrixOld<T>::rsMatrixData<T>::freeMemory()
{
  delete[] m;
  delete[] mFlat;
}
*/

//=================================================================================================

// construction/destruction:

template<class T>
rsMatrixOld<T>::rsMatrixOld()
{
  data = new rsMatrixData<T>(0, 0);
  data->numReferences = 1;
}
// untested

template<class T>
rsMatrixOld<T>::rsMatrixOld(int numRows, int numColumns, bool initElementsWithZeros)
{
  data = new rsMatrixData<T>(numRows, numColumns);
  data->numReferences = 1;
  if( initElementsWithZeros == true )
    initWithZeros();
}
// untested

template<class T>
rsMatrixOld<T>::rsMatrixOld(int numRows, int numColumns, T **values)
{
  data = new rsMatrixData<T>(numRows, numColumns);
  data->numReferences = 1;
  for(int r = 0; r < numRows; r++)
  {
    for(int c = 0; c < numColumns; c++)
      data->m[r][c] = values[r][c];
  }
  // maybe use memcpy instead of nested loop
}
// untested

template<class T>
rsMatrixOld<T>::rsMatrixOld(const rsMatrixOld<T>& other)
{
  data = other.data;
  data->numReferences += 1;
}
// untested

template<class T>
rsMatrixOld<T>::~rsMatrixOld()
{
  data->numReferences -= 1;
  if( data->numReferences == 0 )
    delete data;
}
// untested


//-------------------------------------------------------------------------------------------------
// setup:


//-------------------------------------------------------------------------------------------------
// inquiry:



/*
//-------------------------------------------------------------------------------------------------
// TNT-based computations:

// copies the content of a rsMatrixOld into a corresponding object of class TNT::Array2D, the 
// latter of which must have the correct size (number of rows and columns) beforehand
void copyRosicMatrixToTntArray2D(const rsMatrixOld<T> &rm, TNT::Array2D<double> &tm)
{
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++)
      tm[r][c] = rm.m[r][c];
  }
}
// not tested

// copies the content of a TNT::Array2D into a corresponding object of class rsMatrixOld<T>, the 
// size of the rsMatrixOld<T> will be adjusted if it doesn't match
void copyTntArray2DToRosicMatrix(const TNT::Array2D<double> &tm, rsMatrixOld<T> &rm)
{
  rm.setSize(tm.dim1(), tm.dim2());
  for(int r=0; r<rm.numRows; r++)
  {
    for(int c=0; c<rm.numColumns; c++) 
      rm.m[r][c] = tm[r][c];
  }
}
// not tested

void rsMatrixOld<T>::getSingularValueDecomposition(rsMatrixOld<T> *U, rsMatrixOld<T> *S, rsMatrixOld<T> *V)
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
void rsMatrixOld<T>::makeDeepCopy()
{
  rsMatrixData<T>* copy = data->getDeepCopy();
  data->numReferences -= 1;
  if( data->numReferences == 0 )
    delete data;
  data = copy;
  data->numReferences = 1;
}

template<class T>
RS_INLINE void rsMatrixOld<T>::makeDeepCopyIfNecessary()
{
  if( data->numReferences > 1 )
    makeDeepCopy();
}

template<class T>
void rsMatrixOld<T>::updateDataPointer(rsMatrixData<T> *newData)
{
  newData->numReferences += 1;
  data   ->numReferences -= 1;
  if( data->numReferences == 0 )  
    delete data;
  data = newData;

  /*
  if( data != newData )
  {
    data->numReferences -= 1;
    if( data->numReferences == 0 )
      delete data;
    data = newData;
    data->numReferences += 1;
  }
  */
}

template<class T>
void rsMatrixOld<T>::print()
{
  printf("%s %d %s %d %s", "rsMatrixOld - rows: ", data->numRows, " columns:", data->numColumns, "\n");
  for(int r = 0; r < data->numRows; r++)
  {
    for(int c = 0; c < data->numColumns; c++)
      printf("%.4f %s", data->m[r][c], "  ");
    printf("%s", "\n");
  }
}
