template<class T>
rsMatrixNew<T>::rsMatrixNew(int numRows, int numColumns)
{
  setSize(numRows, numColumns);
}

template<class T>
rsMatrixNew<T>::rsMatrixNew(int numRows, int numColumns, const std::vector<T>& newData)
  : data(newData)
{
  rsAssert(numRows*numColumns == newData.size());
  this->numRows = numRows;
  this->numCols = numColumns;

  // factor out: updateDataPointer:
  if(data.size() > 0)  
    this->d = &data[0];
  else
    this->d = nullptr;
}

template<class T>
void rsMatrixNew<T>::setSize(int numRows, int numColumns)
{
  this->numRows = numRows;
  this->numCols = numColumns;
  data.resize(this->numRows * this->numCols);

  // factor out: updateDataPointer:
  if(data.size() > 0)  
    this->d = &data[0];
  else
    this->d = nullptr;
  // optionally initialize with zeros
}

// fillWithValue
