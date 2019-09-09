template<class T>
rsMatrixNew<T>::rsMatrixNew(int numRows, int numColumns)
{
  setSize(numRows, numColumns);
}

template<class T>
void rsMatrixNew<T>::setSize(int numRows, int numColumns)
{
  this->numRows = numRows;
  this->numCols = numColumns;
  data.resize(this->numRows * this->numCols);
  if(data.size() > 0)
    this->d = &data[0];
  else
    this->d = nullptr;
}

// fillWithValue
