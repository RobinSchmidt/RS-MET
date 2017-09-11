
template<class T>
rsMatrix<T>::rsMatrix(size_t numRows, size_t numColumns)
{
  N = numRows; 
  M = numColumns;
  data.resize(N*M);
  d = &data[0];
}




