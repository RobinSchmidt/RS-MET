template<class T> bool rsGreater(const T& a, const T& b) { return a > b; }
template<class T> bool rsLess(const T& a, const T& b)    { return a < b; }


template<class T>
rsBuffer<T>::rsBuffer(size_t capacity)
{
  size_t c = RAPT::rsNextPowerOfTwo(capacity);
  data.resize(c);
  mask = c-1;
}

template<class T>
void rsBuffer<T>::initBufferValues(T value)
{
  RAPT::rsArrayTools::fillWithValue(&data[0], (int)data.size(), value);
}



template<class T>
void rsRingBuffer<T>::reset()
{
  this->initBufferValues(0);
  rightIndex = 0;
  updateLeftIndex();
}




/*
template<class T>
void rsDoubleEndedQueue<T>::clear()
{
  initBufferValues(0);
  head = 1;
  tail = 0;
}
*/

/*
template<class T>
rsDoubleEndedQueue<T>::rsDoubleEndedQueue(size_t capacity)
{
  //size_t c = RAPT::rsNextPowerOfTwo(capacity);
  data.resize(capacity);
  //mask = c-1;
}
*/
