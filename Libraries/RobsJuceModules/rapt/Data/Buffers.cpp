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
  RAPT::rsArray::fillWithValue(&data[0], (int)data.size(), value);
}