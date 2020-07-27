template<class T> bool rsGreater(const T& a, const T& b) { return a > b; }
template<class T> bool rsLess(const T& a, const T& b)    { return a < b; }
// todo: remove these - or merge with the definition rsDefaultLess ins SortAndSearch and move to
// Basics


template<class T>
rsBuffer<T>::rsBuffer(size_t capacity)
{
  setCapacity(capacity);
}

template<class T>
void rsBuffer<T>::setCapacity(size_t newCapacity)
{
  size_t c = RAPT::rsNextPowerOfTwo(newCapacity);
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

  writeIndex = 0;
  // maybe zero is not best because in getSample, we increment first which means that the first
  // retrieved value after reset will be the value data[1]. it may be more convenient, when the
  // first returned value is data[0]...although, client code should actually not care at all - if 
  // it does, this points to a design flaw. for debugging, it would probably be most convenient, if
  // the first input sample gets written in data[0], so maybe we should init writeIndex to 
  // length-1...baaah - this pre-increment in getSample creates more problems than it solves!
  // i should really switch to post-increment - but that needs serious unit testing - it may be
  // that some code relies on the pre-increment (which really shouldn't be the case!). check the 
  // rsDoubleEndedQueue and rsMinMaxFilter

  adjustReadIndex();
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
