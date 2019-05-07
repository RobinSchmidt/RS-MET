

template<class T>
rsMovingMaximumFilter<T>::rsMovingMaximumFilter(size_t maxLength) 
  : delayLine(maxLength), maxDeque(maxLength)
{
  //greater = rsGreater;
  //greater = rsGreater<const T&, const T&>;
  //greater = std::function<bool(const T&, const T&)>(&rsGreater);
}