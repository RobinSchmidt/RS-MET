template<class T>
size_t rsInterpolatingFunction<T>::addDataPoint(T x, T y)
{
  // maybe find insertion point and use vector::insert

  xValues.push_back(x);
  yValues.push_back(y);

  // keep arrays sorted according to ascending x:
  size_t i = xValues.size()-1;
  while(i > 0 && isNextValueLess(i-1))
  {
    rsSwap(xValues[i-1], xValues[i]);
    rsSwap(yValues[i-1], yValues[i]);
    i--;
  }
  return i;
}

template<class T>
void rsInterpolatingFunction<T>::removeDataPoint(size_t index)
{
  rsRemove(xValues, index);
  rsRemove(yValues, index);
}

template<class T>
size_t rsInterpolatingFunction<T>::moveDataPoint(size_t i, T newX, T newY)
{
  xValues[i] = newX;
  yValues[i] = newY;

  // keep arrays sorted according to ascending x:
  while(i > 0 && isNextValueLess(i-1))
  {
    rsSwap(xValues[i-1], xValues[i]);
    rsSwap(yValues[i-1], yValues[i]);
    i--;
  }
  while(i < xValues.size()-1 && isNextValueLess(i))
  {
    rsSwap(xValues[i+1], xValues[i]);
    rsSwap(yValues[i+1], yValues[i]);
    i++;
  }
  return i; // is this correct?
}