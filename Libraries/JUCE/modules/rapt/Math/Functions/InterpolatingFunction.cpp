template<class T>
void rsInterpolatingFunction<T>::addDataPoint(T x, T y)
{
  xValues.push_back(x);
  yValues.push_back(y);

  // keep arrays sorted according to ascending x:
  size_t i = xValues.size()-1;
  while(i > 0 && xValues[i] < xValues[i-1])
  {
    rsSwap(xValues[i-1], xValues[i]);
    rsSwap(yValues[i-1], yValues[i]);
    i--;
  }
}

template<class T>
void rsInterpolatingFunction<T>::removeDataPoint(size_t index)
{

}

template<class T>
void rsInterpolatingFunction<T>::moveDataPoint(size_t i, T newX, T newY)
{
  xValues[i] = newX;
  yValues[i] = newY;

  // keep arrays sorted according to ascending x:
  while(i > 0 && xValues[i] < xValues[i-1])
  {
    rsSwap(xValues[i-1], xValues[i]);
    rsSwap(yValues[i-1], yValues[i]);
    i--;
  }
  while(i < xValues.size()-1 && xValues[i] > xValues[i+1])
  {
    rsSwap(xValues[i+1], xValues[i]);
    rsSwap(yValues[i+1], yValues[i]);
    i++;
  }
}