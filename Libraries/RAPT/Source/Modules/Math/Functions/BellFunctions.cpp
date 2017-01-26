//-------------------------------------------------------------------------------------------------
// positive range bell functions:

template<class T>
T rsPositiveBellFunctions<T>::linear(T x)
{
  if(x > 1)
    return 0;
  else
    return T(1) - x;
}

template<class T>
T rsPositiveBellFunctions<T>::cubic(T x)
{
  if(x > 1)
    return 0;
  else
    return 1 + (2*x - 3) * x*x;
}

template<class T>
T rsPositiveBellFunctions<T>::quintic(T x)
{
  if(x > 1)
    return 0;
  else
  {
    T x2 = x*x;
    return 1 + (-10 + 15*x - 6*x2) * x*x2; // 1 - 10*x^3 + 15*x^4 - 6*x^5
  }
}

template<class T>
T rsPositiveBellFunctions<T>::heptic(T x)
{
  if(x > 1)
    return 0;
  else
  {
    T x2 = x*x;
    return 1 + (-35 + 84*x - 70*x2 + 20*x2*x) * x2*x2; // 1 - 35*x^4 + 84*x^5 - 70*x^6 + 20*x^7
  }
}

//-------------------------------------------------------------------------------------------------
// class rsParametricBellFunction:

template<class T>
rsParametricBellFunction<T>::rsParametricBellFunction()
{
  bell   = &rsPositiveBellFunctions<T>::cubic;
  center = 0; 
  setWidth(2);
  setFlatTopWidth(0);
}

template<class T>
void rsParametricBellFunction<T>::setCenter(T newCenter)
{
  center = newCenter; 
}

template<class T>
void rsParametricBellFunction<T>::setWidth(T newWidth)
{
  a = 2 / newWidth;
}

template<class T>
void rsParametricBellFunction<T>::setFlatTopWidth(T newWidth)
{
  flat = newWidth;
  if(flat != 1)
    b = 1 / (1 - flat);
  else
    b = 1;
}

template<class T>
void rsParametricBellFunction<T>::setPrototypeBell(T (*newFunction)(T))
{
  bell = newFunction;
}
