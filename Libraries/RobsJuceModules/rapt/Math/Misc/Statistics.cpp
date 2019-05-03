template<class T>
void rsStatistics::linearRegression(int N, T* x, T* y, T& a, T& b)
{
  T xm = rsArray::mean(x, N);
  T ym = rsArray::mean(y, N);
  T xx = rsArray::sumOfSquares(x, N);
  T xy = rsArray::sumOfProducts(x, y, N);
  a = (xy - N*xm*ym) / (xx - N*xm*xm);
  b = ym - a*xm;
}

template<class T>
T rsStatistics::proportionalRegression(int N, T* x, T* y)
{
  T xx = rsArray::sumOfSquares(x, N);
  T xy = rsArray::sumOfProducts(x, y, N);
  return xy / xx;
}
