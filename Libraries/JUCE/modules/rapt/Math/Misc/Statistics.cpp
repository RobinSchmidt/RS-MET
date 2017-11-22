template<class T>
void rsStatistics::linearRegression(int N, T* x, T* y, T& a, T& b)
{
  T xm = rsArray::rsMean(x, N);
  T ym = rsArray::rsMean(y, N);
  T xx = rsArray::rsSumOfSquares(x, N);
  T xy = rsArray::rsSumOfProducts(x, y, N);
  a = (xy - N*xm*ym) / (xx - N*xm*xm);
  b = ym - a*xm;
}

template<class T>
T rsStatistics::proportionalRegression(int N, T* x, T* y)
{
  T xx = rsArray::rsSumOfSquares(x, N);
  T xy = rsArray::rsSumOfProducts(x, y, N);
  return xy / xx;
}
