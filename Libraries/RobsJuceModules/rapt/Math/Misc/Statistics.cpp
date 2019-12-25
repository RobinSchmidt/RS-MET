template<class T>
void rsStatistics::linearRegression(int N, const T* x, const T* y, T& a, T& b)
{
  T xm = rsArrayTools::mean(x, N);
  T ym = rsArrayTools::mean(y, N);
  T xx = rsArrayTools::sumOfSquares(x, N);
  T xy = rsArrayTools::sumOfProducts(x, y, N);
  a = (xy - N*xm*ym) / (xx - N*xm*xm);
  b = ym - a*xm;
}

template<class T>
T rsStatistics::proportionalRegression(int N, const T* x, const T* y)
{
  T xx = rsArrayTools::sumOfSquares(x, N);
  T xy = rsArrayTools::sumOfProducts(x, y, N);
  return xy / xx;
}
