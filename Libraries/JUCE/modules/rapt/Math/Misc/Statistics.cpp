template<class T>
void rsStatistics::linearRegression(int N, T* x, T* y, T& a, T& b)
{
  T xm = ArrayTools::rsMean(x, N);
  T ym = ArrayTools::rsMean(y, N);
  T xx = ArrayTools::rsSumOfSquares(x, N);
  T xy = ArrayTools::rsSumOfProducts(x, y, N);
  a = (xy - N*xm*ym) / (xx - N*xm*xm);
  b = ym - a*xm;
}

template<class T>
T rsStatistics::proportionalRegression(int N, T* x, T* y)
{
  T xx = ArrayTools::rsSumOfSquares(x, N);
  T xy = ArrayTools::rsSumOfProducts(x, y, N);
  return xy / xx;
}
