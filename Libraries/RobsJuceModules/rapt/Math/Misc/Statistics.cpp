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
// todo: use pointers for output variables

template<class T>
T rsStatistics::proportionalRegression(int N, const T* x, const T* y)
{
  T xx = rsArrayTools::sumOfSquares(x, N);
  T xy = rsArrayTools::sumOfProducts(x, y, N);
  return xy / xx;
}

// todo: 
// -implement multiple linear regression (multiple x-arrays, one y-array)
//  linearRegression(int numInputs, int numDataPoints, const T** x, const T* y, T* a, T* b)
//  x is numInputs x numDataPoints matrix
// -based on that, implement univariate polynomial regression - it uses, 1,x,x^2,x^3,... as the 
//  x-arrays in multiple linear regression - but that should go to CurveFitting.h/cpp