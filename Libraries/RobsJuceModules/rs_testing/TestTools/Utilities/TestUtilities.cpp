#include "TestUtilities.h"

/*
bool detectMemoryLeaks()
{
  #ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
  #else
  return false;
  #endif
}
*/

std::vector<double> rsLinearRangeVector(int N, double min, double max)
{
  std::vector<double> v(N);
  RAPT::rsArray::fillWithRangeLinear(&v[0], N, min, max);
  return v;
}

std::vector<double> rsExponentialRangeVector(int N, double min, double max)
{
  std::vector<double> v(N);
  RAPT::rsArray::fillWithRangeExponential(&v[0], N, min, max);
  return v;
}

std::vector<double> rsRandomVector(int N, double min, double max, int seed)
{
  std::vector<double> v(N);
  RAPT::rsArray::fillWithRandomValues(&v[0], N, min, max, seed);
  return v;
}

std::vector<double> rsApplyFunction(const std::vector<double>& v, double p, 
  double (*f) (double, double))
{
  std::vector<double> r(v.size());
  for(int i = 0; i < r.size(); i++)
    r[i] = f(v[i], p);
  return r;
}


#if defined(_MSC_VER)
// works only on MSVC and we need this function only in the performance tests, which we do only
// with the MSVC compiler anyway:
std::string toString(int n)
{
  return std::to_string((_Longlong)n);
}
#endif

/*
double rsSquare(double x)
{
  return x*x;
}
*/
