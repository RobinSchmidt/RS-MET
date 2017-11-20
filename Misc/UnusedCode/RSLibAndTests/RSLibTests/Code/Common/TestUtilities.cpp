#include "TestUtilities.h"

bool detectMemoryLeaks()
{
  #ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
  #else
  return false;
  #endif
}

rsVectorDbl rsLinearRangeVector(int N, double min, double max)
{
  double *values = new double[N];
  rsFillWithRangeLinear(values, N, min, max);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
}

rsVectorDbl rsExponentialRangeVector(int N, double min, double max)
{
  double *values = new double[N];
  rsFillWithRangeExponential(values, N, min, max);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
}

rsVectorDbl rsRandomVector(int N, double min, double max, int seed)
{
  double *values = new double[N];
  rsFillWithRandomValues(values, N, min, max, seed);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
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
