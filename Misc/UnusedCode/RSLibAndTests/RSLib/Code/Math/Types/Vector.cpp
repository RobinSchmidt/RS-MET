/*
using namespace RSLib;

rsVectorDbl RSLib::rsLinearRangeVector(int N, double min, double max)
{
  double *values = new double[N];
  rsFillWithRangeLinear(values, N, min, max);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
}

rsVectorDbl RSLib::rsExponentialRangeVector(int N, double min, double max)
{
  double *values = new double[N];
  rsFillWithRangeExponential(values, N, min, max);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
}

rsVectorDbl RSLib::rsRandomVector(int N, double min, double max, int seed)
{
  double *values = new double[N];
  rsFillWithRandomValues(values, N, min, max, seed);
  rsVectorDbl v(N, values);
  delete[] values;
  return v;
}
*/
