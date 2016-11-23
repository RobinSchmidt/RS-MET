#include "MathDemos.h"

// to be removed, when available in RAPT:
template <class T>
void rsFillWithValue(T *buffer, int length, T value)
{
  for(int i = 0; i < length; i++)
    buffer[i] = value;
}
template <class T>
void rsFillWithRangeLinear(T *buffer, int length, T min, T max)
{
  if( min == max )
    rsFillWithValue(buffer, length, min);
  else
  {
    double factor = (max-min) / (double) (length-1);
    for(int i = 0; i < length; i++)
      buffer[i] = (T) (factor * T(i) + min);
  }
}
// ------------------------------------------

void sigmoids()
{
  // Plots various sigmoid saturation functions.

  typedef double Real;  // define the real number datatype to use
  int N = 501;          // number of datapoints
  Real xMin = -5;
  Real xMax = +5;

  // create x-axis:
  vector<Real> x(N);   // x-axis values
  rsFillWithRangeLinear(&x[0], N, xMin, xMax);

  // plot functions:
  GNUPlotter plt;
  typedef RAPT::rsNormalizedSigmoids Sigmoids;
  plt.addDataFunctions(N, &x[0], &Sigmoids::clip);
  plt.addDataFunctions(N, &x[0], &Sigmoids::atan);
  plt.addDataFunctions(N, &x[0], &Sigmoids::tanh);
  plt.addDataFunctions(N, &x[0], &Sigmoids::rational);
  plt.addDataFunctions(N, &x[0], &Sigmoids::cubic);
  plt.addDataFunctions(N, &x[0], &Sigmoids::cubicRational);
  plt.addDataFunctions(N, &x[0], &Sigmoids::quartic);
  plt.addDataFunctions(N, &x[0], &Sigmoids::hexic);
  plt.addDataFunctions(N, &x[0], &Sigmoids::softClipHexic);
  plt.plot();
}