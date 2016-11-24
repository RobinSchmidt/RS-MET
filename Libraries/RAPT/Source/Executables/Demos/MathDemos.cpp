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

void parametricBell()
{
  // Plots parametric bell functions (of various shapes) with adjustable center, width and flat-top 
  // width.

  static const int N = 1000;
  double center =  10.0;  // center of the bell
  double width  =  4.0;   // width
  double flat   =  0.5;   // relative length of flat top zone (between 0 and 1) 

  // create and set up parametric bell object:
  rsParametricBellFunction bell;
  bell.setCenter(center);
  bell.setWidth(width);
  bell.setFlatTopWidth(flat);

  // create x-axis and allocate y arrays:
  double xMin =  center - 1.2 * width/2;
  double xMax =  center + 1.2 * width/2;
  double x[N];
  rsFillWithRangeLinear(x, N, xMin, xMax);
  double yl[N], yc[N], yq[N], yh[N]; // linear, cubic, quintic, heptic
  int n;

  // create the family of curves (we look at different shapes for the prototype bell):
  bell.setPrototypeBell(&rsPositiveBellFunctions::linear);
  for(n = 0; n < N; n++)
    yl[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions::cubic);
  for(n = 0; n < N; n++)
    yc[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions::quintic);
  for(n = 0; n < N; n++)
    yq[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions::heptic);
  for(n = 0; n < N; n++)
    yh[n] = bell.getValue(x[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, x, yl, yc, yq, yh);
  plt.plot();
}

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
  typedef RAPT::NormalizedSigmoids<Real> Sigmoids;
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

