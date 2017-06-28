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
  float center = 10;    // center of the bell
  float width  = 4;     // width
  float flat   = 0.5;   // relative length of flat top zone (between 0 and 1) 

  // create and set up parametric bell object:
  rsParametricBellFunction<float> bell;
  bell.setCenter(center);
  bell.setWidth(width);
  bell.setFlatTopWidth(flat);

  // create x-axis and allocate y arrays:
  float xMin =  center - 1.2f * width/2;
  float xMax =  center + 1.2f * width/2;
  float x[N];
  rsFillWithRangeLinear(x, N, xMin, xMax);
  float yl[N], yc[N], yq[N], yh[N]; // linear, cubic, quintic, heptic
  int n;

  // create the family of curves (we look at different shapes for the prototype bell):
  bell.setPrototypeBell(&rsPositiveBellFunctions<float>::linear);
  for(n = 0; n < N; n++)
    yl[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<float>::cubic);
  for(n = 0; n < N; n++)
    yc[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<float>::quintic);
  for(n = 0; n < N; n++)
    yq[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<float>::heptic);
  for(n = 0; n < N; n++)
    yh[n] = bell.getValue(x[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, x, yl, yc, yq, yh);
  plt.plot();
}

void scaledAndShiftedSigmoid()
{
  // Plots a scaled and shifted sigmoid function

  static const int N = 500;   // number of datapoints
  float center = 10;          // center of the sigmoid
  float width  = 4;           // width

  // set up sigmoid:
  ScaledAndShiftedSigmoid<float> sigmoid;
  sigmoid.setCenter(center);
  sigmoid.setWidth(width);
  //sigmoid.setPrototypeSigmoid(&RAPT::NormalizedSigmoids<float>::softClipHexic);

  // create data:
  float x[N], y[N];
  float xMin =  center - 3 * width/2;
  float xMax =  center + 3 * width/2;
  rsFillWithRangeLinear(x, N, xMin, xMax);
  for(int n = 0; n < N; n++)
    y[n] = sigmoid.getValue(x[n]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.plot();
}

void sigmoids()
{
  // Plots various sigmoid saturation functions.

  int N = 501;          // number of datapoints
  float xMin = -5;
  float xMax = +5;

  // create x-axis:
  vector<float> x(N);   // x-axis values
  rsFillWithRangeLinear(&x[0], N, xMin, xMax);

  // plot functions:
  GNUPlotter plt;
  typedef RAPT::NormalizedSigmoids<float> Sigmoids;
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

