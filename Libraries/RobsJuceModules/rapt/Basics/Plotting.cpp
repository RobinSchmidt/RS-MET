#include "GNUPlotter.cpp"
#include "../rapt.h"

template<class T>
void rsStemPlot(std::vector<T> v)
{
  std::vector<T> x(v.size());
  RAPT::rsArrayTools::fillWithIndex(&x[0], (int) x.size());
  rsStemPlot((int) x.size(), &x[0], &v[0]);
}
template void rsStemPlot(std::vector<double> v);

template<class T>
void rsPlotMarkers(T* markers, int numMarkers)
{
  std::vector<T> zeros(numMarkers);    // y values for plotting (all zero)
  RAPT::rsArrayTools::fillWithZeros(&zeros[0], numMarkers);
  GNUPlotter plt;
  plt.addDataArrays(numMarkers,   markers, &zeros[0]);
  plt.setGraphStyles("points");
  plt.plot();
}
template void rsPlotMarkers(double* markers, int numMarkers);

template<class T>
void rsPlotSignalWithMarkers(T* signal, int signalLength, T* markers, int numMarkers)
{
  std::vector<T> zeros(numMarkers);    // y values for plotting (all zero)
  RAPT::rsArrayTools::fillWithZeros(&zeros[0], numMarkers);
  GNUPlotter plt;
  plt.addDataArrays(signalLength, signal);
  plt.addDataArrays(numMarkers,   markers, &zeros[0]);
  plt.setGraphStyles("lines", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
}
template void rsPlotSignalWithMarkers(
  double* signal, int signalLength, double* markers, int numMarkers);


template<class T>
void rsPlotArraysXYWithMarks(const T* x, const T* y, int N, const std::vector<int>& marks)
{
  std::vector<T> xm = RAPT::rsSelect(x, marks);
  std::vector<T> ym = RAPT::rsSelect(y, marks);
  int M = (int) marks.size();
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(M, &xm[0], &ym[0]);
  plt.setGraphStyles("lines", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
}
template void rsPlotArraysXYWithMarks(const double* x, const double* y, int N, 
  const std::vector<int>& marks);



template<class T>
void rsPlotDecibels(int N, T* x, T *mag)
{
  T* dB = new T[N];
  for(int i = 0; i < N; i++)
    dB[i] = RAPT::rsAmpToDbWithCheck(mag[i], 0.00000001);
  GNUPlotter plt;
  plt.addDataArrays(N, x, dB);
  plt.plot();
  delete[] dB;
}
template void rsPlotDecibels(int N, double* x, double *mag);

