#pragma once

/** Plotting functions - if the macro RS_PLOTTING is defined in you project, they will compile to
actual invocations of GNUPlotter and plot stuff. If it is not defined, they will compile to empty
dummy no-op functions that will get optimized away. With this mechanism, we can simply inject calls
to functions like rsPlotVector() anywhere in RAPT code for debugging purposes which will get 
optimized out in cases when they are not needed. In TestsRosicAndRapt.jucer it is defined, so 
plotting functions wil actually invoke the plotter in this project. */


#include "GNUPlotter.h"


template<class T>
inline void rsPlotArray(T* x, int N)
{
  GNUPlotter plt;
  plt.plotArrays(N, x);
}

template<class T>
inline void rsPlotArrays(int N, T* a1, T* a2 = nullptr, T* a3 = nullptr, T* a4 = nullptr,
  T* a5 = nullptr)
{
  GNUPlotter plt;
  plt.plotArrays(N, a1, a2, a3, a4, a5);
}
// maybe allow for more than 5

template<class T>
inline void rsPlotVector(std::vector<T> v)
{
  rsPlotArray(&v[0], (int) v.size());
  //GNUPlotter plt;
  //plt.plotArrays((int) v.size(), &v[0]);
}

/** Plots a bunch of vectors. */
template<class T>
inline void rsPlotVectors(
  std::vector<T> v1, 
  std::vector<T> v2 = std::vector<T>(),
  std::vector<T> v3 = std::vector<T>(),
  std::vector<T> v4 = std::vector<T>(),
  std::vector<T> v5 = std::vector<T>()  )
{
  // make a function that can take more vectors...maybe a vector of vectors?
  GNUPlotter plt;
  if(v1.size() > 0) plt.addDataArrays((int) v1.size(), &v1[0]);
  if(v2.size() > 0) plt.addDataArrays((int) v2.size(), &v2[0]);
  if(v3.size() > 0) plt.addDataArrays((int) v3.size(), &v3[0]);
  if(v4.size() > 0) plt.addDataArrays((int) v4.size(), &v4[0]);
  if(v5.size() > 0) plt.addDataArrays((int) v5.size(), &v5[0]);
  plt.plot();
}





/** Plots a whole bunch of vectors which are themselves put together into a vector of vectors. */
//template<class T>
//inline void rsPlotVectors(std::vector<std::vector<T>> v)
//{
//  GNUPlotter plt;
//  for(size_t i = 0; i < v.size; i++)
//    plt.addDataArrays((int) v[i].size(), &v[i][0]);
//  plt.plot;
//}

template<class T>
inline void rsPlotSignalWithMarkers(T* signal, int signalLength, T* markers, int numMarkers)
{
  std::vector<T> zeros(numMarkers);    // y values for plotting (all zero)
  RAPT::rsArray::fillWithZeros(&zeros[0], numMarkers);
  GNUPlotter plt;
  plt.addDataArrays(signalLength, signal);
  plt.addDataArrays(numMarkers,   markers, &zeros[0]);
  plt.setGraphStyles("lines", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
}

template<class T>
inline void rsPlotSpectrum(std::vector<T> fftMagnitudes, T sampleRate, 
  T floorDb = -std::numeric_limits<T>::infinity())
{
  int N = (int)fftMagnitudes.size();
  std::vector<T> f(N), db(N);
  for(int k = 0; k < N; k++) {
    f[k] = k * sampleRate / (2*N); // 2 bcs we assume that we get an array of only positive freq bins
    db[k] = std::max(floorDb, 20*log10(fftMagnitudes[k]));
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &f[0], &db[0]);
  plt.plot();
}
