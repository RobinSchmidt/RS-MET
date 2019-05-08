#pragma once

/** Plotting functions - if the macro RS_PLOTTING is defined in you project, they will compile to
actual invocations of GNUPlotter and plot stuff. If it is not defined, they will compile to empty
dummy no-op functions that will get optimized away. With this mechanism, we can simply inject calls
to functions like rsPlotVector() anywhere in RAPT code for debugging purposes which will get 
optimized out in cases when they are not needed. In TestsRosicAndRapt.jucer it is defined, so 
plotting functions wil actually invoke the plotter in this project. */
// ...hmm...this is not actually true - but i think, i should make it true....


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

template<class T>
inline void rsStemPlot(int N, T *x, T *y)
{
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(N, x, y); // can probably be done without adding the data twice
  plt.setGraphStyles("impulses", "points pt 7 ps 1.2");
  plt.plot();
}

template<class T>
inline void rsStemPlot(std::vector<T> v)
{
  std::vector<T> x(v.size());
  RAPT::rsArray::fillWithIndex(&x[0], (int) x.size());
  rsStemPlot((int) x.size(), &x[0], &v[0]);
}


/** Plots a bunch of vectors. */
template<class T>
inline void rsPlotVectors(
  std::vector<T> v0, 
  std::vector<T> v1 = std::vector<T>(),
  std::vector<T> v2 = std::vector<T>(),
  std::vector<T> v3 = std::vector<T>(),
  std::vector<T> v4 = std::vector<T>(),  
  std::vector<T> v5 = std::vector<T>(),
  std::vector<T> v6 = std::vector<T>(),
  std::vector<T> v7 = std::vector<T>(),
  std::vector<T> v8 = std::vector<T>(), 
  std::vector<T> v9 = std::vector<T>()
  )
{
  // make a function that can take more vectors...maybe a vector of vectors?
  GNUPlotter plt;
  if(v0.size() > 0) plt.addDataArrays((int) v0.size(), &v0[0]);
  if(v1.size() > 0) plt.addDataArrays((int) v1.size(), &v1[0]);
  if(v2.size() > 0) plt.addDataArrays((int) v2.size(), &v2[0]);
  if(v3.size() > 0) plt.addDataArrays((int) v3.size(), &v3[0]);
  if(v4.size() > 0) plt.addDataArrays((int) v4.size(), &v4[0]);
  if(v5.size() > 0) plt.addDataArrays((int) v5.size(), &v5[0]);
  if(v6.size() > 0) plt.addDataArrays((int) v6.size(), &v6[0]);
  if(v7.size() > 0) plt.addDataArrays((int) v7.size(), &v7[0]);
  if(v8.size() > 0) plt.addDataArrays((int) v8.size(), &v8[0]);
  if(v9.size() > 0) plt.addDataArrays((int) v9.size(), &v9[0]);
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

/** Plots markers on the x-axis at the values given in the array. */
template<class T>
inline void rsPlotMarkers(T* markers, int numMarkers)
{
  std::vector<T> zeros(numMarkers);    // y values for plotting (all zero)
  RAPT::rsArray::fillWithZeros(&zeros[0], numMarkers);
  GNUPlotter plt;
  plt.addDataArrays(numMarkers,   markers, &zeros[0]);
  plt.setGraphStyles("points");
  plt.plot();
}

/** Plots a signal together with a set of markers (useful, for example, to mark time-instants of 
signal features such as zero-crossings) */
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

/** Plots the given spectral magnitude values as decibels. You may pass a sampleRate and 
floor/threshold for the dB values. If you pass no sample-rate (or zero), the frequency axis will 
show the bin-index, otherwise the physical frequency in Hz. The spectrum may also optionally be 
normalized such that the maximum value is at 0 dB. */
template<class T>
inline void rsPlotSpectrum(std::vector<T> fftMagnitudes, T sampleRate = T(0), 
  T floorDb = T(-200), bool normalize = false)
{
  int N = (int)fftMagnitudes.size();

  T scl = T(1);
  if(normalize == true) {
    T maxVal = T(0);
    for(int k = 0; k < N; k++)
      if(fftMagnitudes[k] > maxVal)
        maxVal = fftMagnitudes[k];
    scl = T(1) / maxVal;
  }

  std::vector<T> f(N), db(N);
  for(int k = 0; k < N; k++) {
    if(sampleRate > T(0))
      f[k] = k * sampleRate / (2*N); // 2 bcs we assume that we get an array of only positive freq bins
    else
      f[k] = T(k);
    db[k] = std::max(floorDb, 20*log10(scl*fftMagnitudes[k]));
  }

  GNUPlotter plt;
  plt.addDataArrays(N, &f[0], &db[0]);
  plt.plot();
}

// somewhat redundant with rsPlotSpectrum...
template<class T>
inline void rsPlotDecibels(int N, T* x, T *mag)
{
  T* dB = new T[N];
  for(int i = 0; i < N; i++)
    dB[i] = rsAmpToDbWithCheck(mag[i], 0.00000001);
  GNUPlotter plt;
  plt.addDataArrays(N, x, dB);
  plt.plot();
  delete[] dB;
}