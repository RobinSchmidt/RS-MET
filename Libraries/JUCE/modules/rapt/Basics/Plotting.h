#pragma once

/** Plotting functions - if the macro RS_PLOTTING is defined in you project, they will compile to
actual invocations of GNUPlotter and plot stuff. If it is not defined, they will compile to empty
dummy no-op functions that will get optimized away. With this mechanism, we can simply inject calls
to functions like rsPlotVector() anywhere in RAPT code for debugging purposes which will get 
optimized out in cases when they are not needed. In TestsRosicAndRapt.jucer it is defined, so 
plotting functions wil actually invoke the plotter in this project. */


#include "GNUPlotter.h"

template<class T>
inline void rsPlotVector(std::vector<T> v)
{
  GNUPlotter plt;
  plt.plotArrays((int) v.size(), &v[0]);
}

template<class T>
inline void plotSignalWithMarkers(T* signal, int signalLength, T* markers, int numMarkers)
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
