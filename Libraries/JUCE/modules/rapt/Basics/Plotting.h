#pragma once

/** Plotting functions - if the macro RS_PLOTTING is defined in you project, they will compile to
actual invocations of GNUPlotter and plot stuff. If it is not defined, they will compile to empty
dummy no-op functions that will get optimized away. With this mechanism, we can simply inject calls
to functions like rsPlotVector() anywhere in RAPT code for debugging purposes which will get 
optimized out in cases when they are not needed. In TestsRosicAndRapt.jucer it is defined, so 
plotting functions wil actually invoke the plotter in this project. */


#include "../../../Tests/TestsRosicAndRapt/Source/Shared/Plotting/GNUPlotter.h"
// maybe the GNUPlotter code should sit somewhere inside the rapt directory....

template<class T>
inline void rsPlotVector(std::vector<T> v)
{
  GNUPlotter plt;
  plt.plotArrays((int) v.size(), &v[0]);
}


#ifdef RS_PLOTTING





#else

// dummy plotting functions to satisfy the compiler when no plotting is desired:

template<class T> inline void rsPlotVector(std::vector<T> v) {}

//inline void rsPlotVector(std::vector<double> v) {}

#endif
