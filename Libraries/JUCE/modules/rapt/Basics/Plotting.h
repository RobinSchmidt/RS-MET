#pragma once


// the plotting stuff doesn't work because it's all wrapped into the rapt namespace when we 
// declare it here maybe we need to make an extra file for all the plotting stuff and include
// it directly in rapt.h (instead of in Basics.h) in order to have the stuff in the global
// namespace

/** Plotting functions - if the macro RS_PLOTTING is defined in you project, they will compile to
actual invocations of GNUPlotter and plot stuff. If it is not defined, they will compile to empty
dummy functions. With this mechanism, we can simple inject calls to functions like plotVector()
anywhere in RAPT code for debugging purposes which will get optimized out in cases when they are 
not needed. */


// define RS_PLOTTING somewhere in your project and add GNUPlotter.cpp to it for compilation, if 
// you want to use plotting facilities - otherwise, calls to plotting functions will just be a 
// no-op by calling dummy functions:
#ifdef RS_PLOTTING
#include "../../../Tests/TestsRosicAndRapt/Source/Shared/Plotting/GNUPlotter.h"



template<class T>
inline void rsPlotVector(std::vector<T> v)
{
  int N = (int) v.size();
  double *y = new T[N];
  for(int n = 0; n < N; n++)
    y[n] = v[n];
  GNUPlotter plt;
  plt.plotArrays(N, y);
  delete[] y;
}


#else

// dummy plotting functions to satisfy the compiler when no plotting is desired:

inline void rsPlotVector(std::vector<double> v) {}


#endif



// old:
//// uncomment if you want to plot from rapt code (should be done only temporarily for debugging
//// sessions in the test project - trying to actually plot stuff will produce linker errors in other
//// projects):
//#define RS_DEBUG_PLOTTING
//#ifdef RS_DEBUG_PLOTTING
//#include "../../../Tests/TestsRosicAndRapt/Source/Shared/Plotting/GNUPlotter.h"
////#include "../../../Tests/TestsRosicAndRapt/Source/Shared/Plotting/Plotting.h"
//#endif
//// todo: make it work also in other projects (maybe i need to drag the GNUPlotCPP code into the
//// rapt module and conditionally compile it ...but then i will need to uncomment the define
//// wheneve i want to plot from the test project)
//// maybe define it in the jucer file for the test project - all other projects may then uncomment
//// it here but the test project doens't have to




