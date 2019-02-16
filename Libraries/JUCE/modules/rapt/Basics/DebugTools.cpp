#ifdef RS_PLOTTING
#include "../../../Tests/TestsRosicAndRapt/Source/Shared/Plotting/GNUPlotter.h"

template<class T>
void rsPlotVector(std::vector<T> v)
{
  int N = (int) v.size();
  double *y = new T[N];
  for(int n = 0; n < N; n++)
    y[n] = v[n];
  ::GNUPlotter plt;
  plt.plotArrays(N, y);
  delete[] y;
}


#else

// dummy function definitions:
template<class T> void rsPlotVector(std::vector<T> v) {}

#endif

template void rsPlotVector(std::vector<double> v);